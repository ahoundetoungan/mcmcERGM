#' @title Simulate the posterior distribution of the parameters of an exponential random symmetric graph model
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doRNG "%dorng%"
#' @importFrom doParallel registerDoParallel
#' @importFrom ddpcr quiet
#' @importFrom utils head tail
#' @export
mcmcSymNet <- function(network, 
                       formula.u, 
                       formula.w, 
                       f_to_var.u, 
                       f_to_var.w,
                       starting    = NULL,
                       nblock      = 2,
                       ncores      = 1, 
                       Etheta      = NULL,
                       invVtheta   = NULL,
                       Sigma       = NULL, 
                       simutheta   = 1e3,
                       simunet     = 1e3,
                       mcmc.ctr    = list(target  = NULL, kappa = 0.6, jmin = 1e-6, jmax = 3), 
                       data        = NULL){
  stopifnot(is.list(network))
  if(missing(formula.u) & missing(formula.w)){
    stop("both formula.u and formula.w cannot be missing")
  }
  if(missing(formula.u)){
    if(!missing(f_to_var.u)){stop("f_to_var.u is defined without formula.u")}
  } else{
    if(missing(f_to_var.u)){stop("formula.u is defined without f_to_var.u")}
    stopifnot(tolower(unlist(f_to_var.u)) %in% c("sum", "prod", "same", "adiff"))
    f_to_var.u  <- ff_to_var(f_to_var.u)
  }
  if(missing(formula.w)){
    if(!missing(f_to_var.w)){stop("f_to_var.w is defined without formula.w")}
  } else{
    if(missing(f_to_var.w)){stop("formula.w is defined without f_to_var.w")}
    stopifnot(tolower(unlist(f_to_var.w)) %in% c("sum", "prod", "same", "adiff"))
    f_to_var.w  <- ff_to_var(f_to_var.w)
  }
  
  M           <- length(network)
  nvec        <- sapply(network, nrow)
  nveccum     <- c(0, cumsum(nvec))
  
  # data
  Xu          <- NULL
  Xw          <- NULL
  inter.u     <- FALSE
  inter.w     <- FALSE
  name.u      <- NULL
  name.w      <- NULL
  vname.u     <- NULL
  vname.w     <- NULL
  Ku          <- NULL
  Kw          <- NULL
  nvaru       <- 0
  nvarw       <- 0
  if(!missing(formula.u)){
    Xu        <- formula.to.data(formula.u, data)
    inter.u   <- Xu$intercept
    name.u    <- Xu$names
    Xu        <- Xu$X
    Ku        <- ncol(Xu)
    if(length(f_to_var.u) != Ku) stop("length(f_to_var.u) does not match with formula.u")
    nvaru     <- length(unlist(f_to_var.u))
    Xu        <- lapply(1:M, function(m) Xu[(nveccum[m] + 1):nveccum[m + 1],, drop = FALSE])
    vname.u   <- fvarnames(f_to_var.u, name.u, Ku, inter.u, "D")
  }
  if(!missing(formula.w)){
    Xw        <- formula.to.data(formula.w, data)
    inter.w   <- Xw$intercept
    name.w    <- Xw$names
    Xw        <- Xw$X
    Kw        <- ncol(Xw)
    if(length(f_to_var.w) != Kw) stop("length(f_to_var.w) does not match with formula.w")
    nvarw     <- length(unlist(f_to_var.w))
    Xw        <- lapply(1:M, function(m) Xw[(nveccum[m] + 1):nveccum[m + 1],, drop = FALSE])
    vname.w   <- fvarnames(f_to_var.w, name.w, Kw, inter.w, "I")
  }
  
  # Construct cluster
  cl          <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  if(!is.null(Xu)){
    Xu        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xu[[m]], f_to_var.u, nvaru, Ku)}
  }
  if(!is.null(Xw)){
    Xw        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xw[[m]], f_to_var.w, nvarw, Kw)}
  }
  
  # number of parameters
  nparu       <- length(vname.u)
  nparw       <- length(vname.w)
  npar        <- nparu + nparw
  
  # starting
  if(is.null(starting)){
    starting  <- rep(0, npar)
  } else{
    if(length(starting) != npar){
      stop("length(starting) not equal to the number of parameters to be estimated")
    }
  }
  theta       <- starting
  
  # Mcmc control
  target      <- mcmc.ctr$target
  kappa       <- mcmc.ctr$kappa
  jmin        <- mcmc.ctr$jmin
  jmax        <- mcmc.ctr$jmax
  
  if(is.null(target)){
    target    <- ifelse(npar == 1, 0.44, 
                        ifelse(npar == 2, 0.35, 
                               ifelse(npar == 3, 0.3125, 
                                      ifelse(npar == 4, 0.275, 
                                             ifelse(npar == 5, 0.25, 0.234)))))
  }
  if(is.null(kappa)){
    kappa     <- 0.6
  }
  if(is.null(jmin)){
    jmin      <- 1e-6
  }
  if(is.null(jmax)){
    jmax      <- 3
  }
  
  # prior
  # expectation
  if(is.null(Etheta)){
    Etheta    <- rep(0, npar)
  } else{
    if(length(Etheta) != npar){
      stop("length(Etheta) not equal to the number of parameters to be estimated")
    }
  }
  # inverse of the variance
  if(is.null(invVtheta)){
    invVtheta <- matrix(0, npar, npar)
  } else{
    if(length(invVtheta) == 1){invVtheta <- diag(npar)*c(invVtheta)}
    if((nrow(invVtheta) != npar) | (ncol(invVtheta) != npar)){
      stop("dim(invVtheta) not match to the number of parameters to be estimated")
    }
  }
  
  #Sigma
  if(is.null(Sigma)){
    Sigma     <- diag(npar)
  } else{
    if(length(Sigma) == 1){Sigma <- diag(npar)*c(Sigma)}
    if((nrow(Sigma) != npar) | (ncol(Sigma) != npar)){
      stop("dim(Sigma) not match to the number of parameters to be estimated")
    }
  }
  
  # Utility
  uu          <- rep(list(matrix(0, 1, 1)), times = M)
  uw          <- rep(list(matrix(0, 1, 1)), times = M)
  uu          <- futil(M, Xu, uu, theta[1:nparu], nparu, nvec, inter.u)
  uw          <- futil(M, Xw, uw, theta[(nparu + 1):npar], nparw, nvec, inter.w)
  
  # The potential function value at the starting point
  pot_ath     <- fpotensym(M, network, uu, uw, nparu, nparw, nvec)
  
  # MCMC
  jstheta     <- 1  #jumping scale for the adaptive MCMC
  atheta      <- 0  #number of times the draw is accepted
  uup         <- uu
  uwp         <- uw
  combr       <- ffindcom(nblock) 
  ncombr      <- 2^nblock
  ztncombr    <- 0:(ncombr - 1)
  idrows      <- fIDsym(M, nvec)
  idcols      <- idrows$idcols
  ident       <- idrows$ident
  idrows      <- idrows$idrows
  quiet(gc())
  
  posterior   <- as.data.frame(matrix(0, simutheta, npar)); colnames(posterior) <- c(vname.u, vname.w)
  spot        <- c()
  for (s in 1:simutheta) {
    # simulate theta prime
    cat("********************* iteration: ", s, "/", simutheta, "\n", sep = "")
    if(nparu > 0){
      cat("direct links\n")
      cat(head(theta, nparu), "\n")
    }
    if(nparw > 0){
      cat("indirect links\n")
      cat(tail(theta, nparw), "\n")
    }
    cat("potential\n")
    cat(pot_ath, "\n")
    
    thetap    <- fsimtheta(theta, Sigma, jstheta, npar)  #**
    
    # utility at thetap
    uup       <- futil(M, Xu, uup, thetap[1:nparu], nparu, nvec, inter.u) #**
    uwp       <- futil(M, Xw, uwp, thetap[(nparu + 1):npar], nparw, nvec, inter.w) #**
    
    # potential function at a (ie, network) and thetap
    pot_athp  <- fpotensym(M, network, uup, uwp, nparu, nparw, nvec) #**
    
    # Gibbs to simulate ap
    networkp  <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
      fGibbsym(network[[m]], nblock, ncombr, combr, idrows[[m]], idcols[[m]], ident[[m]], ztncombr, 
               uup[[m]], uwp[[m]], nparu, nparw, nvec[[m]], simunet)}
    cat("Gibbs executed\n")
    
    # potential function at ap 
    pot_apth  <- fpotensym(M, networkp, uu, uw, nparu, nparw, nvec)
    pot_apthp <- fpotensym(M, networkp, uup, uwp, nparu, nparw, nvec)
    
    # acceptance rate of theta
    lalpha    <- pot_apth - pot_ath + pot_athp - pot_apthp  + propdnorm(thetap, Etheta, invVtheta) - propdnorm(theta, Etheta, invVtheta) 
    lalpha    <- min(0, lalpha)
    if(runif(1) < exp(lalpha)){
      theta   <- thetap
      uu      <- uup
      uw      <- uwp
      pot_ath <- pot_athp
      atheta  <- atheta + 1
      cat("proposal accepted -- acceptance rate: ", round(100*atheta/s, 1), "%\n", sep = "")
    } else{
      cat("proposal rejected -- acceptance rate: ", round(100*atheta/s, 1), "%\n", sep = "")
    }
    jstheta       <- fupdate_js(jstheta, atheta, s, target, kappa, jmin, jmax)
    posterior[s,] <- c(theta)
    spot[s]       <- pot_ath
  }
  out  <- list(posterior = posterior, potential = spot); class(out) <- "mcmcSymNet"
  out
}


#' @title Simulate the posterior distribution of the parameters of an exponential random symmetric graph model
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doRNG "%dorng%"
#' @importFrom doParallel registerDoParallel
#' @importFrom ddpcr quiet
#' @importFrom utils head tail
#' @export
mcmcDirNet <- function(network, 
                       formula.u, 
                       formula.v,
                       formula.w, 
                       f_to_var.u, 
                       f_to_var.v,
                       f_to_var.w,
                       starting    = NULL,
                       nblock      = 2,
                       ncores      = 1, 
                       Etheta      = NULL,
                       invVtheta   = NULL,
                       Sigma       = NULL, 
                       simutheta   = 1e3,
                       simunet     = 1e3,
                       mcmc.ctr    = list(target  = NULL, kappa = 0.6, jmin = 1e-6, jmax = 3), 
                       data        = NULL){
  stopifnot(is.list(network))
  if(missing(formula.u) & missing(formula.v) & missing(formula.w)){
    stop("formula.u, formula.v, and formula.w cannot be missing")
  }
  
  if(missing(formula.u)){
    if(!missing(f_to_var.u)){stop("f_to_var.u is defined without formula.u")}
  } else{
    if(missing(f_to_var.u)){stop("formula.u is defined without f_to_var.u")}
    stopifnot(tolower(unlist(f_to_var.u)) %in% c("i", "j", "sum", "prod", "same", "adiff", "lower", "greater"))
    f_to_var.u  <- ff_to_var(f_to_var.u)
  }
  if(missing(formula.v)){
    if(!missing(f_to_var.v)){stop("f_to_var.v is defined without formula.u")}
  } else{
    if(missing(f_to_var.v)){stop("formula.v is defined without f_to_var.v")}
    stopifnot(tolower(unlist(f_to_var.v)) %in% c("sum", "prod", "same", "adiff"))
    f_to_var.v  <- ff_to_var(f_to_var.v)
  }
  if(missing(formula.w)){
    if(!missing(f_to_var.w)){stop("f_to_var.w is defined without formula.w")}
  } else{
    if(missing(f_to_var.w)){stop("formula.w is defined without f_to_var.w")}
    stopifnot(tolower(unlist(f_to_var.w)) %in% c("sum", "prod", "same", "adiff"))
    f_to_var.w  <- ff_to_var(f_to_var.w)
  }
  
  M           <- length(network)
  nvec        <- sapply(network, nrow)
  nveccum     <- c(0, cumsum(nvec))
  
  # data
  Xu          <- NULL
  Xv          <- NULL
  Xw          <- NULL
  inter.u     <- FALSE
  inter.v     <- FALSE
  inter.w     <- FALSE
  name.u      <- NULL
  name.v      <- NULL
  name.w      <- NULL
  vname.u     <- NULL
  vname.v     <- NULL
  vname.w     <- NULL
  Ku          <- NULL
  Kv          <- NULL
  Kw          <- NULL
  nvaru       <- 0
  nvarv       <- 0
  nvarw       <- 0
  if(!missing(formula.u)){
    Xu        <- formula.to.data(formula.u, data)
    inter.u   <- Xu$intercept
    name.u    <- Xu$names
    Xu        <- Xu$X
    Ku        <- ncol(Xu)
    if(length(f_to_var.u) != Ku) stop("length(f_to_var.u) does not match with formula.u")
    nvaru     <- length(unlist(f_to_var.u))
    Xu        <- lapply(1:M, function(m) Xu[(nveccum[m] + 1):nveccum[m + 1],, drop = FALSE])
    vname.u   <- fvarnames(f_to_var.u, name.u, Ku, inter.u, "D")
  }
  if(!missing(formula.v)){
    Xv        <- formula.to.data(formula.v, data)
    inter.v   <- Xv$intercept
    name.v    <- Xv$names
    Xv        <- Xv$X
    Kv        <- ncol(Xv)
    if(length(f_to_var.v) != Kv) stop("length(f_to_var.v) does not match with formula.v")
    nvarv     <- length(unlist(f_to_var.v))
    Xv        <- lapply(1:M, function(m) Xv[(nveccum[m] + 1):nveccum[m + 1],, drop = FALSE])
    vname.v   <- fvarnames(f_to_var.v, name.v, Kv, inter.v, "M")
  }
  if(!missing(formula.w)){
    Xw        <- formula.to.data(formula.w, data)
    inter.w   <- Xw$intercept
    name.w    <- Xw$names
    Xw        <- Xw$X
    Kw        <- ncol(Xw)
    if(length(f_to_var.w) != Kw) stop("length(f_to_var.w) does not match with formula.w")
    nvarw     <- length(unlist(f_to_var.w))
    Xw        <- lapply(1:M, function(m) Xw[(nveccum[m] + 1):nveccum[m + 1],, drop = FALSE])
    vname.w   <- fvarnames(f_to_var.w, name.w, Kw, inter.w, "I")
  }
  
  # Construct cluster
  cl          <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  if(!is.null(Xu)){
    Xu        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xu[[m]], f_to_var.u, nvaru, Ku)}
  }
  if(!is.null(Xv)){
    Xv        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xv[[m]], f_to_var.v, nvarv, Kv)}
  }
  if(!is.null(Xw)){
    Xw        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xw[[m]], f_to_var.w, nvarw, Kw)}
  }
  
  # number of parameters
  nparu       <- length(vname.u)
  nparv       <- length(vname.v)
  nparw       <- length(vname.w)
  npar        <- nparu + nparv + nparw
  
  # starting
  if(is.null(starting)){
    starting  <- rep(0, npar)
  } else{
    if(length(starting) != npar){
      stop("length(starting) not equal to the number of parameters to be estimated")
    }
  }
  theta       <- starting
  
  # Mcmc control
  target      <- mcmc.ctr$target
  kappa       <- mcmc.ctr$kappa
  jmin        <- mcmc.ctr$jmin
  jmax        <- mcmc.ctr$jmax
  
  if(is.null(target)){
    target    <- ifelse(npar == 1, 0.44, 
                        ifelse(npar == 2, 0.35, 
                               ifelse(npar == 3, 0.3125, 
                                      ifelse(npar == 4, 0.275, 
                                             ifelse(npar == 5, 0.25, 0.234)))))
  }
  if(is.null(kappa)){
    kappa     <- 0.6
  }
  if(is.null(jmin)){
    jmin      <- 1e-6
  }
  if(is.null(jmax)){
    jmax      <- 3
  }
  
  # prior
  # expectation
  if(is.null(Etheta)){
    Etheta    <- rep(0, npar)
  } else{
    if(length(Etheta) != npar){
      stop("length(Etheta) not equal to the number of parameters to be estimated")
    }
  }
  # inverse of the variance
  if(is.null(invVtheta)){
    invVtheta <- matrix(1/100, npar, npar)
  } else{
    if(length(invVtheta) == 1){invVtheta <- diag(npar)*c(invVtheta)}
    if((nrow(invVtheta) != npar) | (ncol(invVtheta) != npar)){
      stop("dim(invVtheta) not match to the number of parameters to be estimated")
    }
  }
  
  #Sigma
  if(is.null(Sigma)){
    Sigma     <- diag(npar)
  } else{
    if(length(Sigma) == 1){Sigma <- diag(npar)*c(Sigma)}
    if((nrow(Sigma) != npar) | (ncol(Sigma) != npar)){
      stop("dim(Sigma) not match to the number of parameters to be estimated")
    }
  }
  
  # Utility
  uu          <- rep(list(matrix(0, 1, 1)), times = M)
  uv          <- rep(list(matrix(0, 1, 1)), times = M)
  uw          <- rep(list(matrix(0, 1, 1)), times = M)
  uu          <- futil(M, Xu, uu, theta[1:nparu], nparu, nvec, inter.u)
  uv          <- futil(M, Xv, uv, theta[(nparu + 1):(nparu + nparv)], nparv, nvec, inter.v)
  uw          <- futil(M, Xw, uw, theta[(nparu + nparv + 1):npar], nparw, nvec, inter.w)
  
  # The potential function value at the starting point
  pot_ath     <- fpotendir(M, network, uu, uv, uw, nparu, nparv, nparw, nvec)
  
  # MCMC
  jstheta     <- 1  #jumping scale for the adaptive MCMC
  atheta      <- 0  #number of times the draw is accepted
  uup         <- uu
  uvp         <- uv
  uwp         <- uw
  combr       <- ffindcom(nblock) 
  ncombr      <- 2^nblock
  ztncombr    <- 0:(ncombr - 1)
  idrows      <- fIDdir(M, nvec)
  idcols      <- idrows$idcols
  ident       <- idrows$ident
  idrows      <- idrows$idrows
  quiet(gc())
  
  posterior   <- as.data.frame(matrix(0, simutheta, npar)); colnames(posterior) <- c(vname.u, vname.v, vname.w)
  spot        <- c()
  for (s in 1:simutheta) {
    # simulate theta prime
    cat("********************* iteration: ", s, "/", simutheta, "\n", sep = "")
    if(nparu > 0){
      cat("direct links\n")
      cat(head(theta, nparu), "\n")
    }
    if(nparv > 0){
      cat("mutual links\n")
      cat(theta[(nparu + 1):(nparu + nparv)], "\n")
    }
    if(nparw > 0){
      cat("indirect links\n")
      cat(tail(theta, nparw), "\n")
    }
    cat("potential\n")
    cat(pot_ath, "\n")
    
    thetap    <- fsimtheta(theta, Sigma, jstheta, npar)  #**
    # cat("thetap", thetap, "\n")
    # cat("sum(uvp)", sum(uvp[[1]]), "\n")
    # utility at thetap
    uup       <- futil(M, Xu, uup, thetap[1:nparu], nparu, nvec, inter.u) #**
    uvp       <- futil(M, Xv, uvp, thetap[(nparu + 1):(nparu + nparv)], nparv, nvec, inter.v) #**
    uwp       <- futil(M, Xw, uwp, thetap[(nparu + nparv + 1):npar], nparw, nvec, inter.w) #**
    
    # potential function at a (ie, network) and thetap
    pot_athp  <- fpotendir(M, network, uup, uvp, uwp, nparu, nparv, nparw, nvec) #**
    # cat("sum(uvp)", sum(uvp[[1]]), "\n")
    # cat("pot_athp ", pot_athp, "\n")
    # Gibbs to simulate ap
    networkp  <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
      fGibbdir(network[[m]], nblock, ncombr, combr, idrows[[m]], idcols[[m]], ident[[m]], ztncombr, 
               uup[[m]], uvp[[m]], uwp[[m]], nparu, nparv, nparw, nvec[[m]], simunet)}
    cat("Gibbs executed\n")
    
    # potential function at ap 
    pot_apth  <- fpotendir(M, networkp, uu, uv, uw, nparu, nparv, nparw, nvec)
    pot_apthp <- fpotendir(M, networkp, uup, uvp, uwp, nparu, nparv, nparw, nvec)
    
    # acceptance rate of theta
    lalpha    <- pot_apth - pot_ath + pot_athp - pot_apthp  + propdnorm(thetap, Etheta, invVtheta) - propdnorm(theta, Etheta, invVtheta) 
    lalpha    <- min(0, lalpha)
    if(runif(1) < exp(lalpha)){
      theta   <- thetap
      uu      <- uup
      uv      <- uvp
      uw      <- uwp
      pot_ath <- pot_athp
      atheta  <- atheta + 1
      cat("proposal accepted -- acceptance rate: ", round(100*atheta/s, 1), "%\n", sep = "")
    } else{
      cat("proposal rejected -- acceptance rate: ", round(100*atheta/s, 1), "%\n", sep = "")
    }
    jstheta       <- fupdate_js(jstheta, atheta, s, target, kappa, jmin, jmax)
    posterior[s,] <- c(theta)
    spot[s]       <- pot_ath
  }
  out  <- list(posterior = posterior, potential = spot); class(out) <- "mcmcDirNet"
  out
}

# function that computes the utilities
futil        <- function(M, X, u, theta, npar, nvec, inter){
  if(!is.null(X)){
    # cat("npar ", npar, "\n")
    # cat("theta ", theta, "\n")
    # cat("** ", sum(X[[1]]), "\n")
    # cat("** ", sum(u[[1]]), "\n")
    u        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
      futility(X[[m]], theta, npar, nvec[m], inter)}}
  # cat("** ", sum(u[[1]]), "\n")
  u
}

# fonction that computes the potential function value
fpotensym     <- function(M, network, uu, uw, nparu, nparw, nvec){
  sum(unlist(foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
    fQrsym(network[[m]], uu[[m]], uw[[m]], nparu, nparw, nvec[m])}))
}

fpotendir     <- function(M, network, uu, uv, uw, nparu, nparv, nparw, nvec){
  sum(unlist(foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
    fQrdir(network[[m]], uu[[m]], uv[[m]], uw[[m]], nparu, nparv, nparw, nvec[m])}))
}

# formula to variable
ff_to_var  <- function(f_to_var){
  lab      <- c("i", "j", "sum", "prod", "same", "adiff", "lower", "greater")
  code     <- 1:length(lab)
  lapply(f_to_var, function(x) sapply(tolower(x), function(y) code[lab == y]))
}

# variable names
fvarnames  <- function(ftovar, Xnames, K, intercept, prefix){
  len      <- sapply(ftovar, length)
  out      <- unlist(lapply(1:K, function(x) rep(Xnames[x], len[x])))
  lab      <- c("i", "j", "sum", "prod", "same", "adiff", "lower", "greater")
  code     <- 1:length(lab)
  out      <- paste0(out, ".", sapply(unlist(ftovar), function(x) lab[code == x]))
  if(intercept){
    out    <- c("(Intercept)", out)
  }
  out      <- paste0(prefix, ".", out)
}

# update jumping sclace
fupdate_js <- function(jscal, accept, iteration, target, kappa, jmin, jmax){
  out      <- jscal + (accept/iteration - target)/(iteration^kappa)
  if(out < jmin) out <- jmin
  if(out > jmax) out <- jmax
  out
}


