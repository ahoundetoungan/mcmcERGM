#' @title Simulate the posterior distribution of the parameters of an exponential random symmetric graph model
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doRNG "%dorng%"
#' @importFrom doParallel registerDoParallel
#' @importFrom ddpcr quiet
#' @importFrom utils head tail
#' @importFrom stats cov
#' @importFrom matrixcalc is.singular.matrix
#' @export
mcmcSymNet <- function(network, 
                       formula.v, 
                       formula.w, 
                       f_to_var.v, 
                       f_to_var.w,
                       heterogeneity,
                       starting    = NULL,
                       nblock      = 2,
                       ncores      = 1, 
                       prior       = list(),
                       Sigmatheta  = NULL,
                       Sigmahete   = NULL,
                       simutheta   = 1e3,
                       simunet     = 1e3,
                       mcmc.ctr    = list(target  = NULL, kappa = 0.6, jmin = 1e-6, jmax = 3), 
                       data        = NULL){
  stopifnot(is.list(network))
  if(missing(formula.v) & missing(formula.w)){
    stop("both formula.v and formula.w cannot be missing")
  }
  if(missing(formula.v)){
    if(!missing(f_to_var.v)){stop("f_to_var.v is defined without formula.v")}
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
  Xv          <- NULL
  Xw          <- NULL
  inter.v     <- FALSE
  inter.w     <- FALSE
  name.v      <- NULL
  name.w      <- NULL
  vname.v     <- NULL
  vname.w     <- NULL
  Kv          <- NULL
  Kw          <- NULL
  nvaru       <- 0
  nvarw       <- 0
  if(!missing(formula.v)){
    Xv        <- formula.to.data(formula.v, data)
    inter.v   <- Xv$intercept
    name.v    <- Xv$names
    Xv        <- Xv$X
    Kv        <- ncol(Xv)
    if(length(f_to_var.v) != Kv) stop("The length of f_to_var.v does not match formula.v")
    nvaru     <- length(unlist(f_to_var.v))
    Xv        <- lapply(1:M, function(m) Xv[(nveccum[m] + 1):nveccum[m + 1],, drop = FALSE])
    vname.v   <- fvarnames(f_to_var.v, name.v, Kv, inter.v, "M")
  }
  if(!missing(formula.w)){
    Xw        <- formula.to.data(formula.w, data)
    inter.w   <- Xw$intercept
    name.w    <- Xw$names
    Xw        <- Xw$X
    Kw        <- ncol(Xw)
    if(length(f_to_var.w) != Kw) stop("The length of f_to_var.w does not match formula.w")
    nvarw     <- length(unlist(f_to_var.w))
    Xw        <- lapply(1:M, function(m) Xw[(nveccum[m] + 1):nveccum[m + 1],, drop = FALSE])
    vname.w   <- fvarnames(f_to_var.w, name.w, Kw, inter.w, "I")
  }
  
  # Construct cluster
  cl          <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  if(!is.null(Xv)){
    Xv        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xv[[m]], f_to_var.v, nvaru, Kv)}
  }
  if(!is.null(Xw)){
    Xw        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xw[[m]], f_to_var.w, nvarw, Kw)}
  }
  
  # number of parameters
  nparv       <- length(vname.v)
  nparw       <- length(vname.w)
  npar        <- nparv + nparw
  
  # mu and nu
  het         <- FALSE
  mu          <- rep(0, M)
  nu          <- rep(0, M)
  in.mu       <- FALSE #Tells if mu should be inferred
  in.nu       <- FALSE
  if (!missing(heterogeneity)){
    stopifnot(inherits(heterogeneity, "list"))
    if (!all(tolower(names(heterogeneity)) %in% c("mu", "nu"))) stop("'heterogeneity' must be a list containing only 'mu' and 'nu'")
    names(heterogeneity) <- tolower(names(heterogeneity))
    in.mu <- !is.null(heterogeneity$mu)
    in.nu <- !is.null(heterogeneity$nu)
    if (in.mu) {
      mu  <- heterogeneity$mu
    }
    if (in.nu) {
      nu  <- heterogeneity$nu
    }
    het   <- in.mu | in.nu
  }
  if (het) {
    if (length(mu) == 1) {
      mu  <- rep(mu, M)
    } else if (length(mu) != M) {
      stop("'mu' must be either a scalar or a vector of length M")
    }
    if (length(nu) == 1) {
      nu  <- rep(nu, M)
    } else if (length(nu) != M) {
      stop("'nu' must be either a scalar or a vector of length M")
    }
  } 
  khete       <- in.mu + in.nu
  hetval      <- t(cbind(mu, nu))
  if (in.mu && nparv == 0) {
    stop("'mu' heterogeneity cannot be estimated if 'formula.v' is not defined")
  }
  if (in.nu && nparw == 0) {
    stop("'nu' heterogeneity cannot be estimated if 'formula.w' is not defined")
  }
  
  # starting
  if(is.null(starting)){
    starting  <- rep(0, npar)
  } else{
    if(length(starting) != npar){
      stop("The length of 'starting' must match the number of parameters to be estimated")
    }
  }
  theta       <- starting
  
  # Mcmc control
  target      <- mcmc.ctr$target
  kappa       <- mcmc.ctr$kappa
  jmin        <- mcmc.ctr$jmin
  jmax        <- mcmc.ctr$jmax
  
  if(is.null(target)){
    target    <- ifelse((npar + khete*M) == 1, 0.44, 
                        ifelse((npar + khete*M) == 2, 0.35, 
                               ifelse((npar + khete*M) == 3, 0.3125, 
                                      ifelse((npar + khete*M) == 4, 0.275, 
                                             ifelse((npar + khete*M) == 5, 0.25, 0.234)))))
  }
  if(is.null(kappa)){
    kappa     <- 0.6
  }
  if(is.null(jmin)){
    jmin      <- rep(1e-6, 1 + het)
  }
  if(length(jmin) == 1){
    jmin      <- rep(jmin, 1 + het)
  }
  if(is.null(jmax)){
    jmax      <- rep(3, 1 + het)
  }
  if(length(jmax) == 1){
    jmax      <- rep(jmax, 1 + het)
  }
  if (het) {
    if(length(jmin) != 2 | length(jmax) != 2) {
      stop("The minimal and the maximal jumping scales are either scalars (if there is no heterogeneity) or a 2 dimentional vector (if there is heterogeneity)")
    }
  } else {
    if(length(jmin) != 1 | length(jmax) != 1) {
      stop("The minimal and the maximal jumping scales are either scalars (if there is no heterogeneity) or a 2 dimentional vector (if there is heterogeneity)")
    }
  }

  # theta
  # expectation
  Etheta      <- prior$Etheta 
  if(is.null(Etheta)){
    Etheta    <- rep(0, npar)
  } else{
    if(length(Etheta) != npar){
      stop("The length of 'Etheta' must match the number of parameters to be estimated")
    }
  }
  
  # inverse of the variance
  invVtheta   <- prior$invVtheta
  if(is.null(invVtheta)){
    invVtheta <- matrix(0, npar, npar)
  } else{
    if(length(invVtheta) == 1){invVtheta <- diag(npar)*c(invVtheta)}
    if((nrow(invVtheta) != npar) | (ncol(invVtheta) != npar)){
      stop("The dimensions of 'invVtheta' must match the number of parameters to be estimated")
    }
  }
  
  # Variance of proposal
  if(is.null(Sigmatheta)){
    Sigmatheta  <- diag(npar)
  } else{
    if(length(Sigmatheta) == 1){Sigmatheta <- diag(npar)*c(Sigmatheta)}
    if((nrow(Sigmatheta) != npar) | (ncol(Sigmatheta) != npar)){
      stop("The dimensions of 'Sigmatheta' must match the number of parameters to be estimated")
    }
  }
  
  # heterogeneity
  Omega       <- NULL
  scOmega     <- prior$scOmega
  dfOmega     <- prior$dfOmega
  in.het      <- c(in.mu, in.nu)
  int.upd     <- NULL
  if (het) {
    int.upd   <- which(c(vname.v, vname.w) %in% c("M.(Intercept)", "I.(Intercept)")[in.het]) - 1
    if (length(int.upd) > 0){
      hetval[in.het,]  <- frecentering(theta, hetval[in.het, , drop = FALSE], int.upd)
    }
    # variance
    Omega     <- cov(t(hetval[in.het, , drop = FALSE])) 
    if (is.singular.matrix(Omega)) {
      Omega   <- diag(khete)
    }
    
    # Scale parameter Omega
    if(is.null(scOmega)){
      scOmega <- matrix(0, khete, khete)
    } else{
      if(length(scOmega) == 1){scOmega <- diag(khete)*c(scOmega)}
      if((nrow(scOmega) != khete) | (ncol(scOmega) != khete)){
        stop("The dimensions of 'scOmega' must match the number of heterogeneity parameters")
      }
    }
    
    # df of Omega
    if(is.null(dfOmega)){
      dfOmega <- khete
    } else{
      if(dfOmega <= 0){
        stop("Negative degree of freedom")
      }
    }
    
    #  Variance of proposal
    if(is.null(Sigmahete)){
      Sigmahete   <- diag(khete)
    } else{
      if(length(Sigmahete) == 1){Sigmahete <- diag(khete)*c(Sigmahete)}
      if((nrow(Sigmahete) != khete) | (ncol(Sigmahete) != khete)){
        stop("The dimensions of 'Sigmahete' dmust match the number of heterogeneity parameters")
      }
    }
  }
  
  # Utility
  uv          <- rep(list(matrix(0, 1, 1)), times = M)
  uw          <- rep(list(matrix(0, 1, 1)), times = M)
  uv          <- futil(M, Xv, uv, theta[1:nparv], hetval[1,], nparv, nvec, inter.v)
  uw          <- futil(M, Xw, uw, theta[(nparv + 1):npar], hetval[2,], nparw, nvec, inter.w)
  
  # The potential function value at the starting point
  pot_ath     <- fpotensym_eachm(M, network, uv, uw, nparv, nparw, nvec)
  
  # MCMC
  jstheta     <- 1  #jumping scale for the adaptive MCMC
  jshete      <- rep(1, M)
  atheta      <- 0  #number of times the draw is accepted
  ahete       <- rep(0, M)
  uvp         <- uv
  uwp         <- uw
  combr       <- ffindcom(nblock) 
  ncombr      <- 2^nblock
  ztncombr    <- 0:(ncombr - 1)
  idrows      <- fIDsym(M, nvec)
  idcols      <- idrows$idcols
  ident       <- idrows$ident
  idrows      <- idrows$idrows
  quiet(gc())
  
  posterior   <- as.data.frame(matrix(0, simutheta, npar)); colnames(posterior) <- c(vname.v, vname.w)
  simhete     <- NULL
  simOmega    <- NULL
  spot        <- c()
  sJS         <- as.data.frame(matrix(0, simutheta, 1 + het*M)) 
  if (het) {
    simhete   <- as.data.frame(matrix(0, simutheta, M*khete)) 
    simOmega  <- as.data.frame(matrix(0, simutheta, khete^2)) 
    nacol     <- c("mu", "nu")[in.het] 
    colnames(simhete)  <- paste0(rep(nacol, each = M), rep(1:M, khete))
    colnames(simOmega) <- paste0("s_", sapply(nacol, function(s1) sapply(nacol, function(s2) paste0(s1, s2))))
    colnames(sJS)      <- c("theta", paste0("heter.", 1:M))
  } 
  hetvalp     <- hetval
  
  for (s in 1:simutheta) {
    # simulate theta prime
    cat("********************* Iteration: ", s, "/", simutheta, "\n", sep = "")
    cat("Updating theta\n")
    cat("Jumping Scale: ", jstheta, "\n", sep = "")
    if(nparv > 0){
      cat("mutual links\n")
      cat(head(theta, nparv), "\n")
    }
    if(nparw > 0){
      cat("indirect links\n")
      cat(tail(theta, nparw), "\n")
    }
    
    thetap    <- fsimtheta(theta, Sigmatheta, jstheta, npar)  #**
    
    # utility at thetap
    uvp       <- futil(M, Xv, uvp, thetap[1:nparv], hetval[1,], nparv, nvec, inter.v) #**
    uwp       <- futil(M, Xw, uwp, thetap[(nparv + 1):npar], hetval[2,], nparw, nvec, inter.w) #**
    
    # potential function at a (ie, network) and thetap
    pot_athp  <- fpotensym_eachm(M, network, uvp, uwp, nparv, nparw, nvec) #**
    
    # Gibbs to simulate ap
    networkp  <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
      fGibbsym(network[[m]], nblock, ncombr, combr, idrows[[m]], idcols[[m]], ident[[m]], ztncombr, 
               uvp[[m]], uwp[[m]], nparv, nparw, nvec[[m]], simunet)}
    cat("Gibbs executed\n")
    
    # potential function at ap 
    pot_apth  <- fpotensym_eachm(M, networkp, uv, uw, nparv, nparw, nvec)
    pot_apthp <- fpotensym_eachm(M, networkp, uvp, uwp, nparv, nparw, nvec)
    
    # acceptance rate of theta
    lalpha      <- sum(pot_apth - pot_ath + pot_athp - pot_apthp)  + propdnorm(thetap, Etheta, invVtheta) - propdnorm(theta, Etheta, invVtheta) 
    if(runif(1) < exp(lalpha)){
      theta   <- thetap
      uv      <- uvp
      uw      <- uwp
      pot_ath <- pot_athp
      atheta  <- atheta + 1
      cat("Simulated theta accepted -- acceptance rate: ", round(100*atheta/s, 1), "%\n", sep = "")
    } else{
      cat("Simulated theta rejected -- acceptance rate: ", round(100*atheta/s, 1), "%\n", sep = "")
    }
    jstheta       <- fupdate_jstheta(jstheta, atheta, s, target, kappa, jmin[1], jmax[1])
    
    # Update heterogeneity (if any)
    if (het) {
      cat("\nUpdating heterogeneity\n")
      cat("Jumping Scale\n ", "min: ", min(jshete), " -- median: ", median(jshete), " -- max: ",  
          max(jshete), "\n", sep = "")
      hetvalp[in.het,] <- fsimhete(hetval[in.het, , drop = FALSE], Sigmahete, jshete, khete, M)
      
      # utility at hetep
      if (in.mu) {
        uvp     <- futil(M, Xv, uvp, theta[1:nparv], hetvalp[1,], nparv, nvec, inter.v) #**
      }
      if (in.nu) {
        uwp     <- futil(M, Xw, uwp, theta[(nparv + 1):npar], hetvalp[2,], nparw, nvec, inter.w) #**
      }
      
      # potential function at a (ie, network) and hetep
      pot_athp  <- fpotensym_eachm(M, network, uvp, uwp, nparv, nparw, nvec) #**
      
      # Gibbs to simulate ap
      networkp  <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
        fGibbsym(network[[m]], nblock, ncombr, combr, idrows[[m]], idcols[[m]], ident[[m]], ztncombr, 
                 uvp[[m]], uwp[[m]], nparv, nparw, nvec[[m]], simunet)}
      cat("Gibbs executed\n")
      
      # potential function at ap 
      pot_apth  <- fpotensym_eachm(M, networkp, uv, uw, nparv, nparw, nvec)
      pot_apthp <- fpotensym_eachm(M, networkp, uvp, uwp, nparv, nparw, nvec)
      
      # acceptance rate of theta
      lalpha      <- pot_apth - pot_ath + pot_athp - pot_apthp + propdnorm_eachm(hetvalp[in.het, , drop = FALSE], Omega) - propdnorm_eachm(hetval[in.het, , drop = FALSE], Omega)
      tp          <- runif(M) < exp(lalpha)
      hetval[,tp] <- hetvalp[,tp]
      if (length(int.upd) > 0){
        hetval[in.het,] <- frecentering(theta, hetval[in.het, , drop = FALSE], int.upd)
      }
      if (in.mu) {
        uv[tp]    <- uvp[tp]
      }
      if (in.nu) {
        uw[tp]    <- uwp[tp]
      }
      pot_ath[tp] <- pot_athp[tp]
      ahete[tp]   <- ahete[tp] + 1
      cat("Share of accepted heterogeneity: ", round(mean(tp)*100, 1), "%\n", sep = "")
      tp          <- 100*ahete/s
      cat("Acceptance rate\n min: ", round(min(tp), 1), "% -- median: ", round(median(tp), 1), "% -- max: ",  round(max(tp), 1), 
          "%\n", sep = "")
      jshete      <- fupdate_jshete(jshete, ahete, s, target, kappa, jmin[2], jmax[2])
      simhete[s,] <- c(t(hetval[in.het, , drop = FALSE]))
      sJS[s, 2:(M + 1)] <- jshete
      
      # Update Omega
      Omega        <- updateOmega(dfOmega, scOmega, M, hetval[in.het, , drop = FALSE])
      simOmega[s,] <- c(Omega)
    }
    posterior[s,] <- c(theta)
    sJS[s, 1]     <- jstheta
    spot[s]       <- sum(pot_ath)
    cat("potential: ", spot[s], "\n\n", sep = "")
  }
  ARate   <- c(atheta, ahete)/simutheta; names(ARate) <- c("theta", paste0("heter", 1:M))
  if (!het) {
    sJS   <- unlist(sJS); names(sJS) <- NULL
    ARate <- ARate[1]
  }
  out     <- list(posterior = list(theta         = posterior, 
                                   heterogeneity = simhete, 
                                   Omega         = simOmega),
                  potential  = spot,
                  jumpscale  = sJS,
                  acceptrate = ARate,
                  simutheta  = simutheta,
                  simunet    = simunet)
  class(out) <- "mcmcSymNet"
  out
}



#' @title Simulate the posterior distribution of the parameters of an exponential random symmetric graph model
#' @export
mcmcDirNet <- function(network, 
                       formula.u, 
                       formula.v,
                       formula.w, 
                       f_to_var.u, 
                       f_to_var.v,
                       f_to_var.w,
                       heterogeneity,
                       starting    = NULL,
                       nblock      = 2,
                       ncores      = 1, 
                       prior       = list(),
                       Sigmatheta  = NULL,
                       Sigmahete   = NULL,
                       simutheta   = 1e3,
                       simunet     = 1e3,
                       mcmc.ctr    = list(target  = NULL, kappa = 0.6, jmin = 1e-6, jmax = 3), 
                       data        = NULL){
  stopifnot(is.list(network))
  
  if (missing(formula.u) && missing(formula.v) && missing(formula.w)) {
    stop("formula.u, formula.v, and formula.w cannot all be missing")
  }
  
  if (missing(formula.u)) {
    if (!missing(f_to_var.u)) {
      stop("f_to_var.u is defined without formula.u")
    }
  } else {
    if (missing(f_to_var.u)) {
      stop("formula.u is defined without f_to_var.u")
    }
    stopifnot(tolower(unlist(f_to_var.u)) %in% c("i", "j", "sum", "prod", "same", "adiff", "lower", "greater"))
    f_to_var.u <- ff_to_var(f_to_var.u)
  }
  
  if (missing(formula.v)) {
    if (!missing(f_to_var.v)) {
      stop("f_to_var.v is defined without formula.v")
    }
  } else {
    if (missing(f_to_var.v)) {
      stop("formula.v is defined without f_to_var.v")
    }
    stopifnot(tolower(unlist(f_to_var.v)) %in% c("sum", "prod", "same", "adiff"))
    f_to_var.v <- ff_to_var(f_to_var.v)
  }
  
  if (missing(formula.w)) {
    if (!missing(f_to_var.w)) {
      stop("f_to_var.w is defined without formula.w")
    }
  } else {
    if (missing(f_to_var.w)) {
      stop("formula.w is defined without f_to_var.w")
    }
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
    if(length(f_to_var.u) != Ku) stop("The length of f_to_var.u does not match formula.u")
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
    if(length(f_to_var.v) != Kv) stop("The length of f_to_var.v does not match formula.v")
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
    if(length(f_to_var.w) != Kw) stop("The length of f_to_var.w does not match formula.w")
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
  
  # mu, nu, zeta
  het         <- FALSE
  ze          <- rep(0, M)
  mu          <- rep(0, M)
  nu          <- rep(0, M)
  in.ze       <- FALSE #Tells if ze should be inferred
  in.mu       <- FALSE 
  in.nu       <- FALSE
  if (!missing(heterogeneity)){
    stopifnot(inherits(heterogeneity, "list"))
    if (!all(tolower(names(heterogeneity)) %in% c("mu", "nu", "zeta"))) stop("'heterogeneity' must be a list containing only 'zeta', 'mu', and 'nu'")
    names(heterogeneity) <- tolower(names(heterogeneity))
    in.ze <- !is.null(heterogeneity$zeta)
    in.mu <- !is.null(heterogeneity$mu)
    in.nu <- !is.null(heterogeneity$nu)
    if (in.ze) {
      ze  <- heterogeneity$zeta
    }
    if (in.mu) {
      mu  <- heterogeneity$mu
    }
    if (in.nu) {
      nu  <- heterogeneity$nu
    }
    het   <- in.mu | in.nu | in.ze
  }
  if (het) {
    if (length(ze) == 1) {
      ze  <- rep(ze, M)
    } else if (length(ze) != M) {
      stop("'zeta' must be either a scalar or a vector of length M")
    }
    if (length(mu) == 1) {
      mu  <- rep(mu, M)
    } else if (length(mu) != M) {
      stop("'mu' must be either a scalar or a vector of length M")
    }
    if (length(nu) == 1) {
      nu  <- rep(nu, M)
    } else if (length(nu) != M) {
      stop("'nu' must be either a scalar or a vector of length M")
    }
  } 
  khete       <- in.ze + in.mu + in.nu
  hetval      <- t(cbind(ze, mu, nu))
  if (in.ze && nparu == 0) {
    stop("'zeta' heterogeneity cannot be estimated if 'formula.u' is not defined")
  }
  if (in.mu && nparv == 0) {
    stop("'mu' heterogeneity cannot be estimated if 'formula.v' is not defined")
  }
  if (in.nu && nparw == 0) {
    stop("'nu' heterogeneity cannot be estimated if 'formula.w' is not defined")
  }
  
  # starting
  if(is.null(starting)){
    starting  <- rep(0, npar)
  } else{
    if(length(starting) != npar){
      stop("The length of 'starting' must match the number of parameters to be estimated")
    }
  }
  theta       <- starting
  
  # Mcmc control
  target      <- mcmc.ctr$target
  kappa       <- mcmc.ctr$kappa
  jmin        <- mcmc.ctr$jmin
  jmax        <- mcmc.ctr$jmax
  
  if(is.null(target)){
    target    <- ifelse((npar + khete*M) == 1, 0.44, 
                        ifelse((npar + khete*M) == 2, 0.35, 
                               ifelse((npar + khete*M) == 3, 0.3125, 
                                      ifelse((npar + khete*M) == 4, 0.275, 
                                             ifelse((npar + khete*M) == 5, 0.25, 0.234)))))
  }
  if(is.null(kappa)){
    kappa     <- 0.6
  }
  if(is.null(jmin)){
    jmin      <- rep(1e-6, 1 + het)
  }
  if(length(jmin) == 1){
    jmin      <- rep(jmin, 1 + het)
  }
  if(is.null(jmax)){
    jmax      <- rep(3, 1 + het)
  }
  if(length(jmax) == 1){
    jmax      <- rep(jmax, 1 + het)
  }
  if (het) {
    if(length(jmin) != 2 | length(jmax) != 2) {
      stop("The minimal and the maximal jumping scales are either scalars (if there is no heterogeneity) or a 2 dimentional vector (if there is heterogeneity)")
    }
  } else {
    if(length(jmin) != 1 | length(jmax) != 1) {
      stop("The minimal and the maximal jumping scales are either scalars (if there is no heterogeneity) or a 2 dimentional vector (if there is heterogeneity)")
    }
  }
  
  # theta
  # expectation
  Etheta      <- prior$Etheta 
  if(is.null(Etheta)){
    Etheta    <- rep(0, npar)
  } else{
    if(length(Etheta) != npar){
      stop("The length of 'Etheta' must match the number of parameters to be estimated")
    }
  }
  
  # inverse of the variance
  invVtheta   <- prior$invVtheta
  if(is.null(invVtheta)){
    invVtheta <- matrix(0, npar, npar)
  } else{
    if(length(invVtheta) == 1){invVtheta <- diag(npar)*c(invVtheta)}
    if((nrow(invVtheta) != npar) | (ncol(invVtheta) != npar)){
      stop("The dimensions of 'invVtheta' must match the number of parameters to be estimated")
    }
  }
  
  # Variance of proposal
  if(is.null(Sigmatheta)){
    Sigmatheta  <- diag(npar)
  } else{
    if(length(Sigmatheta) == 1){Sigmatheta <- diag(npar)*c(Sigmatheta)}
    if((nrow(Sigmatheta) != npar) | (ncol(Sigmatheta) != npar)){
      stop("The dimensions of 'Sigmatheta' must match the number of parameters to be estimated")
    }
  }
  
  # heterogeneity
  Omega       <- NULL
  scOmega     <- prior$scOmega
  dfOmega     <- prior$dfOmega
  in.het      <- c(in.ze, in.mu, in.nu)
  int.upd     <- NULL
  if (het) {
    int.upd   <- which(c(vname.u, vname.v, vname.w) %in% c("D.(Intercept)", "M.(Intercept)", "I.(Intercept)")[in.het]) - 1
    if (length(int.upd) > 0){
      hetval[in.het,]  <- frecentering(theta, hetval[in.het, , drop = FALSE], int.upd)
    }
    # variance
    Omega     <- cov(t(hetval[in.het, , drop = FALSE])) 
    if (is.singular.matrix(Omega)) {
      Omega   <- diag(khete)
    }
    
    # Scale parameter Omega
    if (is.null(scOmega)) {
      scOmega   <- matrix(0, khete, khete)
    } else{
      if(length(scOmega) == 1){scOmega <- diag(khete)*c(scOmega)}
      if((nrow(scOmega) != khete) | (ncol(scOmega) != khete)){
        stop("The dimensions of 'scOmega' must match the number of heterogeneity parameters")
      }
    }
    
    # df of Omega
    if(is.null(dfOmega)){
      dfOmega   <- khete
    } else{
      if(dfOmega <= 0){
        stop("Negative degree of freedom")
      }
    }
    
    #  Variance of proposal
    if(is.null(Sigmahete)){
      Sigmahete   <- diag(khete)
    } else{
      if(length(Sigmahete) == 1){Sigmahete <- diag(khete)*c(Sigmahete)}
      if((nrow(Sigmahete) != khete) | (ncol(Sigmahete) != khete)){
        stop("The dimensions of 'Sigmahete' dmust match the number of heterogeneity parameters")
      }
    }
  }
  
  # Utility
  uu          <- rep(list(matrix(0, 1, 1)), times = M)
  uv          <- rep(list(matrix(0, 1, 1)), times = M)
  uw          <- rep(list(matrix(0, 1, 1)), times = M)
  uu          <- futil(M, Xu, uu, head(theta, nparu), hetval[1,], nparu, nvec, inter.u)
  uv          <- futil(M, Xv, uv, theta[(nparu + 1):(nparu + nparv)], hetval[2,], nparv, nvec, inter.v)
  uw          <- futil(M, Xw, uw, tail(theta, nparw), hetval[3,], nparw, nvec, inter.w)
  
  # The potential function value at the starting point
  pot_ath     <- fpotendir_eachm(M, network, uu, uv, uw, nparu, nparv, nparw, nvec)
  
  # MCMC
  jstheta     <- 1  #jumping scale for the adaptive MCMC
  jshete      <- rep(1, M)
  atheta      <- 0  #number of times the draw is accepted
  ahete       <- rep(0, M)
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
  simhete     <- NULL
  simOmega    <- NULL
  spot        <- c()
  sJS         <- as.data.frame(matrix(0, simutheta, 1 + het*M)) 
  if (het) {
    simhete   <- as.data.frame(matrix(0, simutheta, M*khete)) 
    simOmega  <- as.data.frame(matrix(0, simutheta, khete^2)) 
    nacol     <- c("zeta", "mu", "nu")[in.het] 
    colnames(simhete)  <- paste0(rep(nacol, each = M), rep(1:M, khete))
    colnames(simOmega) <- paste0("s_", sapply(nacol, function(s1) sapply(nacol, function(s2) paste0(s1, s2))))
    colnames(sJS)      <- c("theta", paste0("heter.", 1:M))
  } 
  hetvalp     <- hetval
  
  for (s in 1:simutheta) {
    # simulate theta prime
    cat("********************* Iteration: ", s, "/", simutheta, "\n", sep = "")
    cat("Updating theta\n")
    cat("Jumping Scale: ", jstheta, "\n", sep = "")
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
    
    thetap    <- fsimtheta(theta, Sigmatheta, jstheta, npar)  #**

    # utility at thetap
    uup       <- futil(M, Xu, uup, head(thetap, nparu), hetval[1,], nparu, nvec, inter.u) #**
    uvp       <- futil(M, Xv, uvp, thetap[(nparu + 1):(nparu + nparv)], hetval[2,], nparv, nvec, inter.v) #**
    uwp       <- futil(M, Xw, uwp, tail(thetap, nparw), hetval[3,], nparw, nvec, inter.w) #**
    
    # potential function at a (ie, network) and thetap
    pot_athp  <- fpotendir_eachm(M, network, uup, uvp, uwp, nparu, nparv, nparw, nvec) #**

    # Gibbs to simulate ap
    networkp  <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
      fGibbdir(network[[m]], nblock, ncombr, combr, idrows[[m]], idcols[[m]], ident[[m]], ztncombr, 
               uup[[m]], uvp[[m]], uwp[[m]], nparu, nparv, nparw, nvec[[m]], simunet)}
    cat("Gibbs executed\n")
    
    # potential function at ap 
    pot_apth  <- fpotendir_eachm(M, networkp, uu, uv, uw, nparu, nparv, nparw, nvec)
    pot_apthp <- fpotendir_eachm(M, networkp, uup, uvp, uwp, nparu, nparv, nparw, nvec)
    
    # acceptance rate of theta
    lalpha      <- sum(pot_apth - pot_ath + pot_athp - pot_apthp) + propdnorm(thetap, Etheta, invVtheta) - propdnorm(theta, Etheta, invVtheta) 
    if(runif(1) < exp(lalpha)){
      theta   <- thetap
      uu      <- uup
      uv      <- uvp
      uw      <- uwp
      pot_ath <- pot_athp
      atheta  <- atheta + 1
      cat("Simulated theta accepted -- acceptance rate: ", round(100*atheta/s, 1), "%\n", sep = "")
    } else{
      cat("Simulated theta rejected -- acceptance rate: ", round(100*atheta/s, 1), "%\n", sep = "")
    }
    jstheta       <- fupdate_jstheta(jstheta, atheta, s, target, kappa, jmin[1], jmax[1])
    
    # Update heterogeneity (if any)
    if (het) {
      cat("\nUpdating heterogeneity\n")
      cat("Jumping Scale\n ", "min: ", min(jshete), " -- median: ", median(jshete), " -- max: ",  
          max(jshete), "\n", sep = "")
      hetvalp[in.het,] <- fsimhete(hetval[in.het, , drop = FALSE], Sigmahete, jshete, khete, M)

      # utility at hetep
      if (in.ze) {
        uup     <- futil(M, Xu, uup, head(theta, nparu), hetvalp[1,], nparu, nvec, inter.u) 
      }
      if (in.mu) {
        uvp     <- futil(M, Xv, uvp, theta[(nparu + 1):(nparu + nparv)], hetvalp[2,], nparv, nvec, inter.v) 
      }
      if (in.nu) {
        uwp     <- futil(M, Xw, uwp, tail(theta, nparw), hetvalp[3,], nparw, nvec, inter.w) 
      }
      
      # potential function at a (ie, network) and hetep
      pot_athp  <- fpotendir_eachm(M, network, uup, uvp, uwp, nparu, nparv, nparw, nvec)
      
      # Gibbs to simulate ap
      networkp  <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
        fGibbdir(network[[m]], nblock, ncombr, combr, idrows[[m]], idcols[[m]], ident[[m]], ztncombr, 
                 uup[[m]], uvp[[m]], uwp[[m]], nparu, nparv, nparw, nvec[[m]], simunet)}
      cat("Gibbs executed\n")
      
      # potential function at ap 
      pot_apth  <- fpotendir_eachm(M, networkp, uu, uv, uw, nparu, nparv, nparw, nvec)
      pot_apthp <- fpotendir_eachm(M, networkp, uup, uvp, uwp, nparu, nparv, nparw, nvec)
      
      # acceptance rate of theta
      lalpha      <- pot_apth - pot_ath + pot_athp - pot_apthp + propdnorm_eachm(hetvalp[in.het, , drop = FALSE], Omega) - propdnorm_eachm(hetval[in.het, , drop = FALSE], Omega)
      tp          <- runif(M) < exp(lalpha)
      hetval[,tp] <- hetvalp[,tp]
      if (length(int.upd) > 0){
        hetval[in.het,] <- frecentering(theta, hetval[in.het, , drop = FALSE], int.upd)
      }
      if (in.ze) {
        uu[tp]    <- uup[tp]
      }
      if (in.mu) {
        uv[tp]    <- uvp[tp]
      }
      if (in.nu) {
        uw[tp]    <- uwp[tp]
      }
      pot_ath[tp] <- pot_athp[tp]
      ahete[tp]   <- ahete[tp] + 1
      cat("Share of accepted heterogeneity: ", round(mean(tp)*100, 1), "%\n", sep = "")
      tp          <- 100*ahete/s
      cat("Acceptance rate\n min: ", round(min(tp), 1), "% -- median: ", round(median(tp), 1), "% -- max: ",  
          round(max(tp), 1), "%\n", sep = "")
      jshete      <- fupdate_jshete(jshete, ahete, s, target, kappa, jmin[2], jmax[2])
      simhete[s,] <- c(t(hetval[in.het, , drop = FALSE]))
      sJS[s, 2:(M + 1)] <- jshete
      
      # Update Omega
      Omega        <- updateOmega(dfOmega, scOmega, M, hetval[in.het, , drop = FALSE])
      simOmega[s,] <- c(Omega)
    }
    posterior[s,] <- c(theta)
    sJS[s, 1]     <- jstheta
    spot[s]       <- sum(pot_ath)
    cat("potential: ", spot[s], "\n\n", sep = "")
  }
  ARate   <- c(atheta, ahete)/simutheta; names(ARate) <- c("theta", paste0("heter", 1:M))
  if (!het) {
    sJS   <- unlist(sJS); names(sJS) <- NULL
    ARate <- ARate[1]
  }
  out     <- list(posterior = list(theta         = posterior, 
                                   heterogeneity = simhete, 
                                   Omega         = simOmega),
                  potential  = spot,
                  jumpscale  = sJS,
                  acceptrate = ARate,
                  simutheta  = simutheta,
                  simunet    = simunet)
  
  
  class(out) <- "mcmcDirNet"
  out
}

# function that computes the utilities
futil        <- function(M, X, u, theta, hetval = rep(0, M), npar, nvec, inter){
  if(!is.null(X)){
    # cat("npar ", npar, "\n")
    # cat("theta ", theta, "\n")
    # cat("** ", sum(X[[1]]), "\n")
    # cat("** ", sum(u[[1]]), "\n")
    u        <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
      futility(X = X[[m]], theta = theta, hetval = hetval[m], npar = npar, n = nvec[m], intercept = inter)}}
  # cat("** ", sum(u[[1]]), "\n")
  u
}

# fonction that computes the potential function value
fpotensym     <- function(M, network, uv, uw, nparu, nparw, nvec){
  sum(foreach(m = 1:M, .packages  = "mcmcERGM", .combine = "c") %dorng% {
    fQrsym(network[[m]], uv[[m]], uw[[m]], nparu, nparw, nvec[m])})
}

fpotensym_eachm <- function(M, network, uv, uw, nparu, nparw, nvec){
  foreach(m = 1:M, .packages  = "mcmcERGM", .combine = "cbind") %dorng% {
    fQrsym(network[[m]], uv[[m]], uw[[m]], nparu, nparw, nvec[m])}
}

fpotendir     <- function(M, network, uu, uv, uw, nparu, nparv, nparw, nvec){
  sum(foreach(m = 1:M, .packages  = "mcmcERGM", .combine = "c") %dorng% {
    fQrdir(network[[m]], uu[[m]], uv[[m]], uw[[m]], nparu, nparv, nparw, nvec[m])})
}

fpotendir_eachm <- function(M, network, uu, uv, uw, nparu, nparv, nparw, nvec){
  foreach(m = 1:M, .packages  = "mcmcERGM", .combine = "cbind") %dorng% {
    fQrdir(network[[m]], uu[[m]], uv[[m]], uw[[m]], nparu, nparv, nparw, nvec[m])}
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