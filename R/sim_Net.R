#' @title Simulate networks using an exponential random symmetric graph model
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doRNG "%dorng%"
#' @importFrom doParallel registerDoParallel
#' @importFrom ddpcr quiet
#' @importFrom utils head tail
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom Matrix bdiag
#' @export
simSymNet <- function(formula.u, 
                      formula.w, 
                      f_to_var.u, 
                      f_to_var.w,
                      theta,
                      heterogeneity,
                      size        = NULL,
                      init.net    = NULL, 
                      nblock      = 2,
                      ncores      = 1, 
                      simunet     = 1e3,
                      fgraph      = NULL,
                      fsubgraph   = NULL,
                      data        = NULL){
  if(inherits(init.net, c("data.frame", "matrix"))){
    init.net    <- list(as.matrix(init.net))
  }
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
  if(!is.null(fgraph)){
    stopifnot(inherits(fgraph, "list"))
  }
  if(!is.null(fsubgraph)){
    stopifnot(inherits(fsubgraph, "list"))
  }
  
  # size
  if(is.null(size)){
    if(is.null(init.net)){
      stop("both size and init.net are missing")
    } else{
      size    <- sapply(init.net, nrow)
    }
  } 
  M           <- length(size)
  nvec        <- size
  nveccum     <- c(0, cumsum(nvec))
  if(is.null(init.net)){
    init.net  <- lapply(nvec, function(x) matrix(0, x, x))
  } else{
    stopifnot(is.list(init.net))
  }
  nfun        <- length(fgraph) + length(fsubgraph)
  namefun     <- c(names(fgraph), names(fsubgraph))
  
  
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
    if(length(f_to_var.u) != Ku) stop("The length of f_to_var.u does not match formula.u")
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
    if(length(f_to_var.w) != Kw) stop("The length of f_to_var.w does not match formula.w")
    nvarw     <- length(unlist(f_to_var.w))
    Xw        <- lapply(1:M, function(m) Xw[(nveccum[m] + 1):nveccum[m + 1],, drop = FALSE])
    vname.w   <- fvarnames(f_to_var.w, name.w, Kw, inter.w, "I")
  }
  
  # Construct cluster
  cl          <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  
  # number of parameters
  nparu       <- length(vname.u)
  nparw       <- length(vname.w)
  npar        <- nparu + nparw
  
  #X
  if(!is.null(Xu)){
    Xu           <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xu[[m]], f_to_var.u, nvaru, Ku)}
    attr(Xu, "doRNG_version") <- NULL
    attr(Xu, "rng")           <- NULL
  }
  if(!is.null(Xw)){
    Xw           <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xw[[m]], f_to_var.w, nvarw, Kw)}
    attr(Xw, "doRNG_version") <- NULL
    attr(Xw, "rng")           <- NULL
  }
  
  # theta
  # print(theta)
  # print(npar)
  if(length(theta) != npar){
    stop("length(theta) not equal to the number of parameters to be estimated")
  }
  
  # ze and nu
  het         <- FALSE
  ze          <- rep(0, M)
  nu          <- rep(0, M)
  in.ze       <- FALSE 
  in.nu       <- FALSE
  if (!missing(heterogeneity)) {
    stopifnot(inherits(heterogeneity, "list"))
    if (!all(tolower(names(heterogeneity)) %in% c("zeta", "nu"))) stop("'heterogeneity' must be a list containing only 'zeta' and 'nu'")
    names(heterogeneity) <- tolower(names(heterogeneity))
    in.ze <- !is.null(heterogeneity$ze)
    in.nu <- !is.null(heterogeneity$nu)
    if (in.ze) {
      ze  <- heterogeneity$ze
    }
    if (in.nu) {
      nu  <- heterogeneity$nu
    }
    het   <- in.ze | in.nu
  }
  if (het) {
    if (length(ze) == 1) {
      ze  <- rep(ze, M)
    } else if (length(ze) != M) {
      stop("'zeta' must be either a scalar or a vector of length M")
    }
    if (length(nu) == 1) {
      nu  <- rep(nu, M)
    } else if (length(nu) != M) {
      stop("'nu' must be either a scalar or a vector of length M")
    }
  } 
  if (in.ze && nparu == 0) {
    stop("'zeta' heterogeneity cannot be estimated if 'formula.u' is not defined")
  }
  if (in.nu && nparw == 0) {
    stop("'nu' heterogeneity cannot be estimated if 'formula.w' is not defined")
  }
  
  # Utility
  uu          <- rep(list(matrix(0, 1, 1)), times = M)
  uw          <- rep(list(matrix(0, 1, 1)), times = M)
  uu          <- futil(M = M, X = Xu, u = uu, theta = theta[1:nparu], hetval = ze, npar = nparu, nvec = nvec, inter = inter.u)
  uw          <- futil(M = M, X = Xw, u = uw, theta = theta[(nparu + 1):npar], hetval = nu, npar = nparw, nvec = nvec, inter = inter.w)
  
  # Simulations
  combr       <- ffindcom(nblock) 
  ncombr      <- 2^nblock
  ztncombr    <- 0:(ncombr - 1)
  idrows      <- fIDsym(M, nvec)
  idcols      <- idrows$idcols
  ident       <- idrows$ident
  idrows      <- idrows$idrows
  quiet(gc())
  out         <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
    fGibbsym2(init.net[[m]], nblock, ncombr, combr, idrows[[m]], idcols[[m]], ident[[m]], 
              ztncombr, uu[[m]], uw[[m]], nparu, nparw, nvec[[m]], simunet)}
  
  network     <- lapply(1:M, function(m) out[[m]]$net)
  degree      <- do.call(rbind, lapply(1:M, function(m) out[[m]]$deg))
  potent      <- Reduce("+", lapply(1:M, function(m) c(out[[m]]$pot)))
  
  # export data
  if(inter.u){
    vname.u   <- vname.u[-1]
  }
  if(inter.w){
    vname.w   <- vname.w[-1]
  }

  if(!is.null(Xu)){
    Xu          <- as.data.frame(lapply(1:nvaru, function(x) mat.to.vec(lapply(1:M, function(m) Xu[[m]][,,x]))))
    colnames(Xu)<- vname.u}
  if(!is.null(Xw)){
    Xw          <- as.data.frame(lapply(1:nvarw, function(x) mat.to.vec(lapply(1:M, function(m) Xw[[m]][,,x]))))
    colnames(Xw)<- vname.w}
  
  # output
  out         <- list(network = network, data = list(Xu = Xu, Xw = Xw), degree = degree, potential = potent)
  class(out)  <- "simSymNet"
  out
}



#' @title Simulate networks using an exponential random directed graph model
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doRNG "%dorng%"
#' @importFrom doParallel registerDoParallel
#' @importFrom ddpcr quiet
#' @importFrom utils head tail
#' @export
simDirNet <- function(formula.u, 
                      formula.v,
                      formula.w, 
                      f_to_var.u, 
                      f_to_var.v,
                      f_to_var.w,
                      theta,
                      heterogeneity,
                      size        = NULL,
                      init.net    = NULL, 
                      nblock      = 2,
                      ncores      = 1, 
                      simunet     = 1e3,
                      fgraph      = NULL,
                      fsubgraph   = NULL,
                      data        = NULL){
  if(inherits(init.net, c("data.frame", "matrix"))){
    init.net    <- list(as.matrix(init.net))
  }
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
    if(!missing(f_to_var.v)){stop("f_to_var.v is defined without formula.v")}
  } else{
    if(missing(f_to_var.v)){stop("formula.v is defined without f_to_var.u")}
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
  if(!is.null(fgraph)){
    stopifnot(inherits(fgraph, "list"))
  }
  if(!is.null(fsubgraph)){
    stopifnot(inherits(fsubgraph, "list"))
  }
  # size
  if(is.null(size)){
    if(is.null(init.net)){
      stop("both size and init.net are missing")
    } else{
      size    <- sapply(init.net, nrow)
    }
  } 
  M           <- length(size)
  nvec        <- size
  nveccum     <- c(0, cumsum(nvec))
  if(is.null(init.net)){
    init.net  <- lapply(nvec, function(x) matrix(0, x, x))
  } else{
    stopifnot(is.list(init.net))
  }
  nfun        <- length(fgraph) + length(fsubgraph)
  namefun     <- c(names(fgraph), names(fsubgraph))
  
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
  
  # number of parameters
  nparu       <- length(vname.u)
  nparv       <- length(vname.v)
  nparw       <- length(vname.w)
  npar        <- nparu + nparv + nparw
  
  #X
  if(!is.null(Xu)){
    Xu           <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xu[[m]], f_to_var.u, nvaru, Ku)}
    attr(Xu, "doRNG_version") <- NULL
    attr(Xu, "rng")           <- NULL
  }
  
  if(!is.null(Xv)){
    Xv           <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xv[[m]], f_to_var.v, nvarv, Kv)}
    attr(Xv, "doRNG_version") <- NULL
    attr(Xv, "rng")           <- NULL
  }
  
  if(!is.null(Xw)){
    Xw           <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {fdatar(Xw[[m]], f_to_var.w, nvarw, Kw)}
    attr(Xu, "doRNG_version") <- NULL
    attr(Xu, "rng")           <- NULL
  }
  
  # theta
  # print(theta)
  # print(npar)
  # print(nparu)
  # print(nparv)
  # print(nparw)
  if(length(theta) != npar){
    stop("length(theta) not equal to the number of parameters to be estimated")
  }
  
  # mu, nu, zeta
  het         <- FALSE
  ze          <- rep(0, M)
  mu          <- rep(0, M)
  nu          <- rep(0, M)
  in.ze       <- FALSE 
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
  if (in.ze && nparu == 0) {
    stop("'zeta' heterogeneity cannot be estimated if 'formula.u' is not defined")
  }
  if (in.mu && nparv == 0) {
    stop("'mu' heterogeneity cannot be estimated if 'formula.v' is not defined")
  }
  if (in.nu && nparw == 0) {
    stop("'nu' heterogeneity cannot be estimated if 'formula.w' is not defined")
  }
  
  
  # Utility
  uu          <- rep(list(matrix(0, 1, 1)), times = M)
  uv          <- rep(list(matrix(0, 1, 1)), times = M)
  uw          <- rep(list(matrix(0, 1, 1)), times = M)

  uu          <- futil(M = M, X = Xu, u = uu, theta = head(theta, nparu), hetval = ze, npar = nparu, nvec = nvec, inter = inter.u)
  uv          <- futil(M = M, X = Xv, u = uv, theta = theta[(nparu + 1):(nparu + nparv)], hetval = mu, npar = nparv, nvec = nvec, inter = inter.v)
  uw          <- futil(M = M, X = Xw, u = uw, theta = tail(theta, nparw), hetval = nu, npar = nparw, nvec = nvec, inter = inter.w)
  
  # Simulations
  combr       <- ffindcom(nblock) 
  ncombr      <- 2^nblock
  ztncombr    <- 0:(ncombr - 1)
  idrows      <- fIDdir(M, nvec)
  idcols      <- idrows$idcols
  ident       <- idrows$ident
  idrows      <- idrows$idrows
  quiet(gc())
  out         <- foreach(m = 1:M, .packages  = "mcmcERGM") %dorng% {
    fGibbdir2(init.net[[m]], nblock, ncombr, combr, idrows[[m]], idcols[[m]], ident[[m]], 
              ztncombr, uu[[m]], uv[[m]], uw[[m]], nparu, nparv, nparw, nvec[[m]], simunet)}
  
  network     <- lapply(1:M, function(m) out[[m]]$net)
  degree      <- do.call(rbind, lapply(1:M, function(m) out[[m]]$deg))
  potent      <- Reduce("+", lapply(1:M, function(m) c(out[[m]]$pot)))
  
  # Export data
  if(inter.u){
    vname.u   <- vname.u[-1]
  }
  if(inter.v){
    vname.v   <- vname.v[-1]
  }
  if(inter.w){
    vname.w   <- vname.w[-1]
  }
  
  if(!is.null(Xu)){
    Xu          <- as.data.frame(lapply(1:nvaru, function(x) mat.to.vec(lapply(1:M, function(m) Xu[[m]][,,x]))))
    colnames(Xu)<- vname.u}
  if(!is.null(Xv)){
    Xv          <- as.data.frame(lapply(1:nvarv, function(x) mat.to.vec(lapply(1:M, function(m) Xv[[m]][,,x]))))
    colnames(Xv)<- vname.v}
  if(!is.null(Xw)){
    Xw          <- as.data.frame(lapply(1:nvarw, function(x) mat.to.vec(lapply(1:M, function(m) Xw[[m]][,,x]))))
    colnames(Xw)<- vname.w}
  
  # output
  out         <- list(network = network, data = list(Xu = Xu, Xv = Xv, Xw = Xw), degree = degree, potential = potent)
  class(out)  <- "simDirNet"
  out
}
