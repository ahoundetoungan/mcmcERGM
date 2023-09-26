vec.to.mat <- function(u, N, normalise = FALSE, byrow = FALSE) {
  M        <- length(N)
  stopifnot(length(u) == sum(N*(N - 1)))
  out      <- NULL
  if (normalise) {
    out    <- frVtoMnorm(u, N, M)
  } else {
    out    <- frVtoM(u, N, M)
  }
  
  if(byrow) {
    out    <- lapply(out, t)
  }
  
  out
}


mat.to.vec <- function(W, ceiled = FALSE, byrow = FALSE) {
  if (!is.list(W)) {
    if (is.matrix(W)) {
      W    <- list(W)
    } else {
      stop("W is neither a matrix nor a list")
    }
  }
  
  M        <- length(W)
  N        <- unlist(lapply(W, nrow))

  out      <- W
  if(byrow) {
    out    <- lapply(W, t)
  }
  if (ceiled) {
    out    <- frMceiltoV(out, N, M)
  } else {
    out    <- frMtoV(out, N, M)
  }
  
  out
}