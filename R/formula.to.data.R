#' @importFrom Formula as.Formula
#' @importFrom formula.tools env
#' @importFrom stats model.frame
#' @importFrom stats terms
#' @importFrom stats update
#' @importFrom stats model.response
#' @importFrom stats model.matrix
#' @importFrom stats delete.response
#' @importFrom Rcpp evalCpp
formula.to.data <- function(formula, data) {
  if (missing(data)) {
    data           <- env(formula)
  } else if(is.null(data)){
    data           <- env(formula)
  }
  
  ## Extract data from the formula
  formula          <- as.Formula(formula)
  stopifnot(length(formula)[1] == 0L, length(formula)[2] %in% 1:2)
  
  ## call model.frame()
  mf               <- model.frame(formula, data = data)
  ## extract terms, model matrices
  mtX              <- terms(formula, data = data, rhs = 1)
  X                <- model.matrix(mtX, mf)
  n.X              <- colnames(X)
  intercept        <- "(Intercept)" %in% n.X
  if(intercept){
    X              <- X[,n.X != "(Intercept)", drop = FALSE]
    n.X            <- colnames(X)
    if(ncol(X) == 0){
      X            <- NULL
      n.X          <- NULL
      intercept    <- FALSE
    }
  }
  list(X = X, names = n.X, intercept = intercept)
}
