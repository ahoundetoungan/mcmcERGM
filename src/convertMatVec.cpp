// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppEigen.h>

using namespace Rcpp;
//using namespace arma;
using namespace std;
//using namespace Eigen;
//using namespace Numer;


//typedef Eigen::Map<Eigen::MatrixXd> MapMat;
//typedef Eigen::Map<Eigen::VectorXd> MapVec;

// Create a vector from a given list a square matrixes
// The size of the length of the vector is the sum(N), where N is the vector of matrice sizes
// The elements in the generated vector are taken from column-wise (ie. the first column is filled up before filling the second column)
// and from the first matrix of the list to the last matrix of the list.
// [[Rcpp::export]]
List frVtoM(const Eigen::VectorXd& u,
            const Rcpp::IntegerVector& N,
            const double& M) {
  List out(M);
  
  int r                                = 0;
  int n;
  
  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    
    n                                  = Nm - 1;
    
    Eigen::MatrixXd outm(Eigen::MatrixXd::Zero(Nm, Nm));
    outm.block(1, 0, n, 1)             = u.segment(r, n);
    
    r                                 += n;
    for(int i(1); i < n; ++i) {
      outm.block(0, i, i, 1)          = u.segment(r, i);
      outm.block(i + 1, i, n - i, 1)  = u.segment(r + i, n - i);
      r                              += n;
    }
    
    outm.block(0, n, n, 1)            = u.segment(r, n);
    r                                += n;
    
    out[m]                            = outm;
  }
  return out;
}


// Same function but the returned matrix are normalized
// [[Rcpp::export]]
List frVtoMnorm(const arma::vec& u,
                const IntegerVector& N,
                const double& M) {
  List out(M);
  
  int r2                               = -1;
  int r1;
  
  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    
    r2                                += Nm - 1;
    r1                                 = r2 - Nm + 2;
    
    arma::mat outm(Nm, Nm, arma::fill::zeros);
    outm.submat(1, 0, Nm - 1, 0)       = u.subvec(r1, r2);
    
    for(int i(1); i < (Nm - 1); ++i) {
      r2                              += Nm - 1;
      r1                               = r2 - Nm + 2;
      outm.submat(0, i, i - 1, i)      = u.subvec(r1, r1 + i - 1);
      outm.submat(i + 1, i, Nm - 1, i) = u.subvec(r1 + i, r2);
    }
    
    r2                                += Nm - 1;
    r1                                 = r2 - Nm + 2;
    outm.submat(0, Nm - 1, Nm - 2, Nm - 1) = u.subvec(r1, r2);
    
    outm                                   = arma::normalise(outm, 1, 1);
    
    out[m]                                 = outm;
  }
  return out;
}


// Create a vector from a given list a square matrixes
// The size of the length of the vector is the sum(N), where N is the vector of matrice sizes
// The elements in the generated vector are taken from column-wise (ie. the first column is filled up before filling the second column)
// and from the first matrix of the list to the last matrix of the list.
// [[Rcpp::export]]
Eigen::VectorXd frMtoV(List& u,
                       const Rcpp::IntegerVector& N,
                       const double& M) {
  int sN                               = sum(N*N - N);
  Eigen::VectorXd out(sN);
  
  int r                                = 0;
  int n;
  
  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    Eigen::MatrixXd um                 = u[m];
    //um                                 = um.array().ceil();
    
    n                                  = Nm - 1;
    
    out.segment(r, n)                  = um.block(1, 0, n, 1);
    r                                 += n;
    for(int i(1); i < n; ++i) {
      out.segment(r, i)                = um.block(0, i, i, 1);
      out.segment(r + i, n - i)        = um.block(i + 1, i, n - i, 1);
      r                               += n;
    }
    
    out.segment(r, n)                  = um.block(0, n, n, 1);
    r                                 += n;
  }
  return out;
}

// same function but the matrixes are ceiled first
// [[Rcpp::export]]
Eigen::VectorXd frMceiltoV(List& u,
                           const Rcpp::IntegerVector& N,
                           const double& M) {
  int sN                               = sum(N*N - N);
  Eigen::VectorXd out(sN);
  
  int r                                = 0;
  int n;
  
  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    Eigen::MatrixXd um                 = u[m];
    um                                 = um.array().ceil();
    
    n                                  = Nm - 1;
    
    out.segment(r, n)                  = um.block(1, 0, n, 1);
    r                                 += n;
    for(int i(1); i < n; ++i) {
      out.segment(r, i)                = um.block(0, i, i, 1);
      out.segment(r + i, n - i)        = um.block(i + 1, i, n - i, 1);
      r                               += n;
    }
    
    out.segment(r, n)                  = um.block(0, n, n, 1);
    r                                 += n;
  }
  return out;
}