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
