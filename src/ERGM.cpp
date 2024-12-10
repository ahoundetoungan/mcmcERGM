// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <wishart.h>
#define NDEBUG 1

using namespace Rcpp;
using namespace arma;
using namespace std;

// this function simulate the inverse Wishart distribution
//[[Rcpp::export]]
arma::mat updateOmega(const double& dfa,
                      const arma::mat& Sca,
                      const int& M,
                      const arma::mat& hete) {
  // arma::mat tp(hete.each_col() - arma::mean(hete, 1));
  return riwish(dfa + M, Sca + hete*hete.t());
}

// this function compute the variable part the dnorm
//[[Rcpp::export]]
double propdnorm(const arma::vec& x, const arma::vec& mu, const arma::mat& invV){
  arma::vec tmp = (x - mu);
  return -0.5*arma::sum(tmp%(invV*tmp));
}

// when there are many points
//[[Rcpp::export]]
arma::rowvec propdnorm_eachm(const arma::mat& x, const arma::mat& V){
  return -0.5*arma::sum(x%arma::solve(V, x), 0);
}


// when there are many points with different jumping scale
//[[Rcpp::export]]
arma::rowvec propdproposal(const arma::mat& x, const arma::mat& V, const arma::rowvec& js){
  return -0.5*arma::sum(x%arma::solve(V, x.each_row()%pow(js, 2)), 0);
}



// this function simulate theta prime
//[[Rcpp::export]]
arma::vec fsimtheta(const arma::vec& mu, const arma::mat& Sigma, const double& js, const int& npar){
  arma::vec x = arma::randn(npar, 1);
  return js*arma::chol(Sigma).t()*x + mu;
}

//[[Rcpp::export]]
arma::mat fsimhete(const arma::mat& mu, const arma::mat& Sigma, const arma::rowvec& js, const int& khete, const int& M){
  arma::mat x(khete, M, arma::fill::randn);
  return arma::chol(Sigma).t()*(x.each_row()%js) + mu;
}

// this function finds the combinations of n entries associated with x
// example, if n = 3
// 1 ---> 0 0 0  
// 2 ---> 1 0 0
// 3 ---> 0 1 0
// 4 ---> 1 1 0
// 5 ---> 0 0 1  
// 6 ---> 1 0 1
// 7 ---> 0 1 1
// 8 ---> 1 1 1
//[[Rcpp::export]]
arma::mat ffindcom(const int& n){
  int ncombr              = pow(2, n);
  arma::mat out(ncombr, n, arma::fill::zeros);
  for(int x(1); x <= ncombr; ++ x){
    int i                 = n;
    int y                 = x;
    while(i > 0){
      int tm              = pow(2, i - 1);
      if(y > tm){
        out(x - 1, i - 1) = 1;
        y                -= tm;
      }
      -- i;
    }
  }
  return out;
}

// this function organizes the data in the suitable shape
// We first compute intermediate function
// Use i characteristics
arma::mat fdatai(const arma::vec& Xk, const int& n){return arma::repmat(Xk, 1, n);}

// Use j characteristics
arma::mat fdataj(const arma::vec& Xk, const int& n){return arma::repmat(Xk.t(), n, 1);}

// sum
arma::mat fdatasum(const arma::vec& Xk, const int& n){
  arma::mat out(n, n, arma::fill::zeros);
  for(int i(0); i < (n - 1); ++ i){
    out.submat(i + 1, i, n - 1, i) = (Xk.subvec(i + 1, n - 1) + Xk(i));
  }
  return out + out.t();
}

// prod
arma::mat fdataprod(const arma::vec& Xk, const int& n){
  arma::mat out(n, n, arma::fill::zeros);
  for(int i(0); i < (n - 1); ++ i){
    out.submat(i + 1, i, n - 1, i) = (Xk.subvec(i + 1, n - 1)*Xk(i));
  }
  return out + out.t();
}

// Use Xi = Xj
arma::umat fdatasame(const arma::vec& Xk, const int& n){
  arma::umat out(n, n, arma::fill::zeros);
  for(int i(0); i < (n - 1); ++ i){
    out.submat(i + 1, i, n - 1, i) = (Xk.subvec(i + 1, n - 1) == Xk(i));
  }
  return out + out.t();
}

// abs(Xi - Xj)
arma::mat fdatadiff(const arma::vec& Xk, const int& n){
  arma::mat out(n, n, arma::fill::zeros);
  for(int i(0); i < (n - 1); ++ i){
    out.submat(i + 1, i, n - 1, i) = abs(Xk.subvec(i + 1, n - 1) - Xk(i));
  }
  return out + out.t();
}

// Xi < Xj
arma::umat fdatalower(const arma::vec& Xk, const int& n){
  arma::umat out(n, n, arma::fill::zeros);
  for(int i(0); i < n; ++ i){
    out.col(i) = (Xk < Xk(i));
  }
  return out;
}

// Xi > Xj
arma::umat fdatagreater(const arma::vec& Xk, const int& n){
  arma::umat out(n, n, arma::fill::zeros);
  for(int i(0); i < n; ++ i){
    out.col(i) = (Xk > Xk(i));
  }
  return out;
}


//[[Rcpp::export]]
arma::cube fdatar(const arma::mat X, List ftovar, const int& nvar, const int& K){
  int n(X.n_rows);
  arma::cube out(n, n, nvar);
  int s(0);
  for(int k(0); k < K; ++ k){
    arma::uvec ftv = ftovar[k];
    int L          = ftv.n_elem;
    for(int l(0); l < L; ++ l){
      switch(ftv(l)) {
      case 1: //Xi
        out.slice(s) = fdatai(X.col(k), n);
        ++ s;
        break;
      case 2: //Xj
        out.slice(s) = fdataj(X.col(k), n);
        ++ s;
        break;
      case 3: //Xi + Xj
        out.slice(s) = fdatasum(X.col(k), n);
        ++ s;
        break;
      case 4: //Xi*Xj
        out.slice(s) = fdataprod(X.col(k), n);
        ++ s;
        break;
      case 5: //Same
        out.slice(s) = arma::conv_to<arma::mat>::from(fdatasame(X.col(k), n));
        ++ s;
        break;
      case 6: //Adiff
        out.slice(s) = fdatadiff(X.col(k), n);
        ++ s;
        break;
      case 7: //1{Xi < Xj}
        out.slice(s) = arma::conv_to<arma::mat>::from(fdatalower(X.col(k), n));
        ++ s;
        break;
      case 8: //1{Xi > Xj}
        out.slice(s) = arma::conv_to<arma::mat>::from(fdatagreater(X.col(k), n));
        ++ s;
      }
    }
  }
  return out;
}

// this function computes the utilities
//[[Rcpp::export]]
arma::mat futility(const arma::cube& X, 
                   const arma::vec& theta,
                   const double& hetval,
                   const int& npar,
                   const int& n,
                   const bool& intercept){
  arma::mat out(n, n, arma::fill::zeros);
  // cout<<"intercept "<<intercept<<endl;
  for(int k(intercept); k < npar; ++ k){
    // cout<<"k - intercept "<<k - intercept<<endl;
    out += theta(k)*X.slice(k - intercept);
  } 
  out   += hetval;
  if(intercept){
    out += theta(0);
  }
  return out;
}

// this function computes the potential function for the case of symmetric networks
//[[Rcpp::export]]
double fQrsym(const arma::mat& ar, 
              const arma::mat& vr, 
              const arma::mat& wr, 
              const int& npu,
              const int& npw,
              const int& nr){
  double out   = 0;
  if(npw > 0){
    arma::vec tmp(1, arma::fill::zeros);
    for(int i(0); i < (nr - 1); ++ i){
      for(int k(i + 1); k < nr; ++ k){
        tmp   += (wr(i,k)*ar.row(i)*ar.col(k));
      }
    }
    out        = arma::accu(tmp);
  }
  
  if(npu > 0){
    out       += arma::accu(ar%vr);
  }
  
  return out;
}

// this function computes the potential function for the case of symmetric networks
//[[Rcpp::export]]
double fQrdir(const arma::mat& ar, 
              const arma::mat& ur, 
              const arma::mat& vr,
              const arma::mat& wr, 
              const int& npu,
              const int& npv,
              const int& npw,
              const int& nr){
  double out   = 0;
  if(npw > 0){
    arma::vec tmp(1, arma::fill::zeros);
    for(int i(0); i < nr; ++ i){
      for(int k(0); k < nr; ++ k){
        if(k != i){
          tmp += (wr(i,k)*ar.row(i)*ar.col(k));
        }
      }
    }
    out        = arma::accu(tmp);
  }
  
  if(npu > 0){
    out       += arma::accu(ar%ur);
  }
  
  if(npv > 0){
    // cout<<accu(ar)<<endl;
    // cout<<accu(vr)<<endl;
    out       += (0.5*arma::accu(ar%vr%ar.t()));
  }
  // cout<<out<<endl;
  
  return out;
}

// this function runs one iteration of the Gibbs on the network. 
// It updates a random samples of nblock entries
void Gibbsymi(arma::mat& arp,  // will be updated
              double& potr,   // will be updated
              const int& nblock,
              const int& ncombr, // number of combinations
              const arma::mat& combr,   
              const arma::uvec& idrows,
              const arma::uvec& idcols,
              const arma::uvec& ident,
              const arma::uvec& ztncombr, //0 to ncombr - 1
              const arma::mat& ur, 
              const arma::mat& wr, 
              const int& npu,
              const int& npw,
              const int& nr){
  arma::uvec selent = Rcpp::RcppArmadillo::sample(ident, nblock, false);
  arma::vec potrp(ncombr); //to store the potential function values
  arma::mat arp2    = arp;      //new network
  // cout<<selent.t()<<endl;
  // we only change selected entries
  for(int s(0); s < ncombr; ++ s){
    for(int m(0); m < nblock; ++ m){
      int i1        = idrows(selent(m));
      int i2        = idcols(selent(m));
      arp2(i1, i2)  = combr(s, m);
      arp2(i2, i1)  = combr(s, m);
    }
    // cout<<arp2<<endl;
    potrp(s)      = fQrsym(arp2, ur, wr, npu, npw, nr);
  }
  // probabilities
  arma::vec pro   = potrp;
  double mpro     = max(pro);
  pro             = exp(pro - mpro);
  pro            /= arma::accu(pro);
  // select one combination
  // cout<<pro.t()<<endl;
  int ss          = sum(Rcpp::RcppArmadillo::sample_main(ztncombr, 1, true, pro));
  // update the network
  for(int m(0); m < nblock; ++ m){
    int i1        = idrows(selent(m));
    int i2        = idcols(selent(m));
    arp(i1, i2)   = combr(ss, m);
    arp(i2, i1)   = combr(ss, m);
  }
  // update the value of the potential function
  potr            = potrp(ss);
}

// same function for the case of directed network
void Gibbdiri(arma::mat& arp,  // will be updated
              double& potr,   // will be updated
              const int& nblock,
              const int& ncombr,
              const arma::mat& combr,   
              const arma::uvec& idrows,
              const arma::uvec& idcols,
              const arma::uvec& ident,
              const arma::uvec& ztncombr, //0 to ncombr - 1
              const arma::mat& ur, 
              const arma::mat& vr, 
              const arma::mat& wr, 
              const int& npu,
              const int& npv,
              const int& npw,
              const int& nr){
  arma::uvec selent = Rcpp::RcppArmadillo::sample(ident, nblock, false);
  arma::vec potrp(ncombr);      //to store the potential function values
  arma::mat arp2    = arp;      //new network
  // cout<<selent.t()<<endl;
  // we only change selected entries
  for(int s(0); s < ncombr; ++ s){
    for(int m(0); m < nblock; ++ m){
      int i1        = idrows(selent(m));
      int i2        = idcols(selent(m));
      arp2(i1, i2)  = combr(s, m);
      // arp2(i2, i1)  = combr(s, m);
    }
    // cout<<arp2<<endl;
    potrp(s)      = fQrdir(arp2, ur, vr, wr, npu, npv, npw, nr);
  }
  // probabilities
  arma::vec pro   = potrp;
  double mpro     = max(pro);
  pro             = exp(pro - mpro);
  pro            /= arma::accu(pro);
  // select one combination
  // cout<<pro.t()<<endl;
  int ss          = sum(Rcpp::RcppArmadillo::sample_main(ztncombr, 1, true, pro));
  // update the network
  for(int m(0); m < nblock; ++ m){
    int i1        = idrows(selent(m));
    int i2        = idcols(selent(m));
    arp(i1, i2)   = combr(ss, m);
    // arp(i2, i1)   = combr(ss, m);
  }
  // update the value of the potential function
  potr            = potrp(ss);
}

// this function runs the Gibbs on the network fo simulate a^{prime}
//[[Rcpp::export]]
arma::mat fGibbsym(const arma::mat& ar,  // will be updated
                   const int& nblock,
                   const int& ncombr,
                   const arma::mat& combr,   
                   const arma::uvec& idrows,
                   const arma::uvec& idcols,
                   const arma::uvec& ident,
                   const arma::uvec& ztncombr, //0 to ncombr - 1
                   const arma::mat& ur, 
                   const arma::mat& wr, 
                   const int& npu,
                   const int& npw,
                   const int& nr,
                   const int& R){
  double potr   = 0;
  arma::mat arp = ar;
  for(int t(0); t < R; ++ t){
    Gibbsymi(arp, potr, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr);
  }
  return arp;
}

// same function for the case of directed network
//[[Rcpp::export]]
arma::mat fGibbdir(const arma::mat& ar,  // will be updated
                   const int& nblock,
                   const int& ncombr,
                   const arma::mat& combr,   
                   const arma::uvec& idrows,
                   const arma::uvec& idcols,
                   const arma::uvec& ident,
                   const arma::uvec& ztncombr, //0 to ncombr - 1
                   const arma::mat& ur, 
                   const arma::mat& vr, 
                   const arma::mat& wr, 
                   const int& npu,
                   const int& npv,
                   const int& npw,
                   const int& nr,
                   const int& R){
  double potr   = 0;
  arma::mat arp = ar;
  for(int t(0); t < R; ++ t){
    Gibbdiri(arp, potr, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, vr, wr, npu, npv, npw, nr);
  }
  return arp;
}

// this function runs the Gibbs on the network fo simulate a^{prime} and exports the degree
//[[Rcpp::export]]
List fGibbsym2(const arma::mat& ar,  // will be updated
               const int& nblock,
               const int& ncombr,
               const arma::mat& combr,   
               const arma::uvec& idrows,
               const arma::uvec& idcols,
               const arma::uvec& ident,
               const arma::uvec& ztncombr, //0 to ncombr - 1
               const arma::mat& ur, 
               const arma::mat& wr, 
               const int& npu,
               const int& npw,
               const int& nr,
               const int& R){
  double potr   = 0;
  arma::mat arp = ar;
  arma::mat deg(nr, R);
  arma::vec outpot(R);
  for(int t(0); t < R; ++ t){
    Gibbsymi(arp, potr, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, wr, npu, npw, nr);
    deg.col(t)  = arma::sum(arp, 1);
    outpot(t)   = potr;
  }
  return List::create(Named("net") = arp, Named("deg") = deg, Named("pot") = outpot);
}

// same function for the case of directed network
//[[Rcpp::export]]
List fGibbdir2(const arma::mat& ar,  // will be updated
               const int& nblock,
               const int& ncombr,
               const arma::mat& combr,   
               const arma::uvec& idrows,
               const arma::uvec& idcols,
               const arma::uvec& ident,
               const arma::uvec& ztncombr, //0 to ncombr - 1
               const arma::mat& ur, 
               const arma::mat& vr, 
               const arma::mat& wr, 
               const int& npu,
               const int& npv,
               const int& npw,
               const int& nr,
               const int& R){
  double potr   = 0;
  arma::mat arp = ar;
  arma::mat deg(nr, R);
  arma::vec outpot(R);
  for(int t(0); t < R; ++ t){
    Gibbdiri(arp, potr, nblock, ncombr, combr, idrows, idcols, ident, ztncombr, ur, vr, wr, npu, npv, npw, nr);
    deg.col(t)  = arma::sum(arp, 1);
    outpot(t)   = potr;
  }
  return List::create(Named("net") = arp, Named("deg") = deg, Named("pot") = outpot);
}

// this function computes some identifiers for the MCMC
//[[Rcpp::export]]
List fIDsym(const int& M, const arma::vec& nvec){
  List idrows(M), idcols(M), ident(M);
  for(int m(0); m < M; ++ m){
    int nm    = nvec(m);
    int tmp   = nm*(nm - 1)/2;
    int i1    = 0;
    arma::vec idr(tmp), idc(tmp);
    for(int i(0); i < (nm - 1); ++ i){
      int i2  = i1 + nm - i - 2;
      idr.subvec(i1, i2) = i*arma::ones(nm - i - 1);
      idc.subvec(i1, i2) = arma::regspace(i + 1, nm - 1);
      i1      = i2 + 1;
    }
    arma::vec ide        = arma::regspace(0, tmp - 1);
    idrows[m] = idr;
    idcols[m] = idc;
    ident[m]  = ide;
  }
  return List::create(Named("idrows") = idrows, Named("idcols") = idcols, Named("ident") = ident);
}

// same function for the case of directed network
//[[Rcpp::export]]
List fIDdir(const int& M, const arma::vec& nvec){
  List idrows(M), idcols(M), ident(M);
  for(int m(0); m < M; ++ m){
    int nm    = nvec(m);
    int tmp   = nm*(nm - 1);
    int i1    = 0;
    arma::vec idr(tmp), idc(tmp);
    for(int i(0); i < nm; ++ i){
      int i2  = i1 + nm - 2;
      idr.subvec(i1, i2) = i*arma::ones(nm - 1);
      arma::vec tmp2 = arma::regspace(0, nm - 1); tmp2.shed_row(i);
      idc.subvec(i1, i2) = tmp2;
      i1      = i2 + 1;
    }
    arma::vec ide        = arma::regspace(0, tmp - 1);
    idrows[m] = idr;
    idcols[m] = idc;
    ident[m]  = ide;
  }
  return List::create(Named("idrows") = idrows, Named("idcols") = idcols, Named("ident") = ident);
}

// This function updates jumping scales
//[[Rcpp::export]]
double fupdate_jstheta(const double& jscal, 
                  const double& accept, 
                  const int& iteration, 
                  const double& target, 
                  const double& kappa, 
                  const double& jmin, 
                  const double& jmax){
  double out(jscal + (accept/iteration - target)/pow(iteration, kappa));
  if(out < jmin) return jmin;
  if(out > jmax) return jmax;
  return out;
}

//[[Rcpp::export]]
arma::rowvec fupdate_jshete(const arma::rowvec& jscal, 
                  const arma::rowvec& accept, 
                  const int& iteration, 
                  const double& target, 
                  const double& kappa, 
                  const double& jmin, 
                  const double& jmax){
  arma::rowvec out(jscal + (accept/iteration - target)/pow(iteration, kappa));
  return out.clamp(jmin, jmax);
}

// this function recenters heterogeneity to ensure identification
//[[Rcpp::export]]
arma::mat frecentering(arma::vec& theta,
                       const arma::mat& hetval,
                       const arma::uvec intindex) {
  arma::vec tp(arma::mean(hetval, 1));
  theta.elem(intindex) += tp;
  return hetval.each_col() - tp;
}