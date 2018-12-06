#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double C_poi(arma::vec cof,arma::mat x,arma::vec y) {
 
 arma::vec fit = x*cof;
 arma::vec lam = exp(fit);
 arma::vec fac = exp(-lam);
 arma::vec pros= y;

 int N = y.n_rows;
 for(int i = 0 ; i < N ; i++) {
  fac(i)=fac(i)*std::pow(lam(i),y(i));
  pros(i)= std::tgamma(y(i)+1);
 }

 pros = arma::trunc_log(fac/pros);
 double pl = arma::accu(pros);
 return pl;
}

// [[Rcpp::export]]
arma::vec C_factorials(arma::vec k){
   k.transform([](double val){ 
    return(std::tgamma(val+1)); 
    });
  return k;
  }
  
  
// [[Rcpp::export]]
arma::vec C_pows(arma::vec a, arma::vec b){
  arma::mat C = arma::join_rows(a,b);
  C.each_row([](arma::rowvec&x){
    x(0)=std::pow(x(0),x(1));
    return(x);
    });
  
  return C.col(0);
}
