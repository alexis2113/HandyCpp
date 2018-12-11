#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double glsp(arma::vec cof,arma::mat x,arma::vec y) {
 
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
arma::vec fc(arma::vec k){
   k.transform([](double val){ 
    return(std::tgamma(val+1)); 
    });
  return k;
  }
// [[Rcpp::export]]
arma::vec pows(arma::vec a, arma::vec b){
  arma::mat C = arma::join_rows(a,b);
  C.each_row([](arma::rowvec&x){
    x(0)=std::pow(x(0),x(1));
    return(x);
    });
  
  return C.col(0);
}


arma::vec rFunction(arma::mat x, Function f,arma::vec& args) {
  NumericVector res = f(x,args);
  return res;
}

arma::vec rFunction(arma::mat x, Function f,arma::vec& args);

// [[Rcpp::export]]
arma::mat C_derv(Function f,arma::mat X,arma::vec& args){
   int K = X.n_cols;
   int N = X.n_rows;
   X.insert_cols(0,4); 
   arma::mat yy(N,K+4);
   yy.zeros();
   arma::vec fp = {1,-8,8,-1};
   for(int i = 3; i < K+4; i++){
    
     X.col(0) =  X.col(i).head(N)-2;
     
     yy.col(i) = rFunction(X,f,args);
    
     X.col(i).head(N) += 1;
     yy.col(i+1) = rFunction(X,f,args);
     
     
     
     X.col(i).head(N) += 2;
     yy.col(i+2) = rFunction(X,f,args);
     
     X.col(i).head(N) += 1;
     yy.col(i+3) = rFunction(X,f,args);
     
    
     
     yy.col(i) = yy.cols(i,i+3)*fp;
    
     X.col(i) -= 2;
     
   }
  yy = yy.cols(0,K-1);
  yy=yy/12;
  
  return yy;
 
}



/*** R


*/
