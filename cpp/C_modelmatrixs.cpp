#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat C_modelmatrixs(const StringVector& x ) {
 // both FACTOR OR CHARACTER can work
 StringVector levs = Rcpp::sort_unique(x);
 
 arma::ivec B = Rcpp::match(x, levs);//transform to numeric
 arma::ivec BS = unique(B);
 int K = BS.size();//number of groups
 int N = B.n_elem;//numbers of samples
 arma::ivec nums(K);
 arma::mat XX(N,K,arma::fill::zeros);
 XX.col(0) = arma::ones<arma::mat>(N,1);
 
 int olds = 0;
 for(arma::uword i=0; i<BS.n_elem; ++i) {
     int nex = arma::accu(B.elem( arma::find(B==BS(i)) ))/BS(i);
    
     XX( olds, i, arma::size(nex, 1)) = arma::ones<arma::vec>(nex);
     olds=olds+nex;
   }
 return XX;
}





