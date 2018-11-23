#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

//construct design matrix,given y and each group's size
//Inputs: 
//q response variables
//q1: size of first group
//q3: size of third group
//
// this function can not handel ties,therefore may 
// differ from R's results
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double kw(const arma::colvec q,const int& q1,const int& q3 ){
    int n = q.n_rows;
    //calculating ranks of each observation
   
    arma::uvec obex = arma::sort_index(q,"ascend");
    arma::uvec rvex = arma::sort_index(obex,"ascend");
    arma::vec ivec = arma::linspace<arma::vec>(n,0,n+1);
    ivec = ivec(rvex);//convert into scalar type
    
    double avk = mean(ivec);
    int q2 = n-q3;
    
    double avk2 = sum(ivec.head(q2));
    
    double avk1 = sum(ivec.head(q1));
    q2=q2-q1;
    
    avk2 = (avk2-avk1)/q2;
    avk1 = avk1/q1;
    
    double avk3 = mean(ivec.tail(q3));
    //std between groups
    double h1 = q1*(avk1-avk)*(avk1-avk)+ q3*(avk3-avk)*(avk3-avk)+ q2*(avk2-avk)*(avk2-avk);
    //std among groups
    double h2=arma::var(ivec);
    double hs = h1/h2;
    return hs;
    
}   



double kw(const arma::colvec q,const int& q1,const int& q3 );

//construct design matrix,given each group's size
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ModelX(const int a,const int b,const int c){ 
    int N = a+b+c;
    arma::mat A(N,3,arma::fill::zeros);
    
    A.col(0) = arma::ones<arma::mat>(N,1);
    A( a, 1, arma::size(b, 1)) = arma::ones<arma::vec>(b);
    A( a+b, 2,arma::size(c, 1)) = arma::ones<arma::vec>(c);
    return A;
    
}

arma::mat ModelX(const int a,const int b,const int c);

//
//
//A more general function
//
//combine all functions together
// work with 3 groups for now
//simple calculations of anova and kw test
//
//
//
//Inputs:
//A : size of first group
//B: size of second group
//C: size of thord group
//
//vector a,b,c: samples set for each groups
//length must longer than the total sizes of variables
//Sim : number of simulations
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List  fastf(const int A,const int B,const int C, arma::colvec a, arma::colvec b, arma::colvec& c, int Sim=10){
        RNGScope scope;	
        arma:: mat X = ModelX(A,B,C);
        const int n = X.n_rows, k = X.n_cols;
        const int q3 = C;
        const int q1 = A;
        int df2 = n - k;
        int df1 = k-1;
        arma::colvec y(n,arma::fill::zeros);
        arma::colvec temp(size(b),arma::fill::zeros);
       
        arma::colvec fs(Sim,arma::fill::zeros);
        arma::colvec hhs(Sim,arma::fill::zeros);
        for(int i = 0 ; i < Sim ; i++) {
                // Sample initial data
                temp = shuffle(b);
                y = temp.head(n);//filling 2st group
              
                temp = shuffle(a);
                
                y.head(q1) = temp.head(q1);//filling 1st group
             
                temp = shuffle(c);
                
                y.tail(q3) = temp.tail(q3);//filling 3rd group
                
                //kruskal.test
                
                double hs = kw(y,q1,q3);
                
                // linear regression
                // remain singular system solution
                arma::colvec coef = arma::solve(X, y,arma::solve_opts::allow_ugly);
                arma::colvec fitted = X*coef - mean(y);
                arma::colvec res  = y - X*coef;
                //mean sum of square
                double mtss = std::inner_product(fitted.begin(),
                                                 fitted.end(), fitted.begin(), 0.0)/df1;
               
                // std.sum of square
                double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/df2;
               // f-statistics
                double fv = mtss/s2;
  
                fs(i) = fv ;
                hhs(i) = hs;
                
        }  
        
        return Rcpp::List::create(Rcpp::Named("F") =  fs,
                             
                                  Rcpp::Named("kws") = hhs
                                
                                 );       
        
}





// Random-Sampling function
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma:: vec Csample(const int N, const int b, const arma::vec& yr){
        RNGScope scope;
        arma::uvec Iid(N);
        arma::colvec y(N);
        Iid= arma::randi<arma::uvec>(N,arma::distr_param(1, b));
        y=yr(Iid);
        return y;
}




//C-lapply

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List clapply(List input, Function f) {
        int n = input.size();
        List out(n);
        
        for(int i = 0; i < n; i++) {
                out[i] = f(input[i]);
        }
        
        return out;
}



//combine all functions together
//can only work with 3 groups for now
//only for simple calculations of anova and kw test
//
//Inputs:
//A : size of first group
//B: size of second group
//C: size of thord group
//
//vector a,b,c: samples pool for each groups
//Sim : number of simulations
//out: numbers of outliers
//
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List  fastfo(const int A,const int B,const int C, arma::colvec a, 
             arma::colvec b, arma::colvec& c, int Sim=10,int out=10){
  RNGScope scope;	
  arma:: mat X = ModelX(A,B,C);
  const int n = X.n_rows, k = X.n_cols;
  const int q3 = C;
  const int q1 = A;
  int df2 = n - k;
  int df1 = k-1;
  //allocate memories
  arma::colvec y(n,arma::fill::zeros);
  
  arma::colvec temp(size(b),arma::fill::zeros);
  //outliers locations
  arma::uvec outliers(out,arma::fill::zeros);
  
  arma::colvec fs(Sim,arma::fill::zeros);
  arma::colvec hhs(Sim,arma::fill::zeros);
  for(int i = 0 ; i < Sim ; i++) {
    // Sample initial data
    temp = shuffle(b);
    y = temp.head(n);//filling 2st group
    
    temp = shuffle(a);
    
    y.head(q1) = temp.head(q1);//filling 1st group
    
    temp = shuffle(c);
    
    y.tail(q3) = temp.tail(q3);//filling 3rd group
    
    outliers = arma::randi<arma::uvec>(out,arma::distr_param(1,n-2));
   
    
    y(outliers)= arma::randi<arma::vec>(out,arma::distr_param(0,800));
    
   
    //kruskal.test
    
    double hs = kw(y,q1,q3);
    
    // linear regression
    // remain singular system solution
    arma::colvec coef = arma::solve(X, y,arma::solve_opts::allow_ugly);
    arma::colvec fitted = X*coef - mean(y);
    arma::colvec res  = y - X*coef;
    //mean sum of square
    double mtss = std::inner_product(fitted.begin(),
                                     fitted.end(), fitted.begin(), 0.0)/df1;
    
    // std.sum of square
    double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/df2;
    // f-statistics
    double fv = mtss/s2;
    
    fs(i) = fv ;
    hhs(i) = hs;
    
  }  
  
  return Rcpp::List::create(Rcpp::Named("F") =  fs,
                            
                            Rcpp::Named("kws") = hhs);
}


