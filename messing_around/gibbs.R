






gibbs_sampleing <- function(X,N,nvars=1,maxit=2000,nth=1,k,y,print=0,burnin=0) {
  #solve matrix of xi with k variables
  solutions<-solve(t(X)%*%X)
  smatrix<-solutions*(N/(N+1))
  covx=chol(smatrix)
  
  
  beta<-rep(0,nvars) 
  z<-qnorm((rank(as.numeric(y),ties.method = "random")/(N+1)),mean = 0,sd=1)
  g<-rep(NA,k-1)
  glist<-matrix(NA,(maxit/nth),k-1)
  betas<-matrix(NA,(maxit/nth),nvars) 
  
  mixing<-0
  mu<-rep(0,k-1) 
  sigma<-rep(100,k-1)
  
  for(i in 1:maxit) 
  {
        #get upper and lower bounds
        for(j in 1:(k-1)) 
        {
          #Reset g 
          lower<-max(z[as.factor(y)==j])
          upper<-min(z[as.factor(y)==j+1])
          u<-runif(1, min=pnorm(((lower-mu[j])/sigma[j])), max=pnorm((upper-mu[j])/sigma[j]) )
          g[j] <-  mu[j] + sigma[j]*qnorm(u)
        }
        
        #recalculate betas 
        means<- smatrix%*%( t(X)%*%z )
        beta<- covx%*%rnorm(nvars)+means
        
        # Sample u from ∼ uniform distributions
        ez <- X %*% t(beta)
        lower <- c(-1000000000,g)[match(as.numeric(y)-1,0:k)]
        upper <- c(g,1000000000)[as.numeric(y)]  
        u <- runif(N, pnorm(lower-ez),pnorm(upper-ez))
        z <-  ez + qnorm(u)
        
        
        #
        mixs <- rnorm(1,0,N^(-1/3))  
        zp <- z+ mixs ; gp <-g+mixs
        lhr <-  sum(dnorm(zp,mean=ez,sd=1,log=T) - dnorm(z,mean=ez,sd=1,log=T),na.rm = TRUE ) + 
          sum(dnorm(gp,mu,sigma,log=T) - dnorm(g,mu,sigma,log=T) ,na.rm = TRUE)
        if(log(runif(1))<lhr) { z<-zp ; g<-gp ; mixing<-mixing+1 }
        
        if(i%%nth==0) 
        { 
          
          betas[i/nth,] <-  beta
          zs<- z
          glist[i/nth,] <- g
        }
        
        if(print > 0 & i %% print == 0){
          cat(i/maxit,mixing/i,"\n")
        }
  } 
  #------output list-----
  posterior = as.data.frame(cbind(betas[burnin:], glist[burnin: ] ))
   for(i in 1:ncol(posterior)){
     if(i<=nvars){
       names(posterior)[i] = paste("Beta",i,sep="")
     }else{
       names(posterior)[i] =paste("cuts",i-nvars,sep="")
       }
   } 
  # transform output into a mcmc object in order to use other package 
  # to perform diagnostic test
  outputs <- mcmc(posterior)
  
  return(outputs)
}






