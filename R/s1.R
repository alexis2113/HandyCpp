teamEM<-function(data,epsilon=1e-08,maxit=1000,...){
  #
  #
  #
  #
  #
  #
  #
  
  #-----------------input-checks----------------
  #only 3 levels of factors is allowed
  #everyone can try to complete this part
  if(!is.data.frame(data)|!is.vector(data$Length)){
    stop("Input data must be a dataframe with column named Length")
  }
  
  if(all(!is.numeric(data$Length))|length(data$Length)<3){
    stop("cannot find valid numeric input")
  }
  
  if(!is.vector(data$Age)){data$Age<-na}
  
  if(!is.numeric(epsilon)|epsilon>=1|epsilon<=0){
    stop("invalid epsilon value")
  }
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if(!is.wholenumber(maxit)|maxit<=0){
    stop("invalid maxit value")
  }
  #------------------------Initialisation---------------
  
  require(dplyr)
  
  train.set<-subset(data,is.na(Age))
  test.set<-na.omit(data)
  tl<-length(data$Length)
  center<-kmeans(train.set$Length,3,iter.max = 10)$centers
  center<-sort(center)
  train.set$Age<-kmeans(train.set$Length,centers = center,iter.max = 10)$cluster  
  
  if(!is.null(test.set)){
  train.set<-rbind(test.set,train.set)
  }
  
  inits<-train.set%>%
    group_by(Age)%>%
    summarise(mu=mean(Length),sigma=sd(Length),lambda=NROW(Length)/tl)%>%
    rename(Parameter=Age)%>%
    arrange(mu)
  
  inits$Parameter<-c("Age1","Age2","Age3")
  
  
  #-----preparation before looping--------
  

  
  fl<-train.set$Length
 
  mu1=inits$mu[1]
  mu2=inits$mu[2]
  mu3=inits$mu[3]
  sd1=inits$sigma[1]
  sd2=inits$sigma[2]
  sd3=inits$sigma[3]
  np1=inits$lambda[1]
  np2=inits$lambda[2]
  np3=inits$lambda[3]
  
  psum<-dnorm(fl,mean = mu1,sd=sd1)*np1+dnorm(fl,mean = mu2,sd=sd2)*np2+dnorm(fl,mean = mu3,sd=sd3)*np3
  lpsum<-sum(log(psum),na.rm = TRUE)
  
  logmax<-0
  logmax[1]<-0
  logmax[2]<-lpsum
  
  
  estimates<-data.frame(Parameter=c("Age1","Age2","Age3"),
                        mu=c(mu1,mu2,mu3),
                        sigma=c(sd1,sd2,sd3),
                        lambda=c(np1,np2,np3))
  
  #-----------start to loop-----
  k=1
  
  while(abs(logmax[k+1]-logmax[k])>epsilon){
    
    if(k-1<=maxit){
      converged<-TRUE
    }else{
      warning("cannot converged in max iteration")
      break}
    #------------------------Expectation---------------
    ge1<-dnorm(fl,mean = mu1,sd=sd1)*np1
    ge2<-dnorm(fl,mean = mu2,sd=sd2)*np2
    ge3<-dnorm(fl,mean = mu3,sd=sd3)*np3
    
    ed.sum<-ge1[!is.na(ge1)]+ge2[!is.na(ge2)]+ge3[!is.na(ge3)]
    
    g1<-ge1/ed.sum
    g2<-ge2/ed.sum
    g3<-ge3/ed.sum
    
    posterior=data.frame(observations=train.set$Length,Age1=g1,Age2=g2,Age3=g3)
    
    sum1<-sum(g1,na.rm = TRUE)
    sum2<-sum(g2,na.rm = TRUE)
    sum3<-sum(g3,na.rm = TRUE)
   #####################
    mu1<-sum(g1*fl,na.rm = TRUE)/sum1
    mu2<-sum(g2*fl,na.rm = TRUE)/sum2
    mu3<-sum(g3*fl,na.rm = TRUE)/sum3
  
    sd1<-sqrt(sum(g1*((fl-mu1)^2),na.rm = TRUE))/sqrt(sum1)
    
    sd2<-sqrt(sum(g2*((fl-mu2)^2),na.rm = TRUE))/sqrt(sum2)
    sd3<-sqrt(sum(g3*((fl-mu3)^2),na.rm = TRUE))/sqrt(sum3)
  
    np1<-sum1/tl
    np2<-sum2/tl
    np3<-sum3/tl
  
    #------------------------Maximization--------------
    psum<-dnorm(fl,mu1,sd1)*np1+dnorm(fl, mu2,sd2)*np2+dnorm(fl,mu3,sd3)*np3
    lpsum<-sum(log(psum,base = exp(1)),na.rm = TRUE)
    
    k=k+1
    logmax[k+1]<-lpsum
    
    if(is.finite(lpsum)&all(!is.na(c(mu1,mu2,mu3,np1,np2,np3,sd1,sd2,sd3)))){
      #if nothing is na then return estimates
      estimates$mu<-c(mu1,mu2,mu3)
      estimates$sigma<-c(sd1,sd2,sd3)
      estimates$lambda<-c(np1,np2,np3)
      
      next}else{
        warning("iteration stop because of numeric error")
        break}
  
  }
  #--------------------Returns-------------------
  
  likelihood<-logmax[2:k]
  message(paste(as.character(k-1),"iteration"))
  return(list(estimates=estimates,posterior=posterior,inits=inits,converged=converged,likelihood=likelihood))
  
}  
  











  