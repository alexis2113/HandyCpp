em.initialization <- function(dataframe, varCol, catCol,...){ 
  
  # Create vectors of values for each of the known classes. 
  colnames(dataframe)[varCol]<-"variables"
  colnames(dataframe)[catCol]<-"category"
  
  test_set <- subset(dataframe,!is.na(category))#withlabels
  train_set<- subset(dataframe,is.na(category))#without labels
  
  if(is.data.frame(test_set)&!all(is.na(test_set))){
    
    means<-tapply(test_set$variables, test_set$category, mean)
    sds<- tapply(test_set$variables, test_set$category, sd)
  }
  
  
  p_matrix <- matrix(train_set$variables, nrow = NROW(train_set$variables), ncol = k)
  
  for (i in 1:k) {
    p_matrix[,i]<- dnorm(train_set$variables,means[i],sds[i])
  }
  
  p_matrix<-as.data.frame(p_matrix)
  #assign lables 
  train_set["category"] <- colnames(p_matrix)[max.col(p_matrix,ties.method="first")]
  m <- regexpr("[0-9]+", train_set$category, perl=TRUE)
  train_set["category"] <- regmatches(train_set$category, m)
  
  dataframe<-rbind(train_set,test_set)
  
  #recalculate parameters
  
  mu<-tapply(dataframe$variables, dataframe$category, mean)
  sigma<- tapply(dataframe$variables, dataframe$category, sd)
  lambda<-tapply(dataframe$variables, dataframe$category, function(x){NROW(x)/NROW(dataframe)})
                 
  return(list(mu = mu,sigma = sigma,lambda = lambda ) )                
}                                 

