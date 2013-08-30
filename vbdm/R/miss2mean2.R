miss2mean2 <- function(x){
  #x[is.na(x)] <- -1;
  #w0 <- x==0;
  #w1 <- x==1;
  #w2 <- x==2;
  #wnk <- which(w0+w1+w2==0)
  #wk <- which(w0+w1+w2==1);
  wnk <- which(is.na(x));
  if(length(wnk)==0){
    #stop("No non-missing genotypes.\n");
    return(x);
  } else{
    #mu <- mean(x[wk]);
    mu <- mean(x,na.rm=T);
    x[wnk]<-mu;
    #print(mu);
    return(x);
  }
}