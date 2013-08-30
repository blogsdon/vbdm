compute_F_old <- function(sigma_full,sigma_reduced,n,df1,df2){
  
  F_num <- (sigma_reduced*n-sigma_full*n)/df1;
  F_denom <- (sigma_full*n)/df2;
  #llf <- -0.5*n*(log(2*pi*sigma_full)+1);
  #llr <- -0.5*n*(log(2*pi*sigma_reduced)+1);
  #llr <- lbr;
  #llf <- lbf;
  
  return(F_num/F_denom);
  
  
}