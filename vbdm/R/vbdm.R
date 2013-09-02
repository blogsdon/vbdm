vbdm <- function(y,
                 G,
                 X=NULL,
                 thres=0.05,
                 genotypes=TRUE,
                 include.mean=TRUE,
                 minor.allele=TRUE,
                 impute="MEAN",
                 eps=1e-4,
                 scaling=1,
                 nperm=NULL,
                 maxit=1000,
                 regress=1,
                 bootstrap=FALSE){

  #data dimension check
  sizeG <- dim(G);
  sizey <- length(y);
  if(!is.null(X)){
    if(!is.matrix(X)){
      stop('covariates are not in matrix form')
    }
  }
  
  if(!is.null(X)){
    sizeX <- dim(X);
    if((sizeX[1]!=sizeG[1])||(sizeX[1]!=sizey)){
      stop('sample size mismatch between covariates and phenotype or genotype')
    }
  }
  if(sizey!=sizeG[1]){
    stop('sample size mismatch between phenotype and genotype')
  }
  
  #genotype matrix check
  if(!is.matrix(G)){
    stop('genotypes are not in matrix form')
  }
    ##add mean parameter if necessary
	if(is.null(X)){
		X <- as.matrix(rep(1,nrow(G)));
	}else{
		#X <- t(X);
		if(include.mean){
			X <- cbind(rep(1,nrow(X)),X);
		}
	}
  #if G is additive encoding of genotypes perform necessary checks and
  #data reformatting
  if(genotypes){
    #genotype check
    if((sum(G<0,na.rm=TRUE)>0)||(sum(G>2,na.rm=TRUE)>0)){
      stop('genotypes with values less than 0 or greater than 2 present')
    }
    cs <- colMeans(G,na.rm=T)/2;
	  cvec <- apply(cbind(cs,1-cs),1,min);
	  keep <- which((cvec<thres)*(cvec>0)==1);
	  cvec2 <- apply(cbind(cs,1-cs),1,which.min);
	  G <- as.matrix(G[,keep]);
	  if(length(keep)==0){
  		stop("No rare variants left.\n");
  	}
    
    #flip allele encoding so minor alleles have dosage 1.
    if(minor.allele){
	    flip <- which(cvec2[keep]==2);
	    if(length(flip)>0){
    		G[,flip] <- 2-G[,flip];
  	  }
    }
   
    if(impute=="MAJOR"){
      G[is.na(G)]<-0;
    } else if (impute=="MEAN"){
      
      G<-apply(G,2,function(x){y<-is.na(x);if(sum(y)>0){x[y]<-mean(x,na.rm=T)};return(x);});
      
    } else {
      if(sum(is.na(G)>0)){
        stop('Missing data present and no imputation method specified')
      }
    }
  }else{
    if(impute=="MAJOR"){
      stop("major allele imputation specified, but design matrix is not additive encoding of genotypes")   
    } else if (impute=="MEAN"){
      G<-apply(G,2,function(x){y<-is.na(x);if(sum(y)>0){x[y]<-mean(x,na.rm=T)};return(x);});
    } else{
      if(sum(is.na(G)>0)){
        stop('Missing data present and no imputation method specified')
      }
    } 
    keep <- 1:sizeG[2];
  }

  if(sum(is.na(X))>0){
    stop('Missing data present in covariates.');
  }
  
	n <- sizey;
	m <- sizeG[2];
	p <- ncol(X);
	Xhat <- t(solve(t(X)%*%X)%*%t(X));
	var_y <- var(y);

	pvec_res <- rep(0,m);
	gamma_res <- rep(0,p);
	theta_res <- 0;
	sigma_res <- 0;
	prob_res <- 0;
	lb_res <- 0;
	test_null <- 0;
	result<-.C("run_pathmix_wrapper",
		as.double(eps),
		as.integer(maxit),
		as.integer(regress),
		as.integer(scaling),
		as.integer(test_null),
		as.double(G),
		as.double(X),
		as.double(Xhat),
		as.double(y),
		as.double(var_y),
		as.integer(n),
		as.integer(m),
		as.integer(p),
		as.double(pvec_res),
		as.double(gamma_res),
		as.double(theta_res),
		as.double(sigma_res),
		as.double(prob_res),
		as.double(lb_res));		

	model <- vector("list",10);
	names(model)<-c("y",
		"G",
		"X",
		"pvec",
		"prob",
		"theta",
		"gamma",
		"sigma",
		"lb",
		"keep");
	model$y <- y;
	model$G <- G;
	model$X <- X;
	model$pvec <- result[[14]];
	model$gamma <- result[[15]];
	model$theta <- result[[16]];
	model$sigma <- result[[17]];
	model$prob <- result[[18]];
	model$lb <- result[[19]];
	model$keep <- keep;
	test_null <- 1;
	result<-.C("run_pathmix_wrapper",
		as.double(eps),
		as.integer(maxit),
		as.integer(regress),
		as.integer(scaling),
		as.integer(test_null),
		as.double(G),
		as.double(X),
		as.double(Xhat),
		as.double(y),
		as.double(var_y),
		as.integer(n),
		as.integer(m),
		as.integer(p),
		as.double(pvec_res),
		as.double(gamma_res),
		as.double(theta_res),
		as.double(sigma_res),
		as.double(prob_res),
		as.double(lb_res));
	sigma_r <- result[[17]];
	lbr <- result[[19]];
	#print(sigma_r);
	model$F <- compute_F(model$lb,lbr,n,1,n-p-1);
	F_dist <- c();
	if(!is.null(nperm)){
		#for (k in 1:nperm){
		k <- 1;
			test_null <- 1;
			result<-.C("run_pathmix_wrapper",
				as.double(eps),
				as.integer(maxit),
				as.integer(regress),
				as.integer(scaling),
				as.integer(test_null),
				as.double(G),
				as.double(X),
				as.double(Xhat),
				as.double(y),
				as.double(var_y),
				as.integer(n),
				as.integer(m),
				as.integer(p),
				as.double(pvec_res),
				as.double(gamma_res),
				as.double(theta_res),
				as.double(sigma_res),
				as.double(prob_res),
				as.double(lb_res));
			sigma_r <- result[[17]];
			lb_r <- result[[19]];
		while(k < nperm && length(which(F_dist > model$F))<10){
			k <- k+1;
			#ynew <- sample(y,n,replace=bootstrap);
			ordn <- sample(1:n,n,replace=bootstrap);
			ynew <- y[ordn];
			Xnew <- X[ordn,];
			Xhatnew <- Xhat[ordn,];
			test_null <- 0;
			result<-.C("run_pathmix_wrapper",
				as.double(eps),
				as.integer(maxit),
				as.integer(regress),
				as.integer(scaling),
				as.integer(test_null),
				as.double(G),
				as.double(Xnew),
				as.double(Xhatnew),
				as.double(ynew),
				as.double(var_y),
				as.integer(n),
				as.integer(m),
				as.integer(p),
				as.double(pvec_res),
				as.double(gamma_res),
				as.double(theta_res),
				as.double(sigma_res),
				as.double(prob_res),
				as.double(lb_res));	
			sigma_f <- result[[17]];
			lb_f <- result[[19]];

			F_statistic <- compute_F(lb_f,lb_r,n,1,n-p-1);
			F_dist <- c(F_dist,F_statistic);
		}		
		model$Fd <- F_dist;
	}
  
	model2 <- list();
	#print(model$F);
	model2$F <- model$F
	model2$Fd <- model$Fd
	model2$X <- X;
	model2$y <- y;
	model2$p.value <- pchisq(model$F,1,lower.tail=F);
	model2$cumul.mac <- sum(colSums(G));
	#model2$post.cumul.mac <- model$pvec%*%colSums(G);
	model2$pvec <- model$pvec;
	model2$G <- G;

	model2$theta <- model$theta;
	return(model2);

}
#set.seed(11);
#n <- 1000;
#m <- 10;

#y <- rnorm(n);

#X <- matrix(rbinom(n*m,2,.01),n,m);
#y <- X%*%rbinom(m,1,.5)+rnorm(n);
#y <- rnorm(n);
#Z <- matrix(rnorm(n*3),n,3);
#ZZ <- cbind(rep(1,n),Z);
#AB <- solve(t(ZZ)%*%ZZ)%*%t(ZZ);
#print(system.time(res_1 <- pathmix(y=y,G=X,X=Z,sided=1,eps=1e-6,nperm=1000)));
#print(system.time(res_2 <- pathmix(y=y,G=X,X=Z,sided=2,eps=1e-6,nperm=1000)));

#source('~/Documents/spike_slab_gibbs/vb_rv.R');
#source('~/Documents/spike_slab_gibbs/vb_rv2.R');


#print(system.time(res1 <- run_vb_rv(y=y,G=X,X=Z,tol=1e-6)));
#print(system.time(res2 <- run_vb_rv2(y=y,G=X,X=Z,tol=1e-6)));
