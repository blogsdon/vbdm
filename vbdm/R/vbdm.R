

#dyn.load('/home/blogsdon/Documents/pathmix/pathmix.so');
#dyn.load('/home/blogsdon/vbdm/C_implementation/vbdm_v1/pathmix.so')
dyn.load('C_implementation/vbdm_v1/pathmix.so')

compute_F_old <- function(sigma_full,sigma_reduced,n,df1,df2){

	F_num <- (sigma_reduced*n-sigma_full*n)/df1;
	F_denom <- (sigma_full*n)/df2;
	#llf <- -0.5*n*(log(2*pi*sigma_full)+1);
	#llr <- -0.5*n*(log(2*pi*sigma_reduced)+1);
	#llr <- lbr;
	#llf <- lbf;
	
	return(F_num/F_denom);


}


compute_F <- function(lbf,lbr,n,df1,df2){

	#F_num <- (sigma_reduced*n-sigma_full*n)/df1;
	#F_denom <- (sigma_full*n)/df2;
	#llf <- -0.5*n*(log(2*pi*sigma_full)+1);
	#llr <- -0.5*n*(log(2*pi*sigma_reduced)+1);
	llr <- lbr;
	llf <- lbf;
	
	#return(F_num/F_denom);
	return(-2*llr+2*llf);


}

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

pathmix <- function(y,G,X=NULL,include.mean=TRUE,eps=1e-4,thres=0.05,scaling=1,sided=1,nperm=NULL,maxit=1000,regress=1,bootstrap=FALSE,
miss.code=9){

	#preprocess covariates, thres is 
	G <- t(G);
	wmc <- which(G==miss.code);
	if(length(wmc)>0){
		G[wmc] <- NA;
	}
	if(is.null(X)){
		X <- as.matrix(rep(1,nrow(G)));
	}else{
		#X <- t(X);
		if(include.mean){
			X <- cbind(rep(1,nrow(X)),X);
		}
	}
	#preprocess G
	cs <- colMeans(G,na.rm=T)/2;
	cvec <- apply(cbind(cs,1-cs),1,min);
	keep <- which((cvec<thres)*(cvec>0)==1);
	cvec2 <- apply(cbind(cs,1-cs),1,which.min);
	G <- as.matrix(G[,keep]);
	if(length(keep)==0){
		stop("No rare variants left.\n");
	}
	flip <- which(cvec2[keep]==2);
	if(length(flip)>0){
		G[,flip] <- 2-G[,flip];
	}

	#cs <- colSums(G,na.rm=T);
	#wsi <- which(cs==1);
	#if(length(wsi)>1){
	#	si_v <- rowSums(G[,wsi],na.rm=T);
	#	G <- G[,-wsi];
	#	G <- cbind(G,si_v);
	#}
		

	G<-apply(G,2,miss2mean2);

	#nav <- which(is.na(G));	
	#if(length(nav)>0){
	#	G[nav]<-0;
	#}
	n <- length(y);
	m <- ncol(G);
	p <- ncol(X);
	Xhat <- t(solve(t(X)%*%X)%*%t(X));
	var_y <- var(y);
	#print(Xhat);
	#G <- scale(G);
	if(sided==1){
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
			as.integer(sided),
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
		model$pvec <- result[[15]];
		model$gamma <- result[[16]];
		model$theta <- result[[17]];
		model$sigma <- result[[18]];
		model$prob <- result[[19]];
		model$lb <- result[[20]];
		model$keep <- keep;
		test_null <- 1;
		result<-.C("run_pathmix_wrapper",
			as.double(eps),
			as.integer(maxit),
			as.integer(regress),
			as.integer(scaling),
			as.integer(sided),
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
		sigma_r <- result[[18]];
		lbr <- result[[20]];
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
					as.integer(sided),
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
				sigma_r <- result[[18]];
				lb_r <- result[[20]];
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
					as.integer(sided),
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
				sigma_f <- result[[18]];
				lb_f <- result[[20]];

				F_statistic <- compute_F(lb_f,lb_r,n,1,n-p-1);
				F_dist <- c(F_dist,F_statistic);
			}		
			model$Fd <- F_dist;
		}
		
		

	} else if(sided==2){


		pvec_res <- rep(0,2*m);
		gamma_res <- rep(0,p);
		theta_res <- c(0,0);
		sigma_res <- 0;
		prob_res <- c(0,0);
		lb_res <- 0;
		test_null <- 0;
		result<-.C("run_pathmix_wrapper",
			as.double(eps),
			as.integer(maxit),
			as.integer(regress),
			as.integer(scaling),
			as.integer(sided),
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


		model <- vector("list",13);
		names(model)<-c("y",
			"G",
			"X",
			"pvec_p",
			"pvec_n",
			"prob_p",
			"prob_n",
			"theta",
			"psi",
			"gamma",
			"sigma",
			"lb",
			"keep");

		model$y <- y;
		model$G <- G;
		model$X <- X;
		model$pvec_p <- result[[15]][1:m];
		model$pvec_n <- result[[15]][(m+1):(2*m)];
		model$gamma <- result[[16]];
		model$theta <- result[[17]][1];
		model$psi <- result[[17]][2];
		model$sigma <- result[[18]];
		model$prob_p <- result[[19]][1];
		model$prob_n <- result[[19]][2];
		model$lb <- result[[20]];
		model$keep <- keep;
		test_null <- 1;
		result<-.C("run_pathmix_wrapper",
			as.double(eps),
			as.integer(maxit),
			as.integer(regress),
			as.integer(scaling),
			as.integer(sided),
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
		sigma_r <- result[[18]];
		lb_r <- result[[20]];
		model$F <- compute_F(model$lb,lb_r,n,2,n-p-2);
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
					as.integer(sided),
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
				sigma_r <- result[[18]];
				lb_r <- result[[20]];

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
					as.integer(sided),
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
				sigma_f <- result[[18]];
				lb_f <- result[[20]];

				F_statistic <- compute_F(lb_f,lb_r,n,2,n-p-2);
				F_dist <- c(F_dist,F_statistic);
			}		
			model$Fd <- F_dist;
		}

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
	if(sided==2){
		model2$pvec <- cbind(model$pvec_p,model$pvec_n);
		model2$psi <- model$psi
		model2$prob <- c(model$prob_p,model$prob_n)
}
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
