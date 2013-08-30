//variational Bayes rare variant mixture test (pathmix) R package C code
//Copyright 2013 Benjamin A Logsdon
#include "vbdm.h"

//get a column of genotype matrix G
inline double * gc(struct model_struct * model, int j){
	return (&(model->data.G[j]))->col;
}


//get a column of covariate matrix X
inline double * xc(struct model_struct * model, int j){
	return (&(model->data.X[j]))->col;
}

//get a column of the Xhat matrix
inline double * hc(struct model_struct * model, int j){
	return (&(model->data.Xhat[j]))->col;
}


inline void ddot_w(int n,double *vect1,double *vect2,double * result){
	const int incxy = 1;
	(*result)=F77_NAME(ddot)(&n,vect1,&incxy,vect2,&incxy);
}


inline void daxpy_w(int n,double *x,double *y,double alpha){
	//y<- ax+y;
	const int incxy =1;
	F77_NAME(daxpy)(&n,&alpha,x,&incxy,y,&incxy);
}

inline void dnrm2_w(int n,double *x,double *result){
	const int incxy=1;
	(*result)=F77_NAME(dnrm2)(&n,x,&incxy);
}

inline void dscal_w(int n,double *x, double alpha){
	const int incxy=1;
	F77_NAME(dscal)(&n,&alpha,x,&incxy);
}


void scale_vector(double * vec,double * ones,int n){
	//mean zero, variance 1
	double mean,sd;
	double nd = (double) n;
	ddot_w(n,vec,ones,&mean);
	mean = mean/nd;
	daxpy_w(n,ones,vec,-mean);
	dnrm2_w(n,vec,&sd);
	dscal_w(n,vec,sqrt(nd-1)/(sd));
}

inline double compute_ssq(double *vec,int n){
	double a;
	ddot_w(n,vec,vec,&a);
	return a;
}

void process_data(struct model_struct * model){
	int j;
	double nd = ((double) model->data.n);

	switch(model->control_param.scaleType){

		case STANDARD:
			//Rprintf("Scaling...\n");
			
			for(j=0;j<model->data.m;j++){
				//if(j>0){
				scale_vector(gc(model,j),model->data.one_vec,model->data.n);
				//}
				model->data.g_sum_sq[j] = nd-1;
			}
			break;

		case NOSCALE:
			//Rprintf("Sum of squares pre-compute...\n");
			for(j=0;j<model->data.m;j++){
				model->data.g_sum_sq[j]=compute_ssq(gc(model,j),model->data.n);
				//Rprintf("gssq[%d]: %g\n",j,model->data.g_sum_sq[j]);
			}
			break;
	}


}


void initialize_model(double * eps, 
			int * maxit, 
			int * regress, 
			int * scale, 
			int * sided, 
			int * test_null,
			double * G, 
			double * X,
			double * Xhat,
			double * y, 
			double * var_y, 
			int * n, 
			int * m, 
			int * p,
			struct model_struct * model){
	int k,l;
	model->control_param.eps = (*eps);
	model->control_param.maxit = (*maxit);
	
	if((*regress)==1){
		model->control_param.regressType = LINEAR;
	} else{
		model->control_param.regressType = LOGISTIC;
	}

	if((*scale)==1){
		model->control_param.scaleType = STANDARD;
	} else{
		model->control_param.scaleType = NOSCALE;
	}

	if((*sided)==1){
		model->control_param.testType = ONESIDED;
	} else{
		model->control_param.testType = TWOSIDED;
	}

	model->control_param.test_null = (*test_null);

	model->data.G = (struct matrix_v *) malloc(sizeof(struct matrix_v)*(*m));
	for(k=0;k<(*m);k++){
		(&(model->data.G[k]))->col = (double *) malloc(sizeof(double)*(*n));
	}
	for(k=0;k<(*m);k++){
		for(l=0;l<(*n);l++){
			(&(model->data.G[k]))->col[l] = G[k*(*n)+l];
		}
	}

	//Rprintf("here1\n");
	model->data.X= (struct matrix_v *) malloc(sizeof(struct matrix_v)*(*p));
	for(k=0;k<(*p);k++){
		(&(model->data.X[k]))->col = (double *) malloc(sizeof(double)*(*n));
	}

	for(k=0;k<(*p);k++){
		for(l=0;l<(*n);l++){
			(&(model->data.X[k]))->col[l] = X[k*(*n)+l];
		}
	}
	//Rprintf("here2\n");
	model->data.Xhat= (struct matrix_v *) malloc(sizeof(struct matrix_v)*(*p));
	for(k=0;k<(*p);k++){
		(&(model->data.Xhat[k]))->col = (double *) malloc(sizeof(double)*(*n));
	}

	for(k=0;k<(*p);k++){
		for(l=0;l<(*n);l++){
			(&(model->data.Xhat[k]))->col[l] = Xhat[k*(*n)+l];
		}
	}
	//Rprintf("here3\n");
	model->data.y = y;

	model->data.var_y = (*var_y);

	model->data.n = (*n);
	model->data.m = (*m);
	model->data.p = (*p);
	model->data.g_sum_sq = (double *) malloc(sizeof(double)*(*m));

	model->data.one_vec = (double *) malloc(sizeof(double)*(*n));
	for(k=0;k<(*n);k++){
		model->data.one_vec[k]= 1.0;
	}
	//Rprintf("here4\n");
	process_data(model);

	model->model_param.pvec_p = (double *) malloc(sizeof(double)*(*m));
	model->model_param.pvec_n = (double *) malloc(sizeof(double)*(*m));

	//Rprintf("here5\n");
	switch(model->control_param.testType){
		case ONESIDED:
			//Rprintf("onesided\n");
			for(k=0;k<(*m);k++){	
				model->model_param.pvec_p[k] = 0.5;
				model->model_param.pvec_n[k] = 0.5;
			}
			break;
		case TWOSIDED:
			//Rprintf("twosided\n");
			for(k=0;k<(*m);k++){	
				model->model_param.pvec_p[k] = 1.0/3.0;
				model->model_param.pvec_n[k] = 1.0/3.0;
			}		
			break;
	}
	//Rprintf("here6\n");
	model->model_param.gamma = (double *) malloc(sizeof(double)*(*p));
	for(k=0;k<(*p);k++){
		model->model_param.gamma[k] = 0.0;
	}
	//Rprintf("here7\n");
	model->model_param.theta = (double *) malloc(sizeof(double)*2);
	for(k=0;k<2;k++){
		model->model_param.theta[k]= 0.0;
	}
	//Rprintf("here8\n");
	model->model_param.entropy = (double *) malloc(sizeof(double)*3);
	for(k=0;k<3;k++){
		model->model_param.entropy[k]=0.0;
	}

	model->model_param.sigma = (*var_y);

	model->model_param.prob = (double *) malloc(sizeof(double)*2);
	switch(model->control_param.testType){
		case ONESIDED:
			for(k=0;k<2;k++){
				model->model_param.prob[k]= 0.5;
			}
			break;
		case TWOSIDED:
			for(k=0;k<2;k++){
				model->model_param.prob[k]= 1.0/3.0;
			}
			break;
	}
	//Rprintf("here9\n");
	model->model_param.resid_vec = (double *) malloc(sizeof(double)*(*n));
	model->model_param.Gp = (double *) malloc(sizeof(double)*(*n));
	model->model_param.Gn = (double *) malloc(sizeof(double)*(*n));
	for(k=0;k<(*n);k++){
		model->model_param.resid_vec[k] = y[k];
		model->model_param.Gp[k] = 0;
		model->model_param.Gn[k] = 0;
	}
	//Rprintf("here10, x1: %g, %g, %g\n",gc(model,0)[0],gc(model,0)[1],gc(model,0)[2]);
	//Rprintf("here10, pvec_p[0]: %g\n",model->model_param.pvec_p[0]);
	for(k=0;k<(*m);k++){
		daxpy_w((*n),gc(model,k),model->model_param.Gp,model->model_param.pvec_p[k]);
		daxpy_w((*n),gc(model,k),model->model_param.Gn,model->model_param.pvec_n[k]);
	}
	//Rprintf("here11\n");	
	model->model_param.lb = -1e100;
	model->model_param.psum_p = 0.0;
	model->model_param.psum_n = 0.0;
	model->model_param.vsum_p = 0.0;
	model->model_param.vsum_n = 0.0;
	model->model_param.covsum = 0.0;

}


void free_model(struct model_struct * model){

	int k;
	for(k=0;k<(model->data.m);k++){
		free((&(model->data.G[k]))->col);
		
	}
	free(model->data.G);
	for(k=0;k<(model->data.p);k++){
		free((&(model->data.X[k]))->col);
		
	}
	free(model->data.X);

	for(k=0;k<(model->data.p);k++){
		free((&(model->data.Xhat[k]))->col);
		
	}
	free(model->data.Xhat);


	free(model->data.g_sum_sq);
	free(model->data.one_vec);

	free(model->model_param.pvec_p);
	free(model->model_param.pvec_n);
	free(model->model_param.gamma);
	free(model->model_param.theta);
	free(model->model_param.prob);
	free(model->model_param.resid_vec);
	free(model->model_param.Gp);
	free(model->model_param.Gn);

}

void update_p(struct model_struct * model){
	int k;
	double pold,vec1,a1,a2,a3,a4,a5,pnew;
	double pold_p,pold_n,a1p,a1n,a2p,a2n,a3p,a3n,POS,NEG,pnew_p,pnew_n;
	double md = (double) model->data.m;
	switch(model->control_param.testType){
		case ONESIDED:
			for(k=0;k<model->data.m;k++){			
				pold = model->model_param.pvec_p[k];
				ddot_w(model->data.n,model->model_param.resid_vec,gc(model,k),&vec1);
				vec1 = vec1 + model->data.g_sum_sq[k]*pold*model->model_param.theta[0];
				a1 = pow(model->model_param.theta[0],2)*model->data.g_sum_sq[k];
				a2 = -(2*model->model_param.theta[0]*vec1);
				a3 = -(digamma((model->model_param.prob[0])*(md+2))-digamma(md+2));
				a4 = digamma((1-model->model_param.prob[0])*(md+2))-digamma(md+2);
				a5 = (1/(2*model->model_param.sigma))*(a1+a2) + a3 + a4;
				pnew = 1/(1+exp(a5));
				model->model_param.pvec_p[k] = pnew;
				daxpy_w(model->data.n,gc(model,k),model->model_param.resid_vec,model->model_param.theta[0]*(pold-pnew));
				daxpy_w(model->data.n,gc(model,k),model->model_param.Gp,pnew-pold);
				model->model_param.psum_p = model->model_param.psum_p + pnew;
				model->model_param.vsum_p = model->model_param.vsum_p + model->data.g_sum_sq[k]*(pnew-pow(pnew,2));
				if(pnew==1){
					//model->model_param.entropy[0] = model->model_param.entropy[0]-pnew*log(pnew);
					//model->model_param.entropy[1] = model->model_param.entropy[1]-(1-pnew)*log(1-pnew);
				}else if(pnew==0){
					//model->model_param.entropy[0] = model->model_param.entropy[0]-pnew*log(pnew);
					//model->model_param.entropy[1] = model->model_param.entropy[1]-(1-pnew)*log(1-pnew);
				} else {
					model->model_param.entropy[0] = model->model_param.entropy[0]-pnew*log(pnew);
					model->model_param.entropy[1] = model->model_param.entropy[1]-(1-pnew)*log(1-pnew);
				}
			}
			model->model_param.prob[0] = (model->model_param.psum_p+1)/(md+2);

		break;

		case TWOSIDED:

			for(k=0;k<model->data.m;k++){
				pold_p = model->model_param.pvec_p[k];
				pold_n = model->model_param.pvec_n[k];
				//vec1 <- t(model$resid)%*%model$G[,j]+model$gssq[j]*(pold_p*model$theta) - model$gssq[j]*(pold_n*model$psi);
				ddot_w(model->data.n,model->model_param.resid_vec,gc(model,k),&vec1);
				vec1 = vec1 + model->data.g_sum_sq[k]*pold_p*model->model_param.theta[0];
				vec1 = vec1 - model->data.g_sum_sq[k]*pold_n*model->model_param.theta[1];
				a1p = pow(model->model_param.theta[0],2)*model->data.g_sum_sq[k];
				a1n = pow(model->model_param.theta[1],2)*model->data.g_sum_sq[k];
				a2p = -(2*model->model_param.theta[0]*vec1);
				a2n = (2*model->model_param.theta[1]*vec1);
				a3p = -(digamma((model->model_param.prob[0])*(md+3))-digamma(md+3));
				a3n = -(digamma((model->model_param.prob[1])*(md+3))-digamma(md+3));
				a4 = digamma((1-model->model_param.prob[0]-model->model_param.prob[1])*(md+3))-digamma(md+3);
				POS = (1/(2*model->model_param.sigma))*(a1p+a2p) + a3p;
				NEG = (1/(2*model->model_param.sigma))*(a1n+a2n) + a3n;
				pnew_p = 1/(1+exp(POS-NEG)+exp(POS+a4));
				pnew_n = 1/(1+exp(NEG-POS)+exp(NEG+a4));
				

				model->model_param.pvec_p[k] = pnew_p;
				model->model_param.pvec_n[k] = pnew_n;
				daxpy_w(model->data.n,gc(model,k),model->model_param.resid_vec,model->model_param.theta[0]*(pold_p-pnew_p));
				daxpy_w(model->data.n,gc(model,k),model->model_param.resid_vec,model->model_param.theta[1]*(pnew_n-pold_n));
				daxpy_w(model->data.n,gc(model,k),model->model_param.Gp,pnew_p-pold_p);
				daxpy_w(model->data.n,gc(model,k),model->model_param.Gn,pnew_n-pold_n);

				model->model_param.psum_p = model->model_param.psum_p + pnew_p;
				model->model_param.psum_n = model->model_param.psum_n + pnew_n;

				model->model_param.vsum_p = model->model_param.vsum_p + model->data.g_sum_sq[k]*(pnew_p-pow(pnew_p,2));
				model->model_param.vsum_n = model->model_param.vsum_n + model->data.g_sum_sq[k]*(pnew_n-pow(pnew_n,2));

				model->model_param.covsum = model->model_param.covsum + model->data.g_sum_sq[k]*(pnew_p*pnew_n);
				model->model_param.entropy[0] = model->model_param.entropy[0]-pnew_p*log(pnew_p);
				model->model_param.entropy[1] = model->model_param.entropy[1]-pnew_n*log(pnew_n);
				model->model_param.entropy[2] = model->model_param.entropy[2]-(1-pnew_p-pnew_n)*log(1-pnew_p-pnew_n);
				//if(k==0){
				//	Rprintf("%g,%g,%g,%g,%g,%g,%g\n",pold_n,pnew_n,a1n,a2n,a3n,NEG,a4);
				//}
			}
			model->model_param.prob[0] = (model->model_param.psum_p+1)/(md+3);
			model->model_param.prob[1] = (model->model_param.psum_n+1)/(md+3);			
		break;

	}
}

void update_theta_gamma(struct model_struct * model){


	double theta_old, theta_new, const1,const2;
	double psi_old,psi_new;
	int p = model->data.p;
	int k;
	double gamma_old[p];
	double gamma_new;
	for (k=0;k<p;k++){
		gamma_old[k] = model->model_param.gamma[k];
	}
	
	switch(model->control_param.testType){
		case ONESIDED:
			//Rprintf("test_null: %d\n",model->control_param.test_null);
			//update theta
			if(model->control_param.test_null==1){
				theta_new = 0.0;
			}else{
				theta_old = model->model_param.theta[0];
				ddot_w(model->data.n,model->model_param.resid_vec,model->model_param.Gp,&theta_new);
				ddot_w(model->data.n,model->model_param.Gp,model->model_param.Gp,&const1);

				theta_new = theta_new + const1*theta_old;
				const2 = const1+model->model_param.vsum_p;
				theta_new = theta_new/const2;
			
				model->model_param.theta[0] = theta_new;
				daxpy_w(model->data.n,model->model_param.Gp,model->model_param.resid_vec,theta_old-theta_new);
			}

			//update gamma
			for (k=0;k<p;k++){
				ddot_w(model->data.n,model->model_param.resid_vec,hc(model,k),&gamma_new);
				ddot_w(model->data.n,hc(model,k),xc(model,k),&const1);
				//Rprintf("gamma_new: %g, const1: %g\n",gamma_new,const1);
				gamma_new = gamma_new + const1*gamma_old[k];
				//Rprintf("gamma_new: %g, const1: %g gamma_old[%d]: %g\n",gamma_new,const1,k,gamma_old[k]);
				//Rprintf("hc[0]: %g, hc[1]: %g, hc[2]: %g\n",hc(model,k)[0],hc(model,k)[1],hc(model,k)[2]);
				model->model_param.gamma[k]=gamma_new;
				daxpy_w(model->data.n,xc(model,k),model->model_param.resid_vec,gamma_old[k]-gamma_new);
			}
		break;
	
		case TWOSIDED:
			//update theta
			if(model->control_param.test_null==1){
				theta_new = 0.0;
			}else{
				theta_old = model->model_param.theta[0];
				ddot_w(model->data.n,model->model_param.resid_vec,model->model_param.Gp,&theta_new);
				ddot_w(model->data.n,model->model_param.Gp,model->model_param.Gp,&const1);
				theta_new = theta_new + const1*theta_old;
				const2 = const1+model->model_param.vsum_p;
				theta_new = theta_new/const2;
			
				model->model_param.theta[0] = theta_new;
	
				daxpy_w(model->data.n,model->model_param.Gp,model->model_param.resid_vec,theta_old-theta_new);
			}

			//update psi
			if(model->control_param.test_null==1){
				psi_new = 0.0;
			}else{				
				psi_old = model->model_param.theta[1];
				ddot_w(model->data.n,model->model_param.resid_vec,model->model_param.Gn,&psi_new);
				ddot_w(model->data.n,model->model_param.Gn,model->model_param.Gn,&const1);
				psi_new = const1*psi_old -psi_new;
				const2 = const1+model->model_param.vsum_n;
				psi_new = psi_new/const2;
				model->model_param.theta[1] = psi_new;
				daxpy_w(model->data.n,model->model_param.Gn,model->model_param.resid_vec,psi_new-psi_old);
			}


			//update gamma
			for (k=0;k<p;k++){
				ddot_w(model->data.n,model->model_param.resid_vec,hc(model,k),&gamma_new);
				ddot_w(model->data.n,hc(model,k),xc(model,k),&const1);
				gamma_new = gamma_new + const1*gamma_old[k];
				model->model_param.gamma[k]=gamma_new;
				daxpy_w(model->data.n,xc(model,k),model->model_param.resid_vec,gamma_old[k]-gamma_new);
			}


		break;


	}

}

void update_sigma(struct model_struct * model){
	double sigma;
	double nd = (double) model->data.n;
	switch(model->control_param.testType){
		case ONESIDED:
			ddot_w(model->data.n,model->model_param.resid_vec,model->model_param.resid_vec,&sigma);
			sigma = sigma+model->model_param.vsum_p*pow(model->model_param.theta[0],2);
			sigma = sigma/nd;
			model->model_param.sigma = sigma;
		break;
	
		case TWOSIDED:

			ddot_w(model->data.n,model->model_param.resid_vec,model->model_param.resid_vec,&sigma);
			sigma= sigma+model->model_param.vsum_p*pow(model->model_param.theta[0],2);
			sigma= sigma+model->model_param.vsum_n*pow(model->model_param.theta[1],2);
			sigma= sigma-2*model->model_param.covsum*model->model_param.theta[0]*model->model_param.theta[1];
			sigma = sigma/nd;
			model->model_param.sigma = sigma;
		break;
	}
}


void update_lb(struct model_struct * model){

	double lb,alpha1,beta1,gamma1,a0;
	double nd = (double) model->data.n;
	double md = (double) model->data.m;
	switch(model->control_param.testType){
		case ONESIDED:
			lb = -0.5*(nd*(log(2*M_PI*model->model_param.sigma)+1));
			//Rprintf("lb ll: %g\n",lb);
			lb = lb + (digamma((model->model_param.prob[0])*(md+2))-digamma(md+2))*(model->model_param.psum_p+1);
			lb = lb + (digamma((1-model->model_param.prob[0])*(md+2))-digamma(md+2))*(md-model->model_param.psum_p+1);
			//Rprintf("lb elp: %g\n",lb);
			lb = lb + model->model_param.entropy[0];
			lb = lb + model->model_param.entropy[1];
			//Rprintf("lb entropy: %g\n",lb);
			alpha1 = model->model_param.psum_p +1;
			beta1 = md - model->model_param.psum_p+1;
			
			lb = lb + lbeta(alpha1,beta1) - (alpha1-1)*digamma(alpha1)-(beta1-1)*digamma(beta1)+(alpha1+beta1-2)*digamma(md+2);
			//Rprintf("lb entropy beta: %g\n",lb);
			model->model_param.lb = lb;

		break;
	
		case TWOSIDED:
			lb = -0.5*(nd*(log(2*M_PI*model->model_param.sigma)+1));
			lb = lb + (digamma((model->model_param.prob[0])*(md+3))-digamma(md+3))*(model->model_param.psum_p+1);
			lb = lb + (digamma((model->model_param.prob[1])*(md+3))-digamma(md+3))*(model->model_param.psum_n+1);
			lb = lb + (digamma((1-model->model_param.prob[0]-model->model_param.prob[1])*(md+3))-digamma(md+3))*(md-model->model_param.psum_p-model->model_param.psum_n+1);
			lb = lb + model->model_param.entropy[0];
			lb = lb + model->model_param.entropy[1];
			lb = lb + model->model_param.entropy[2];
			alpha1 = model->model_param.psum_p+1;
			beta1 = model->model_param.psum_n+1;
			gamma1 = md-model->model_param.psum_p-model->model_param.psum_n+1;
			a0 = alpha1+beta1+gamma1;
			lb = lb + lgamma(alpha1)+lgamma(beta1)+lgamma(gamma1)-lgamma(a0);
			lb = lb + digamma(a0)*(a0-3)- (alpha1-1)*digamma(alpha1)-(beta1-1)*digamma(beta1)-(gamma1-1)*digamma(gamma1);
			model->model_param.lb = lb;

		break;
	}
}


void collapse_results(struct model_struct * model,
		double * pvec_res,
		double * gamma_res,
		double * theta_res,
		double * sigma_res,
		double * prob_res,
		double * lb_res){

	int k;
	int n = model->data.n;
	int m = model->data.m;
	int p = model->data.p;
	
	switch(model->control_param.testType){
		case ONESIDED:
			for(k=0;k<m;k++){
				pvec_res[k] = model->model_param.pvec_p[k];
			}
			for(k=0;k<p;k++){
				gamma_res[k] = model->model_param.gamma[k];
			}
			theta_res[0] = model->model_param.theta[0];
			sigma_res[0] = model->model_param.sigma;
			prob_res[0] = model->model_param.prob[0];
			lb_res[0] = model->model_param.lb;
		break;
		case TWOSIDED:
			for(k=0;k<m;k++){
				pvec_res[k] = model->model_param.pvec_p[k];
			}
			for(k=m;k<(2*m);k++){
				pvec_res[k] = model->model_param.pvec_n[k-m];
			}
			for(k=0;k<p;k++){
				gamma_res[k] = model->model_param.gamma[k];
			}
			theta_res[0] = model->model_param.theta[0];
			theta_res[1] = model->model_param.theta[1];
			sigma_res[0] = model->model_param.sigma;
			prob_res[0] = model->model_param.prob[0];
			prob_res[1] = model->model_param.prob[1];
			lb_res[0] = model->model_param.lb;
		break;
	}

}




void run_pathmix(struct model_struct * model){
	double tol=1;
	double lb_old;
	int count = 0;
	//Rprintf("tol: %g, eps: %g\n",fabs(tol),model->control_param.eps);
	while(fabs(tol)>model->control_param.eps && count < model->control_param.maxit){
		lb_old = model->model_param.lb;

		model->model_param.psum_p = 0.0;
		model->model_param.psum_n = 0.0;
		model->model_param.vsum_p = 0.0;
		model->model_param.vsum_n = 0.0;
		model->model_param.covsum = 0.0;
		model->model_param.entropy[0] = 0.0;
		model->model_param.entropy[1] = 0.0;
		model->model_param.entropy[2] = 0.0;

		//if(model->control_param.test_null==0){
			update_p(model);
		//}
		update_theta_gamma(model);
		//if(count==1&&model->control_param.test_null!=1){
		//	model->model_param.theta[0] = rnorm(0,1);
			//Rprintf("theta: %g\n",model->model_param.theta[0]);
		//}
		update_sigma(model);
		update_lb(model);
		//Rprintf("LOWER BOUND: %g\n",model->model_param.lb);
		tol = lb_old - model->model_param.lb;
		count = count+1;
	}
	//Rprintf("LOWER BOUND: %g\n",model->model_param.lb);
}

void run_pathmix_wrapper(double * eps,
			int * maxit,
			int * regress,
			int * scale,
			int * sided,
			int * test_null,
			double * G,
			double * X,
			double * Xhat,
			double * y,
			double * var_y,
			int * n,
			int * m,
			int * p,
			double * pvec_res,
			double * gamma_res,
			double * theta_res,
			double * sigma_res,
			double * prob_res,
			double * lb_res){

	struct model_struct model;
	//Rprintf("Initializing model...\n");
	initialize_model(eps,maxit,regress,scale,sided,test_null,G,X,Xhat,y,var_y,n,m,p,&model);
	//Rprintf("Model initialized, running model...\n");
	run_pathmix(&model);
	//Rprintf("Model run, collapsing results...\n");
	collapse_results(&model,pvec_res,gamma_res,theta_res,sigma_res,prob_res,lb_res);
	//Rprintf("Results collapsed, freeing memory...\n");
	free_model(&model);
	//Rprintf("Memory freed\n");
}


