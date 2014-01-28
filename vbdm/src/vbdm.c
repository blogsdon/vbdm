//variational Bayes rare variant mixture test (vbdm) R package C code
//Copyright 2013 Benjamin A Logsdon
#include "vbdm.h"

//get a column of genotype matrix G
double * gc(struct model_struct * model, int j){
	return (&(model->data.G[j]))->col;
}


//get a column of covariate matrix X
double * xc(struct model_struct * model, int j){
	return (&(model->data.X[j]))->col;
}

//get a column of the Xhat matrix
double * hc(struct model_struct * model, int j){
	return (&(model->data.Xhat[j]))->col;
}

double * xcf(struct model_struct * model, int j){
  return (&(model->data.X_fixed[j]))->col;
}

//get a column of the Xhat matrix
double * hcf(struct model_struct * model, int j){
	return (&(model->data.Xhat_fixed[j]))->col;
}


void ddot_w(int n,double *vect1,double *vect2,double * result){
	const int incxy = 1;
	(*result)=F77_NAME(ddot)(&n,vect1,&incxy,vect2,&incxy);
}


void daxpy_w(int n,double *x,double *y,double alpha){
	//y<- ax+y;
	const int incxy =1;
	F77_NAME(daxpy)(&n,&alpha,x,&incxy,y,&incxy);
}

void dnrm2_w(int n,double *x,double *result){
	const int incxy=1;
	(*result)=F77_NAME(dnrm2)(&n,x,&incxy);
}

void dscal_w(int n,double *x, double alpha){
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

double compute_ssq(double *vec,int n){
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
			double * G, 
			double * X,
			double * Xhat,
			double * y, 
			double * var_y,
      double * hyper,
			int * n, 
			int * m, 
			int * p,
      int * nperm,
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

	model->control_param.test_null = 0;
  model->control_param.nperm = (*nperm);
  //Rprintf("allocating G\n");
	model->data.G = (struct matrix_v *) malloc(sizeof(struct matrix_v)*(*m));
	for(k=0;k<(*m);k++){
		(&(model->data.G[k]))->col = (double *) malloc(sizeof(double)*(*n));
	}
  //Rprintf("allocating G 1 %d %d\n",*m,*n);
	for(k=0;k<(*m);k++){
		for(l=0;l<(*n);l++){
			(&(model->data.G[k]))->col[l] = G[k*(*n)+l];
		}
	}
  //Rprintf("allocating G 2\n");

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
  
  model->data.X_fixed= (struct matrix_v *) malloc(sizeof(struct matrix_v)*(*p));
	for(k=0;k<(*p);k++){
		(&(model->data.X_fixed[k]))->col = (double *) malloc(sizeof(double)*(*n));
	}

	for(k=0;k<(*p);k++){
		for(l=0;l<(*n);l++){
			(&(model->data.X_fixed[k]))->col[l] = X[k*(*n)+l];
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
  
  model->data.Xhat_fixed= (struct matrix_v *) malloc(sizeof(struct matrix_v)*(*p));
	for(k=0;k<(*p);k++){
		(&(model->data.Xhat_fixed[k]))->col = (double *) malloc(sizeof(double)*(*n));
	}

	for(k=0;k<(*p);k++){
		for(l=0;l<(*n);l++){
			(&(model->data.Xhat_fixed[k]))->col[l] = Xhat[k*(*n)+l];
		}
	}
  
	//Rprintf("here3\n");
	model->data.y = y;
  model->data.y_fixed = (double *) malloc(sizeof(double)*(*n));
  for(k=0;k<(*n);k++){
    model->data.y_fixed[k] = y[k];
  }
  
  model->data.perm = (int *) malloc(sizeof(int)*(*n));
  for(k=0;k<(*n);k++){
    model->data.perm[k] = 0;
  }
  
  model->data.ans = (int *) malloc(sizeof(int)*(*n));
  for(k=0;k<(*n);k++){
    model->data.ans[k] = 0;
  }
  
  double nd = (double) (*n);
  model->data.probv = (double *) malloc(sizeof(double)*(*n));
  for(k=0;k<(*n);k++){
    model->data.probv[k] = 1.0/nd;
  }

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
  
  model->model_param.hyper = hyper;

	model->model_param.pvec = (double *) malloc(sizeof(double)*(*m));

	for(k=0;k<(*m);k++){	
	  model->model_param.pvec[k] = 0.5;
	}


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
  
  if((*nperm)>0){
    model->model_param.lb_perm = (double *) malloc(sizeof(double)*(*nperm));
    model->model_param.lb_perm_null = (double *) malloc(sizeof(double)*(*nperm));
    for(k=0;k<(*nperm);k++){
      model->model_param.lb_perm[k] = 0.0;
      model->model_param.lb_perm_null[k] = 0.0;
    }
  }

	model->model_param.prob = (double *) malloc(sizeof(double)*1);
	for(k=0;k<1;k++){
    model->model_param.prob[k]= 0.5;
  }

  //Rprintf("here9\n");
	model->model_param.resid_vec = (double *) malloc(sizeof(double)*(*n));
	model->model_param.Gp = (double *) malloc(sizeof(double)*(*n));
	for(k=0;k<(*n);k++){
		model->model_param.resid_vec[k] = y[k];
		model->model_param.Gp[k] = 0;
	}
	//Rprintf("here10, x1: %g, %g, %g\n",gc(model,0)[0],gc(model,0)[1],gc(model,0)[2]);
	//Rprintf("here10, pvec[0]: %g\n",model->model_param.pvec[0]);
	for(k=0;k<(*m);k++){
		daxpy_w((*n),gc(model,k),model->model_param.Gp,model->model_param.pvec[k]);
	}
	//Rprintf("here11\n");	
	model->model_param.lb = -1e100;
	model->model_param.psum = 0.0;
	model->model_param.vsum = 0.0;

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
  
  for(k=0;k<(model->data.p);k++){
		free((&(model->data.X_fixed[k]))->col);
		
	}
	free(model->data.X_fixed);

	for(k=0;k<(model->data.p);k++){
		free((&(model->data.Xhat_fixed[k]))->col);
		
	}
	free(model->data.Xhat_fixed);  
  
  
  if(model->control_param.nperm>0){
    free(model->model_param.lb_perm);
    free(model->model_param.lb_perm_null);
  }


	free(model->data.g_sum_sq);
	free(model->data.one_vec);
  free(model->data.perm);
  free(model->data.probv);
  free(model->data.y_fixed);
  free(model->data.ans);

	free(model->model_param.pvec);
  free(model->model_param.entropy);
	free(model->model_param.gamma);
	free(model->model_param.theta);
	free(model->model_param.prob);
	free(model->model_param.resid_vec);
	free(model->model_param.Gp);

}

void reset_response(struct model_struct * model){
  int k;
  for(k=0;k<model->data.n;k++){
    model->data.y[k] = model->data.y_fixed[k];
  }
  
}

void generatePermutation(struct model_struct * model){
  int k;
  for (k=0;k<model->data.n;k++){
    model->data.ans[k] = k;
    model->data.probv[k] = unif_rand();
    if(k==0){
    //Rprintf("random sample: %g\n",model->data.probv[k]);
    }
  }
  rsort_with_index(model->data.probv,model->data.ans,model->data.n);
  //Rprintf("permutation: %d %d %d\n",model->data.ans[0],model->data.ans[1],model->data.ans[2]);
  
}

void permutey(struct model_struct * model){
  int k;
  for(k=0;k<model->data.n;k++){
    model->data.y[k] = model->data.y_fixed[model->data.ans[k]];
  }
}

void reset_xhat(struct model_struct * model){
  int j,k;
  for(j=0;j<model->data.p;j++){
    for(k=0;k<model->data.n;k++){
      (hc(model,j))[k] = (hcf(model,j))[k];
      //permute Xhat
    }
  }  
}

void permutexhat(struct model_struct * model){
  int j,k;
  for(j=0;j<model->data.p;j++){
    for(k=0;k<model->data.n;k++){
      (hc(model,j))[k] = (hcf(model,j))[model->data.ans[k]];
      //permute Xhat
    }
  }
}

void reset_x(struct model_struct * model){
  int j,k;
  for(j=0;j<model->data.p;j++){
    for(k=0;k<model->data.n;k++){
      (xc(model,j))[k] = (xcf(model,j))[k];
      //permute Xhat
    }
  } 
}

void permutex(struct model_struct * model){
  int j,k;
  for(j=0;j<model->data.p;j++){
    for(k=0;k<model->data.n;k++){
      (xc(model,j))[k] = (xcf(model,j))[model->data.ans[k]];
      //permute Xhat
    }
  }
}

void reset_model(struct model_struct * model){
  int k;
	for(k=0;k<(model->data.m);k++){	
	  model->model_param.pvec[k] = 0.5;
	}
	for(k=0;k<(model->data.p);k++){
		model->model_param.gamma[k] = 0.0;
	}

	for(k=0;k<2;k++){
		model->model_param.theta[k]= 0.0;
	}
	//Rprintf("here8\n");
  for(k=0;k<3;k++){
		model->model_param.entropy[k]=0.0;
	}

	model->model_param.sigma = model->data.var_y;

	for(k=0;k<1;k++){
    model->model_param.prob[k]= 0.5;
  }

  //Rprintf("here9\n");
	for(k=0;k<model->data.n;k++){
		model->model_param.resid_vec[k] = model->data.y[k];
		model->model_param.Gp[k] = 0;
	}
	//Rprintf("here10, x1: %g, %g, %g\n",gc(model,0)[0],gc(model,0)[1],gc(model,0)[2]);
	//Rprintf("here10, pvec[0]: %g\n",model->model_param.pvec[0]);
	for(k=0;k<(model->data.m);k++){
		daxpy_w((model->data.n),gc(model,k),model->model_param.Gp,model->model_param.pvec[k]);
	}
	//Rprintf("here11\n");	
	model->model_param.lb = -1e100;
	model->model_param.psum = 0.0;
	model->model_param.vsum = 0.0;
  
}


void update_p(struct model_struct * model){
	int k;
	double pold,vec1,a1,a2,a3,a4,a5,pnew;
	double halpha,hbeta,albe;
  
	double md = (double) model->data.m;
  halpha = model->model_param.hyper[0];
  hbeta = model->model_param.hyper[1];  
  albe = halpha+hbeta;
  
	for(k=0;k<model->data.m;k++){			
	  pold = model->model_param.pvec[k];
		ddot_w(model->data.n,model->model_param.resid_vec,gc(model,k),&vec1);
		vec1 = vec1 + model->data.g_sum_sq[k]*pold*model->model_param.theta[0];
		a1 = pow(model->model_param.theta[0],2)*model->data.g_sum_sq[k];
		a2 = -(2*model->model_param.theta[0]*vec1);
		a3 = -(digamma((model->model_param.prob[0])*(md+albe))-digamma(md+albe));
		a4 = digamma((1-model->model_param.prob[0])*(md+albe))-digamma(md+albe);
		a5 = (1/(2*model->model_param.sigma))*(a1+a2) + a3 + a4;
		pnew = 1/(1+exp(a5));
		model->model_param.pvec[k] = pnew;
		daxpy_w(model->data.n,gc(model,k),model->model_param.resid_vec,model->model_param.theta[0]*(pold-pnew));
		daxpy_w(model->data.n,gc(model,k),model->model_param.Gp,pnew-pold);
		model->model_param.psum = model->model_param.psum + pnew;
		model->model_param.vsum = model->model_param.vsum + model->data.g_sum_sq[k]*(pnew-pow(pnew,2));
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
  model->model_param.prob[0] = (model->model_param.psum+halpha)/(md+albe);
}

void update_theta_gamma(struct model_struct * model){


	double theta_old, theta_new, const1,const2;
	int p = model->data.p;
	int k;
	double gamma_old[p];
	double gamma_new;
	for (k=0;k<p;k++){
		gamma_old[k] = model->model_param.gamma[k];
	}
  //Rprintf("test_null: %d\n",model->control_param.test_null);
	//update theta
	if(model->control_param.test_null==1){
	  theta_new = 0.0;
	}else{
	  theta_old = model->model_param.theta[0];
		ddot_w(model->data.n,model->model_param.resid_vec,model->model_param.Gp,&theta_new);
		ddot_w(model->data.n,model->model_param.Gp,model->model_param.Gp,&const1);

		theta_new = theta_new + const1*theta_old;
		const2 = const1+model->model_param.vsum;
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

}

void update_sigma(struct model_struct * model){
	double sigma;
	double nd = (double) model->data.n;
	ddot_w(model->data.n,model->model_param.resid_vec,model->model_param.resid_vec,&sigma);
	sigma = sigma+model->model_param.vsum*pow(model->model_param.theta[0],2);
	sigma = sigma/nd;
	model->model_param.sigma = sigma;
}


void update_lb(struct model_struct * model){

	double lb,alpha1,beta1,halpha,hbeta,albe;
	double nd = (double) model->data.n;
	double md = (double) model->data.m;
  halpha = model->model_param.hyper[0];
  hbeta = model->model_param.hyper[1];
  albe = halpha+hbeta;
  
	lb = -0.5*(nd*(log(2*M_PI*model->model_param.sigma)+1));
	//Rprintf("lb ll: %g\n",lb);
	lb = lb + (digamma((model->model_param.prob[0])*(md+albe))-digamma(md+albe))*(model->model_param.psum+1);
	lb = lb + (digamma((1-model->model_param.prob[0])*(md+albe))-digamma(md+albe))*(md-model->model_param.psum+1);
	//Rprintf("lb elp: %g\n",lb);
	lb = lb + model->model_param.entropy[0];
	lb = lb + model->model_param.entropy[1];
	//Rprintf("lb entropy: %g\n",lb);
	alpha1 = model->model_param.psum +halpha;
	beta1 = md - model->model_param.psum+hbeta;
			
	lb = lb + lbeta(alpha1,beta1) - (alpha1-1)*digamma(alpha1)-(beta1-1)*digamma(beta1)+(alpha1+beta1-2)*digamma(md+2);
	//Rprintf("lb entropy beta: %g\n",lb);
	model->model_param.lb = lb;
}


void collapse_results(struct model_struct * model,
		double * pvec_res,
		double * gamma_res,
		double * theta_res,
		double * sigma_res,
		double * prob_res,
		double * lb_res,
    double * lb_null_res){

	int k;
	//int n = model->data.n;
	int m = model->data.m;
	int p = model->data.p;
	for(k=0;k<m;k++){
		pvec_res[k] = model->model_param.pvec[k];
	}
	for(k=0;k<p;k++){
		gamma_res[k] = model->model_param.gamma[k];
	}
	theta_res[0] = model->model_param.theta[0];
	sigma_res[0] = model->model_param.sigma;
	prob_res[0] = model->model_param.prob[0];
	lb_res[0] = model->model_param.lb;
  lb_null_res[0] = model->model_param.lb_null;
  if(model->control_param.nperm>0){
    for(k=0;k<model->control_param.nperm;k++){
      lb_res[k+1] = model->model_param.lb_perm[k];
      lb_null_res[k+1] = model->model_param.lb_perm_null[k];
      
    }
  }
}

void run_vbdm(struct model_struct * model){
	double tol=1;
	double lb_old;
	int count = 0;
  int i;
	//Rprintf("tol: %g, eps: %g\n",fabs(tol),model->control_param.eps);
  //run null model.
  //grab lower bound.
  model->control_param.test_null = 1;
  while(fabs(tol)>model->control_param.eps && count < model->control_param.maxit){
    //Rprintf("fitting null model %d\n",count);
		lb_old = model->model_param.lb;
		model->model_param.psum = 0.0;
		model->model_param.vsum = 0.0;
		model->model_param.entropy[0] = 0.0;
		model->model_param.entropy[1] = 0.0;
		model->model_param.entropy[2] = 0.0;
		update_p(model);
		update_theta_gamma(model);
		update_sigma(model);
		update_lb(model);
		tol = lb_old - model->model_param.lb;
		count = count+1;
	}
  //Rprintf("null model lb: %g, iter: %d, theta: %g\n",model->model_param.lb,count,model->model_param.theta[0]);
  
  count = 0;
  tol = 1;
  model->model_param.lb_null = model->model_param.lb;
  reset_model(model);
  model->control_param.test_null = 0;

  //un null model
  //if permutations, run permutations
  //int i;
  GetRNGstate();
  if(model->control_param.nperm>0){
    for(i=0;i<model->control_param.nperm;i++){
      generatePermutation(model);
      permutey(model);
      permutex(model);
      permutexhat(model);
      reset_model(model);
      //Rprintf("y[0]: %g, y[1]: %g\n",model->data.y[0],model->data.y[1]);
      model->control_param.test_null = 1;
      while(fabs(tol)>model->control_param.eps && count < model->control_param.maxit){
        //Rprintf("fitting null model %d\n",count);
  	    lb_old = model->model_param.lb;
		    model->model_param.psum = 0.0;
		    model->model_param.vsum = 0.0;
		    model->model_param.entropy[0] = 0.0;
		    model->model_param.entropy[1] = 0.0;
		    model->model_param.entropy[2] = 0.0;
		    update_p(model);
		    update_theta_gamma(model);
		    update_sigma(model);
		    update_lb(model);
		    tol = lb_old - model->model_param.lb;
		    count = count+1;
	    }
      count=0;
      tol=1;
      model->model_param.lb_perm_null[i] = model->model_param.lb;
      reset_model(model);
      model->control_param.test_null= 0;
      
      while(fabs(tol)>model->control_param.eps && count < model->control_param.maxit){
  		  lb_old = model->model_param.lb;
		    model->model_param.psum = 0.0;
  		  model->model_param.vsum = 0.0;
	  	  model->model_param.entropy[0] = 0.0;
  		  model->model_param.entropy[1] = 0.0;
		    model->model_param.entropy[2] = 0.0;
		    update_p(model);
		    update_theta_gamma(model);
		    update_sigma(model);
		    update_lb(model);
        tol = lb_old - model->model_param.lb;
		    count = count+1;
	    }
      //Rprintf("perm alternative lb: %g, iter: %d, theta: %g\n",model->model_param.lb,count,model->model_param.theta[0]);

      model->model_param.lb_perm[i]= model->model_param.lb;
      
      tol=1;
      count=0;
    }
  }
  PutRNGstate();
  reset_response(model);
  reset_x(model);
  reset_xhat(model);
  reset_model(model);

  tol=1;
  count=0;
  model->control_param.test_null = 0;
  //reset parameters
  //run alternative model
	while(fabs(tol)>model->control_param.eps && count < model->control_param.maxit){
    //Rprintf("fitting alternative model %d\n",count);
		lb_old = model->model_param.lb;
		model->model_param.psum = 0.0;
		model->model_param.vsum = 0.0;
		model->model_param.entropy[0] = 0.0;
		model->model_param.entropy[1] = 0.0;
		model->model_param.entropy[2] = 0.0;
		update_p(model);
		update_theta_gamma(model);
		update_sigma(model);
		update_lb(model);
		tol = lb_old - model->model_param.lb;
		count = count+1;
	}
  //Rprintf("true alternative lb: %g, iter: %d, theta: %g\n",model->model_param.lb,count,model->model_param.theta[0]);
}

void run_vbdm_wrapper(double * eps,
			int * maxit,
			int * regress,
			int * scale,
			double * G,
			double * X,
			double * Xhat,
			double * y,
			double * var_y,
      double * hyper,
			int * n,
			int * m,
			int * p,
      int * nperm,
			double * pvec_res,
			double * gamma_res,
			double * theta_res,
			double * sigma_res,
			double * prob_res,
			double * lb_res,
      double * lb_null_res){

	struct model_struct model;
	//Rprintf("Initializing model...\n");
	initialize_model(eps,maxit,regress,scale,G,X,Xhat,y,var_y,hyper,n,m,p,nperm,&model);
	//Rprintf("Model initialized, running model...\n");
	run_vbdm(&model);
	//Rprintf("Model run, collapsing results...\n");
	collapse_results(&model,pvec_res,gamma_res,theta_res,sigma_res,prob_res,lb_res,lb_null_res);
	//Rprintf("Results collapsed, freeing memory...\n");
	free_model(&model);
	//Rprintf("Memory freed\n");
}


