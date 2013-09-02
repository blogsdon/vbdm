//variational Bayes rare variant mixture test (vbdm) R package C library declarations
//Copyright 2013, Benjamin Logsdon


#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

typedef enum {LINEAR, LOGISTIC} regressionType;
typedef enum {STANDARD, NOSCALE} scalingType;


struct matrix_v {

	//x_col: a given vector of length n
	double * col;

};


struct control_param_struct {

	//convergence tolerance:
	double eps;

	//maximum number of iterations:
	double maxit;

	//type of regression
	regressionType regressType;

	//type of scaling
	scalingType scaleType;
  
  //number of permutations;
  int nperm;

	//whether or not to test reduced model
	int test_null;

};

struct model_param_struct {

	//posterior probabilities
	double * pvec;

	//fixed covariate effects
	double * gamma;

	//rare variant effect(s)
	double * theta;

	//error variance parameter
	double sigma;

	//mixing probability parameter(s)
	double * prob;

	//residual vector
	double * resid_vec;

	//positive reweighted vector
	double * Gp;

	//entropy of approximating distribution
	double * entropy;
  
	//sum of pvec;
	double psum;
	
	//positive variance correction
	double vsum;
  
  //hyper parameters for beta prior over probability parameter
  double * hyper;
  
  //lower bound of null model
  double lb_null;
  
  //lower bound of permutations
  double * lb_perm;

	//lower bound
	double lb;

};


struct data_struct { 

	//Genotype matrix
	struct matrix_v * G;

	//Covariate matrix
	struct matrix_v * X;

	//Covariate hat matrix
	struct matrix_v * Xhat;

	//response vector
	double * y;
  
  //fixed response vector
  double * y_fixed;
  
  //perm vector for resampling function
  int * perm;
  
  //ans vector for resampling function
  int * ans;

	//variance of response vector
	double var_y;

	//number of observations
	int n;

	//number of rare variants
	int m;

	//number of covariates
	int p;

	//sum of G squares
	double * g_sum_sq;

	//vector of 1s
	double * one_vec;
  
  //vector of 1/ns for sampling w/o replacement
  double * probv;

};


//declare full model structure
struct model_struct { 
	//control parameters
	struct control_param_struct control_param;

	//data set
	struct data_struct data;

	//model parameters
	struct model_param_struct model_param;
};

inline double * gc(struct model_struct * model, int j);

inline double * xc(struct model_struct * model, int j);

inline double * hc(struct model_struct * model, int j);

inline void ddot_w(int n,double *vect1,double *vect2,double * result);

inline void daxpy_w(int n,double *x,double *y,double alpha);

inline void dnrm2_w(int n,double *x,double *result);

inline void dscal_w(int n,double *x, double alpha);

void scale_vector(double * vec,double * ones,int n);

inline double compute_ssq(double *vec,int n);

void process_data(struct model_struct * model);

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
			struct model_struct * model);

void free_model(struct model_struct * model);

void update_p(struct model_struct * model);

void update_theta_gamma(struct model_struct * model);

void update_sigma(struct model_struct * model);

void update_lb(struct model_struct * model);

void run_vbdm(struct model_struct * model);

void reset_response(struct model_struct * model);

void permutey(struct model_struct * model);

void collapse_results(struct model_struct * model,
		double * pvec_res,
		double * gamma_res,
		double * theta_res,
		double * sigma_res,
		double * prob_res,
		double * lb_res,
    double * lb_null_res);
    
void reset_model(struct model_struct * model);

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
      double * lb_null_res);

