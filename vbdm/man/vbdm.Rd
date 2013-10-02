\name{vbdm}
\alias{vbdm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fit a discrete mixture model}
\description{Fits a discrete mixture model for rare variant association analysis.  Uses an approximate variational Bayes coordinate ascent algorithm for a computationally efficient solution.}
\usage{
vbdm(y, G, X=NULL, thres=0.05, genotypes=TRUE,
     include.mean=TRUE, minor.allele=TRUE, impute="MEAN",
     eps=1e-4, scaling=TRUE, nperm=0, maxit=1000, hyper=c(2,2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of continuous phenotypes.}
  \item{G}{A matrix of genotypes or variables of interest.  }
  \item{X}{An optional matrix of covariates.}
  \item{thres}{If the matrix is of genotypes, then this specifies a minor allele frequency threshold.  Variants with a MAF greater than this threshold are excluded from the analysis.}
  \item{genotypes}{This specifies whether or not to treat G as a matrix of genotypes.  If it is treated as genotypes then it will be filtered based on \code{thres}, and there are more options for missing data imputations.  The default genotype encoding is additive (e.g. genotypes are encoded as 0,1,2).  Also if \code{G} is a genotype matrix vbdm will flip the encoding such that the homozygous major allele genotype is encoded as 0, the heterozygote as 1, and the homozygous minor allele genotype as 2 unless \code{minor.allele=FALSE}}
  \item{include.mean}{This specifies whether to add an interecept term to the model.  If no covariates are provided it is automatically added, but if there are covariates provided it can be optional.}
  \item{minor.allele}{When \code{minor.allele=TRUE} and \code{genotypes=TRUE} the genotypes are flipped so that the major allele genotype is encoded as 0.}
  \item{impute}{If there is missing data in \code{G} this specifies the method with which to impute the missing data. There are two options \code{impute="MEAN"} which sets any missing genotype to the expected dosage given the MAF, or \code{impute="MAJOR"} which sets any missing genoypte to the homozygous genotype of the major allele.  If the matrix is not treated as a genotype matrix (e.g. \code{genotype=FALSE}), then only \code{impute="MEAN"} will work.  Also, missing data is not allowed in the covariates X.}
  \item{eps}{The tolerance for convergence of the coordinate ascent algorithm based on the change in the lower bound of the log marginal likeilhood.}
  \item{scaling}{Whether or not to scale the genotypes to have mean 0 and variance 1.}
  \item{nperm}{Optional parameter defining the number of null permutations of the vbdm likelihood ratio test.  This can be used to generate an exact p-value}
  \item{maxit}{The maximum number of iterations allowed for the vbdm algorithm.}
  \item{hyper}{The hyperparameters for the prior defined over the mixing probability parameter.  The first hyperparameter is the alpha parameter, and the second is the beta parameter.}
  
  
}

\value{

  \item{y}{The phenotype vector passed to vbdm.}
  \item{G}{The genotype matrix passed to vbdm.  Note that any variables that were dropped will be dropped from this matrix.}
  \item{X}{The covariate matrix passed to vbdm.  Will include intercept term if it was added earlier.}
  \item{keep}{A vector of indices of the kept variables in G (if any were excluded based on \code{thres})}
  \item{pvec}{The vector of estimated posterior probabilities for each variable in G.}
  \item{gamma}{A vector of additive covariate effect estimates.}
  \item{theta}{The estimated effect of the variables in G.}
  \item{sigma}{The estimated error variance.}
  \item{prob}{The estimated mixing parameter.}
  \item{lb}{The lower bound of the marginal log likelihood.}
  \item{lbnull}{The lower bound of the marginal log likelihood under the null model.}
  \item{lrt}{The approximate likelihood ratio test based on the lower bounds.}
  \item{p.value}{A p-value computed based on \code{lrt} with the assumption that \code{lrt~chi^2_1}}
  \item{lbperm}{If \code{nperm>0}, the lower bound of the fitted null permutations.}
  \item{lrtperm}{If \code{nperm>0}, the likelihood ratio test of the fitted null permutations.}
  \item{p.value.perm}{If \code{nperm>0}, the empirical p-value based on the fitted null permutations.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Benjamin A. Logsdon (blogsdon@uw.edu)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vbdmR},\link{burdenPlot}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#generate some test data
library(vbdm)
set.seed(3)
n <- 1000
m <- 20
G <- matrix(rbinom(n*m,2,.01),n,m);
beta1 <- rbinom(m,1,.2)
y <- G\%*\%beta1+rnorm(n,0,1.3)

#with scaling:
res <- vbdm(y=y,G=G);
T5 <- summary(lm(y~rowSums(scale(G))))$coef[2,4];
cat('vbdm p-value:',res$p.value,'\nT5 p-value:',T5,'\n')
#vbdm p-value: 0.001345869 
#T5 p-value: 0.9481797 

#without scaling:
res <- vbdm(y=y,G=G,scaling=FALSE)
T5 <- summary(lm(y~rowSums(G)))$coef[2,4];
cat('vbdm p-value:',res$p.value,'\nT5 p-value:',T5,'\n')
#vbdm p-value: 0.0005315836 
#T5 p-value: 0.904476 

#run 10000 permutations
set.seed(2)
res <- vbdm(y=y,G=G,scaling=FALSE,nperm=1e4);
cat('vbdm approximate p-value:',res$p.value,'\nvbdm permutation p-value:',res$p.value.perm,'\n');
#vbdm approximate p-value: 0.0005315836 
#vbdm permutation p-value: 8e-04 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{vbdm}
\keyword{association}
\keyword{genetic}
\keyword{rare}
\keyword{variational}% __ONLY ONE__ keyword per line
