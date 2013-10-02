\name{vbdmR}
\alias{vbdmR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fit a discrete mixture model (R implementation)}
\description{Fits a discrete mixture model for rare variant association analysis.  Uses an approximate variational Bayes coordinate ascent algorithm for a computationally efficient solution.  This is the slow but well documented R implementation.}
\usage{
vbdmR(y, G, X=NULL, tol=1e-4, thres=0.05, scaling=TRUE, 
      hyper=c(2,2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of continuous phenotypes.}
  \item{G}{A matrix of genotypes or variables of interest.  }
  \item{X}{An optional matrix of covariates.}
  \item{tol}{The tolerance for convergence based on the change in the lower bound on the marginal log likelihood in the vbdm algorithm.}
  \item{thres}{If the matrix is of genotypes, then this specifies a minor allele frequency threshold.  Variants with a MAF greater than this threshold are excluded from the analysis.}
  \item{scaling}{Whether or not to scale the genotypes to have mean 0 and variance 1.}
  \item{hyper}{The hyperparameters for the prior defined over the mixing probability parameter.  The first hyperparameter is the alpha parameter, and the second is the beta parameter.}  
}
\details{
This function contains the much slower, but well documented R implementation of the vbdm algorithm.  This function does not have all of the sanity checks that \code{vbdm} has, and should therefore only be used for diagnostic purposes.
}
\value{
  \item{y}{The phenotype vector passed to vbdmR.}
  \item{G}{The genotype matrix passed to vbdmR.  Note that any variables that were dropped will be dropped from this matrix.}
  \item{X}{The covariate matrix passed to vbdmR.  Will include intercept term if it was added earlier.}
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
}

\author{
Benjamin A. Logsdon (blogsdon@uw.edu)
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vbdm}}, \code{\link{burdenPlot}}
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

#compare implementations
res1 <- vbdm(y=y,G=G);
res2 <- vbdmR(y=y,G=G);
T5 <- summary(lm(y~rowSums(scale(G))))$coef[2,4];
cat('vbdm p-value:',res1$p.value,
  '\nvbdmR p-value:',res2$p.value,
  '\nT5 p-value:',T5,'\n')
#vbdm p-value: 0.001345869 
#vbdmR p-value: 0.001345869 
#T5 p-value: 0.9481797 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{vbdm}
\keyword{association}
\keyword{genetic}
\keyword{rare}
\keyword{variational}% __ONLY ONE__ keyword per line

