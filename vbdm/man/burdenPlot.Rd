\name{burdenPlot}
\alias{burdenPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plotting function for rare variant vbdm test}
\description{This function produces a plot depicting the phenotype distribution, the distribution of carriers of the rare alleles, and also can include the posterior probaiblity of association as estimated by the vbdm algorithm.}
\usage{
burdenPlot(y, G, annotation = rep('missense',ncol(G)), title='',
          order='mean', legend='keep', type='lines',
          post=NULL, name.snp=NULL)

}
\arguments{
  \item{y}{A vector of continuous phenotypes.}
  \item{G}{A matrix of genotypes that are encoded as 0, 1, or 2.}
  \item{annotation}{A vector of SNP annotations for the columns of G.  The default is \code{'missense'}.  Possible annotations are \code{'missense'}, \code{'nonsense'},\code{'splice'}, and \code{'synonymous'}.  Default is for all variants to be missense.}
  \item{title}{An optional title for the plot.}
  \item{order}{How the variants should be ordered in the bottom quadrant of the plot.  Possible options include \code{order='mean'} which orders based on the mean phenotypic value of carriers of the rare variants, \code{order='MAF'} which orders the variants based on minor allele frequency, \code{order='MAF.mean'} which orders variants first by MAF, then by mean phenotypic value, \code{order='anno'} which orders by annotation first, then by mean phenotypic value, and \code{order=''} which removes any reordering.}
  \item{legend}{If \code{legend='keep'} then a legend is depicted with variant annotations.}
  \item{type}{If \code{type='lines'} then the range of phenotype values for carriers are shown with a horizontal lines, with phenotypic values at the vertical tick marks.  If \code{type='points'} then the phenotypic values are shown with points.}
  \item{post}{Optional vector of posterior probabilities from \code{\link{vbdm}} result.}
  \item{name.snp}{Optional vector of snp names.}

}
\references{Logsdon, B.A., et al. (2014)
  \emph{A Variational Bayes Discrete Mixture Test for
  Rare Variant Association.},
  \emph{Genetic Epidemiology, Vol. 38(1), 21-30 2014}
}
\author{
Paul L. Auer (paul.wl.auer@gmail.com), Benjamin A. Logsdon (blogsdon@uw.edu)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vbdm}}, \code{\link{vbdmR}}
}
\examples{
#generate some test data
library(vbdm)
set.seed(1)
n <- 1000
m <- 30
G <- matrix(rbinom(n*m,2,.01),n,m);
beta1 <- rbinom(m,1,.2)
y <- G\%*\%beta1+rnorm(n,0,2)
res <- vbdm(y=y,G=G,scaling=FALSE);
bp<-burdenPlot(y=y,G=G,post=res$pvec,name.snp=1:30);




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{vbdm}
\keyword{association}
\keyword{genetic}
\keyword{rare}
\keyword{variational}% __ONLY ONE__ keyword per line
\keyword{plot}
