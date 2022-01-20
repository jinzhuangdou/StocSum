%% Question? do I need to delete the unnecessary arguments in this function? Z=NULL, group.idx=NULL...
\name{StocSum.LDSC.glmmkin2randomvec}
\alias{StocSum.LDSC.glmmkin2randomvec}
\title{
Generate random vectors from multivariate normal distribution. 
}
\description{
Generate random vectors from multivariate normal distribution with mean 0 and covariance matrix P, in which \eqn P = I-1inv((t(1)1))t(1). 
}
\usage{
glmmkin2randomvec(obj, Z = NULL, N.randomvec = 1000, 
  group.idx=NULL, cluster.idx=NULL, robust = FALSE)
}
\arguments{
  \item{obj}{
a class glmmkin object, returned by fitting the null GLMM using \code{glmmkin}.
}
  \item{Z}{
Z = NULL. %%a list of matrices, the number of which should match the number of variance components in \code{obj}. The metrices is the upper triangular factor of the Choleski decomposition of the relationship matrices in kins (i.e., the arguments in \code{glmmkin}) (default = NULL).
}
  \item{N.randomvec}{
the number of replicates to simulate a length N random vector rb from a multivariate normal distribution with mean 0 and covariance matrix P. All the N.randomvec random vectors will be combined an N x N.randomvec matrix R. (default = 1000).
}
  \item{group.idx}{
group.idx = NULL. %%an optional categorical variable indicating the groups used in a heteroscedastic linear mixed model (allowing residual variances in different groups to be different) (defult = NULL).
}
  \item{cluster.idx}{
cluster.idx = NULL. 
}
  \item{robust}{
robust = FALSE.
}
}

\value{
\code{glmmkin2randomvec} returns a data list with the following components:
  \item{random.vectors}{
random vectors generated from multivariate normal distribution with mean 0 and covariance matrix P, which is from \code{glmmkin}.

}
  \item{theta}{
obj$theta set to be 1. %%a vector or a list of variance component parameter estimates. This option is inherited from \code{glmmkin}. 
}
  \item{scaled.residuals}{
a vector or a matrix for the scaled residuals, calculated as the original residuals divided by the dispersion parameter (in heteroscedastic linear mixed models, corresponding residual variance estimates by each group). This opition is inherited from \code{glmmkin}.
} 
  \item{id_include}{
a vector indicating the id of rows in data with nonmissing outcome and covariates, thus are included in the model fit. This opition is inherited from \code{glmmkin}.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Han Chen, Nannan Wang
}
\seealso{
\code{\link{glmmkin}}
}
\examples{
\donttest{
library(GMMAT)
data(example)
attach(example)
GRM.file <- system.file(“extdata”, “GRM.txt.bz2”, package = “StocSum”)
GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = “id”, family = binomial(link = “logit”))
obj <- glmmkin2randomvec(nullmod)}
}
\keyword{ random vector }
