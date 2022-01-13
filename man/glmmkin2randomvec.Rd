\name{glmmkin2randomvec}
\alias{glmmkin2randomvec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
glmmkin2randomvec
}
\description{
Generate the random vectors whose covariance approximate to the projection matrix of null model fitted from \code{glmmkin} function in \code{GMMAT} package.
}
\usage{
glmmkin2randomvec(obj, Z = NULL, N.randomvec = 1000, 
  group.idx=NULL, cluster.idx=NULL, robust = FALSE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{obj}{
a class glmmkin, returned by fitting the null GLMM using \code{glmmkin}.
}
  \item{Z}{
a list of matrices, the number of which should match the number of variance components in \code{obj}.
}
  \item{N.randomvec}{
the number of replicates (default = 1000).
}
  \item{group.idx}{
xxx.
}
  \item{cluster.idx}{
xxx.
}
  \item{robust}{
xxx. check the robust parameter. The code will not run if set robust=TRUE.
}

}

\value{
\code{glmmkin2randomvec} returns a data list with the following components:
  \item{random.vectors}{
random vectors whose covariance approximate to the projection matrix of null model fitted from \code{glmmkin}.

}
  \item{theta}{
a vector or a list of variance component parameter estimates. This option is inherited from \code{glmmkin}. 
}
  \item{scaled.residuals}{
a vector or a matrix for the scaled residuals, calculated as the original residuals divided by the dispersion parameter (in heteroscedastic linear mixed models, corresponding residual variance estimates by each group). This opition is inherited from \code{glmmkin}.
} 
  \item{id_include}{
a vector indicating the id of rows in data with nonmissing outcome and covariates, thus are included in the model fit.This opition is inherited from \code{glmmkin}.
}
}


\references{
%% ~put references to the literature/web site here ~
}
\author{
Nannan Wang, Han Chen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{glmmkin}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ random vector }
