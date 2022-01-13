\name{StocSum.stat}
\alias{StocSum.stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
StocSum.stat
}
\description{
Calculate the compact and informatics statics by leveraging the stochastic approximation.
}
\usage{
StocSum.stat(null.obj, geno.file, meta.file.prefix, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, missing.method = "impute2mean", nperbatch = 10000, ncores = 1)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{null.obj}{
a class glmmkin, returned by fitting the null GLMM using \code{glmmkin}.
}
  \item{geno.file}{
the .gds file name or an object of class SeqVarGDSClass for the full genotypes. The \code{sample.id} in \code{geno.file} should overlap \code{id_include} in \code{null.obj}. It is recommended that \code{sample.id} in \code{geno.file} include the full samples (at least all samples as specified in \code{id_include} of \code{null.obj}). It is not necessary for the user to take a subset of \code{geno.file} before running the analysis. If \code{geno.file} is an object of class SeqVarGDSClass, the .gds file will be closed upon successful completion of the function.
}
  \item{meta.file.prefix}{
prefix of intermediate files (*.sample.1 and *.resample.1) required in a meta-analysis.
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
}
  \item{missing.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
}
  \item{missing.method}{
method of handling missing genotypes. Either "impute2mean" or "impute2zero" (default = "impute2mean").
}
  \item{nperbatch}{
an integer for how many SNPs should be tested in a batch (default = 10000).
The computational time can increase dramatically if this value is either small
or large. The optimal value for best performance depends on the user’s system.
}
  \item{ncores}{
a positive integer indicating the number of cores to be used in parallel computing
(default = 1).
}


}

\value{
NULL. 
\code{StocSum.stat} will store the compact and informatics statics in two files with the prefix specified by the user.
}
  
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nannan Wang, Han Chen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{glmmkin}, \code{StocSum.pval}, \code{StocSum.svt}, \code{StocSum.meta}
}

\examples{
\donttest{
library(GMMAT)
data(example)
attach(example)

# run glmmkinrandomvec to generate the null model object
GRM.file <- system.file(“extdata”, “GRM.txt.bz2”, package = “StocSum”)
GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = “id”, family = binomial(link = “logit”))
obj <- glmmkin2randomvec(nullmod)}
out.prefix <- "test"
gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")

# run StocSum.stat 
StocSum.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compact and informatics statics }
