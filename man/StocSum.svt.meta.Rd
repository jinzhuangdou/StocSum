\name{StocSum.svt.meta}
\alias{StocSum.svt.meta}
\title{
Performing meta-analysis for GLMM bashed socre test results
}
\description{
Use output files from \code{StocSum.stat} to peform meta-analysis.
}
\usage{
StocSum.svt.meta(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), outfile.prefix, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, auto.flip = FALSE, nperbatch = 10000)
}
\arguments{
  \item{meta.file.prefix}{
a character vector for prefix of intermediate files (*.sample.* and *.resample.*) which stores summary statistics calculated from \code{StocSum.stat}
}
  \item{n.files}{
an integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix)))
}
  \item{outfile.prefix}{
a character vector for prefix of intermediate files (*.sample.* and *.resample.*) which stores summary statistics of the meta-analysis
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
}
  \item{missing.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
}
  \item{auto.flip}{
a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
}
  \item{nperbatch}{
an integer for how many SNPs should be tested in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
}
}
\value{
NULL. 
\code{StocSum.svt.meta} will store the summary statistics and the statistics which to be used to calculate LD statistics in two files with the prefix specified by the user.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Han Chen, Nannan Wang
}
\seealso{
\code{glmmkin2randomvec}, \code{StocSum.pval},\code{StocSum.stat}, \code{StocSum.svt}
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
StocSum.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1

# run StocSum.svt
out<-StocSum.svt(out.prefix, n.file=ncores, MAF.range=c(0,0.5), miss.cutoff = 1, auto.flip=F)
}
}
\keyword{ meta-analysis }
\keyword{ score test }
