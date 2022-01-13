\name{StocSum.svt}
\alias{StocSum.svtt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
StocSum.stat
}
\description{
Single variant tests using the compact and informatics statics.
}
\usage{
StocSum.svt(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), 
MAF.range = c(1e-7, 0.5), miss.cutoff = 1, auto.flip = FALSE, nperbatch = 10000)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{meta.file.prefix---> name is not good, meta is not necessary}{
prefix of intermediate files (*.sample.1 and *.resample.1) required in the single variant test. 
}
  \item{n.files--->check}{
an integer vector with the same length as meta.files.prefix, indicating how many
sets of intermediate files (.score.* and .var.*) are expected from each cohort,
usually as the result of multi-threading in creating the intermediate files (default
= rep(1, length(meta.files.prefix)))
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
}
  \item{missing.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
}
  \item{auto.flip}{
a logical switch for whether to enable automatic allele flipping if a variant with
alleles ref/alt is not found at a position, but a variant at the same position with
alleles alt/ref is found (default = FALSE). Use with caution for whole genome
sequence data, as both ref/alt and alt/ref variants at the same position are not
uncommon, and they are likely two different variants, rather than allele flipping.
}
  \item{nperbatch}{
an integer for how many SNPs should be tested in a batch (default = 10000).
The computational time can increase dramatically if this value is either small
or large. The optimal value for best performance depends on the user’s system.
}
}

\value{

\code{StocSum.stat} returns a data frame with the following components:
  \item{SNP}{
SNP name.
}
  \item{chr}{
chromosome name.
}
  \item{pos}{
the genome location of SNP.
}
  \item{ref}{
allele of reference.
}
  \item{alt}{
alternative allele.
}
  \item{N}{
total sample size.
}
 \item{missrate}{
missing rate of variants.
}
  \item{altfreq}{
alternative allele frequency.
}
  \item{SCORE}{
the summary score of the alternaive allele.
}
  \item{VAL}{
the variance of the summary score.
}
  \item{PVAL}{
the p-value of the suammry score.
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
\code{glmmkin}, \code{StocSum.pval},\code{StocSum.stat}, \code{StocSum.meta}
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
)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compact and informatics statics }
