\name{StocSum.LDSC.ldscWinFastParrel}
\alias{StocSum.LDSC.ldscWinFastParrel}
\title{
Calculate LDscore
}
\description{
Use summary statistics from StocSum.stat to estimate LDscore.
}
\usage{
StocSum.LDSC.ldscWinFastParrel(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), N.randomvec=1000, MAF.range = c(1e-3, 0.5), MAF.weights.beta = c(0.5, 0.5), miss.cutoff = 1, wind.b=1000000, use.minor.allele = FALSE, auto.flip = FALSE, nperbatch = 100000, ncores = 1)
}
\arguments{
  \item{meta.file.prefix}{
a character vector for prefix of intermediate files (*.sample.* and *.resample.*). 
}
  \item{n.files}{
an integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix)))
}
 \item{N.randomvec}{
the number of replicates to simulate the random vectors in \code{StocSum.LDSC.glmmkin2randomvec} (default = 1000).
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-3, 0.5)).
}
  \item{MAF.weight.beta}{
MAF.weight.beta = c(0.5, 0.5). %%a numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies. This internal minor allele frequency weight is multiplied by the external weight given by the group.file. To turn off internal minor allele frequency weight and only use the external weight given by the group.file, use c(1, 1) to assign flat weights (default = c(1, 25)). Applied to the combined samples.
}
  \item{missing.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
}
 \item{wind.b}{
a positive integer define the window size in bases (default = 1000000).
}
  \item{use.minor.allele}{
a logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SKAT results, but Burden (as well as SKAT-O and hybrid test to combine the burden test and SKAT) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction. Use with caution, as major/minor alleles may flip in different cohorts. In that case, minor allele will be determined based on the allele frequency in the combined samples.
}
  \item{auto.flip}{
a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
}
\item{nperbatch}{
an integer for how many SNPs to be included in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
}
  \item{ncores}{
a positive integer indicating the number of cores to be used in parallel computing
(default = 1).
}
}
\value{
\code{StocSum.LDSC.ldscWinFastParrel} returns a data frame with the following components:
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
  \qitem{N}{
total sample size.
}
  \item{altfreq}{
alternative allele frequency.
}
  \item{idx}{
indicating SNPs
}
 \item{pos.end}{
To which SNP (position) is included to esitmate LD score for the current SNP.
}
 \item{pos.start}{
From which SNP (position) is included to estimate LD score for the current SNP. 
}
  \item{LDscore}{
LD score for each SNP.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Han Chen, Nannan Wang
}
\seealso{
\code{StocSum.LDSC.glmmkin2randomvec}, \code{StocSum.LDSC}
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

# run StocSum.pval
out<-StocSum.pval(out.prefix, group.file = group.file, MAF.range=c(0,0.5), miss.cutoff = 1)
}
}
\keyword{ LD score }
