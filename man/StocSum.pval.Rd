\name{StocSum.pval}
\alias{StocSum.pval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
StocSum.pval
}
\description{
variant set tests using the compact and informatics statics.
}
\usage{
StocSum.pval(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), cohort.group.idx = NULL, group.file, group.file.sep = "\t", MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25), miss.cutoff = 1, method = "davies", tests = "E", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE, auto.flip = FALSE)
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
 \item{cohort.group.idx--->check this one}{
a vector with the same length as meta.files.prefix, indicating which cohorts
should be grouped together in the meta-analysis assuming homogeneous genetic
effects. For example, c("a","b","a","a","b") means cohorts 1, 3, 4 are assumed
to have homogeneous genetic effects, and cohorts 2, 5 are in another group with
homogeneous genetic effects (but possibly heterogeneous with group "a"). If
NULL, all cohorts are in the same group (default = NULL).
}
 \item{group.file}{
a plain text file with 6 columns defining the test units. There should be no
headers in the file, and the columns are group name, chromosome, position,
reference allele, alternative allele and weight, respectively.
}
 \item{group.file.sep}{
the delimiter in group.file (default = "\t").
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
}
  \item{MAF.weight.beta}{
a numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies. This internal minor allele frequency weight is multiplied by the external weight given by the group.file. To turn off
internal minor allele frequency weight and only use the external weight given
by the group.file, use c(1, 1) to assign flat weights (default = c(1, 25)). Applied
to the combined samples
}
  \item{missing.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
}
  \item{method}{
a method to compute p-values for SKAT-type test statistics (default = "davies").
"davies" represents an exact method that computes a p-value by inverting the
characteristic function of the mixture chisq distribution, with an accuracy of
1e-6. When "davies" p-value is less than 1e-5, it defaults to method "kuonen".
"kuonen" represents a saddlepoint approximation method that computes the tail
probabilities of the mixture chisq distribution. When "kuonen" fails to compute
a p-value, it defaults to method "liu". "liu" is a moment-matching approximation
method for the mixture chisq distribution.
}
  \item{tests}{
a character vector indicating which SMMAT tests should be performed ("B" for
the burden test, "S" for SKAT, "O" for SKAT-O and "E" for the efficient hybrid
test of the burden test and SKAT). The burden test and SKAT are automatically
included when performing "O", and the burden test is automatically included
when performing "E" (default = "E").
}
  \item{rho}{
a numeric vector defining the search grid used in SMMAT-O for SKAT-O (see
the SKAT-O paper for details). Not used for SMMAT-B for the burden test,
SMMAT-S for SKAT or SMMAT-E for the efficient hybrid test of the burden
test and SKAT (default = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1)).
}
  \item{use.minor.allele}{
a logical switch for whether to use the minor allele (instead of the alt allele) as
the coding allele (default = FALSE). It does not change SMMAT-S results, but
SMMAT-B (as well as SMMAT-O and SMMAT-E) will be affected. Along with
the MAF filter, this option is useful for combining rare mutations, assuming rare
allele effects are in the same direction. Use with caution, as major/minor alleles
may flip in different cohorts. In that case, minor allele will be determined based
on the allele frequency in the combined samples.
}
  \item{auto.flip}{
a logical switch for whether to enable automatic allele flipping if a variant with
alleles ref/alt is not found at a position, but a variant at the same position with
alleles alt/ref is found (default = FALSE). Use with caution for whole genome
sequence data, as both ref/alt and alt/ref variants at the same position are not
uncommon, and they are likely two different variants, rather than allele flipping.
}
}

\value{

\code{StocSum.pval} returns a data frame with the following components:
  \item{groups}{
name of the test unit group.
}
  \item{n.variants}{
number of variants in the test unit group that pass the missing rate and allele frequency filters.
}
  \item{B.score}{
burden test score statistic (avaiable if test is Burden, SKATO, or SMMAT).
}
  \item{B.var}{
variance of burden test score statistic (avaiable if test is Burden, SKATO, or SMMAT).
}
  \item{B.pval}{
burden test p-value (avaiable if test is Burden, SKATO, or SMMAT).
}
  \item{S.pval}{
SKAT p-value (avaiable if test is SKATO or SKAT).
}
  \item{O.pval}{
SKAT-O p-value (avaiable if test is SKATO).
}
 \item{O.minp}{
minimum p-value in the SKAT-O search grid (avaiable if test is SKATO).
}
  \item{O.minp.rho}{
rho value at the minimum p-value in the SKAT-O search grid (avaiable if test is SKATO).
}
  \item{O.E.pval}{
efficient hybrid test of the burden test and SKAT p-value (avaiable if test is SMMAT).
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
\code{glmmkin}, \code{StocSum.stat},\code{StocSum.svt}, \code{StocSum.meta}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compact and informatics statics }
