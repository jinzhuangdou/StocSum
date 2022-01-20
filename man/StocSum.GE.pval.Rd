\name{StocSum.GE.pval}
\alias{StocSum.GE.pval}
\title{
Mixed model association tests for gene-environment interactions
}
\description{
Use the summary statistic from \code{StocSum.GE.stat} to perform variant set-based main effect tests, gene-environment interaction tests, and joint tests for association with genotypes in a GDS file (.gds). 7 user-defined tests are included: Main effect variance component test (MV), Main effect hybrid test of burden and variance component test using Fisher's method (MF), Interaction variance component test (IV), Interaction hybrid test of burden and variance component test using Fisher's method (IF), Joint variance component test (JV), Joint hybrid test of burden and variance component test using Fisher's method (JF), and Joint hybrid test of burden and variance component test using double Fisher's procedures (JD). 
}
\usage{
StocSum.GE.pval(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), n.pheno = 1, cohort.group.idx = NULL, group.file, group.file.sep = "\t", interatcion, interaction.covariates = NULL, MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(0.5, 0.5), miss.cutoff = 1, method = "davies", tests = "JF", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE, auto.flip = FALSE)

}
\arguments{
  \item{meta.file.prefix}{
a character vector for prefix of intermediate files (*.sample.* and *.resample.*). 
}
  \item{n.files}{
an integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix)))
}
 \item{n.pheno}{
an integer indicating the number of phenotypes in multiple phenotype analysis (for single phenotype analysis, \code{n.pheno = 1}) (default = 1).
}
 \item{cohort.group.idx}{
a vector with the same length as meta.files.prefix, indicating which cohorts should be grouped together in the meta-analysis assuming homogeneous genetic effects. For example, c("a","b","a","a","b") means cohorts 1, 3, 4 are assumed to have homogeneous genetic effects, and cohorts 2, 5 are in another group with homogeneous genetic effects (but possibly heterogeneous with group "a"). If NULL, all cohorts are in the same group (default = NULL).
}
 \item{group.file}{
a plain text file with 6 columns defining the test units. There should be no headers in the file, and the columns are group name, chromosome, position, reference allele, alternative allele and weight, respectively.
}
 \item{group.file.sep}{
the delimiter in group.file (default = "\\t").
}
 \item{interaction}{
a numeric or a character vector indicating the environmental factors. If a numberic vector, it represents which indices in the order of covariates are the environmental factors; if a character vector, it represents the variable names of the environmental factors.
}
 \item{interaction.covariates}{
a numeric or a character vector indicating the interaction covariates. If a numeric vector, it represents which indices in the order of covariates are the interaction covariates; if a character vector, it represents the variable names of the interaction covariates.
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
}
  \item{MAF.weight.beta}{
a numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies. This internal minor allele frequency weight is multiplied by the external weight given by the group.file. To turn off internal minor allele frequency weight and only use the external weight given by the group.file, use c(1, 1) to assign flat weights (default = c(0.5, 0.5)). Applied to the combined samples.
}
  \item{missing.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
}
  \item{method}{
a method to compute p-values for SKAT-type test statistics (default = "davies"). "davies" represents an exact method that computes a p-value by inverting the characteristic function of the mixture chisq distribution, with an accuracy of 1e-6. When "davies" p-value is less than 1e-5, it defaults to method "kuonen". "kuonen" represents a saddlepoint approximation method that computes the tail probabilities of the mixture chisq distribution. When "kuonen" fails to compute a p-value, it defaults to method "liu". "liu" is a moment-matching approximation method for the mixture chisq distribution.
}
  \item{tests}{
a character vector indicating which tests should be performed ("MV" for the main effect variance component test, "MF" for the main effect combined test of the burden and variance component tests using Fisher's method, "IV" for the interaction variance component test, "IF" for the interaction combined test of the burden and variance component tests using Fisher's method, "JV" for the joint variance component test for main effect and interaction, "JF" for the joint combined test of the burden and variance component tests for main effect and interaction using Fisher's method, or "JD" for the joint combined test of the burden and variance component tests for main effect and interaction using double Fisher's method.). The "MV" and "IV" test are automatically included when performing "JV", and the "MF" and "IF" test are automatically included when performing "JF" or "JD" (default = "JF").
}
  \item{rho}{
a numeric vector defining the search grid used in SKAT-O (see the SKAT-O paper for details). Not used for the burden test, SKAT or the efficient hybrid test of the burden test and SKAT (default = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1)).
}
  \item{use.minor.allele}{
a logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SKAT results, but Burden (as well as SKAT-O and hybrid test to combine the burden test and SKAT) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction. Use with caution, as major/minor alleles may flip in different cohorts. In that case, minor allele will be determined based on the allele frequency in the combined samples.
}
  \item{auto.flip}{
a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
}
}
\value{
\code{StocSum.GE.pval} returns a data frame with the following components:
  \item{groups}{
name of the test unit group.
}
  \item{n.variants}{
number of variants in the test unit group that pass the missing rate and allele frequency filters.
}
  \item{MV.pval}{
MV test p-value.
}
  \item{MF.pval}{
MF test p-value.
}
  \item{IV.pval}{
IV test p-value.
}
  \item{IF.pval}{
IF test p-value.
}
  \item{JV.pval}{
JV test p-value.
}
  \item{JF.pval}{
JF test p-value.
}
  \item{JD.pval}{
JD test p-value.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Han Chen, Nannan Wang
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{glmmkin2randomvec}, \code{StocSum.GE.stat}
}
\examples{
\donttest{
library(GMMAT)
data(example)
attach(example)
GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
obj <- StocSum.GE.glmmkin2randomvec(nullmod)
out.prefix <- "test.GE"
gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
interaction <- c("sex")
StocSum.GE.stat(obj, interaction = interaction, geno.file = gdsfile, meta.file.prefix = out.prefix)
out <- StocSum.GE.pval(out.prefix, group.file = group.file, tests=c("JV", "JF", "JD"))
}
}
\keyword{ variant set-based test }
\keyword{ generalized linear mixed model }
\keyword{ gene-environment interaction }
\keyword{ joint test }
