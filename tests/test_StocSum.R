context("StocSum")

test_that("cross-sectional id le 400 binomial", {
  gdsfile <- system.file("extdata", "geno.gds", package = "StocSum")
  group.file <- system.file("extdata", "SetID.withweights.txt", package = "StocSum")
  data(example)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  pheno <- rbind(example$pheno, example$pheno[1:100, ])
  pheno$id <- 1:500
  pheno$disease[sample(1:500,20)] <- NA
  pheno$age[sample(1:500,20)] <- NA
  pheno$sex[sample(1:500,20)] <- NA
  pheno <- pheno[sample(1:500,450), ]
  pheno <- pheno[pheno$id <= 400, ]
  kins <- example$GRM
  nullmod <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"))
  if(is.null(nullmod$P)) {
    kinship <- kins[match(nullmod$id_include, rownames(kins)), match(nullmod$id_include, colnames(kins))]
    kinship.chol <- chol(kinship)
    obj <- glmmkin2randomvec(nullmod, Z = list(t(kinship.chol)), group.idx = NULL)#change group.idx
  } else {
    obj <- glmmkin2randomvec(nullmod, Z = NULL, group.idx = NULL)
  }
  out.prefix <- "test"
  # out1 <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
  StocSum.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix, MAF.range=c(0,0.5), miss.cutoff = 1)
  out<-StocSum.pval(out.prefix, group.file = group.file, MAF.range=c(0,0.5), miss.cutoff = 1)
    
  # expect_equal(signif(range(out$B.pval)), signif(c(0.1540811, 0.9500619)))
  # expect_equal(signif(range(out$E.pval)), signif(c(0.01675254, 0.76516832)))   
})

