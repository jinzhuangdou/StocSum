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
  #variant set test
  if(is.null(nullmod$P)) {
    kinship <- kins[match(nullmod$id_include, rownames(kins)), match(nullmod$id_include, colnames(kins))]
    kinship.chol <- chol(kinship)
    obj <- glmmkin2randomvec(nullmod, Z = list(t(kinship.chol)), group.idx = NULL)#change group.idx
  } else {
    obj <- glmmkin2randomvec(nullmod, Z = NULL, group.idx = NULL)
  }
  out.prefix <- "test.vst"
  StocSum.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix, MAF.range=c(0,0.5), miss.cutoff = 1)
  #single variant test
  out1 <- StocSum.svt(out.prefix, n.files = ncores, MAF.range=c(1e-7,0.5), miss.cutoff = 0.1, auto.flip=F)
  #variant set test
  out2 <- StocSum.pval(out.prefix, group.file = group.file, MAF.range=c(0,0.5), miss.cutoff = 1)
  # #meta analysis
  # out3 <- StocSum.pval(meta.files.prefix = out.prefix, group.file=File.SetID)
  
  # expect_equal(signif(range(out$B.pval)), signif(c(0.1540811, 0.9500619)))
  # expect_equal(signif(range(out$E.pval)), signif(c(0.01675254, 0.76516832)))   
})

