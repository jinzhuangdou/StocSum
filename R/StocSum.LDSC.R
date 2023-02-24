#' Create random vectors for a glmmkin object for LD Score regresion
#' @description Generate random vectors from multivariate normal distribution with mean 0 and covariance matrix P, in which \eqn P = I-1inv((t(1)1))t(1).
#' @param obj The glmmkin object.
#' @param N.randomvec The number of random vectors to generate (default = 1000).
#' @return A list of class glmmkin.randomvec
#' \item{theta}{ set to be 1.}
#' \item{scaled.residuals}{inherited from the glmmkin object. A vector or a matrix for the scaled residuals, calculated as the original residuals divided by the dispersion parameter (in heteroscedastic linear mixed models, corresponding residual variance estimates by each group).}
#' \item{random.vectors}{A random matrix with dimensions equal to the sample size multiplied by \code{N.randomvec}. as.vector(r) is from multivate normal distribution with mean 0 and covariance matrix \eqn P = I-1inv((t(1)1))t(1).}
#' \item{r}{A random matrix with dimensions equal to the sample size multiplied by \code{N.randomvec}. as.vector(r) is from multivate normal distribution with mean 0 and covariance matrix I.}
#' \item{id_include}{inherited from the glmmkin object. A vector indivating the samples included in model fit.}
#' @reference 
#' @author Han Chen, Nannan Wang
#' @examples
#' \donttest{
#' library(GMMAT)
#' data(example)
#' attach(example)
#' seed <- 12345
#' set.seed(seed)
#' GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
#' GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
#' nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
#' obj <- LDSC.glmmkin2randomvec(nullmod)
#' }
#' @keywords random vector
#' @export
LDSC.glmmkin2randomvec <- function(obj, N.randomvec = 1000) {
    if(class(obj) != "glmmkin") stop("Error: \"obj\" must be a class glmmkin object.")
    N <- length(obj$id_include)
    random.vectors <- matrix(rnorm(N*N.randomvec), nrow = N, ncol = N.randomvec)
    r<-random.vectors
    obj$P<-NULL
    obj$Sigma_i<-diag(N)
    obj$Sigma_iX<-matrix(1,N,1)
    obj$cov<-as.matrix(1/N)
    obj$theta<-1
    obj$n.groups<-1
    Z <- NULL
    group.idx = NULL
    if(is.null(group.idx)) group.idx <- rep(1, N)
    random.vectors <- sqrt(obj$theta[group.idx]) * random.vectors
    random.vectors <- crossprod(obj$Sigma_i, random.vectors) - tcrossprod(obj$Sigma_iX, tcrossprod(crossprod(random.vectors, obj$Sigma_iX), obj$cov))
    out <- list(theta = obj$theta, scaled.residuals = obj$scaled.residuals, random.vectors = as.matrix(random.vectors),r=as.matrix(r), id_include = obj$id_include)
    class(out) <- "glmmkin.randomvec"  #LDSC.randomvec
    return(out)
}


####LDSC.stat is same as G.stat
#' Calculate summary statistics and stochastic statistics
#' @description Calculate summary statistics and stochastic statistics.
#' @param null.obj A class glmmkin.randomvec object, returned by generating random vectors using \code{glmmkin2randomvec}.
#' @param geno.file The .gds file name or an object of class SeqVarGDSClass for the full genotypes. The \code{sample.id} in \code{geno.file} should overlap \code{id_include} in \code{null.obj}. It is recommended that \code{sample.id} in \code{geno.file} include the full samples (at least all samples as specified in \code{id_include} of \code{null.obj}). It is not necessary for the user to take a subset of \code{geno.file} before running the analysis. If \code{geno.file} is an object of class SeqVarGDSClass, the .gds file will be closed upon successful completion of the function.
#' @param meta.file.prefix Prefix of intermediate files (*.sample.1 and *.resample.1) required in \code{G.pval}.
#' @param MAF.range A numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param miss.cutoff The maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param missing.method Method of handling missing genotypes. Either "impute2mean" or "impute2zero" (default = "impute2mean").
#' @param nperbatch An integer for how many SNPs to be included in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @param ncores A positive integer indicating the number of cores to be used in parallel computing (default = 1).
#' @return NULL. \code{G.stat} will store the summary statistics and the stochastic statistics in two files with the prefix specified by user.
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{glmmkin2randomvec}, \code{G.pval}
#' @examples
#' \donttest{
#' library(StocSum)
#' data(example)
#' attach(example)
#' seed <- 12345
#' set.seed(seed)
#' GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
#' GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
#' nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
#' if(!is.null(nullmod$P)){
#'   obj <- glmmkin2randomvec(nullmod)
#' }else{
#'   kinship.chol <- chol(GRM)
#'   obj<-glmmkin2randomvec(nullmod, Z = list(t(kinship.chol)))
#' }
#' out.prefix <- "test"
#' gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
#' G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix, MAF.range=c(0,0.5), miss.cutoff = 1)
#' }
#' @keywords summary statistics
#' @export
G.stat <- function(null.obj, geno.file, meta.file.prefix, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, missing.method = "impute2mean", nperbatch = 10000, ncores = 1)
{
  if(Sys.info()["sysname"] == "Windows" && ncores > 1) {
    warning("The package doMC is not available on Windows... Switching to single thread...")
    ncores <- 1
  }
  missing.method <- try(match.arg(missing.method, c("impute2mean", "impute2zero")))
  if(class(missing.method) == "try-error") stop("Error: \"missing.method\" must be \"impute2mean\" or \"impute2zero\".")
  if(class(null.obj) != "glmmkin.randomvec") stop("Error: \"null.obj\" must be a class glmmkin.randomvec object.")
  if(any(duplicated(null.obj$id_include))) {
    J <- sapply(unique(null.obj$id_include), function(x) 1*(null.obj$id_include==x))
    residuals <- crossprod(J, null.obj$scaled.residuals)
    residuals2 <- crossprod(J, null.obj$random.vectors)
    rm(J)
  } else {
    residuals <- null.obj$scaled.residuals
    residuals2 <- null.obj$random.vectors
  }
  if(!grepl("\\.gds$", geno.file)) stop("Error: currently only .gds format is supported in geno.file!")
  gds <- SeqArray::seqOpen(geno.file)
  sample.id <- SeqArray::seqGetData(gds, "sample.id")
  if(any(is.na(match(null.obj$id_include, sample.id)))) warning("Check your data... Some individuals in null.obj$id_include are missing in sample.id of geno.file!")
  sample.id <- sample.id[sample.id %in% null.obj$id_include]
  if(length(sample.id) == 0) stop("Error: null.obj$id_include does not match sample.id in geno.file!")
  match.id <- match(sample.id, unique(null.obj$id_include))
  residuals <- residuals[match.id]
  residuals2 <- residuals2[match.id, , drop = FALSE]
  variant.idx.all <- SeqArray::seqGetData(gds, "variant.id")
  SeqArray::seqClose(gds)
  p.all <- length(variant.idx.all)
  ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
  if(ncores > 1) {
    doMC::registerDoMC(cores = ncores)
    p.percore <- (p.all-1) %/% ncores + 1
    n.p.percore_1 <- p.percore * ncores - p.all
    foreach(b=1:ncores, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
      variant.idx <- if(b <= n.p.percore_1) variant.idx.all[((b-1)*(p.percore-1)+1):(b*(p.percore-1))] else variant.idx.all[(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1-1)*p.percore+1):(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1)*p.percore)]
      p <- length(variant.idx)
      if(.Platform$endian!="little") stop("Error: platform must be little endian.")
      meta.file.sample <- paste0(meta.file.prefix, ".sample.", b)
      meta.file.resample <- paste0(meta.file.prefix, ".resample.", b)
      write.table(t(c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")), meta.file.sample, quote = FALSE, row.names = FALSE, col.names = FALSE)
      meta.file.resample.handle <- file(meta.file.resample, "wb")
      writeBin(as.integer(ncol(residuals2)), meta.file.resample.handle, size = 4)
      gds <- SeqArray::seqOpen(geno.file)
      SeqArray::seqSetFilter(gds, sample.id = sample.id)
      nbatch.flush <- (p-1) %/% nperbatch + 1
      for(i in 1:nbatch.flush) {
        gc()
        tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*nperbatch+1):p] else variant.idx[((i-1)*nperbatch+1):(i*nperbatch)]
        SeqArray::seqSetFilter(gds, variant.id  = tmp.variant.idx, verbose = FALSE)
        miss <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
        freq <- 1 - SeqVarTools::alleleFrequency(gds)
        include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
        if(sum(include) == 0) next
        tmp.variant.idx <- tmp.variant.idx[include]
        tmp.p <- length(tmp.variant.idx)
        SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx)
        SNP <- SeqArray::seqGetData(gds, "annotation/id")
        SNP[SNP == ""] <- NA
        out <- data.frame(SNP = SNP, chr = SeqArray::seqGetData(gds, "chromosome"), pos = SeqArray::seqGetData(gds, "position"))
        rm(SNP)
        alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
        out$ref <- unlist(lapply(alleles.list, function(x) x[1]))
        out$alt <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
        out$missrate <- miss[include]
        out$altfreq <- freq[include]
        rm(alleles.list, include)
        SeqArray::seqSetFilter(gds, variant.id  = tmp.variant.idx, verbose = FALSE)
        geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
        out$N <- nrow(geno) - colSums(is.na(geno))
        if(max(out$missrate)>0) {
          miss.idx <- which(is.na(geno))
          geno[miss.idx] <- if(missing.method=="impute2mean") 2*out$altfreq[ceiling(miss.idx/nrow(geno))] else 0
        }
        out$SCORE <- as.vector(crossprod(geno, residuals))
        write.table(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")], meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
        writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
        rm(out)
      }
      SeqArray::seqClose(gds)
      close(meta.file.resample.handle)
    }
  } else { # use a single core
    variant.idx <- variant.idx.all
    rm(variant.idx.all)
    p <- length(variant.idx)
    if(.Platform$endian!="little") stop("Error: platform must be little endian.")
    meta.file.sample <- paste0(meta.file.prefix, ".sample.1")
    meta.file.resample <- paste0(meta.file.prefix, ".resample.1")
    write.table(t(c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")), meta.file.sample, quote = FALSE, row.names = FALSE, col.names = FALSE)
    meta.file.resample.handle <- file(meta.file.resample, "wb")
    writeBin(as.integer(ncol(residuals2)), meta.file.resample.handle, size = 4)
    gds <- SeqArray::seqOpen(geno.file)
    SeqArray::seqSetFilter(gds, sample.id = sample.id)
    nbatch.flush <- (p-1) %/% nperbatch + 1
    for(i in 1:nbatch.flush) {
      gc()
      tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*nperbatch+1):p] else variant.idx[((i-1)*nperbatch+1):(i*nperbatch)]
      SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
      miss <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
      freq <- 1 - SeqVarTools::alleleFrequency(gds)
      include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
      if(sum(include) == 0) next
      tmp.variant.idx <- tmp.variant.idx[include]
      tmp.p <- length(tmp.variant.idx)
      SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx)
      SNP <- SeqArray::seqGetData(gds, "annotation/id")
      SNP[SNP == ""] <- NA
      out <- data.frame(SNP = SNP, chr = SeqArray::seqGetData(gds, "chromosome"), pos = SeqArray::seqGetData(gds, "position"))
      rm(SNP)
      alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
      out$ref <- unlist(lapply(alleles.list, function(x) x[1]))
      out$alt <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
      out$missrate <- miss[include]
      out$altfreq <- freq[include]
      rm(alleles.list, include)
      SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
      geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
      out$N <- nrow(geno) - colSums(is.na(geno))
      if(max(out$missrate)>0) {
        miss.idx <- which(is.na(geno))
        geno[miss.idx] <- if(missing.method=="impute2mean") 2*out$altfreq[ceiling(miss.idx/nrow(geno))] else 0
      }
      out$SCORE <- as.vector(crossprod(geno, residuals))
      write.table(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")], meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
      writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
      rm(out)
    }
    SeqArray::seqClose(gds)
    close(meta.file.resample.handle)
  } 
  return(invisible(NULL))
}

#' Calculate LDscore
#' @description Use summary statistics and stochastic statistics from G.stat to estimate LDscore.
#' @param meta.file.prefix A character for prefix of intermediate files (*.sample.* and *.resample.*). 
#' @param n.files An integer indicating how many sets of intermediate files (.score.* and .var.*), usually as the result of multi-threading in creating the intermediate files (default = 1)).
#' @param N.randomvec The number of replicates to simulate the random vectors in \code{StocSum.LDSC.glmmkin2randomvec} (default = 1000).
#' @param MAF.range A numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param miss.cutoff The maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param wind.b A positive integer define the window size in bases (default = 1000000).
#' @param use.minor.allele A logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SKAT results, but Burden (as well as SKAT-O and hybrid test to combine the burden test and SKAT) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction. Use with caution, as major/minor alleles may flip in different cohorts. In that case, minor allele will be determined based on the allele frequency in the combined samples.
#' @param auto.flip A logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
#' @param nperbatch An integer for how many SNPs to be included in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @param ncores A positive integer indicating the number of cores to be used in parallel computing (default = 1).
#' @return \code{LDSC.win} returns a data frame with the following components:
#' \item{SNP}{SNP name.}
#' \item{chr}{chromosome name.}
#' \item{pos}{the genome location of SNP.}
#' \item{ref}{allele of reference.}
#' \item{alt}{alternative allele.}
#' \item{N}{total sample size.}
#' \item{altfreq}{alternative allele frequency.}
#' \item{idx}{indicating SNPs.}
#' \item{pos.end}{To which SNP (position) is included to esitmate LD score for the current SNP.}
#' \item{pos.start}{From which SNP (position) is included to estimate LD score for the current SNP.}
#' \item{LDscore}{LD score for each SNP.}
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{LDSC.glmmkin2randomvec}, \code{G.stat}
#' @examples
#' \donttest{
#' library(GMMAT)
#' data(example)
#' attach(example)
#' seed <- 12345
#' set.seed(seed)
#' GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
#' GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
#' nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
#' obj <- LDSC.glmmkin2randomvec(nullmod)
#' out.prefix <- "test"
#' gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
#' G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix, MAF.range=c(0,0.5), miss.cutoff = 1)
#' out<-LDSC.win(out.prefix, use.minor.allele = FALSE, auto.flip = FALSE, wind.b = 1000000, nperbatch = 10000))
#' }
#' @keywords LD Score
#' @export

LDSC.win <- function(meta.files.prefix, n.files = 1, N.randomvec=1000, MAF.range = c(1e-3, 0.5), miss.cutoff = 1, wind.b=1000000, use.minor.allele = FALSE, auto.flip = FALSE, nperbatch = 10000, ncores = 1)
{
  MAF.weights.beta <- c(0.5, 0.5)
  if(.Platform$endian!="little") stop("Error: platform must be little endian.")
  group.info <- NULL
  for(j in 1:n.files) {
    tmp <- fread(paste0(meta.files.prefix, ".sample.", j), header=TRUE, data.table = FALSE)
    if (class(tmp) == "try-error") {
      stop(paste0("Error: cannot read ", meta.files.prefix, ".sample.", j, "!"))
    }
    tmp <- tmp[,c("chr", "pos","N","missrate", "altfreq", "SCORE" )]
    tmp <- tmp[!duplicated(tmp$pos), , drop = FALSE]
    tmp$file <- j
    tmp$variant.idx <- 1:nrow(tmp)
    if(!is.null(group.info)) {
        tmp <- subset(tmp, !pos %in% group.info$pos)
    } 
    if(nrow(tmp) > 0) group.info <- rbind(group.info, tmp)
    rm(tmp)
  }
  group.info <- group.info[order(group.info$chr, group.info$pos), ]
  group.info <- group.info[!duplicated(group.info$pos), , drop=FALSE]
  n.groups.all<-nrow(group.info)
  group.info$idx <- 1:nrow(group.info)
  group.info$pos.end<-group.info$pos+wind.b
  group.info$pos.end[group.info$pos.end>group.info$pos[n.groups.all]]<-group.info$pos[n.groups.all]
  group.info$pos.start<-group.info$pos-wind.b
  group.info$pos.start[group.info$pos.start<group.info$pos[1]]<-group.info$pos[1]
  scores <- cons <- vector("list", 1)
  N.resampling <- rep(0, 1)
  cons <- file(paste0(meta.files.prefix, ".resample.1"), "rb")
  N.resampling <- readBin(cons, what = "integer", n = 1, size = 4)
  current.lines <- current.cons <- rep(1, 1)
  all.out <- NULL
  ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
  if(ncores > 1) {
    doMC::registerDoMC(cores = ncores)
    n.groups.percore <- (n.groups.all-1) %/% ncores + 1
    n.groups.percore_1 <- n.groups.percore * ncores - n.groups.all
    b <- NULL
    all.out <- foreach(b=1:ncores, .combine=rbind, .multicombine = TRUE, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
      idx0 <- if(b <= n.groups.percore_1) group.info$idx[((b-1)*(n.groups.percore-1)+1):(b*(n.groups.percore-1))] else group.info$idx[(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1-1)*n.groups.percore+1):(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1)*n.groups.percore)]
      n.groups <- length(idx0)
      nbatch.flush <- (n.groups-1) %/% nperbatch + 1
      for(i in 1:nbatch.flush) {
        ii<-idx0[((i-1)*nperbatch+1):min((i*nperbatch),n.groups)]
        itmp.idx<-group.info$idx[ii]
        dis<-group.info$pos-group.info[group.info$idx==itmp.idx[1],]$pos.start
        tmp.idx1<-group.info$idx[dis==min(dis[dis>=0])]
        dis<-group.info$pos-group.info[group.info$idx==itmp.idx[length(itmp.idx)],]$pos.end
        tmp.idx2<-group.info$idx[dis==min(dis[dis>=0])]
        tmp.idx<-group.info$idx[tmp.idx1:tmp.idx2]
        rm(dis)
        U.list <- V.list <- vector("list", 1)
        variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
        tmp.scores <- group.info[group.info$idx %in% tmp.idx, , drop = FALSE]
        if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
          U.list <- tmp.scores[tmp.include, , drop = FALSE]
          tmp.V <- matrix(NA, sum(tmp.include), N.resampling)
          for(ij in 1:sum(tmp.include)) {
            if(U.list$file[ij]!=current.cons) {
              close(cons)
              current.cons <- U.list$file[ij]
              cons <- file(paste0(meta.files.prefix, ".resample.", current.cons), "rb")
              tmp.N.resampling <- readBin(cons, what = "integer", n = 1, size = 4)
              if(tmp.N.resampling != N.resampling) stop(paste0("Error: N.resampling in ", meta.files.prefix, ".resample.", current.cons, " does not match that in ",meta.files.prefix, ".resample.1"))
            }
            current.cons <- U.list$file[ij]
            cons <- file(paste0(meta.files.prefix, ".resample.", current.cons), "rb")
            ind <- 4+4*N.resampling*(U.list$variant.idx[ij]-1)
            seek(cons, where = ind, origin = "start", rw = "read")
            tmp.V[ij,] <- readBin(cons, what = "numeric", n = N.resampling, size = 4)
          }
          V.list <- tmp.V / sqrt(N.resampling)
          rm(tmp.V)
          variant.indices <- c(variant.indices, U.list$idx)
          tmp.N <- c(tmp.N, U.list$N)
          tmp.Nmiss <- c(tmp.Nmiss, U.list$N * U.list$missrate/(1-U.list$missrate))
          tmp.AC <- c(tmp.AC, 2*U.list$N*U.list$altfreq)
        }
        if(length(variant.indices) == 0) next
        tmp.variant.indices <- variant.indices
        variant.indices <- sort(unique(variant.indices))
        N <- sapply(variant.indices, function(x) sum(tmp.N[tmp.variant.indices==x]))
        Nmiss <- sapply(variant.indices, function(x) sum(tmp.Nmiss[tmp.variant.indices==x]))
        AF <- sapply(variant.indices, function(x) sum(tmp.AC[tmp.variant.indices==x]))/2/N
        include <- (Nmiss/(N+Nmiss) <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
        rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices)
        if(sum(include) == 0) next
        variant.indices <- variant.indices[include]
        out <- group.info[match(itmp.idx,group.info$idx),]
        n.p <- length(variant.indices)
        V <- matrix(0, n.p, sum(N.resampling))
        if(!is.null(U.list) & !is.null(V.list)) {
          IDX <- match(U.list$idx, variant.indices)
          if(sum(!is.na(IDX)) == 0) next
          IDX2 <- which(!is.na(IDX))
          IDX <- IDX[IDX2]
          V[IDX, 1:N.resampling] <- V.list[IDX2,]
        }
        n.batchi<-length(itmp.idx)
        tmp.weight <- MAF.weights.beta.fun(AF[include], MAF.weights.beta[1], MAF.weights.beta[2])
        if(use.minor.allele) tmp.weight[AF[include] > 0.5] <- -tmp.weight[AF[include] > 0.5]
        V <- V*tmp.weight*(gamma(MAF.weights.beta[1])*gamma(MAF.weights.beta[2])/gamma(MAF.weights.beta[1]+MAF.weights.beta[2])/sqrt(2))
        LDscore<-rep(NA,n.batchi)
        for (j in 1:n.batchi){
          jidx<-itmp.idx[j]
          if (jidx %in% variant.indices) {
            jtmp.idx<-group.info[group.info$pos>=group.info$pos.start[ii[j]] & group.info$pos<=group.info$pos.end[ii[j]],]$idx
            jV<-V[variant.indices==jidx,]
            jV0<-match(jtmp.idx,variant.indices)
            jV0<-jV0[!is.na(jV0)]
            nvariants<-length(jV0)
            V0<-V[jV0,]
            jN<-N[variant.indices==jidx]
            LDscore[j]<- (1+1/(N.randomvec-2)+1/(jN-2))*rowSums(tcrossprod((tcrossprod(jV,V0))/(jN-1)))-nvariants/(N.randomvec-2)-nvariants/(jN-2)
            rm(jtmp.idx,jV0,V0)
          }
        }
        out$LDscore<-LDscore
        all.out <- rbind(all.out, out)
      }
      return(all.out)
    }
  } else { # use a single core
    n.groups <- n.groups.all  #length(variant.id1)
    nbatch.flush <- (n.groups-1) %/% nperbatch + 1
    for(i in 1:nbatch.flush) {
      itmp.idx<-group.info$idx[((i-1)*nperbatch+1):min((i*nperbatch),n.groups)]
      dis<-group.info$pos-group.info[group.info$idx==itmp.idx[1],]$pos.start
      tmp.idx1<-group.info$idx[dis==min(dis[dis>=0])]
      dis<-group.info$pos-group.info[group.info$idx==itmp.idx[length(itmp.idx)],]$pos.end
      tmp.idx2<-group.info$idx[dis==min(dis[dis>=0])]
      tmp.idx<-group.info$idx[tmp.idx1:tmp.idx2]
      rm(dis)
      U.list <- V.list <- vector("list", 1)
      variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
      tmp.scores <- group.info[group.info$idx %in% tmp.idx, , drop = FALSE]
      if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
        U.list <- tmp.scores[tmp.include, , drop = FALSE]
        tmp.V <- matrix(NA, sum(tmp.include), N.resampling)
        for(ij in 1:sum(tmp.include)) {
          if(U.list$file[ij]!=current.cons) {
              close(cons)
              current.cons <- U.list$file[ij]
              cons <- file(paste0(meta.files.prefix, ".resample.", current.cons), "rb")
              tmp.N.resampling <- readBin(cons, what = "integer", n = 1, size = 4)
              if(tmp.N.resampling != N.resampling) stop(paste0("Error: N.resampling in ", meta.files.prefix, ".resample.", current.cons, " does not match that in ",meta.files.prefix, ".resample.1"))
              current.lines <- 1
          }
          if(U.list$variant.idx[ij]!=current.lines) seek(cons, where = 4*N.resampling*(U.list$variant.idx[ij]-current.lines), origin = "current", rw = "read")
          tmp.V[ij,] <- readBin(cons, what = "numeric", n = N.resampling, size = 4)
          current.lines <- U.list$variant.idx[ij]+1
        }
        V.list <- tmp.V / sqrt(N.resampling)
        rm(tmp.V)
        variant.indices <- c(variant.indices, U.list$idx)
        tmp.N <- c(tmp.N, U.list$N)
        tmp.Nmiss <- c(tmp.Nmiss, U.list$N * U.list$missrate/(1-U.list$missrate))
        tmp.AC <- c(tmp.AC, 2*U.list$N*U.list$altfreq)
      }
      if(length(variant.indices) == 0) next
      tmp.variant.indices <- variant.indices
      variant.indices <- sort(unique(variant.indices))
      N <- sapply(variant.indices, function(x) sum(tmp.N[tmp.variant.indices==x]))
      Nmiss <- sapply(variant.indices, function(x) sum(tmp.Nmiss[tmp.variant.indices==x]))
      AF <- sapply(variant.indices, function(x) sum(tmp.AC[tmp.variant.indices==x]))/2/N
      include <- (Nmiss/(N+Nmiss) <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
      rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices)
      if(sum(include) == 0) next
      variant.indices <- variant.indices[include]
      out <- group.info[match(itmp.idx,group.info$idx),]
      n.p <- length(variant.indices)
      V <- matrix(0, n.p, sum(N.resampling))
      if(!is.null(U.list) & !is.null(V.list)) {
        IDX <- match(U.list$idx, variant.indices)
        if(sum(!is.na(IDX)) == 0) next
        IDX2 <- which(!is.na(IDX))
        IDX <- IDX[IDX2]
        V[IDX, 1:N.resampling] <- V.list[IDX2,]
      }
      n.batchi<-length(itmp.idx)
      tmp.weight <- MAF.weights.beta.fun(AF[include], MAF.weights.beta[1], MAF.weights.beta[2])
      if(use.minor.allele) tmp.weight[AF[include] > 0.5] <- -tmp.weight[AF[include] > 0.5]
      V <- V*tmp.weight*(gamma(MAF.weights.beta[1])*gamma(MAF.weights.beta[2])/gamma(MAF.weights.beta[1]+MAF.weights.beta[2])/sqrt(2))
      LDscore<-rep(NA,n.batchi)
      for (j in 1:n.batchi){
        jidx<-itmp.idx[j]
        if (jidx %in% variant.indices) {
          jtmp.idx<-group.info[group.info$pos>=group.info$pos.start[(i-1)*nperbatch+j] & group.info$pos<=group.info$pos.end[(i-1)*nperbatch+j],]$idx
          jV<-V[variant.indices==jidx,]
          jV0<-match(jtmp.idx,variant.indices)
          jV0<-jV0[!is.na(jV0)]
          nvariants<-length(jV0)
          V0<-V[jV0,]
          jN<-N[variant.indices==jidx]
          LDscore[j]<- (1+1/(N.randomvec-2)+1/(jN-2))*rowSums(tcrossprod((tcrossprod(jV,V0))/(jN-1)))-nvariants/(N.randomvec-2)-nvariants/(jN-2)
          rm(jtmp.idx,jV0,V0)
        }
      }
      out$LDscore<-LDscore
      all.out <- rbind(all.out, out)
    }
  }
  close(cons)
  return(all.out)
}


MAF.weights.beta.fun <- function(freq, beta1, beta2) {
  freq[freq > 0.5] <- 1 - freq[freq > 0.5]
  ifelse(freq <= 0, 0, dbeta(freq, beta1, beta2))
}
