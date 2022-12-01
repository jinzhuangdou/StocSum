#' Create random vectors for a glmmkin object
#' @description Generate random vectors from multivariate normal distribution with mean 0 and covariance matrix P.
#' @param obj The glmmkin object.
#' @param Z A list of design matrices for the random effects. The length must match the number of variance components.
#' @param N.randomvec The number of random vectors to generate (default = 1000).
#' @param group.idx A length N index vector showing which observation belongs to which variance group, for heteroscedastic linear mixed models (default = NULL for homoscedastic linear mixed models).
#' @param cluster.idx A length N index vector showing which observation belongs to which cluster (default = NULL for no clusters).
#' @param robust A logical switch: whether robust variance should be used (default = FALSE).
#' @return A list of class glmmkin.randomvec
#' \item{theta}{Variance estimates, inherited from the glmmkin object.}
#' \item{scaled.residuals}{Scaled residuals, inherited from the glmmkin object.}
#' \item{random.vectors}{An N by N.randomvec matrix of the random vectors generated.}
#' \item{X}{model matrix for the fixed effects from the glmmkin object}
#' \item{id_include}{The ID vector of included samples, inherited from the glmmkin object.}
#' @reference 
#' @author Han Chen, Nannan Wang
#' @examples
#' \donttest{
#' library(GMMAT)
#' data(example)
#' attach(example)
#' GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
#' GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
#' nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
#' if(!is.null(nullmod$P)){
#'   obj <- glmmkin2randomvec(nullmod)
#' }else{
#'   kinship.chol <- chol(GRM)
#'   obj<-glmmkin2randomvec(nullmod, Z = list(t(kinship.chol)))
#' }
#' }
#' @keywords random vector
#' @export

glmmkin2randomvec <- function(obj, Z = NULL, N.randomvec = 1000, group.idx=NULL, cluster.idx=NULL, robust = FALSE) 
{
    if(class(obj) != "glmmkin") stop("Error: \"obj\" must be a class glmmkin object.")
    N <- length(obj$id_include)
    random.vectors <- matrix(rnorm(N*N.randomvec),nrow = N,ncol = N.randomvec)
    if(!is.null(obj$P) && !robust) {
        eig <- eigen(obj$P, symmetric = TRUE)
        random.vectors <- tcrossprod(eig$vectors, t(random.vectors * sqrt(pmax(eig$values, 0))))
        rm(eig)
    } else {
        if(obj$n.groups != 1 && (is.null(group.idx) || !all.equal(seq_len(obj$n.groups), sort(unique(group.idx))))) stop("Error: heteroscedastic linear mixed models should include a valid group.idx argument.")
        if(is.null(group.idx)) group.idx <- rep(1, N)
        if(!robust) random.vectors <- sqrt(obj$theta[group.idx]) * random.vectors
#       else random.vectors <- random.vectors * abs(obj$residuals)
        else {
            res <- as.numeric(obj$Y - tcrossprod(obj$X, t(obj$coefficient)))
            if(is.null(cluster.idx)) random.vectors <- random.vectors * res
            else random.vectors <- random.vectors[match(cluster.idx, unique(cluster.idx)),] * res
        }
        if(!is.null(Z)) {
            if(class(Z) != "list") stop("Error: \"Z\" must be a list of matrices.")
            if(length(Z) != length(obj$theta) - obj$n.groups) stop("Error: number of matrices in \"Z\" does not match the number of variance components in \"obj\".")
            for(i in 1:length(Z)) {
                if(nrow(Z[[i]]) != N) stop("Error: \"Z\" matrix ", i, " is not compatible in sample size with \"obj\".")
                p <- ncol(Z[[i]])
                if(obj$theta[i+obj$n.groups] < 0) stop("Error: negative variance component estimates are not allowed.")
                if(obj$theta[i+obj$n.groups] == 0) next
                random.vectors2 <- matrix(rnorm(p*N.randomvec), nrow=N.randomvec, ncol=p)
                random.vectors <- random.vectors + sqrt(obj$theta[i+obj$n.groups]) * tcrossprod(Z[[i]], random.vectors2)
            }
        }
        if(!is.null(obj$P)) random.vectors <- crossprod(obj$P, random.vectors)
        else random.vectors <- crossprod(obj$Sigma_i, random.vectors) - tcrossprod(obj$Sigma_iX, tcrossprod(crossprod(random.vectors, obj$Sigma_iX), obj$cov))
    }
    out <- list(theta = obj$theta, scaled.residuals = obj$scaled.residuals, random.vectors = as.matrix(random.vectors), id_include = obj$id_include, X = obj$X)
    class(out) <- "glmmkin.randomvec"
    return(out)
}


#' Calculate stochastic statistics
#' @description Calculate summary statistics and stochastic statistics.
#' @param null.obj a class glmmkin.randomvec object, returned by generating random vectors using \code{glmmkin2randomvec}.
#' @param geno.file the .gds file name or an object of class SeqVarGDSClass for the full genotypes. The \code{sample.id} in \code{geno.file} should overlap \code{id_include} in \code{null.obj}. It is recommended that \code{sample.id} in \code{geno.file} include the full samples (at least all samples as specified in \code{id_include} of \code{null.obj}). It is not necessary for the user to take a subset of \code{geno.file} before running the analysis. If \code{geno.file} is an object of class SeqVarGDSClass, the .gds file will be closed upon successful completion of the function.
#' @param meta.file.prefix prefix of intermediate files (*.sample.1 and *.resample.1) required in \code{G.pval}.
#' @param MAF.range a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param missing.cutoff the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param missing.method method of handling missing genotypes. Either "impute2mean" or "impute2zero" (default = "impute2mean").
#' @param nperbatch an integer for how many SNPs to be included in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @param ncores a positive integer indicating the number of cores to be used in parallel computing (default = 1).
#' @return NULL. \code{G.stat} will store the summary statistics and the stochastic statistics in two files with the prefix specified by the user.
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{glmmkin2randomvec}, \code{G.pval}
#' @examples
#' \donttest{
#' library(GMMAT)
#' data(example)
#' attach(example)
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
#' G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1)
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


#' Performing meta-analysis for GLMM bashed socre test results
#' @description Use summary statistics and stochastic statistics from \code{G.stat} to peform meta-analysis.
#' @param meta.file.prefix a character vector for prefix of intermediate files (*.sample.* and *.resample.*) which stores summary statistics and stochastic statistics calculated from \code{G.stat}.
#' @param n.files an integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix))).
#' @param outfile.prefix a character vector for prefix of intermediate files (*.sample.* and *.resample.*) which stores summary statistics and stochatic statistics of meta-analysis.
#' @param MAF.range a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param missing.cutoff the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param auto.flip a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
#' @param nperbatch an integer for how many SNPs to be included in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @return NULL. \code{svt.meta} will store the summary statistics and the stochastic in two files with the prefix specified by the user.
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{glmmkin2randomvec}, \code{G.stat}, \code{G.pval}
#' @examples
#' \donttest{
#' library(GMMAT)
#' data(example)
#' attach(example)
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
#' G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1)
#' NOT DONE!!!!!!!!!!!!!!!!
#' }
#' @keywords meta-analysis
#' @export

svt.meta <- function(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), outfile.prefix, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, auto.flip = FALSE, nperbatch = 10000)
{
    if(.Platform$endian!="little") stop("Error: platform must be little endian.")
    n.cohort <- length(meta.files.prefix)
    if(length(n.files) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and n.files do not match.")
    group.info <- NULL
    for(i in 1:n.cohort) {
        for(j in 1:n.files[i]) {
            tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
            if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
            }
            tmp <- tmp[,c("SNP", "chr", "pos", "ref", "alt")]
            tmp$snpid <- paste(tmp$chr, tmp$pos, tmp$ref, tmp$alt, sep = ":")
            tmp <- tmp[!duplicated(tmp$snpid), , drop = FALSE]
            if(auto.flip) {
                snpid2 <- paste(tmp$chr, tmp$pos, tmp$alt, tmp$ref, sep = ":")
                match.snpid <- match(snpid2, tmp$snpid)
                tmp <- tmp[is.na(match.snpid) | match.snpid > 1:length(match.snpid), , drop = FALSE]
                snpid2 <- snpid2[is.na(match.snpid) | match.snpid > 1:length(match.snpid)]
                rm(match.snpid)
            }
            if(!is.null(group.info)) {
                if(auto.flip) {
                    tmp <- subset(tmp, !snpid %in% group.info$snpid & !snpid2 %in% group.info$snpid)
                    rm(snpid2)
                } else tmp <- subset(tmp, !snpid %in% group.info$snpid)
            } else if(auto.flip) {
                rm(snpid2)
            }
            if(nrow(tmp) > 0) group.info <- rbind(group.info, tmp)
            rm(tmp)
        }
    }
    group.info <- group.info[order(group.info$chr, group.info$pos), ]
    p <- nrow(group.info)
    group.info$idx <- 1:p
    variant.id1 <- group.info$snpid
    scores <- cons <- vector("list", n.cohort)
    N.resampling <- rep(0, n.cohort)
    if(auto.flip) {
        cat("Automatic allele flipping enabled...\nVariants matching alt/ref but not ref/alt alleles will also be included, with flipped effects\n")
        variant.id2 <- paste(group.info$chr, group.info$pos, group.info$alt, group.info$ref, sep = ":")
    }
    for(i in 1:n.cohort) {
        tmp.scores <- NULL
        for(j in 1:n.files[i]) {
            tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
            if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
            }
            tmp <- tmp[,c("chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")]
            tmp$variant.idx <- 1:nrow(tmp)
            variant.id <- paste(tmp$chr, tmp$pos, tmp$ref, tmp$alt, sep = ":")
            variant.idx1 <- variant.id %in% variant.id1
            if(auto.flip) {
                variant.idx2 <- variant.id %in% variant.id2
                if(any(variant.idx1 & variant.idx2)) {
                    tmp.dups <- which(variant.idx1 & variant.idx2)
                    cat("The following ambiguous variants were found in",paste0(meta.files.prefix[i], ".sample.", j),":\n")
                    cat("chr:", tmp$chr[tmp.dups], "\n")
                    cat("pos:", tmp$pos[tmp.dups], "\n")
                    cat("ref:", tmp$ref[tmp.dups], "\n")
                    cat("alt:", tmp$alt[tmp.dups], "\n")
                    cat("Warning: both variants with alleles ref/alt and alt/ref were present at the same position and coding should be double checked!\nFor these variants, only those with alleles ref/alt were used in the analysis...\n")
                    variant.idx2[tmp.dups] <- FALSE
                    rm(tmp.dups)
                }
                tmpallele <- tmp$ref[variant.idx2]
                tmp$ref[variant.idx2] <- tmp$alt[variant.idx2]
                tmp$alt[variant.idx2] <- tmpallele
                rm(tmpallele)
                tmp$altfreq[variant.idx2] <- 1-tmp$altfreq[variant.idx2]
                tmp$flip <- 1-2*(variant.idx2)
                tmp <- tmp[variant.idx1 | variant.idx2, ]
            } else {
                tmp$flip <- 1
                tmp <- tmp[variant.idx1, ]
            }
            tmp$file <- j
            tmp.scores <- rbind(tmp.scores, tmp)
            rm(tmp)
        }
        scores[[i]] <- cbind(group.info, tmp.scores[match(variant.id1, paste(tmp.scores$chr, tmp.scores$pos, tmp.scores$ref, tmp.scores$alt, sep = ":")), c("N", "missrate", "altfreq", "SCORE", "file", "variant.idx", "flip")])
        rm(tmp.scores)
        cons[[i]] <- file(paste0(meta.files.prefix[i], ".resample.1"), "rb")
        N.resampling[i] <- readBin(cons[[i]], what = "integer", n = 1, size = 4)
    }
    minN <- min(N.resampling)
    meta.file.sample <- paste0(outfile.prefix, ".sample.1")
    meta.file.resample <- paste0(outfile.prefix, ".resample.1")
    write.table(t(c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")), meta.file.sample, quote = FALSE, row.names = FALSE, col.names = FALSE)
    meta.file.resample.handle <- file(meta.file.resample, "wb")
    writeBin(as.integer(minN), meta.file.resample.handle, size = 4)
    current.lines <- current.cons <- rep(1, n.cohort)
    nbatch.flush <- (p-1) %/% nperbatch + 1
    for(i in 1:nbatch.flush) {
        tmp.idx <- if(i == nbatch.flush) group.info$idx[((i-1)*nperbatch+1):p] else group.info$idx[((i-1)*nperbatch+1):(i*nperbatch)]
        #tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
        U.list <- V.list <- vector("list", n.cohort)
        variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
        for(j in 1:n.cohort) {
            tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
            if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
                U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
                U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$flip
                tmp.V <- matrix(NA, sum(tmp.include), N.resampling[j])
                for(ij in 1:sum(tmp.include)) {
                    if(U.list[[j]]$file[ij]!=current.cons[j]) {
                        close(cons[[j]])
                        current.cons[j] <- U.list[[j]]$file[ij]
                        cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.", current.cons[j]), "rb")
                        tmp.N.resampling <- readBin(cons[[j]], what = "integer", n = 1, size = 4)
                        if(tmp.N.resampling != N.resampling[j]) stop(paste0("Error: N.resampling in ", meta.files.prefix[j], ".resample.", current.cons[j], " does not match that in ",meta.files.prefix[j], ".resample.1"))
                        current.lines[j] <- 1
                    }
                    if(U.list[[j]]$variant.idx[ij]!=current.lines[j]) seek(cons[[j]], where = 4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-current.lines[j]), origin = "current", rw = "read")
                    tmp.V[ij,] <- readBin(cons[[j]], what = "numeric", n = N.resampling[j], size = 4)
                    current.lines[j] <- U.list[[j]]$variant.idx[ij]+1
                }
                V.list[[j]] <- tmp.V * U.list[[j]]$flip
                rm(tmp.V)
                variant.indices <- c(variant.indices, U.list[[j]]$idx)
                tmp.N <- c(tmp.N, U.list[[j]]$N)
                tmp.Nmiss <- c(tmp.Nmiss, U.list[[j]]$N * U.list[[j]]$missrate/(1-U.list[[j]]$missrate))
                tmp.AC <- c(tmp.AC, 2*U.list[[j]]$N*U.list[[j]]$altfreq)
            }
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
        out <- data.frame(SNP = group.info$SNP[variant.indices], chr = group.info$chr[variant.indices], pos = group.info$pos[variant.indices], ref = group.info$ref[variant.indices], alt = group.info$alt[variant.indices], N = N[include], missrate = (Nmiss/(N+Nmiss))[include], altfreq = AF[include])
        n.p <- length(variant.indices)
        U <- rep(0, n.p)
        V <- matrix(0, n.p, minN)
        for(j in 1:n.cohort) {
            if(!is.null(U.list[[j]]) & !is.null(V.list[[j]])) {
                IDX <- match(U.list[[j]]$idx, variant.indices)
                if(sum(!is.na(IDX)) == 0) next
                IDX2 <- which(!is.na(IDX))
                IDX <- IDX[IDX2]
                U[IDX] <- U[IDX]+U.list[[j]]$SCORE[IDX2]
                V[IDX, ] <- V[IDX, , drop = FALSE]+V.list[[j]][IDX2, 1:minN, drop = FALSE]
            }
        }
        out$SCORE <- U
        write.table(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")], meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
        writeBin(as.numeric(t(V)), meta.file.resample.handle, size = 4)
        rm(out)
    }
    close(meta.file.resample.handle)
    for(i in 1:n.cohort) close(cons[[i]])
    return(invisible(NULL))
}


#' Performing GLMM based score tests
#' @description Perform score tests for association with genotypes in a GDS file .gds, with summary statistics calculated from \code{glmmkin2randomvec}.
#' @param meta.file.prefix a character vector for prefix of intermediate files (*.sample.* and *.resample.*). 
#' @param n.files an integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix))).
#' @param MAF.range a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param MAF.weight.beta a numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies. This internal minor allele frequency weight is multiplied by the external weight given by the group.file. To turn off internal minor allele frequency weight and only use the external weight given by the group.file, use c(1, 1) to assign flat weights (default = c(1, 25)). Applied to the combined samples.
#' @param missing.cutoff the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param auto.flip a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
#' @param nperbatch an integer for how many SNPs should be tested in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @return \code{svt.pval} returns a data frame with the following components:
#' \item{SNP}{SNP name.}
#' \item{chr}{chromosome name.}
#' \item{pos}{the genome location of SNP.}
#' \item{ref}{allele of reference.}
#' \item{alt}{alternative allele.}
#' \item{N}{total sample size.}
#' \item{missrate}{missing rate of variants.}
#' \item{altfreq}{alternative allele frequency.}
#' \item{SCORE}{the summary score of the alternaive allele.}
#' \item{VAL}{the variance of the summary score.}
#' \item{PVAL}{the p-value of the suammry score.}
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{glmmkin2randomvec}, \code{G.stat}
#' @examples
#' \donttest{
#' library(GMMAT)
#' data(example)
#' attach(example)
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
#' G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1)
#' out <- svt.pval(out.prefix, MAF.range=c(1e-7, 0.5), miss.cutoff = 1, auto.flip=F)
#' }
#' @keywords score test, generalized linear mixed model
#' @export

svt.pval <- function(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), MAF.range = c(1e-7, 0.5), miss.cutoff = 1, auto.flip = FALSE, nperbatch = 10000)
{
    if(.Platform$endian!="little") stop("Error: platform must be little endian.")
    n.cohort <- length(meta.files.prefix)
    if(length(n.files) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and n.files do not match.")
    group.info <- NULL
    for(i in 1:n.cohort) {
        for(j in 1:n.files[i]) {
            tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
            if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
            }
            tmp <- tmp[,c("SNP", "chr", "pos", "ref", "alt")]
            tmp$snpid <- paste(tmp$chr, tmp$pos, tmp$ref, tmp$alt, sep = ":")
            tmp <- tmp[!duplicated(tmp$snpid), , drop = FALSE]
            if(auto.flip) {
                snpid2 <- paste(tmp$chr, tmp$pos, tmp$alt, tmp$ref, sep = ":")
                match.snpid <- match(snpid2, tmp$snpid)
                tmp <- tmp[is.na(match.snpid) | match.snpid > 1:length(match.snpid), , drop = FALSE]
                snpid2 <- snpid2[is.na(match.snpid) | match.snpid > 1:length(match.snpid)]
                rm(match.snpid)
            }
            if(!is.null(group.info)) {
                if(auto.flip) {
                    tmp <- subset(tmp, !snpid %in% group.info$snpid & !snpid2 %in% group.info$snpid)
                    rm(snpid2)
                } else tmp <- subset(tmp, !snpid %in% group.info$snpid)
            } else if(auto.flip) {
                rm(snpid2)
            }
            if(nrow(tmp) > 0) group.info <- rbind(group.info, tmp)
            rm(tmp)
        }
    }
    group.info <- group.info[order(group.info$chr, group.info$pos), ]
    p <- nrow(group.info)
    group.info$idx <- 1:p
    variant.id1 <- group.info$snpid
    scores <- cons <- vector("list", n.cohort)
    N.resampling <- rep(0, n.cohort)
    if(auto.flip) {
        cat("Automatic allele flipping enabled...\nVariants matching alt/ref but not ref/alt alleles will also be included, with flipped effects\n")
        variant.id2 <- paste(group.info$chr, group.info$pos, group.info$alt, group.info$ref, sep = ":")
    }
    for(i in 1:n.cohort) {
        tmp.scores <- NULL
        for(j in 1:n.files[i]) {
            tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
            if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
            }
            tmp <- tmp[,c("chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")]
            tmp$variant.idx <- 1:nrow(tmp)
            variant.id <- paste(tmp$chr, tmp$pos, tmp$ref, tmp$alt, sep = ":")
            variant.idx1 <- variant.id %in% variant.id1
            if(auto.flip) {
                variant.idx2 <- variant.id %in% variant.id2
                if(any(variant.idx1 & variant.idx2)) {
                    tmp.dups <- which(variant.idx1 & variant.idx2)
                    cat("The following ambiguous variants were found in",paste0(meta.files.prefix[i], ".sample.", j),":\n")
                    cat("chr:", tmp$chr[tmp.dups], "\n")
                    cat("pos:", tmp$pos[tmp.dups], "\n")
                    cat("ref:", tmp$ref[tmp.dups], "\n")
                    cat("alt:", tmp$alt[tmp.dups], "\n")
                    cat("Warning: both variants with alleles ref/alt and alt/ref were present at the same position and coding should be double checked!\nFor these variants, only those with alleles ref/alt were used in the analysis...\n")
                    variant.idx2[tmp.dups] <- FALSE
                    rm(tmp.dups)
                }
                tmpallele <- tmp$ref[variant.idx2]
                tmp$ref[variant.idx2] <- tmp$alt[variant.idx2]
                tmp$alt[variant.idx2] <- tmpallele
                rm(tmpallele)
                tmp$altfreq[variant.idx2] <- 1-tmp$altfreq[variant.idx2]
                tmp$flip <- 1-2*(variant.idx2)
                tmp <- tmp[variant.idx1 | variant.idx2, ]
            } else {
                tmp$flip <- 1
                tmp <- tmp[variant.idx1, ]
            }
            tmp$file <- j
            tmp.scores <- rbind(tmp.scores, tmp)
            rm(tmp)
        }
        scores[[i]] <- cbind(group.info, tmp.scores[match(variant.id1, paste(tmp.scores$chr, tmp.scores$pos, tmp.scores$ref, tmp.scores$alt, sep = ":")), c("N", "missrate", "altfreq", "SCORE", "file", "variant.idx", "flip")])
        rm(tmp.scores)
        cons[[i]] <- file(paste0(meta.files.prefix[i], ".resample.1"), "rb")
        N.resampling[i] <- readBin(cons[[i]], what = "integer", n = 1, size = 4)
    }
    current.lines <- current.cons <- rep(1, n.cohort)
    nbatch.flush <- (p-1) %/% nperbatch + 1
    all.out <- NULL
    for(i in 1:nbatch.flush) {
        tmp.idx <- if(i == nbatch.flush) group.info$idx[((i-1)*nperbatch+1):p] else group.info$idx[((i-1)*nperbatch+1):(i*nperbatch)]
        #tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
        U.list <- V.list <- vector("list", n.cohort)
        variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
        for(j in 1:n.cohort) {
            tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
            if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
                U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
                U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$flip
                tmp.V <- matrix(NA, sum(tmp.include), N.resampling[j])
                for(ij in 1:sum(tmp.include)) {
                    if(U.list[[j]]$file[ij]!=current.cons[j]) {
                        close(cons[[j]])
                        current.cons[j] <- U.list[[j]]$file[ij]
                        cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.", current.cons[j]), "rb")
                        tmp.N.resampling <- readBin(cons[[j]], what = "integer", n = 1, size = 4)
                        if(tmp.N.resampling != N.resampling[j]) stop(paste0("Error: N.resampling in ", meta.files.prefix[j], ".resample.", current.cons[j], " does not match that in ",meta.files.prefix[j], ".resample.1"))
                        current.lines[j] <- 1
                    }
                    if(U.list[[j]]$variant.idx[ij]!=current.lines[j]) seek(cons[[j]], where = 4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-current.lines[j]), origin = "current", rw = "read")
                    tmp.V[ij,] <- readBin(cons[[j]], what = "numeric", n = N.resampling[j], size = 4)
                    current.lines[j] <- U.list[[j]]$variant.idx[ij]+1
                }
                V.list[[j]] <- tmp.V * U.list[[j]]$flip / sqrt(N.resampling[j])
                rm(tmp.V)
                variant.indices <- c(variant.indices, U.list[[j]]$idx)
                tmp.N <- c(tmp.N, U.list[[j]]$N)
                tmp.Nmiss <- c(tmp.Nmiss, U.list[[j]]$N * U.list[[j]]$missrate/(1-U.list[[j]]$missrate))
                tmp.AC <- c(tmp.AC, 2*U.list[[j]]$N*U.list[[j]]$altfreq)
            }
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
        out <- data.frame(SNP = group.info$SNP[variant.indices], chr = group.info$chr[variant.indices], pos = group.info$pos[variant.indices], ref = group.info$ref[variant.indices], alt = group.info$alt[variant.indices], N = N[include], missrate = (Nmiss/(N+Nmiss))[include], altfreq = AF[include])
        n.p <- length(variant.indices)
        U <- rep(0, n.p)
        V <- matrix(0, n.p, sum(N.resampling))
        for(j in 1:n.cohort) {
            if(!is.null(U.list[[j]]) & !is.null(V.list[[j]])) {
                IDX <- match(U.list[[j]]$idx, variant.indices)
                if(sum(!is.na(IDX)) == 0) next
                IDX2 <- which(!is.na(IDX))
                IDX <- IDX[IDX2]
                U[IDX] <- U[IDX]+U.list[[j]]$SCORE[IDX2]
                V[IDX, (sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])] <- V[IDX,(sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])]+V.list[[j]][IDX2,]
            }
        }
        out$SCORE <- U
        out$VAR <- rowSums(V^2)
        out$PVAL <- pchisq(out$SCORE^2/out$VAR, 1, lower = FALSE)
        all.out <- rbind(all.out, out)
    }
    for(i in 1:n.cohort) close(cons[[i]])
    return(all.out)
}


#' Variant Set Mixed Model Association Tests (StocSum)
#' @description Variant Set Mixed Model Association Tests (the burden test, SKAT, SKAT-O, and an efficient hybrid test to combine the burden test and SKAT) for multiple user-defined test units and a null generalized linear mixed model.
#' @param meta.file.prefix a character vector for prefix of intermediate files (*.sample.* and *.resample.*). 
#' @param n.files an integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix))).
#' @param cohort.group.idx a vector with the same length as meta.files.prefix, indicating which cohorts should be grouped together in the meta-analysis assuming homogeneous genetic effects. For example, c("a","b","a","a","b") means cohorts 1, 3, 4 are assumed to have homogeneous genetic effects, and cohorts 2, 5 are in another group with homogeneous genetic effects (but possibly heterogeneous with group "a"). If NULL, all cohorts are in the same group (default = NULL).
#' @param group.file a plain text file with 6 columns defining the test units. There should be no headers in the file, and the columns are group name, chromosome, position, reference allele, alternative allele and weight, respectively.
#' @param group.file.sep the delimiter in group.file (default = "\\t").
#' @param MAF.range a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param MAF.weight.beta a numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies. This internal minor allele frequency weight is multiplied by the external weight given by the group.file. To turn off internal minor allele frequency weight and only use the external weight given by the group.file, use c(1, 1) to assign flat weights (default = c(1, 25)). Applied to the combined samples.
#' @param missing.cutoff the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param method a method to compute p-values for SKAT-type test statistics (default = "davies"). "davies" represents an exact method that computes a p-value by inverting the characteristic function of the mixture chisq distribution, with an accuracy of 1e-6. When "davies" p-value is less than 1e-5, it defaults to method "kuonen". "kuonen" represents a saddlepoint approximation method that computes the tail probabilities of the mixture chisq distribution. When "kuonen" fails to compute a p-value, it defaults to method "liu". "liu" is a moment-matching approximation method for the mixture chisq distribution.
#' @param tests a character vector indicating which tests should be performed ("B" for the burden test, "S" for SKAT, "O" for SKAT-O and "E" for the efficient hybrid test of the burden test and SKAT). The burden test and SKAT are automatically included when performing "O", and the burden test is automatically included when performing "E" (default = "E").
#' @param rho a numeric vector defining the search grid used in SKAT-O (see the SKAT-O paper for details). Not used for the burden test, SKAT or the efficient hybrid test of the burden test and SKAT (default = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1)).
#' @param use.minor.allele a logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SKAT results, but Burden (as well as SKAT-O and hybrid test to combine the burden test and SKAT) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction. Use with caution, as major/minor alleles may flip in different cohorts. In that case, minor allele will be determined based on the allele frequency in the combined samples.
#' @param auto.flip a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
#' @return \code{G.pval} returns a data frame with the following components:
#' \item{groups}{name of the test unit group.}
#' \item{n.variants}{number of variants in the test unit group that pass the missing rate and allele frequency filters.}
#' \item{B.score}{burden test score statistic (avaiable if test is Burden, SKATO, or SMMAT).}
#' \item{B.var}{variance of burden test score statistic (avaiable if test is Burden, SKATO, or SMMAT).}
#' \item{B.pval}{burden test p-value (avaiable if test is Burden, SKATO, or SMMAT).}
#' \item{S.pval}{SKAT p-value (avaiable if test is SKATO or SKAT).}
#' \item{O.pval}{SKAT-O p-value (avaiable if test is SKATO).}
#' \item{E.pval}{efficient hybrid test of the burden test and SKAT p-value (avaiable if test is SMMAT).}
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{glmmkin2randomvec}, \code{G.stat}
#' @examples
#' \donttest{
#' library(GMMAT)
#' data(example)
#' attach(example)
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
#' G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1)
#' group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
#' out2 <- G.pval(out.prefix, group.file = group.file, MAF.range=c(0,0.5), miss.cutoff = 1)
#' }
#' @keywords variant set-based test
#' @export

G.prep <- function(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), cohort.group.idx = NULL, group.file, group.file.sep = "\t", auto.flip = FALSE)
{
    if(.Platform$endian!="little") stop("Error: platform must be little endian.")
    n.cohort <- length(meta.files.prefix)
    if(length(n.files) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and n.files do not match.")
    if(!is.null(cohort.group.idx)) {
        if(length(cohort.group.idx) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and cohort.group.idx do not match.")
        cohort.group.idx <- as.numeric(factor(cohort.group.idx))
        n.cohort.groups <- length(unique(cohort.group.idx))
    }
    # group.info <- try(read.table(group.file, header = FALSE, stringsAsFactors = FALSE, sep = group.file.sep), silent = TRUE)
    group.info<- try(fread(group.file,header=FALSE,stringsAsFactors = FALSE, sep=group.file.sep, data.table = FALSE), silent = TRUE)
    if (class(group.info) == "try-error") {
        stop("Error: cannot read group.file!")
    }
    colnames(group.info) <- c("group", "chr", "pos", "ref", "alt", "weight")
    group.info <- group.info[!duplicated(paste(group.info$group, group.info$chr, group.info$pos, group.info$ref, group.info$alt, sep = ":")), ]
    group.info <- group.info[order(group.info$chr, group.info$pos), ]
    groups <- unique(group.info$group)
    n.groups <- length(groups)
    group.info$group.idx <- match(group.info$group, groups)
    group.info <- group.info[order(group.info$group.idx), ]
    group.info$idx <- 1:nrow(group.info)
    group.idx.end <- findInterval(1:n.groups, group.info$group.idx)
    group.idx.start <- c(1, group.idx.end[-n.groups] + 1)
    variant.id1 <- paste(group.info$chr, group.info$pos, group.info$ref, group.info$alt, sep = ":")
    scores <- cons <- vector("list", n.cohort)
    N.resampling <- rep(0, n.cohort)
    if(auto.flip) {
        cat("Automatic allele flipping enabled...\nVariants matching alt/ref but not ref/alt alleles will also be included, with flipped effects\n")
        variant.id2 <- paste(group.info$chr, group.info$pos, group.info$alt, group.info$ref, sep = ":")
    }
    for(i in 1:n.cohort) {
        tmp.scores <- NULL
        for(j in 1:n.files[i]) {
            # tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
            tmp <- try(fread(paste0(meta.files.prefix[i], ".sample.", j), header=TRUE, data.table = FALSE))
            if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
            }
            tmp <- tmp[,c("chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")]
            tmp$variant.idx <- 1:nrow(tmp)
            variant.id <- paste(tmp$chr, tmp$pos, tmp$ref, tmp$alt, sep = ":")
            variant.idx1 <- variant.id %in% variant.id1
            if(auto.flip) {
                variant.idx2 <- variant.id %in% variant.id2
                if(any(variant.idx1 & variant.idx2)) {
                    tmp.dups <- which(variant.idx1 & variant.idx2)
                    cat("The following ambiguous variants were found in",paste0(meta.files.prefix[i], ".sample.", j),":\n")
                    cat("chr:", tmp$chr[tmp.dups], "\n")
                    cat("pos:", tmp$pos[tmp.dups], "\n")
                    cat("ref:", tmp$ref[tmp.dups], "\n")
                    cat("alt:", tmp$alt[tmp.dups], "\n")
                    cat("Warning: both variants with alleles ref/alt and alt/ref were present at the same position and coding should be double checked!\nFor these variants, only those with alleles ref/alt were used in the analysis...\n")
                    variant.idx2[tmp.dups] <- FALSE
                    rm(tmp.dups)
                }
                tmpallele <- tmp$ref[variant.idx2]
                tmp$ref[variant.idx2] <- tmp$alt[variant.idx2]
                tmp$alt[variant.idx2] <- tmpallele
                rm(tmpallele)
                tmp$altfreq[variant.idx2] <- 1-tmp$altfreq[variant.idx2]
                tmp$flip <- 1-2*(variant.idx2)
                tmp <- tmp[variant.idx1 | variant.idx2, ]
            } else {
                tmp$flip <- 1
                tmp <- tmp[variant.idx1, ]
            }
            if (nrow(tmp)==0) next
            tmp$file <- j
            tmp.scores <- rbind(tmp.scores, tmp)
            rm(tmp)
        }
        scores[[i]] <- cbind(group.info, tmp.scores[match(variant.id1, paste(tmp.scores$chr, tmp.scores$pos, tmp.scores$ref, tmp.scores$alt, sep = ":")), c("N", "missrate", "altfreq", "SCORE", "file", "variant.idx", "flip")])
        rm(tmp.scores)
    }
    if(!is.null(cohort.group.idx)) {
        out <- list(meta.files.prefix = meta.files.prefix, n.cohort = n.cohort, cohort.group.idx = cohort.group.idx, n.cohort.groups = n.cohort.groups, n.groups = n.groups, groups = groups, group.idx.start = group.idx.start, group.idx.end = group.idx.end, scores = scores, auto.filp = auto.flip)
    } else{
        out <- list(meta.files.prefix = meta.files.prefix, n.cohort = n.cohort, cohort.group.idx = cohort.group.idx, n.groups = n.groups, groups = groups, group.idx.start = group.idx.start, group.idx.end = group.idx.end, scores = scores, auto.filp = auto.flip)
    }
    class(out) <- "G.prep"
    return(out)
}

G.pval <- function(G.prep.obj, MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25), miss.cutoff = 1, method = "davies", tests = "E", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE)
{
    if(class(G.prep.obj) != "G.prep") stop("Error: G.prep.obj must be a class G.prep object!")
    if(.Platform$endian!="little") stop("Error: platform must be little endian.")
    meta.files.prefix <- G.prep.obj$meta.files.prefix
    n.groups <- G.prep.obj$n.groups
    n.cohort <- G.prep.obj$n.cohort
    cohort.group.idx <- G.prep.obj$cohort.group.idx 
    if(!is.null(cohort.group.idx)) n.cohort.groups <- G.prep.obj$n.cohort.groups
    n.cohort.groups <- G.prep.obj$n.cohort.groups
    groups <- G.prep.obj$groups
    group.idx.start <- G.prep.obj$group.idx.start
    group.idx.end <- G.prep.obj$group.idx.end
    scores <- G.prep.obj$scores
    # N.resampling <- G.prep.obj$N.resampling
    # cons <- G.prep.obj$cons
    auto.flip <- G.prep.obj$auto.flip
    if(any(!tests %in% c("B", "S", "O", "E"))) stop("Error: \"tests\" should only include \"B\" for the burden test, \"S\" for SKAT, \"O\" for SKAT-O or \"E\" for the efficient hybrid test of the burden test and SKAT.")
    Burden <- "B" %in% tests
    SKAT <- "S" %in% tests
    SKATO <- "O" %in% tests
    SMMAT <- "E" %in% tests
    cons <- vector("list", n.cohort)
    N.resampling <- rep(0, n.cohort)
    # if(auto.flip) {
    #     cat("Automatic allele flipping enabled...\nVariants matching alt/ref but not ref/alt alleles will also be included, with flipped effects\n")
    #     variant.id2 <- paste(group.info$chr, group.info$pos, group.info$alt, group.info$ref, sep = ":")
    # }
    for(i in 1:n.cohort) {
        cons[[i]] <- file(paste0(meta.files.prefix[i], ".resample.1"), "rb")
        N.resampling[i] <- readBin(cons[[i]], what = "integer", n = 1, size = 4)
    }
    n.variants <- rep(0,n.groups)
    if(Burden | SKATO | SMMAT) {
        Burden.score <- rep(NA, n.groups)
        Burden.var <- rep(NA, n.groups)
        Burden.pval <- rep(NA, n.groups)
    }
    if(SKAT | SKATO) SKAT.pval <- rep(NA, n.groups)
    if(SKATO) {
        SKATO.pval <- rep(NA, n.groups)
        SKATO.minp <- rep(NA, n.groups)
        SKATO.minp.rho <- rep(NA, n.groups)
    }
    if(SMMAT) SMMAT.pval <- rep(NA, n.groups)
    current.lines <- current.cons <- rep(1, n.cohort)
    for(i in 1:n.groups) {
        tmp.idx <- group.idx.start[i]:group.idx.end[i]
        #tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
        U.list <- V.list <- vector("list", n.cohort)
        variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
        for(j in 1:n.cohort) {
            tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
            if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
                U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
                U.list[[j]]$weight <- U.list[[j]]$weight * U.list[[j]]$flip
                U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$weight
                tmp.V <- matrix(NA, sum(tmp.include), N.resampling[j])
                for(ij in 1:sum(tmp.include)) {
                    if(U.list[[j]]$file[ij]!=current.cons[j]) {
                        close(cons[[j]])
                        current.cons[j] <- U.list[[j]]$file[ij]
                        cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.", current.cons[j]), "rb")
                        tmp.N.resampling <- readBin(cons[[j]], what = "integer", n = 1, size = 4)
                        if(tmp.N.resampling != N.resampling[j]) stop(paste0("Error: N.resampling in ", meta.files.prefix[j], ".resample.", current.cons[j], " does not match that in ",meta.files.prefix[j], ".resample.1"))
                        current.lines[j] <- 1
                    }
                    if(U.list[[j]]$variant.idx[ij]!=current.lines[j]) seek(cons[[j]], where = 4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-current.lines[j]), origin = "current", rw = "read")
                    tmp.V[ij,] <- readBin(cons[[j]], what = "numeric", n = N.resampling[j], size = 4)
                    current.lines[j] <- U.list[[j]]$variant.idx[ij]+1
                }
                V.list[[j]] <- tmp.V * U.list[[j]]$weight / sqrt(N.resampling[j])
                rm(tmp.V)
                variant.indices <- c(variant.indices, U.list[[j]]$idx)
                tmp.N <- c(tmp.N, U.list[[j]]$N)
                tmp.Nmiss <- c(tmp.Nmiss, U.list[[j]]$N * U.list[[j]]$missrate/(1-U.list[[j]]$missrate))
                tmp.AC <- c(tmp.AC, 2*U.list[[j]]$N*U.list[[j]]$altfreq)
            }
        }
        if(length(variant.indices) == 0) next
        tmp.variant.indices <- variant.indices
        variant.indices <- sort(unique(variant.indices))
        N <- sapply(variant.indices, function(x) sum(tmp.N[tmp.variant.indices==x]))
        Nmiss <- sapply(variant.indices, function(x) sum(tmp.Nmiss[tmp.variant.indices==x]))
        AF <- sapply(variant.indices, function(x) sum(tmp.AC[tmp.variant.indices==x]))/2/N
        include <- (Nmiss/(N+Nmiss) <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
        rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices, N, Nmiss)
        if(sum(include) == 0) next
        variant.indices <- variant.indices[include]
        n.p <- length(variant.indices)
        n.variants[i] <- n.p
        U <- if(!is.null(cohort.group.idx)) rep(0, n.cohort.groups*n.p) else rep(0, n.p)
        V <- if(!is.null(cohort.group.idx)) matrix(0, n.cohort.groups*n.p, sum(N.resampling)) else matrix(0, n.p, sum(N.resampling))
        for(j in 1:n.cohort) {
            if(!is.null(U.list[[j]]) & !is.null(V.list[[j]])) {
                IDX <- match(U.list[[j]]$idx, variant.indices)
                if(sum(!is.na(IDX)) == 0) next
                IDX2 <- which(!is.na(IDX))
                IDX <- IDX[IDX2]
                if(!is.null(cohort.group.idx)) IDX <- IDX+n.p*(cohort.group.idx[j]-1)
                U[IDX] <- U[IDX]+U.list[[j]]$SCORE[IDX2]
                V[IDX, (sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])] <- V[IDX,(sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])]+V.list[[j]][IDX2,]
            }
        }
        tmp.weight <- MAF.weights.beta.fun(AF[include], MAF.weights.beta[1], MAF.weights.beta[2])
        if(use.minor.allele) tmp.weight[AF[include] > 0.5] <- -tmp.weight[AF[include] > 0.5]
        if(!is.null(cohort.group.idx)) tmp.weight <- rep(tmp.weight, n.cohort.groups)
        U <- U*tmp.weight
        V <- V*tmp.weight
        if(max(diff(V)) < sqrt(.Machine$double.eps)) {
            burden.score <- sum(U)
            burden.var <- sum(colSums(V)^2)
            burden.pval <- pchisq(burden.score^2/burden.var, df=1, lower.tail=FALSE)
            if(Burden | SKATO | SMMAT) {
                Burden.score[i] <- burden.score
                Burden.var[i] <- burden.var
                Burden.pval[i] <- burden.pval
            }
            if(SKAT | SKATO) SKAT.pval[i] <- burden.pval
            if(SKATO) {
                SKATO.pval[i] <- burden.pval
                SKATO.minp[i] <- burden.pval
                SKATO.minp.rho[i] <- 1
            }
            if(SMMAT) SMMAT.pval[i] <- burden.pval
        } else {
            if(SKATO) {
                re <- try(.skato_Vpval(U = U, V = V, rho = rho, method = method))
                if(class(re)[1] != "try-error") {
                    Burden.score[i] <- re$Burden.score
                    Burden.var[i] <- re$Burden.var
                    Burden.pval[i] <- re$Burden.pval
                    SKAT.pval[i] <- re$SKAT.pval
                    SKATO.pval[i] <-re$p
                    SKATO.minp[i] <- re$minp
                    SKATO.minp.rho[i] <- re$minp.rho
                }
            } else {
                if(SKAT) SKAT.pval[i] <- tryCatch(.quad_Vpval(U = U, V = V, method = method), error = function(e) { NA })
                if(Burden | SMMAT) {
                    Burden.score[i] <- sum(U)
                    Burden.var[i] <- sum(colSums(V)^2)
                    Burden.pval[i] <- pchisq(Burden.score[i]^2/Burden.var[i], df=1, lower.tail=FALSE)
                }
            }
            if(SMMAT) {
                Vsum <- colSums(V)
                U <- U - crossprod(t(V),Vsum) * Burden.score[i] / Burden.var[i]
                V <- V - crossprod(crossprod(Vsum,t(V)),t(Vsum)) / Burden.var[i]
                if(mean(abs(V)) < sqrt(.Machine$double.eps)) SMMAT.pval[i] <- Burden.pval[i]
                else SMMAT.pval[i] <- tryCatch(pchisq(-2*log(Burden.pval[i])-2*log(.quad_Vpval(U = U, V = V, method = method)), df = 4, lower.tail = FALSE), error = function(e) { Burden.pval[i] })
            }
        }
    }
    for(i in 1:n.cohort) close(cons[[i]])
    out <- data.frame(group=groups, n.variants=n.variants)
    if(Burden | SKATO | SMMAT) {
        out$B.score <- Burden.score
        out$B.var <- Burden.var
        out$B.pval <- Burden.pval
    }
    if(SKAT | SKATO) out$S.pval <- SKAT.pval
    if(SKATO) {
        out$O.pval <- SKATO.pval
        out$O.minp <- SKATO.minp
        out$O.minp.rho <- SKATO.minp.rho
    }
    if(SMMAT) out$E.pval <- SMMAT.pval
    return(out)
}

Cond.svt.pval <- function(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), tagChr, StartPos, EndPos, tagPos, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, auto.flip = FALSE, nperbatch = 10000, tol=1e-5)
{
    if(.Platform$endian!="little") stop("Error: platform must be little endian.")
    n.cohort <- length(meta.files.prefix)
    if(length(n.files) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and n.files do not match.")
    group.info <- NULL
    for(i in 1:n.cohort) {
        for(j in 1:n.files[i]) {
            # tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
            tmp <- try(fread(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, data.table = FALSE))
            if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
            }
            tmp <- tmp[,c("SNP", "chr", "pos", "ref", "alt")]
            tmp$snpid <- paste(tmp$chr, tmp$pos, tmp$ref, tmp$alt, sep = ":")
            tmp <- tmp[!duplicated(tmp$snpid), , drop = FALSE]
            if(auto.flip) {
                snpid2 <- paste(tmp$chr, tmp$pos, tmp$alt, tmp$ref, sep = ":")
                match.snpid <- match(snpid2, tmp$snpid)
                tmp <- tmp[is.na(match.snpid) | match.snpid > 1:length(match.snpid), , drop = FALSE]
                snpid2 <- snpid2[is.na(match.snpid) | match.snpid > 1:length(match.snpid)]
                rm(match.snpid)
            }
            if(!is.null(group.info)) {
                if(auto.flip) {
                    tmp <- subset(tmp, !snpid %in% group.info$snpid & !snpid2 %in% group.info$snpid)
                    rm(snpid2)
                } else tmp <- subset(tmp, !snpid %in% group.info$snpid)
            } else if(auto.flip) {
                rm(snpid2)
            }
            if(nrow(tmp) > 0) group.info <- rbind(group.info, tmp)
            rm(tmp)
        }
    }
    group.info <- group.info[order(group.info$chr, group.info$pos), ]
    group.info <- group.info[(group.info$chr == tagChr) & (group.info$pos > StartPos) & (group.info$pos < EndPos),]
    p <- nrow(group.info)
    group.info$idx <- 1:p
    variant.id1 <- group.info$snpid
    scores <- cons <- vector("list", n.cohort)
    N.resampling <- rep(0, n.cohort)
    if(auto.flip) {
        cat("Automatic allele flipping enabled...\nVariants matching alt/ref but not ref/alt alleles will also be included, with flipped effects\n")
        variant.id2 <- paste(group.info$chr, group.info$pos, group.info$alt, group.info$ref, sep = ":")
    }
    for(i in 1:n.cohort) {
        tmp.scores <- NULL
        for(j in 1:n.files[i]) {
            # tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
            tmp <- try(fread(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, data.table = FALSE))
            if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
            }
            tmp <- tmp[,c("chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")]
            tmp$variant.idx <- 1:nrow(tmp)
            variant.id <- paste(tmp$chr, tmp$pos, tmp$ref, tmp$alt, sep = ":")
            variant.idx1 <- variant.id %in% variant.id1
            if(auto.flip) {
                variant.idx2 <- variant.id %in% variant.id2
                if(any(variant.idx1 & variant.idx2)) {
                    tmp.dups <- which(variant.idx1 & variant.idx2)
                    cat("The following ambiguous variants were found in",paste0(meta.files.prefix[i], ".sample.", j),":\n")
                    cat("chr:", tmp$chr[tmp.dups], "\n")
                    cat("pos:", tmp$pos[tmp.dups], "\n")
                    cat("ref:", tmp$ref[tmp.dups], "\n")
                    cat("alt:", tmp$alt[tmp.dups], "\n")
                    cat("Warning: both variants with alleles ref/alt and alt/ref were present at the same position and coding should be double checked!\nFor these variants, only those with alleles ref/alt were used in the analysis...\n")
                    variant.idx2[tmp.dups] <- FALSE
                    rm(tmp.dups)
                }
                tmpallele <- tmp$ref[variant.idx2]
                tmp$ref[variant.idx2] <- tmp$alt[variant.idx2]
                tmp$alt[variant.idx2] <- tmpallele
                rm(tmpallele)
                tmp$altfreq[variant.idx2] <- 1-tmp$altfreq[variant.idx2]
                tmp$flip <- 1-2*(variant.idx2)
                tmp <- tmp[variant.idx1 | variant.idx2, ]
            } else {
                tmp$flip <- 1
                tmp <- tmp[variant.idx1, ]
            }
            if (nrow(tmp)==0) next
            tmp$file <- j
            tmp.scores <- rbind(tmp.scores, tmp)
            rm(tmp)
        }
        scores[[i]] <- cbind(group.info, tmp.scores[match(variant.id1, paste(tmp.scores$chr, tmp.scores$pos, tmp.scores$ref, tmp.scores$alt, sep = ":")), c("N", "missrate", "altfreq", "SCORE", "file", "variant.idx", "flip")])
        rm(tmp.scores)
        cons[[i]] <- file(paste0(meta.files.prefix[i], ".resample.1"), "rb")
        N.resampling[i] <- readBin(cons[[i]], what = "integer", n = 1, size = 4)
    }
    current.lines <- current.cons <- rep(1, n.cohort)

    tmp.idx <- group.info$idx[group.info$pos %in% tagPos]
    U.list <- V.list <- vector("list", n.cohort)
    variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
    for(j in 1:n.cohort) {
        tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
        if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
            U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
            U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$flip
            tmp.V <- matrix(NA, sum(tmp.include), N.resampling[j])
            for(ij in 1:sum(tmp.include)) {
                if(U.list[[j]]$file[ij]!=current.cons[j]) {
                    close(cons[[j]])
                    current.cons[j] <- U.list[[j]]$file[ij]
                    cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.", current.cons[j]), "rb")
                    tmp.N.resampling <- readBin(cons[[j]], what = "integer", n = 1, size = 4)
                    if(tmp.N.resampling != N.resampling[j]) stop(paste0("Error: N.resampling in ", meta.files.prefix[j], ".resample.", current.cons[j], " does not match that in ",meta.files.prefix[j], ".resample.1"))
                    current.lines[j] <- 1
                }
                if(U.list[[j]]$variant.idx[ij]!=current.lines[j]) seek(cons[[j]], where = 4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-current.lines[j]), origin = "current", rw = "read")
                tmp.V[ij,] <- readBin(cons[[j]], what = "numeric", n = N.resampling[j], size = 4)
                current.lines[j] <- U.list[[j]]$variant.idx[ij]+1
            }
            V.list[[j]] <- tmp.V * U.list[[j]]$flip / sqrt(N.resampling[j])
            rm(tmp.V)
            variant.indices <- c(variant.indices, U.list[[j]]$idx)
            tmp.N <- c(tmp.N, U.list[[j]]$N)
            tmp.Nmiss <- c(tmp.Nmiss, U.list[[j]]$N * U.list[[j]]$missrate/(1-U.list[[j]]$missrate))
            tmp.AC <- c(tmp.AC, 2*U.list[[j]]$N*U.list[[j]]$altfreq)
        }
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
    n.p <- length(variant.indices)
    Ub1 <- rep(0, n.p)
    Vb1 <- matrix(0, n.p, sum(N.resampling))
    for(j in 1:n.cohort) {
        if(!is.null(U.list[[j]]) & !is.null(V.list[[j]])) {
            IDX <- match(U.list[[j]]$idx, variant.indices)
            if(sum(!is.na(IDX)) == 0) next
            IDX2 <- which(!is.na(IDX))
            IDX <- IDX[IDX2]
            Ub1[IDX] <- Ub1[IDX]+U.list[[j]]$SCORE[IDX2]
            Vb1[IDX, (sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])] <- Vb1[IDX,(sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])]+V.list[[j]][IDX2,]
        }
    }
    for(i in 1:n.cohort) close(cons[[i]])

    current.lines <- current.cons <- rep(1, n.cohort)
    nbatch.flush <- (p-1) %/% nperbatch + 1
    all.out <- NULL
    for (j in 1:n.cohort) {
        cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.1"), "rb")
        N.resampling[j] <- readBin(cons[[j]], what = "integer", n = 1, size = 4)
    }
    for(i in 1:nbatch.flush) {
        # print(i)
        tmp.idx <- if(i == nbatch.flush) group.info$idx[((i-1)*nperbatch+1):p] else group.info$idx[((i-1)*nperbatch+1):(i*nperbatch)]
        #tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
        U.list <- V.list <- vector("list", n.cohort)
        variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
        for(j in 1:n.cohort) {
            # cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.1"), "rb")
            tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
            if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
                U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
                U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$flip
                tmp.V <- matrix(NA, sum(tmp.include), N.resampling[j])
                for(ij in 1:sum(tmp.include)) {
                    if(U.list[[j]]$file[ij]!=current.cons[j]) {
                        close(cons[[j]])
                        current.cons[j] <- U.list[[j]]$file[ij]
                        cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.", current.cons[j]), "rb")
                        tmp.N.resampling <- readBin(cons[[j]], what = "integer", n = 1, size = 4)
                        if(tmp.N.resampling != N.resampling[j]) stop(paste0("Error: N.resampling in ", meta.files.prefix[j], ".resample.", current.cons[j], " does not match that in ",meta.files.prefix[j], ".resample.1"))
                        current.lines[j] <- 1
                    }
                    # a <-seek(cons[[j]], where = 0, origin = "current", rw = "read")    # 
                    # print(a) #
                    if(U.list[[j]]$variant.idx[ij]!=current.lines[j]) seek(cons[[j]], where = 4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-current.lines[j]), origin = "current", rw = "read")  
                    tmp.V[ij,] <- readBin(cons[[j]], what = "numeric", n = N.resampling[j], size = 4)
                    # a <-seek(cons[[j]], where = 0, origin = "current", rw = "read")    # 
                    # print(a) # 
                    # print(U.list[[j]]$variant.idx[ij])
                    # pause(2)
                    current.lines[j] <- U.list[[j]]$variant.idx[ij]+1
                }
                V.list[[j]] <- tmp.V * U.list[[j]]$flip / sqrt(N.resampling[j])
                rm(tmp.V)
                variant.indices <- c(variant.indices, U.list[[j]]$idx)
                tmp.N <- c(tmp.N, U.list[[j]]$N)
                tmp.Nmiss <- c(tmp.Nmiss, U.list[[j]]$N * U.list[[j]]$missrate/(1-U.list[[j]]$missrate))
                tmp.AC <- c(tmp.AC, 2*U.list[[j]]$N*U.list[[j]]$altfreq)
            }
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
        out <- data.frame(SNP = group.info$SNP[variant.indices], chr = group.info$chr[variant.indices], pos = group.info$pos[variant.indices], ref = group.info$ref[variant.indices], alt = group.info$alt[variant.indices], N = N[include], missrate = (Nmiss/(N+Nmiss))[include], altfreq = AF[include])
        n.p <- length(variant.indices)
        U <- rep(0, n.p)
        V <- matrix(0, n.p, sum(N.resampling))
        for(j in 1:n.cohort) {
            if(!is.null(U.list[[j]]) & !is.null(V.list[[j]])) {
                IDX <- match(U.list[[j]]$idx, variant.indices)
                if(sum(!is.na(IDX)) == 0) next
                IDX2 <- which(!is.na(IDX))
                IDX <- IDX[IDX2]
                U[IDX] <- U[IDX]+U.list[[j]]$SCORE[IDX2]
                V[IDX, (sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])] <- V[IDX,(sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])]+V.list[[j]][IDX2,]
            }
        }
        Vb1_V<-tcrossprod(Vb1) #*_V variance matrix
        Vb1_V_i <- try(solve(Vb1_V), silent = TRUE)
        Tadj<-tcrossprod(tcrossprod(V,Vb1), Vb1_V_i)
        # Tadj<-(V %*% t(Vb1)) %*% Vb1_V_i
        U.adj<-U- Tadj %*% Ub1
        V.adj<- V - Tadj %*% Vb1
        # out$cor <- cor(V, Vb1)
        out$SCORE <- U.adj
        out$VAR <- rowSums(V.adj^2)
        # out$PVAL <- pchisq(out$SCORE^2/out$VAR, 1, lower = FALSE)
        out$PVAL <- ifelse(out$VAR>tol, pchisq(out$SCORE^2/out$VAR, 1, lower = FALSE), 1)
        cor<-rep(0, length(out$PVAL))
        for(iii in 1:length(out$PVAL)){
            tmpcor<-c()
            for (jjj in seq(1,nrow(Vb1))){
                vb1<-as.numeric(Vb1[jjj,])
                tmpcor<-c(tmpcor,cor(V[iii,],vb1))
            }
            cor[iii]<-max(abs(tmpcor))
        }
        out$cor<-cor
        all.out <- rbind(all.out, out)
    }
    for(i in 1:n.cohort) close(cons[[i]])
    return(all.out)
}

.skato_Vpval <- function(U, V, rho, method = "davies") {
    n.r <- length(rho)
    n.p <- length(U)
    lambdas <- vector("list", n.r)
    pval <- qval <- rep(NA, n.r)
    Q <- (1-rho)*sum(U^2)+rho*sum(U)^2
    Burden.score <- Burden.var <- Burden.pval <- SKAT.pval <- NA
    Vsum <- colSums(V)
    for(i in 1:n.r) {
        if(rho[i]==1) {
            Burden.score <- sum(U)
            Burden.var <- sum(Vsum^2)
            Burden.pval <- pchisq(Burden.score^2/Burden.var, df=1, lower.tail=FALSE)
            lambdas[[i]] <- Burden.var
            pval[i] <- Burden.pval
            next
        }
        if(rho[i]!=0) {
            R.M <- matrix(rho[i], n.p, n.p)
            diag(R.M) <- 1
            R.M.chol <- t(chol(R.M, pivot = TRUE))
            V.temp <- crossprod(R.M.chol, V)
        } else V.temp <- V
    V.temp.svd <- try(svd(V.temp))
        if(class(V.temp.svd)[1] == "try-error") {
            dim1 <- nrow(V.temp)
            dim2 <- ncol(V.temp)
            if(dim1 < dim2) lambda <- zapsmall(eigen(tcrossprod(V.temp), symmetric = T, only.values = T)$values
)
            else lambda <- zapsmall(eigen(crossprod(V.temp), symmetric = T, only.values = T)$values)
        } else lambda <- zapsmall(V.temp.svd$d^2)
        lambdas[[i]] <- lambda[lambda > 0]
        pval[i] <- .Q_pval(Q[i], lambdas[[i]], method = method)
        if(rho[i]==0) SKAT.pval <- pval[i]
    }
    minp <- min(pval)
    for(i in 1:n.r) {
        df <- sum(lambdas[[i]]^2)^2/sum(lambdas[[i]]^4)
        qval[i] <- (qchisq(minp, df, lower.tail = FALSE)-df)/sqrt(2*df)*sqrt(2*sum(lambdas[[i]]^2))+sum(lambdas[[i]])
    }
    ZMZ <- crossprod(crossprod(Vsum,t(V)),t(Vsum))/sum(Vsum^2)
    V.temp <- V - ZMZ
    V.temp.svd <- try(svd(V.temp))
    if(class(V.temp.svd)[1] == "try-error") {
            dim1 <- nrow(V.temp)
            dim2 <- ncol(V.temp)
            if(dim1 < dim2) lambda <- zapsmall(eigen(tcrossprod(V.temp), symmetric = T, only.values = T)$values)
            else lambda <- zapsmall(eigen(crossprod(V.temp), symmetric = T, only.values = T)$values)
        } else lambda <- zapsmall(V.temp.svd$d^2)
    lambda <- lambda[lambda > 0]
    muq <- sum(lambda)
    trace_ZMZ_Vtemp <- if(nrow(V.temp)>ncol(V.temp)) sum(crossprod(ZMZ, V.temp)^2) else sum(tcrossprod(ZMZ)*tcrossprod(V.temp))
    varq <- sum(lambda^2) * 2 + trace_ZMZ_Vtemp * 4
    df <- sum(lambda^2)^2/sum(lambda^4)
    tau <- rho * sum(Vsum^2) + sum(crossprod(Vsum,t(V))^2)/sum(Vsum^2) * (1 - rho)
    re <- tryCatch({
        integrate(function(x){
            t1 <- tau %x% t(x)
            re<-pchisq((apply((qval - t1)/(1-rho),2,min) - muq)/sqrt(varq)*sqrt(2*df) + df, df=df) * dchisq(x,df=1)
            return(re)
        }, lower = 0, upper = 40, subdivisions = 2000, abs.tol = 10^-25)
    }, error=function(e) NA)
    return(list(p = min(1-re[[1]], minp*n.r), minp = minp, minp.rho = rho[which.min(pval)], Burden.score=Burden.score, Burden.var=Burden.var, Burden.pval=Burden.pval, SKAT.pval=SKAT.pval))
}

.Q_pval <- function(Q, lambda, method = "davies") {
    if(method == "davies") {
        tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6))
        pval <- tmp$Qq
    if((tmp$ifault > 0) | (pval <= 1e-5) | (pval >= 1)) method <- "kuonen"
    }
    if(method == "kuonen") {
        pval <- .pKuonen(x = Q, lambda = lambda)
    if(is.na(pval)) method <- "liu"
    }
    if(method == "liu") pval <- CompQuadForm::liu(q = Q, lambda = lambda)
    return(pval)
}

.quad_Vpval <- function(U, V, method = "davies") {
    Q <- sum(U^2)
    V.svd <- try(svd(V))
    if(class(V.svd)[1] == "try-error") {
        dim1 <- nrow(V)
        dim2 <- ncol(V)
        if(dim1 < dim2) lambda <- zapsmall(eigen(tcrossprod(V), symmetric = T, only.values = T)$values)
        else lambda <- zapsmall(eigen(crossprod(V), symmetric = T, only.values = T)$values)
    } else lambda <- zapsmall(V.svd$d^2)
    lambda <- lambda[lambda > 0]
    pval <- .Q_pval(Q, lambda, method = method)
    return(pval)
}

.pKuonen<-function (x, lambda, delta = rep(0, length(lambda)), df = rep(1, length(lambda)))
{
    delta <- delta[lambda != 0]
    df <- df[lambda != 0]
    lambda <- lambda[lambda != 0]
    if(length(lambda) != length(delta)) stop("Error: inconsistent length in lambda and delta!")
    if(length(lambda) != length(df)) stop("Error: inconsistent length in lambda and df!")
    if (length(lambda) == 1) {
        pchisq(x/lambda, df = df, ncp = delta, lower.tail = FALSE)
    }
    d <- max(lambda)
    lambda <- lambda/d
    x <- x/d
    k0 <- function(zeta) {
        -sum(df * log(1 - 2 * zeta * lambda))/2 + sum((delta * lambda *
            zeta)/(1 - 2 * zeta * lambda))
    }
    kprime0 <- function(zeta) {
        sapply(zeta, function(zz) {
            sum(((delta + df) * lambda)/(1 - 2 * zz * lambda) + 2 * (delta *
                zz * lambda^2)/(1 - 2 * zz * lambda)^2)
        })
    }
    kpprime0 <- function(zeta) {
        sum((2 * (2 * delta + df) * lambda^2)/(1 - 2 * zeta * lambda)^2 + 8 *
            delta * zeta * lambda^3/(1 - 2 * zeta * lambda)^3)
    }
    if (any(lambda < 0)) {
        lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
    }
    else if (x > sum((df+delta)*lambda)) {
        lmin <- -0.01
    }
    else {
        lmin <- -length(lambda)*max(df+delta)/(2 * x)
    }
    lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
    hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
        upper = lmax, tol = 1e-08)$root
    w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v <- hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04)
        NA
    else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}

MAF.weights.beta.fun <- function(freq, beta1, beta2) {
    freq[freq > 0.5] <- 1 - freq[freq > 0.5]
    ifelse(freq <= 0, 0, dbeta(freq, beta1, beta2))
}
