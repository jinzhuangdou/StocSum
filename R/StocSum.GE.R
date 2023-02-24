#' Create random vectors for a glmmkin object
#' @description Generate random vectors from multivariate normal distribution with mean 0 and covariance matrix P.
#' @param obj A glmmkin object.
#' @param Z A list of design matrices for the random effects. The length must match the number of variance components.
#' @param N.randomvec The number of random vectors to generate (default = 1000).
#' @param group.idx A length N index vector showing which observation belongs to which variance group, for heteroscedastic linear mixed models (default = NULL for homoscedastic linear mixed models).
#' @param cluster.idx A length N index vector showing which observation belongs to which cluster (default = NULL for no clusters).
#' @param robust A logical switch: whether robust variance should be used (default = FALSE).
#' @return A list of class glmmkin.randomvec
#' \item{theta}{inherited from the glmmkin object. A vector or a list of variance component parameter estimates.}
#' \item{scaled.residuals}{inherited from the glmmkin object. A vector or a matrix for the scaled residuals, calculated as the original residuals divided by the dispersion parameter (in heteroscedastic linear mixed models, corresponding residual variance estimates by each group).}
#' {random.vectors}{A random matrix with dimensions equal to the sample size multiplied by \code{N.randomvec}.}
#' \item{X}{inherited from the glmmkin object. Model matrix for the fixed effects.}
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
#' if(!is.null(nullmod$P)){
#'   obj <- glmmkin2randomvec(nullmod)
#' }else{
#'   kinship.chol <- chol(GRM)
#'   obj<-glmmkin2randomvec(nullmod, Z = list(t(kinship.chol)))
#' }
#' }
#' @keywords random vector
#' @export
glmmkin2randomvec <- function(obj, Z = NULL, N.randomvec = 1000, group.idx=NULL,cluster.idx=NULL,robust = FALSE) 
{
    if(class(obj) != "glmmkin") stop("Error: \"obj\" must be a class glmmkin object.")
    N <- length(obj$id_include)
    random.vectors <- matrix(rnorm(N*N.randomvec),nrow=N,ncol=N.randomvec)
    if(!is.null(obj$P) && !robust) {
        eig <- eigen(obj$P, symmetric = TRUE)
        random.vectors <- tcrossprod(eig$vectors, t(random.vectors * sqrt(pmax(eig$values, 0))))
        rm(eig)
    } else {
        if(obj$n.groups != 1 && (is.null(group.idx) || !all.equal(seq_len(obj$n.groups), sort(unique(group.idx))))) stop("Error: heteroscedastic linear mixed models should include a valid group.idx argument.")
        if(is.null(group.idx)) group.idx <- rep(1, N)
        if(!robust) random.vectors <- sqrt(obj$theta[group.idx]) * random.vectors
        else {
            res <- as.numeric(obj$Y - tcrossprod(obj$X, t(obj$coefficient)))
            if(is.null(cluster.idx)) random.vectors <- random.vectors * res
            else random.vectors <- random.vectors[match(cluster.idx,unique(cluster.idx)),] * res
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
    out <- list(theta = obj$theta, scaled.residuals = obj$scaled.residuals, random.vectors = as.matrix(random.vectors), id_include = obj$id_include, X=obj$X)
    class(out) <- "glmmkin.randomvec"
    return(out)
}


#' Calculate summary statistics and  stochastic statistics for gene-environment tests
#' @description Calculate summary statistics and stochastic statistics in the context of gene-environment tests.
#' @param null.obj A class glmmkin.randomvec object, returned by generating random vectors using \code{glmmkin2randomvec}.
#' @param interaction A numeric or a character vector indicating the environmental factors. If a numberic vector, it represents which indices in the order of covariates are the environmental factors; if a character vector, it represents the variable names of the environmental factors.
#' @param interaction.covariates A numeric or a character vector indicating the interaction covariates. If a numeric vector, it represents which indices in the order of covariates are the interaction covariates; if a character vector, it represents the variable names of the interaction covariates.
#' @param geno.file The .gds file name or an object of class SeqVarGDSClass for the full genotypes. The \code{sample.id} in \code{geno.file} should overlap \code{id_include} in \code{null.obj}. It is recommended that \code{sample.id} in \code{geno.file} include the full samples (at least all samples as specified in \code{id_include} of \code{null.obj}). It is not necessary for the user to take a subset of \code{geno.file} before running the analysis. If \code{geno.file} is an object of class SeqVarGDSClass, the .gds file will be closed upon successful completion of the function.
#' @param meta.file.prefix Prefix of intermediate files (*.sample.1 and *.resample.1).
#' @param MAF.range A numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param miss.cutoff The maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param missing.method Method of handling missing genotypes. Either "impute2mean" or "impute2zero" (default = "impute2mean").
#' @param nperbatch An integer for how many SNPs to be included in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @param ncores A positive integer indicating the number of cores to be used in parallel computing (default = 1).
#' @return NULL. \code{GE.stat} will store the summary statistics and the stochastic statistics in two files with the prefix specified by user.
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{glmmkin2randomvec}, \code{GE.svt.pval}
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
#' out.prefix <- "test.GE"
#' gdsfile <- system.file("extdata", "geno.gds", package = "StocSum")
#' GE.stat(obj, interaction='sex', geno.file = gdsfile, meta.file.prefix = out.prefix)
#' }
#' @keywords summary statistics, gene-environment tests
#' @export
GE.stat <- function(null.obj, interaction, interaction.covariates = NULL, geno.file, meta.file.prefix, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, missing.method = "impute2mean", nperbatch = 10000, ncores = 1)
{
    if(Sys.info()["sysname"] == "Windows" && ncores > 1) {
        warning("The package doMC is not available on Windows... Switching to single thread...")
        ncores <- 1
    }
    missing.method <- try(match.arg(missing.method, c("impute2mean", "impute2zero")))
    if(class(missing.method) == "try-error") stop("Error: \"missing.method\" must be \"impute2mean\" or \"impute2zero\".")
    if(class(null.obj) != "glmmkin.randomvec") stop("Error: \"null.obj\" must be a class glmmkin.randomvec object.")
    if(!class(interaction) %in% c("integer", "numeric", "character")) stop("Error: \"interaction\" should be an integer, numeric, or character vector.")
    # if(any(duplicated(null.obj$id_include))) {
    #   J <- sapply(unique(null.obj$id_include), function(x) 1*(null.obj$id_include==x))
    #   residuals <- crossprod(J, null.obj$scaled.residuals)
    #   residuals2 <- crossprod(J, null.obj$random.vectors)
    #   rm(J)
    # } else {
    #   residuals <- null.obj$scaled.residuals
    #   residuals2 <- null.obj$random.vectors
    # }
    residuals <- null.obj$scaled.residuals
    residuals2 <- null.obj$random.vectors
    qi <- length(interaction.covariates) # number of covariates with interaction effects but we don't test
    ei <- length(interaction) # number of covariates with interaction effects that we want to test
    if(class(interaction)=="character") {
        if (is.null(interaction.covariates)) {
            if (!all(interaction %in% colnames(null.obj$X))) {stop("there are interactions not in column name of covariate matrix.")}
            E <- as.matrix(null.obj$X[,interaction])
        } else {
            if (any(interaction.covariates %in% interaction)) {stop("there are interaction.covariates also specified as interaction.")}
            interaction <- c(interaction.covariates, interaction)
            if (!all(interaction %in% colnames(null.obj$X))) {stop("there are interaction or interaction.covariates not in column name of covariate matrix.")}
            E <- as.matrix(null.obj$X[,interaction])
        }
    } else {
        if (is.null(interaction.covariates)) {
            E <- as.matrix(null.obj$X[,interaction+1])
        } else {
            if (any(interaction.covariates %in% interaction)) {stop("there are interaction.covariates also specified as interaction.")}
            interaction <- c(interaction.covariates, interaction)
            E <- as.matrix(null.obj$X[,interaction+1])
        }
    }
    interaction <- as.character(interaction)
    n.E <- as.numeric(dim(E)[2]) # n.E = qi + ei
    # E <- scale(E, scale = FALSE)
    if(!grepl("\\.gds$", geno.file)) stop("Error: currently only .gds format is supported in geno.file!")
    gds <- SeqArray::seqOpen(geno.file)
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    if(any(is.na(match(null.obj$id_include, sample.id)))) warning("Check your data... Some individuals in null.obj$id_include are missing in sample.id of geno.file!")
    sample.id <- sample.id[sample.id %in% null.obj$id_include]
    if(length(sample.id) == 0) stop("Error: null.obj$id_include does not match sample.id in geno.file!")
    # match.id <- match(sample.id, unique(null.obj$id_include))
    if(any(duplicated(null.obj$id_include))) {
        match.id <- null.obj$id_include %in% sample.id
        null.obj$id_include <- null.obj$id_include[match.id]
        J <- t(sparseMatrix(i=1:length(null.obj$id_include), j=match(null.obj$id_include, unique(null.obj$id_include)[match(sample.id, unique(null.obj$id_include))]), x=1))
    } else match.id <- match(sample.id, null.obj$id_include)
    E <- as.matrix(E[match.id, , drop = FALSE])
    E <- scale(E, scale = FALSE)
    strata <- apply(E, 1, paste, collapse = ":")
    strata <- if(length(unique(strata))>length(strata)/100) NULL else as.numeric(as.factor(strata)) 
    if(!is.null(strata)) strata.list <- lapply(unique(strata), function(x) which(strata==x)) 
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
            write.table(t(c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "G.SCORE",paste("K.SCORE",1:ei,sep=""),"freq.strata.min","freq.strata.max")), meta.file.sample, quote = FALSE, row.names = FALSE, col.names = FALSE)
            meta.file.resample.handle <- file(meta.file.resample, "wb")
            writeBin(as.integer(ncol(residuals2)), meta.file.resample.handle, size = 4)
            gds <- SeqArray::seqOpen(geno.file)
            SeqArray::seqSetFilter(gds, sample.id = sample.id)
            nbatch.flush <- (p-1) %/% nperbatch + 1
            for(i in 1:nbatch.flush) {
                print(i)
                gc()
                tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*nperbatch+1):p] else variant.idx[((i-1)*nperbatch+1):(i*nperbatch)]
                SeqArray::seqSetFilter(gds, variant.id  = tmp.variant.idx, verbose = FALSE)
                miss <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
                freq <- 1 - SeqVarTools::alleleFrequency(gds)
                include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
                n.p<-sum(include)
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
                out$freq.strata.min <- rep(NA,sum(include))
                out$freq.strata.max <- rep(NA,sum(include))
                rm(alleles.list, include)
                SeqArray::seqSetFilter(gds, variant.id  = tmp.variant.idx, verbose = FALSE)
                geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
                out$N <- nrow(geno) - colSums(is.na(geno))
                if(any(duplicated(null.obj$id_include))) geno <- crossprod(J, geno)
                if(!is.null(strata)) { # E is not continuous
                    freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2) # freq.tmp is a matrix, each column is a strata, and each row is a varirant 
                    if (length(dim(freq.tmp)) == 2) freq_strata <- apply(freq.tmp, 1, range) else freq_strata <- as.matrix(range(freq.tmp)) # freq_strata is the range of allele freq across strata.list
                    rm(freq.tmp)
                } else {
                    freq_strata <- matrix(data=NA, nrow=2, ncol=dim(out)[1])
                }
                if(max(out$missrate)>0) {
                    miss.idx <- which(is.na(geno))
                    geno[miss.idx] <- if(missing.method=="impute2mean") 2*out$altfreq[ceiling(miss.idx/nrow(geno))] else 0
                }
                out$SCORE <- as.vector(crossprod(geno, residuals))
                K <- do.call(cbind, sapply((1+qi):ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
                SK <- matrix(crossprod(K,residuals), ncol=ei)
                write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK,as.data.frame(t(freq_strata))), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na="NA")
                writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
                rm(out,K,SK)
            }
            for (ii in seq(1,ei)){
                for(i in 1:nbatch.flush) {
                    gc()
                    tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*nperbatch+1):p] else variant.idx[((i-1)*nperbatch+1):(i*nperbatch)]
                    SeqArray::seqSetFilter(gds, variant.id  = tmp.variant.idx, verbose = FALSE)
                    miss <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
                    freq <- 1 - SeqVarTools::alleleFrequency(gds)
                    include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
                    n.p<-sum(include)
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
                    if(any(duplicated(null.obj$id_include))) geno <- crossprod(J, geno)
                    if(max(out$missrate)>0) {
                        miss.idx <- which(is.na(geno))
                        geno[miss.idx] <- if(missing.method=="impute2mean") 2*out$altfreq[ceiling(miss.idx/nrow(geno))] else 0
                    }
                    K <- do.call(cbind, sapply(qi+ii, function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
                    writeBin(as.numeric(crossprod(residuals2, K[,((ii-1)*n.p+1):(ii*n.p)])), meta.file.resample.handle, size = 4)
                    rm(out,K)
                }
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
        write.table(t(c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "G.SCORE",paste("K.SCORE",1:ei,sep=""),"freq.strata.min","freq.strata.max")), meta.file.sample, quote = FALSE, row.names = FALSE, col.names = FALSE)
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
            n.p<-sum(include)
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
            out$freq.strata.min<-rep(NA,dim(out)[1])
            out$freq.strata.max<-rep(NA,dim(out)[1])
            rm(alleles.list, include)
            SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
            geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
            out$N <- nrow(geno) - colSums(is.na(geno))
            if(any(duplicated(null.obj$id_include))) geno <- crossprod(J, geno)
            if(!is.null(strata)) { # E is not continuous
                freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2) # freq.tmp is a matrix, each column is a strata, and each row is a varirant 
                if (length(dim(freq.tmp)) == 2) freq_strata <- apply(freq.tmp, 1, range) else freq_strata <- as.matrix(range(freq.tmp)) # freq_strata is the range of allele freq across strata.list
                rm(freq.tmp) 
            } else {
                freq_strata <- matrix(data=NA, nrow=2, ncol=dim(out)[1])
            }
            if(max(out$missrate)>0) {
                miss.idx <- which(is.na(geno))
                geno[miss.idx] <- if(missing.method=="impute2mean") 2*out$altfreq[ceiling(miss.idx/nrow(geno))] else 0
            }
            out$SCORE <- as.vector(crossprod(geno, residuals))
            K <- do.call(cbind, sapply((1+qi):ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
            SK <- matrix(crossprod(K,residuals), ncol=ei)
            write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK,as.data.frame(t(freq_strata))), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na="NA")
            writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
            rm(out,K,SK)
        }
        for (ii in seq(1,ei)){
            for(i in 1:nbatch.flush) {
                gc()
                tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*nperbatch+1):p] else variant.idx[((i-1)*nperbatch+1):(i*nperbatch)]
                SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
                miss <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
                freq <- 1 - SeqVarTools::alleleFrequency(gds)
                include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
                n.p<-sum(include)
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
                if(any(duplicated(null.obj$id_include))) geno <- crossprod(J, geno)
                if(max(out$missrate)>0) {
                    miss.idx <- which(is.na(geno))
                    geno[miss.idx] <- if(missing.method=="impute2mean") 2*out$altfreq[ceiling(miss.idx/nrow(geno))] else 0
                }
                K <- do.call(cbind, sapply(qi+ii, function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
                writeBin(as.numeric(crossprod(residuals2, K[,((ii-1)*n.p+1):(ii*n.p)])), meta.file.resample.handle, size = 4)
                rm(out,K)
            }
        }
        SeqArray::seqClose(gds)
        close(meta.file.resample.handle)
    }
    return(invisible(NULL))
}

#' GLMM based single variant tests for gene-environment interactions
#' @description Use the summary statistic and stochastic statistics from \code{GE.stat} to perform single variant main effect tests, gene-environment interaction tests, and joint tests for association with genotypes in a GDS file (.gds). 7 user-defined tests are included: Main effect variance component test (MV), Main effect hybrid test of burden and variance component test using Fisher's method (MF), Interaction variance component test (IV), Interaction hybrid test of burden and variance component test using Fisher's method (IF), Joint variance component test (JV), Joint hybrid test of burden and variance component test using Fisher's method (JF), and Joint hybrid test of burden and variance component test using double Fisher's procedures (JD). 
#' @param meta.file.prefix A character vector for prefix of intermediate files (*.sample.* and *.resample.*). 
#' @param out.file The output file name.
#' @param n.files An integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix))).
#' @param n.pheno An integer indicating the number of phenotypes in multiple phenotype analysis (for single phenotype analysis, \code{n.pheno = 1}) (default = 1).
#' @param interaction A numeric or a character vector indicating the environmental factors. If a numberic vector, it represents which indices in the order of covariates are the environmental factors; if a character vector, it represents the variable names of the environmental factors.
#' @param interaction.covariates A numeric or a character vector indicating the interaction covariates. If a numeric vector, it represents which indices in the order of covariates are the interaction covariates; if a character vector, it represents the variable names of the interaction covariates.
#' @param MAF.range A numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param miss.cutoff The maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param auto.flip A logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
#' @param nperbatch An integer for how many SNPs should be tested in a batch (default = 100). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user's system.
#' @return NULL
#' @reference 
#' @author Han Chen, Nannan Wang
#' @seealso \code{glmmkin2randomvec}, \code{GE.stat}
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
#' if(!is.null(nullmod$P)){
#'   obj <- glmmkin2randomvec(nullmod)
#' }else{
#'   kinship.chol <- chol(GRM)
#'   obj<-glmmkin2randomvec(nullmod, Z = list(t(kinship.chol)))
#' }
#' out.prefix <- "test.GE.svt"
#' gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
#' interaction <- c("sex")
#' GE.stat(obj, interaction = interaction, geno.file = gdsfile, meta.file.prefix = out.prefix,MAF.range=c(0,0.5), miss.cutoff = 1)
#' out.file<-paste0(out.prefix,".out")
#' GE.svt.pval(meta.files.prefix = out.prefix, out.file, n.files = 1, n.pheno = 1, interaction=interaction, MAF.range=c(1e-7,0.5), miss.cutoff = 1, auto.flip=F)
#' }
#' @keywords generalized linear mixed model, gene-environment interaction
#' @export
GE.svt.pval <- function(meta.files.prefix, out.file, n.files = rep(1, length(meta.files.prefix)), n.pheno = 1, interaction, interaction.covariates = NULL, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, auto.flip = FALSE, nperbatch = 100)
{
    if(.Platform$endian!="little") stop("Error: platform must be little endian.")
    n.cohort <- length(meta.files.prefix)
    if(length(n.files) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and n.files do not match.")
    if(!class(interaction) %in% c("integer", "numeric", "character")) stop("Error: \"interaction\" should be an integer, numeric, or character vector.")
    qi <- length(interaction.covariates) # number of covariates with interaction effects but we don't test
    ei <- length(interaction) # number of covariates with interaction effects that we want to test
    if (!is.null(interaction.covariates)) {
        if (any(interaction.covariates %in% interaction)) {stop("there are interaction.covariates also specified as interaction.")}
        interaction <- c(interaction.covariates, interaction)
    }
    interaction <- as.character(interaction)
    n.E = qi + ei
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
    nSNP.resample<-vector("list",n.cohort)
    for(i in 1:n.cohort) {
        tmp.scores <- NULL
        for(j in 1:n.files[i]) {
            tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
            nSNP.resample[[i]][j]<-nrow(tmp)
            if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
            }
            tmp <- tmp[,c("chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "G.SCORE",paste("K.SCORE",1:ei,sep=""),"freq.strata.min","freq.strata.max")]
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
        scores[[i]] <- cbind(group.info, tmp.scores[match(variant.id1, paste(tmp.scores$chr, tmp.scores$pos, tmp.scores$ref, tmp.scores$alt, sep = ":")), c("N", "missrate", "altfreq", "G.SCORE", paste("K.SCORE",1:ei,sep=""), "freq.strata.min","freq.strata.max" ,"file", "variant.idx", "flip")])
        rm(tmp.scores)
        cons[[i]] <- file(paste0(meta.files.prefix[i], ".resample.1"), "rb")
        N.resampling[i] <- readBin(cons[[i]], what = "integer", n = 1, size = 4)
    }
    current.lines <- current.cons <- rep(1, n.cohort)
    nbatch.flush <- (p-1) %/% nperbatch + 1
    all.out <- NULL
    write.table(t(c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "freq.strata.min", "freq.strata.max", "PVAL.MAIN", paste("PVAL.GE",1:ei,sep=""), "PVAL.JOINT")), out.file, quote = FALSE, row.names = FALSE, col.names = FALSE)     
    for(i in 1:nbatch.flush) {
        tmp.idx <- if(i == nbatch.flush) group.info$idx[((i-1)*nperbatch+1):p] else group.info$idx[((i-1)*nperbatch+1):(i*nperbatch)]
        U.list <- V.list <- vector("list", n.cohort)
        variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- tmp.strata.min <- tmp.strata.max <- c()
        for(j in 1:n.cohort) {
            tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
            if(any(tmp.include <- !is.na(tmp.scores$G.SCORE))) {
                U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
                U.list[[j]]$G.SCORE <- U.list[[j]]$G.SCORE * U.list[[j]]$flip
                colName <- colnames( U.list[[j]])
                for (k in seq(1:ei)){
                    index <- which(colName==paste0("K.SCORE",k),arr.ind = TRUE)
                    U.list[[j]][,index] <- U.list[[j]][,index] * U.list[[j]]$flip
                }
                tmp.V <- vector("list", 1+ei)
                for (k in seq(1,1+ei)){
                    tmp.V[[k]] <- matrix(NA, sum(tmp.include), N.resampling[j])
                }
                for(ij in 1:sum(tmp.include)) {
                    if(U.list[[j]]$file[ij]!=current.cons[j]) {
                        close(cons[[j]])
                        current.cons[j] <- U.list[[j]]$file[ij]
                        cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.", current.cons[j]), "rb")
                        tmp.N.resampling <- readBin(cons[[j]], what = "integer", n = 1, size = 4)
                        if(tmp.N.resampling != N.resampling[j]) stop(paste0("Error: N.resampling in ", meta.files.prefix[j], ".resample.", current.cons[j], " does not match that in ",meta.files.prefix[j], ".resample.1"))
                        current.lines[j] <- 1
                    }
                    for (k in seq(1,1+ei)){
                        ind <- 4+4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-1)+(k-1)*4*N.resampling[j]*nSNP.resample[[j]][current.cons[j]]
                        seek(cons[[j]], where = ind, origin = "start", rw = "read")
                        tmp.V[[k]][ij,] <- readBin(cons[[j]], what = "numeric", n = N.resampling[j], size = 4)
                    }
                    current.lines[j] <- U.list[[j]]$variant.idx[ij]+1
                }
                V.list[[j]] <- vector("list", 1+ei)
                for (k in seq(1,1+ei)){
                    V.list[[j]][[k]] <- tmp.V[[k]] * U.list[[j]]$flip / sqrt(N.resampling[j])
                }
                rm(tmp.V)
                variant.indices <- c(variant.indices, U.list[[j]]$idx)
                tmp.N <- c(tmp.N, U.list[[j]]$N)
                tmp.Nmiss <- c(tmp.Nmiss, U.list[[j]]$N * U.list[[j]]$missrate/(1-U.list[[j]]$missrate))
                tmp.AC <- c(tmp.AC, 2*U.list[[j]]$N*U.list[[j]]$altfreq)
                tmp.strata.min<-c(tmp.strata.min, U.list[[j]]$freq.strata.min)
                tmp.strata.max<-c(tmp.strata.max, U.list[[j]]$freq.strata.max)
            }
        }
        if(length(variant.indices) == 0) next
        tmp.variant.indices <- variant.indices
        variant.indices <- sort(unique(variant.indices))
        N <- sapply(variant.indices, function(x) sum(tmp.N[tmp.variant.indices==x]))
        Nmiss <- sapply(variant.indices, function(x) sum(tmp.Nmiss[tmp.variant.indices==x]))
        AF <- sapply(variant.indices, function(x) sum(tmp.AC[tmp.variant.indices==x]))/2/N
        include <- (Nmiss/(N+Nmiss) <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
        if(any(!is.na(tmp.strata.min))) { # E is not continuous
            include <- include & !is.na(tmp.strata.min) & !is.na(tmp.strata.max) & tmp.strata.min >= MAF.range[1] & tmp.strata.max <= 1-MAF.range[1]
        }
        rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices)
        if(sum(include) == 0) next
        variant.indices <- variant.indices[include]
        out <- data.frame(SNP = group.info$SNP[variant.indices], chr = group.info$chr[variant.indices], pos = group.info$pos[variant.indices], ref = group.info$ref[variant.indices], alt = group.info$alt[variant.indices], N = N[include], missrate = (Nmiss/(N+Nmiss))[include], altfreq = AF[include], freq.strata.min=tmp.strata.min[include], freq.strata.max=tmp.strata.max[include])
        n.p <- length(variant.indices)
        U <- matrix(0, n.p, 1+ei)
        V <- vector("list", 1+ei)
        for (k in seq(1,1+ei)){
            V[[k]] <- matrix(0, n.p, sum(N.resampling))
        }
        for(j in 1:n.cohort) {
            if(!is.null(U.list[[j]]) & !is.null(V.list[[j]][1])) {
                IDX <- match(U.list[[j]]$idx, variant.indices)
                if(sum(!is.na(IDX)) == 0) next
                IDX2 <- which(!is.na(IDX))
                IDX <- IDX[IDX2]
                colName <- colnames( U.list[[j]])
                for (k in seq(1,1+ei)){
                    if (k==1) {
                        U[IDX,1] <- U[IDX,1]+U.list[[j]]$G.SCORE[IDX2]
                    } else {
                        index <- which(colName==paste0("K.SCORE",(k-1)),arr.ind = TRUE)
                        U[IDX,k] <- U[IDX,k]+U.list[[j]][IDX2,index]
                    }
                    V[[k]][IDX, (sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])] <- V[[k]][IDX,(sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])]+V.list[[j]][[k]][IDX2,]
                }
            }
        }

        UGC<-as.vector(U[,1:(1+qi)])     #U of geno and interaction.covariates
        VGC_u<-do.call(rbind, lapply(V[1:(1+qi)], unlist))  #rbind(V[[1:(1+qi)]]); 
        c1 <- rep(1:n.p,n.pheno)+rep((0:(n.pheno-1))*n.p*(1+qi), each=n.p) # index for GPG and row.index for GPC
        if(!is.null(interaction.covariates)) {
            c2 <- rep((n.p+1):(n.p*(1+qi)),n.pheno)+rep((0:(n.pheno-1))*n.p*(1+qi), each=n.p*qi) # index for CPC and col.index for GPC
        }
        
        SK <- as.vector(U[,(1+qi+1):(1+n.E)])
        VG_u<-do.call(rbind, lapply(V[1:n.pheno], unlist)) 
        VK_u<-do.call(rbind, lapply(V[(1+qi+1):(1+n.E)], unlist)) # rbind(V[[(1+qi+1):(1+n.E)]]) #
        UG <- UGC[c1]
        VG_u <- VGC_u[c1,]
        if(length(c1)==1){
            VG_u <- as.matrix(VG_u)
            VG_u <- t(VG_u)
        }
        ng <- nrow(VG_u)

        out$SCORE.MAIN <- UG
        out$VAR.MAIN <- rowSums(VG_u^2)
        out$PVAL.MAIN<- pchisq(out$SCORE.MAIN^2/out$VAR.MAIN, 1, lower = FALSE)

        ng1   <- ng+1
        ei1<-ei+1
        ngei1<-ng*ei1
        
        Vgg <- rowSums(VG_u^2)
        VGG<-matrix(0,ng,ng)
        diag(VGG)<-Vgg
        Vkk <- rowSums(VK_u^2)
        VKK<-matrix(0,ng,ng)
        diag(VKK)<-Vkk
        Vgk <- rowSums(VG_u*VK_u)
        VGK<-matrix(0,ng,ng)
        diag(VGK)<-Vgk
        V_GK<- rbind(cbind(VGG,VGK),cbind(t(VGK),VKK))
        U_GK <- c(UG, SK)
        U_GK <- (rep(1, 1+ei) %x% diag(ng)) * as.vector(U_GK)
        V_GK_i <- try(solve(V_GK),silent = TRUE)
        if(class(V_GK_i)[1] == "try-error") V_GK_i <- MASS::ginv(V_GK)
        BETA.GK <- crossprod(V_GK_i, U_GK) #cbind(out$BETA.MAIN,BETA.K) #crossprod(V_GK_i, U_GK)

        IV.K_i<-try(solve(V_GK_i[ng1:ngei1, ng1:ngei1]), silent=TRUE)
        if(class(IV.K_i)[1] == "try-error") IV.K_i <- try(MASS::ginv(V_GK_i[ng1:ngei1, ng1:ngei1]), silent = TRUE)
        STAT.K<-diag(crossprod(BETA.GK[ng1:ngei1,], crossprod(IV.K_i, BETA.GK[ng1:ngei1,])))
        STAT.JOINT <- diag(crossprod(BETA.GK, crossprod(V_GK, BETA.GK)))

        PVAL.K   <- pchisq(STAT.K, df=ei, lower.tail=FALSE)
        PVAL.JOINT <- ifelse(is.na(out$PVAL.MAIN), NA, pchisq(STAT.JOINT, df=1+ei, lower.tail=FALSE))
        
        write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "freq.strata.min", "freq.strata.min", "PVAL.MAIN")], PVAL.K, PVAL.JOINT), out.file, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na="NA")
    }
    for(i in 1:n.cohort) close(cons[[i]])
    # return(all.out)
}


fisher_pval <- function(p) {
    is.valid.p <- !is.na(p) & p > 0 & p <= 1
    if(sum(is.valid.p) == 0) return(NA)
    p <- p[is.valid.p]
    pchisq(-2*sum(log(p)), df = 2*length(p), lower.tail = FALSE)
}

.Q_pval <- function(Q, lambda, method = "davies") {
    if(method == "davies") {
        tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6))
        pval <- tmp$Qq
        if((tmp$ifault > 0) | (pval <= 1e-5) | (pval >= 1)) method <- "kuonen"
    }
    if(method == "kuonen") {
        pval <- pKuonen(x = Q, lambda = lambda)
        if(is.na(pval)) method <- "liu"
    }
    if(method == "liu") pval <- CompQuadForm::liu(q = Q, lambda = lambda)
    return(pval)
}

.quad_pval <- function(U, V, method = "davies") {
    Q <- sum(U^2)
    lambda <- eigen(V, only.values = TRUE, symmetric=TRUE)$values
    lambda <- lambda[lambda > 0]
    pval <- .Q_pval(Q, lambda, method = method)
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
            if(dim1 < dim2) lambda <- zapsmall(eigen(tcrossprod(V.temp), symmetric = T, only.values = T)$values)
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
