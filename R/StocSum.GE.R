# modified from GMMAT v1.1.0
# 20201204: implemented group.idx for heteroscedastic linear mixed models
#       added single variant tests
# 20201209: switched to eigen() when svd() fails

StocSum.GE.glmmkin2randomvec <- function(obj, Z = NULL, N.randomvec = 1000, group.idx=NULL,cluster.idx=NULL,robust = FALSE) {
  if(class(obj) != "glmmkin") stop("Error: \"obj\" must be a class glmmkin object.")
  N <- length(obj$id_include)
  random.vectors <- matrix(rnorm(N*N.randomvec),nrow=N,ncol=N.randomvec)
  # random.vectors<-diag(N)*sqrt(N)
  # N.randomvec<-N
  if(!is.null(obj$P) && !robust) {
    eig <- eigen(obj$P, symmetric = TRUE)
    random.vectors <- tcrossprod(eig$vectors, t(random.vectors * sqrt(pmax(eig$values, 0))))
    rm(eig)
  } else {
    if(obj$n.groups != 1 && (is.null(group.idx) || !all.equal(seq_len(obj$n.groups), sort(unique(group.idx))))) stop("Error: heteroscedastic linear mixed models should include a valid group.idx argument.")
    if(is.null(group.idx)) group.idx <- rep(1, N)
    if(!robust) random.vectors <- sqrt(obj$theta[group.idx]) * random.vectors
    #   else random.vectors <- random.vectors * abs(obj$residuals)
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

StocSum.GE.stat <- function(null.obj, interaction, interaction.covariates = NULL, geno.file, meta.file.prefix, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, missing.method = "impute2mean", nperbatch = 10000, ncores = 1)
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
  # write.table(cbind(sample.id,strata),file=paste0(meta.file.prefix[1], ".strata"), quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")  #sample.id, strata
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
        out$freq.strata.min<-rep(NA,length(include))
        out$freq.strata.max<-rep(NA,length(include))
        rm(alleles.list, include)
        SeqArray::seqSetFilter(gds, variant.id  = tmp.variant.idx, verbose = FALSE)
        geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
        out$N <- nrow(geno) - colSums(is.na(geno))
        if(any(duplicated(null.obj$id_include))) geno <- crossprod(J, geno)
        if(!is.null(strata)) { # E is not continuous
          freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2) # freq.tmp is a matrix, each column is a strata, and each row is a varirant 
          if (length(dim(freq.tmp)) == 2) freq_strata <- apply(freq.tmp, 1, range) else freq_strata <- as.matrix(range(freq.tmp)) # freq_strata is the range of allele freq across strata.list
          rm(freq.tmp)
          # write.table(cbind(out[,c("SNP","chr","pos","ref","alt")], t(freq_strata)),file=paste0(meta.file.prefix[1], ".strata"), quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")  
        }
        if(max(out$missrate)>0) {
          miss.idx <- which(is.na(geno))
          geno[miss.idx] <- if(missing.method=="impute2mean") 2*out$altfreq[ceiling(miss.idx/nrow(geno))] else 0
        }
        out$SCORE <- as.vector(crossprod(geno, residuals))
        K <- do.call(cbind, sapply((1+qi):ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
        SK <- matrix(crossprod(K,residuals), ncol=ei)
        if(!is.null(strata)){
          write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK,as.data.frame(t(freq_strata))), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na="NA")
        } else {
          write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
        }
        # write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK,t(freq_strata)), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
        writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
        # for (i in seq(1,ei)){
        #   writeBin(as.numeric(crossprod(residuals2, K[,((i-1)*n.p+1):(i*n.p)])), meta.file.resample.handle, size = 4)
        # }
        rm(out,K,SK)
      }
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
        # out$SCORE <- as.vector(crossprod(geno, residuals))
        K <- do.call(cbind, sapply((1+qi):ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
        # SK <- matrix(crossprod(K,residuals), ncol=ei)
        # write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
        # writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
        for (i in seq(1,ei)){
          writeBin(as.numeric(crossprod(residuals2, K[,((i-1)*n.p+1):(i*n.p)])), meta.file.resample.handle, size = 4)
        }
        rm(out,K)
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
          # write.table(cbind(out[,c("SNP","chr","pos","ref","alt")], t(freq_strata)),file=paste0(meta.file.prefix[1], ".strata"), quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")  
      }
      if(max(out$missrate)>0) {
        miss.idx <- which(is.na(geno))
        geno[miss.idx] <- if(missing.method=="impute2mean") 2*out$altfreq[ceiling(miss.idx/nrow(geno))] else 0
      }
      out$SCORE <- as.vector(crossprod(geno, residuals))
      K <- do.call(cbind, sapply((1+qi):ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
      SK <- matrix(crossprod(K,residuals), ncol=ei)
      if(!is.null(strata)){
        write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK,as.data.frame(t(freq_strata))), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na="NA")
      } else {
        write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
      }
      # write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK,t(freq_strata)), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
      writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
      # for (ii in seq(1,ei)){
      #   writeBin(as.numeric(crossprod(residuals2, K[,((ii-1)*n.p+1):(ii*n.p)])), meta.file.resample.handle, size = 4)
      # }
      rm(out,K,SK)
    }
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
      # out$SCORE <- as.vector(crossprod(geno, residuals))
      K <- do.call(cbind, sapply((1+qi):ncol(E), function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
      # SK <- matrix(crossprod(K,residuals), ncol=ei)
      # write.table(cbind(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")],SK), meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
      # writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
      for (ii in seq(1,ei)){
        writeBin(as.numeric(crossprod(residuals2, K[,((ii-1)*n.p+1):(ii*n.p)])), meta.file.resample.handle, size = 4)
      }
      rm(out,K)
    }
    SeqArray::seqClose(gds)
    close(meta.file.resample.handle)
  }
  return(invisible(NULL))
}

StocSum.GE.pval <- function(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), n.pheno = 1, cohort.group.idx = NULL, group.file, group.file.sep = "\t", interatcion, interaction.covariates = NULL, MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(0.5, 0.5), miss.cutoff = 1, method = "davies", tests = "JF", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE, auto.flip = FALSE)
{
  if(.Platform$endian!="little") stop("Error: platform must be little endian.")
  n.cohort <- length(meta.files.prefix)
  if(length(n.files) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and n.files do not match.")
  if(!is.null(cohort.group.idx)) {
    if(length(cohort.group.idx) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and cohort.group.idx do not match.")
    cohort.group.idx <- as.numeric(factor(cohort.group.idx))
    n.cohort.groups <- length(unique(cohort.group.idx))
  }
  if(any(!tests %in% c("MV", "MF", "IV", "IF", "JV", "JF", "JD"))) stop("Error: \"tests\" should only include \"MV\" for the main effect variance component test, \"MF\" for the main effect combined test of the burden and variance component tests using Fisher\'s method, \"IV\" for the interaction variance component test, \"IF\" for the interaction combined test of the burden and variance component tests using Fisher\'s method, \"JV\" for the joint variance component test for main effect and interaction, \"JF\" for the joint combined test of the burden and variance component tests for main effect and interaction using Fisher\'s method, or \"JD\" for the joint combined test of the burden and variance component tests for main effect and interaction using double Fisher\'s method.")
  MV <- "MV" %in% tests
  MF <- "MF" %in% tests
  IV <- "IV" %in% tests
  IF <- "IF" %in% tests
  JV <- "JV" %in% tests
  JF <- "JF" %in% tests
  JD <- "JD" %in% tests
  if(!class(interaction) %in% c("integer", "numeric", "character")) stop("Error: \"interaction\" should be an integer, numeric, or character vector.")
  qi <- length(interaction.covariates) # number of covariates with interaction effects but we don't test
  ei <- length(interaction) # number of covariates with interaction effects that we want to test
  if (!is.null(interaction.covariates)) {
    if (any(interaction.covariates %in% interaction)) {stop("there are interaction.covariates also specified as interaction.")}
    interaction <- c(interaction.covariates, interaction)
  }
  interaction <- as.character(interaction)
  n.E = qi + ei
  # strata <- read.table(file=paste0(meta.files.prefix[1], ".strata"), header = F)
  # strata <- if(length(unique(strata))>length(strata)/100) NULL else as.numeric(as.factor(strata)) 
  # if(!is.null(strata)) strata.list <- lapply(unique(strata), function(x) which(strata==x)) 
  group.info <- try(read.table(group.file, header = FALSE, stringsAsFactors = FALSE, sep = group.file.sep), silent = TRUE)
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
    # rm(tmp.scores)
    cons[[i]] <- file(paste0(meta.files.prefix[i], ".resample.1"), "rb")
    N.resampling[i] <- readBin(cons[[i]], what = "integer", n = 1, size = 4)
    # if (i==1){
    #   tmp.strata<-try(read.table(paste0(meta.files.prefix[i], ".strata"), header = TRUE, as.is = TRUE))
    #   if ('try-error' %in% class(tmp.strata)){
    #     strataflag<-0
    #   } else {
    #     strataflag<-1
    #     freq_strata<-cbind(group.info,tmp.strata[match(variant.id1,paste(tmp.strata$chr,tmp.strata$pos,tmp.strata$ref,tmp.strata$alt,sep=":")),])
    #     rm(tmp.strata)
    #   }
    # }
    rm(tmp.scores)
  }
  n.variants <- rep(0,n.groups)
  # miss.min <- rep(NA,n.groups)
  # miss.mean <- rep(NA, n.groups)
  # miss.max <- rep(NA, n.groups)
  # freq.min <- rep(NA, n.groups)
  # freq.mean <- rep(NA, n.groups)
  # freq.max <- rep(NA, n.groups)
  # freq.strata.min <- rep(NA, n.groups)
  # freq.strata.max <- rep(NA, n.groups)
  if(MV | JV) MV.pval <- rep(NA, n.groups)
  if(IV | JV) IV.pval <- rep(NA, n.groups)
  if(JV) JV.pval <- rep(NA, n.groups)
  if(MF | JF | JD) MF.pval <- rep(NA, n.groups)
  if(IF | JF | JD) IF.pval <- rep(NA, n.groups)
  if(JF) JF.pval <- rep(NA, n.groups)
  if(JD) JD.pval <- rep(NA, n.groups)
  current.lines <- current.cons <- rep(1, n.cohort)
  for(i in 1:n.groups) {
    tmp.idx <- group.idx.start[i]:group.idx.end[i]
    #tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
    U.list <- V.list <- vector("list", n.cohort)
    variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- tmp.strata.min <- tmp.strata.max <- c()
    for(j in 1:n.cohort) {
      tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
      if(any(tmp.include <- !is.na(tmp.scores$G.SCORE))) {
        U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
        U.list[[j]]$weight <- U.list[[j]]$weight * U.list[[j]]$flip
        U.list[[j]]$G.SCORE <- U.list[[j]]$G.SCORE * U.list[[j]]$weight
        colName <- colnames( U.list[[j]])
        for (k in seq(1:ei)){
          index <- which(colName==paste0("K.SCORE",k),arr.ind = TRUE)
          U.list[[j]][,index] <- U.list[[j]][,index] * U.list[[j]]$weight
        }
        # U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$weight
        tmp.V <- vector("list", 1+ei)
        for (k in seq(1,1+ei)){
          tmp.V[[k]] <- matrix(NA, sum(tmp.include), N.resampling[j])
        }
        # tmp.V <- matrix(NA, sum(tmp.include), N.resampling[j])
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
          # if(U.list[[j]]$variant.idx[ij]!=current.lines[j]) seek(cons[[j]], where = 4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-current.lines[j]), origin = "current", rw = "read")
          # tmp.V[ij,] <- readBin(cons[[j]], what = "numeric", n = N.resampling[j], size = 4)
          current.lines[j] <- U.list[[j]]$variant.idx[ij]+1
        }
        # V.list[[j]] <- tmp.V * U.list[[j]]$weight / sqrt(N.resampling[j])
        V.list[[j]] <- vector("list", 1+ei)
        for (k in seq(1,1+ei)){
          V.list[[j]][[k]] <- tmp.V[[k]] * U.list[[j]]$weight / sqrt(N.resampling[j])
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
          # ifreq_strata<-freq_strata[variant.indices,]
        include <- include & !is.na(tmp.strata.min) & !is.na(tmp.strata.max) & tmp.strata.min >= MAF.range[1] & tmp.strata.max <= 1-MAF.range[1]
    }
    rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices, N, Nmiss)
    if(sum(include) == 0) next
    variant.indices <- variant.indices[include]
    n.p <- length(variant.indices)
    n.variants[i] <- n.p
    # U <- if(!is.null(cohort.group.idx)) rep(0, n.cohort.groups*n.p) else rep(0, n.p)
    # V <- if(!is.null(cohort.group.idx)) matrix(0, n.cohort.groups*n.p, sum(N.resampling)) else matrix(0, n.p, sum(N.resampling))
    U <- if(!is.null(cohort.group.idx)) matrix(0,n.cohort.groups*n.p,1+ei) else matrix(0, n.p, 1+ei)
    V <- vector("list", 1+ei)
    for (k in seq(1,1+ei)){
      V[[k]] <- if(!is.null(cohort.group.idx)) matrix(0, n.cohort.groups*n.p, sum(N.resampling)) else matrix(0, n.p, sum(N.resampling))
    }
    for(j in 1:n.cohort) {
      if(!is.null(U.list[[j]]) & !is.null(V.list[[j]][1])) {
        IDX <- match(U.list[[j]]$idx, variant.indices)
        if(sum(!is.na(IDX)) == 0) next
        IDX2 <- which(!is.na(IDX))
        IDX <- IDX[IDX2]
        if(!is.null(cohort.group.idx)) IDX <- IDX+n.p*(cohort.group.idx[j]-1)
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
        # U[IDX] <- U[IDX]+U.list[[j]]$SCORE[IDX2]
        # V[IDX, (sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])] <- V[IDX,(sum(N.resampling[1:j])-N.resampling[j]+1):sum(N.resampling[1:j])]+V.list[[j]][IDX2,]
      }
    }
    tmp.weight <- MAF.weights.beta.fun(AF[include], MAF.weights.beta[1], MAF.weights.beta[2])
    if(use.minor.allele) tmp.weight[AF[include] > 0.5] <- -tmp.weight[AF[include] > 0.5]
    if(!is.null(cohort.group.idx)) tmp.weight <- rep(tmp.weight, n.cohort.groups)
    U <- U*tmp.weight
    for (k in seq((1+qi),(1+n.E))){
      V[[k]] <- V[[k]]*tmp.weight
    }
    UGC<-as.vector(U[,1:(1+qi)])     #U of geno and interaction.covariates
    VGC_u<-do.call(rbind, lapply(V[1:(1+qi)], unlist))  #rbind(V[[1:(1+qi)]]); 
    VGC<-tcrossprod(VGC_u)       #covariance matrix V of geno and interaction.covariates
    if(MV | MF | JV | JF | JD) c1 <- rep(1:n.p,n.pheno)+rep((0:(n.pheno-1))*n.p*(1+qi), each=n.p) # index for GPG and row.index for GPC
    if(!is.null(interaction.covariates) && (JV | JF | JD)) {
      c2 <- rep((n.p+1):(n.p*(1+qi)),n.pheno)+rep((0:(n.pheno-1))*n.p*(1+qi), each=n.p*qi) # index for CPC and col.index for GPC
      CPC_i <- try(solve(VGC[c2,c2]), silent = TRUE)
      if(class(CPC_i)[1] == "try-error") CPC_i <- MASS::ginv(VGC[c2,c2])
      U.adj <- UGC[c1] - tcrossprod(tcrossprod(VGC[c1,c2],CPC_i),t(UGC[c2]))
      V.adj <- VGC[c1,c1] - tcrossprod(tcrossprod(VGC[c1,c2],CPC_i),VGC[c1,c2])
    }
    if(IV | IF | JV | JF | JD) {
      SK <- as.vector(U[,(1+qi+1):(1+n.E)])
      VG_u<-do.call(rbind, lapply(V[1:n.pheno], unlist)) 
      # VE_u<-do.call(rbind, lapply(V[(1+n.pheno):(1+n.E)], unlist)) 
      # VEG<-tcrossprod(VE_u,VG_u)    #E=[C,K]
      VK_u<-do.call(rbind, lapply(V[(1+qi+1):(1+n.E)], unlist)) # rbind(V[[(1+qi+1):(1+n.E)]]) #
      VK<-tcrossprod(VK_u)
      VKG<-tcrossprod(VK_u,VG_u)
      VGC_i <- try(solve(VGC), silent = TRUE)
      if(class(VGC_i)[1] == "try-error") VGC_i <- MASS::ginv(VGC)
      IV.U<-SK-tcrossprod(tcrossprod(VKG,VGC_i),t(UGC))
      IV.V<-VK-tcrossprod(tcrossprod(VKG,VGC_i),VKG)
    }
    if(MV | MF | JV | JF | JD) {
      UG <- UGC[c1]
      VG <- VGC[c1,c1]
    }
    if(MV | JV) MV.pval[i] <- tryCatch(.quad_pval(U = UG, V = VG, method = method), error = function(e) { NA })
    if(IV | JV) IV.pval[i] <- tryCatch(.quad_pval(U = IV.U, V = IV.V, method = method), error = function(e) { NA })
    if(JV && is.null(interaction.covariates)) JV.pval[i] <- tryCatch(fisher_pval(c(MV.pval[i], IV.pval[i])), error = function(e) { MV.pval[i] })
    if(JV && !is.null(interaction.covariates)) {
      MV.adj.pval <- tryCatch(.quad_pval(U = U.adj, V = V.adj, method = method), error = function(e) { NA })
      JV.pval[i] <- tryCatch(fisher_pval(c(MV.adj.pval, IV.pval[i])), error = function(e) { MV.adj.pval })
    }
    if(MF | JF | JD) {
      MF.BU <- sum(UG)
      MF.BV <- sum(VG)
      MF.Bp <- pchisq(MF.BU^2/MF.BV,df=1,lower.tail=FALSE)
      V.rowSums <- if(length(VG)>1) rowSums(VG) else VG[1]
      MF.U <- UG - V.rowSums * MF.BU / MF.BV
      MF.V <- VG - tcrossprod(V.rowSums) / MF.BV
      if(MF.BV == 0 | mean(abs(MF.V)) < sqrt(.Machine$double.eps)) MF.p <- NA
      else MF.p <- tryCatch(.quad_pval(U = MF.U, V = MF.V, method = method), error = function(e) { NA })
      MF.pval[i] <- tryCatch(fisher_pval(c(MF.Bp, MF.p)), error = function(e) { MF.Bp })
    }
    if((JF | JD) && !is.null(interaction.covariates)) {
      MF.BU.adj <- sum(U.adj)
      MF.BV.adj <- sum(V.adj)
      MF.Bp.adj <- pchisq(MF.BU.adj^2/MF.BV.adj,df=1,lower.tail=FALSE)
      V.adj.rowSums <- if(length(V.adj)>1) rowSums(V.adj) else V.adj[1]
      MF.U.adj <- U.adj  - V.adj.rowSums * MF.BU.adj / MF.BV.adj
      MF.V.adj  <- V.adj  - tcrossprod(V.adj.rowSums) / MF.BV.adj
      if(MF.BV.adj == 0 | mean(abs(MF.V.adj)) < sqrt(.Machine$double.eps)) MF.adj.p <- NA
      else MF.adj.p <- tryCatch(.quad_pval(U = MF.U.adj, V = MF.V.adj, method = method), error = function(e) { NA })
      MF.adj.pval <- tryCatch(fisher_pval(c(MF.Bp.adj, MF.adj.p)), error = function(e) { MF.Bp.adj })
    }
    if(IF | JF | JD) {
      IF.BU <- sum(IV.U)
      IF.BV <- sum(IV.V)
      IF.Bp <- pchisq(IF.BU^2/IF.BV,df=1,lower.tail=FALSE)
      IV.V.rowSums <- if(length(IV.V)>1) rowSums(IV.V) else IV.V[1]
      IF.U <- IV.U - IV.V.rowSums * IF.BU / IF.BV
      IF.V <- IV.V - tcrossprod(IV.V.rowSums) / IF.BV
      if(IF.BV == 0 | mean(abs(IF.V)) < sqrt(.Machine$double.eps)) IF.p <- NA
      else IF.p <- tryCatch(.quad_pval(U = IF.U, V = IF.V, method = method), error = function(e) { NA })
      IF.pval[i] <- tryCatch(fisher_pval(c(IF.Bp, IF.p)), error = function(e) { IF.Bp })
    }
    if(JF && is.null(interaction.covariates)) JF.pval[i] <- tryCatch(fisher_pval(c(MF.Bp, MF.p, IF.Bp, IF.p)), error = function(e) { MF.Bp })
    if(JF && !is.null(interaction.covariates)) JF.pval[i] <- tryCatch(fisher_pval(c(MF.Bp.adj, MF.adj.p, IF.Bp, IF.p)), error = function(e) { MF.Bp.adj })
    if(JD && is.null(interaction.covariates)) JD.pval[i] <- tryCatch(fisher_pval(c(MF.pval[i], IF.pval[i])), error = function(e) { MF.pval[i] })
    if(JD && !is.null(interaction.covariates)) JD.pval[i] <- tryCatch(fisher_pval(c(MF.adj.pval, IF.pval[i])), error = function(e) { MF.adj.pval })
  }
  for(i in 1:n.cohort) close(cons[[i]])
  out <- data.frame(group=groups, n.variants=n.variants)
  if(MV | JV) out$MV.pval <- MV.pval
  if(MF | JF | JD) out$MF.pval <- MF.pval
  if(IV | JV) out$IV.pval <- IV.pval
  if(IF | JF | JD) out$IF.pval <- IF.pval
  if(JV) out$JV.pval <- JV.pval
  if(JF) out$JF.pval <- JF.pval
  if(JD) out$JD.pval <- JD.pval
  return(out)
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
