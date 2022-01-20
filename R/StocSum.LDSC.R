# modified from GMMAT v1.1.0
# 20201204: implemented group.idx for heteroscedastic linear mixed models
# 	    added single variant tests
# 20201209: switched to eigen() when svd() fails

StocSum.LDSC.glmmkin2randomvec <- function(obj, Z = NULL, N.randomvec = 1000, group.idx = NULL, cluster.idx = NULL, robust = FALSE) {
    if(class(obj) != "glmmkin") stop("Error: \"obj\" must be a class glmmkin object.")
    N <- length(obj$id_include)
    random.vectors <- matrix(rnorm(N*N.randomvec),nrow=N,ncol=N.randomvec)
    r<-random.vectors
    obj$P<-NULL
    obj$Sigma_i<-diag(N)
    obj$Sigma_iX<-matrix(1,N,1)
    obj$cov<-as.matrix(1/N)
    obj$theta<-1
    obj$n.groups<-1
    Z <- NULL
    if(!is.null(obj$P) && !robust) {
        eig <- eigen(obj$P, symmetric = TRUE)
        random.vectors <- tcrossprod(eig$vectors, t(random.vectors * sqrt(pmax(eig$values, 0))))
        rm(eig)
    } else {
        if(obj$n.groups != 1 && (is.null(group.idx) || !all.equal(seq_len(obj$n.groups), sort(unique(group.idx))))) stop("Error: heteroscedastic linear mixed models should include a valid group.idx argument.")
        if(is.null(group.idx)) group.idx <- rep(1, N)
        if(!robust) random.vectors <- sqrt(obj$theta[group.idx]) * random.vectors
        #		else random.vectors <- random.vectors * abs(obj$residuals)
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
        else {
            random.vectors <- crossprod(obj$Sigma_i, random.vectors) - tcrossprod(obj$Sigma_iX, tcrossprod(crossprod(random.vectors, obj$Sigma_iX), obj$cov))
        }
    }
    out <- list(theta = obj$theta, scaled.residuals = obj$scaled.residuals, random.vectors = as.matrix(random.vectors),r=as.matrix(r), id_include = obj$id_include)
    class(out) <- "glmmkin.randomvec"
    return(out)
}

StocSum.LDSC.stat <- function(null.obj, geno.file, meta.file.prefix, MAF.range = c(1e-3, 0.5), miss.cutoff = 1, missing.method = "impute2mean", nperbatch = 10000, ncores = 1)
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
        #adding by NW----start
        # sdgeno<-apply(geno,2,sd)
        # geno<-sweep(geno,2,sdgeno,'/')
        # for(i in seq(1, ncol(geno),1)){
        #   geno[,i] <- geno[,i]/sd(geno[,i])
        #   # geno[,i] <- (geno[,i]-mean(geno[,i]))
        #   # geno[,i] <- geno[,i]/sd(geno[,i])
        # }
        #adding by NW----end
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
      #adding by NW----start
      # sdgeno<-apply(geno,2,sd)
      # geno<-sweep(geno,2,sdgeno,'/')
      # for(i in seq(1, ncol(geno),1)){
      #   geno[,i] <- geno[,i]/sd(geno[,i])
      #   # geno[,i] <- (geno[,i]-mean(geno[,i]))
      #   # geno[,i] <- geno[,i]/sd(geno[,i])
      # }
      #adding by NW----end
      out$SCORE <- as.vector(crossprod(geno, residuals))
      write.table(out[,c("SNP", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")], meta.file.sample, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na=".")
      writeBin(as.numeric(crossprod(residuals2, geno)), meta.file.resample.handle, size = 4)
      rm(out)
    }
    SeqArray::seqClose(gds)
    close(meta.file.resample.handle)
  } 
  # tp<-list()
  # tp$geno<-geno
  # tp$R<-residuals2
  # tp$U<-crossprod(residuals2, geno)
  # tp$matchid<-match.id
  # return(tp)
  return(invisible(NULL))
}

StocSum.LDSC.ldscWinFastParrel <- function(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), N.randomvec=1000, MAF.range = c(1e-3, 0.5), MAF.weights.beta = c(0.5, 0.5), miss.cutoff = 1, wind.b=1000000, use.minor.allele = FALSE, auto.flip = FALSE, nperbatch = 100000, ncores = 1)
{
  cohort.group.idx <- NULL
  if(.Platform$endian!="little") stop("Error: platform must be little endian.")
  n.cohort <- length(meta.files.prefix)
  if(length(n.files) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and n.files do not match.")
  if(!is.null(cohort.group.idx)) {
    if(length(cohort.group.idx) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and cohort.group.idx do not match.")
    cohort.group.idx <- as.numeric(factor(cohort.group.idx))
    n.cohort.groups <- length(unique(cohort.group.idx))
  }
  group.info <- NULL
  for(i in 1:n.cohort) {
    for(j in 1:n.files[i]) {
      tmp <- try(read.table(paste0(meta.files.prefix[i], ".sample.", j), header = TRUE, as.is = TRUE))
      if (class(tmp) == "try-error") {
        stop(paste0("Error: cannot read ", meta.files.prefix[i], ".sample.", j, "!"))
      }
      tmp <- tmp[,c("SNP", "chr", "pos", "ref", "alt","N","altfreq")]
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
  n.groups.all<-nrow(group.info)
  group.info$idx <- 1:nrow(group.info)
  group.info$pos.end<-group.info$pos+wind.b
  group.info$pos.end[group.info$pos.end>group.info$pos[n.groups.all]]<-group.info$pos[n.groups.all]
  group.info$pos.start<-group.info$pos-wind.b
  group.info$pos.start[group.info$pos.start<group.info$pos[1]]<-group.info$pos[1]
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
  # nbatch.flush <- (n.groups.all-1) %/% nperbatch + 1
  all.out <- NULL
  ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
  if(ncores > 1) {
    doMC::registerDoMC(cores = ncores)
    n.groups.percore <- (n.groups.all-1) %/% ncores + 1
    n.groups.percore_1 <- n.groups.percore * ncores - n.groups.all
    b <- NULL
    all.out <- foreach(b=1:ncores, .combine=rbind, .multicombine = TRUE, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
      variant.id0 <- if(b <= n.groups.percore_1) variant.id1[((b-1)*(n.groups.percore-1)+1):(b*(n.groups.percore-1))] else variant.id1[(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1-1)*n.groups.percore+1):(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1)*n.groups.percore)]
      n.groups <- length(variant.id0)
      nbatch.flush <- (n.groups-1) %/% nperbatch + 1
      for(i in 1:nbatch.flush) {
        # itmp.idx<-group.info$idx[((i-1)*nperbatch+1):min((i*nperbatch),n.groups)]
        ii<-match(variant.id0[((i-1)*nperbatch+1):min((i*nperbatch),n.groups)],group.info$snpid)
        itmp.idx<-group.info$idx[ii]
        dis<-group.info$pos-group.info[group.info$idx==itmp.idx[1],]$pos.start
        tmp.idx1<-group.info$idx[dis==min(dis[dis>=0])]
        dis<-group.info$pos-group.info[group.info$idx==itmp.idx[length(itmp.idx)],]$pos.end
        tmp.idx2<-group.info$idx[dis==min(dis[dis>=0])]
        tmp.idx<-group.info$idx[tmp.idx1:tmp.idx2]
        rm(dis)
        U.list <- V.list <- vector("list", n.cohort)
        variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
        for(j in 1:n.cohort) {
          tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
          if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
            U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
            # U.list[[j]]$weight <- U.list[[j]]$weight * U.list[[j]]$flip
            # U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$weight
            U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$flip
            tmp.V <- matrix(NA, sum(tmp.include), N.resampling[j])
            for(ij in 1:sum(tmp.include)) {
              if(U.list[[j]]$file[ij]!=current.cons[j]) {
                close(cons[[j]])
                current.cons[j] <- U.list[[j]]$file[ij]
                cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.", current.cons[j]), "rb")
                tmp.N.resampling <- readBin(cons[[j]], what = "integer", n = 1, size = 4)
                if(tmp.N.resampling != N.resampling[j]) stop(paste0("Error: N.resampling in ", meta.files.prefix[j], ".resample.", current.cons[j], " does not match that in ",meta.files.prefix[j], ".resample.1"))
                # current.lines[j] <- 1
              }
              current.cons[j] <- U.list[[j]]$file[ij]
              cons[[j]] <- file(paste0(meta.files.prefix[j], ".resample.", current.cons[j]), "rb")
              # if(U.list[[j]]$variant.idx[ij]!=current.lines[j]) seek(cons[[j]], where = 4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-current.lines[j]), origin = "current", rw = "read")
              ind <- 4+4*N.resampling[j]*(U.list[[j]]$variant.idx[ij]-1)
              seek(cons[[j]], where = ind, origin = "start", rw = "read")
              tmp.V[ij,] <- readBin(cons[[j]], what = "numeric", n = N.resampling[j], size = 4)
              # current.lines[j] <- U.list[[j]]$variant.idx[ij]+1
            }
            # V.list[[j]] <- tmp.V * U.list[[j]]$weight / sqrt(N.resampling[j])
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
        # rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices, N, Nmiss)
        rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices)
        if(sum(include) == 0) next
        variant.indices <- variant.indices[include]
        # out <- data.frame(SNP = group.info$SNP[variant.indices], chr = group.info$chr[variant.indices], pos = group.info$pos[variant.indices], ref = group.info$ref[variant.indices], alt = group.info$alt[variant.indices], N = N[include], missrate = (Nmiss/(N+Nmiss))[include], altfreq = AF[include])
        out <- group.info[match(itmp.idx,group.info$idx),]
        n.p <- length(variant.indices)
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
        n.batchi<-length(itmp.idx)
        tmp.weight <- MAF.weights.beta.fun(AF[include], MAF.weights.beta[1], MAF.weights.beta[2])
        if(use.minor.allele) tmp.weight[AF[include] > 0.5] <- -tmp.weight[AF[include] > 0.5]
        if(!is.null(cohort.group.idx)) tmp.weight <- rep(tmp.weight, n.cohort.groups)
        U <- U*tmp.weight*(gamma(MAF.weights.beta[1])*gamma(MAF.weights.beta[2])/gamma(MAF.weights.beta[1]+MAF.weights.beta[2])/sqrt(2))
        V <- V*tmp.weight*(gamma(MAF.weights.beta[1])*gamma(MAF.weights.beta[2])/gamma(MAF.weights.beta[1]+MAF.weights.beta[2])/sqrt(2))
        LDscore<-rep(NA,n.batchi)
        # ss<-rep(NA,n.batchi)
        # vv<-rep(NA,n.batchi)
        # LDscore<-sapply(as.list(seq(1,n.batchi)),function(x) L2(x))
        for (j in 1:n.batchi){
          jidx<-itmp.idx[j]
          if (jidx %in% variant.indices) {
            # jtmp.idx<-group.info[group.info$pos>group.info$pos.start[(i-1)*nperbatch+j] & group.info$pos<group.info$pos.end[(i-1)*nperbatch+j],]$idx
            jtmp.idx<-group.info[group.info$pos>=group.info$pos.start[ii[j]] & group.info$pos<=group.info$pos.end[ii[j]],]$idx
            jV<-V[variant.indices==jidx,]
            jV0<-match(jtmp.idx,variant.indices)
            jV0<-jV0[!is.na(jV0)]
            nvariants<-length(jV0)
            V0<-V[jV0,]
            jN<-N[variant.indices==jidx]
            LDscore[j]<- (1+1/(N.randomvec-2)+1/(jN-2))*rowSums(tcrossprod((tcrossprod(jV,V0))/(jN-1)))-nvariants/(N.randomvec-2)-nvariants/(jN-2)
            # ss[j]<-sum(jV)
            # vv[j]<-sum(rowSums(V0))
            rm(jtmp.idx,jV0,V0)
          }
        }
        out$LDscore<-LDscore
        # out$ss<-ss
        # out$vv<-vv
        all.out <- rbind(all.out, out)
      }
      return(all.out)
    }
  } else { # use a single core
    n.groups <- length(variant.id1)
    nbatch.flush <- (n.groups-1) %/% nperbatch + 1
    for(i in 1:nbatch.flush) {
      itmp.idx<-group.info$idx[((i-1)*nperbatch+1):min((i*nperbatch),n.groups)]
      dis<-group.info$pos-group.info[group.info$idx==itmp.idx[1],]$pos.start
      tmp.idx1<-group.info$idx[dis==min(dis[dis>=0])]
      dis<-group.info$pos-group.info[group.info$idx==itmp.idx[length(itmp.idx)],]$pos.end
      tmp.idx2<-group.info$idx[dis==min(dis[dis>=0])]
      tmp.idx<-group.info$idx[tmp.idx1:tmp.idx2]
      rm(dis)
      U.list <- V.list <- vector("list", n.cohort)
      variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
      for(j in 1:n.cohort) {
        tmp.scores <- scores[[j]][tmp.idx, , drop = FALSE]
        if(any(tmp.include <- !is.na(tmp.scores$SCORE))) {
          U.list[[j]] <- tmp.scores[tmp.include, , drop = FALSE]
          # U.list[[j]]$weight <- U.list[[j]]$weight * U.list[[j]]$flip
          # U.list[[j]]$SCORE <- U.list[[j]]$SCORE * U.list[[j]]$weight
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
          # V.list[[j]] <- tmp.V * U.list[[j]]$weight / sqrt(N.resampling[j])
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
      # rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices, N, Nmiss)
      rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices)
      if(sum(include) == 0) next
      variant.indices <- variant.indices[include]
      # out <- data.frame(SNP = group.info$SNP[variant.indices], chr = group.info$chr[variant.indices], pos = group.info$pos[variant.indices], ref = group.info$ref[variant.indices], alt = group.info$alt[variant.indices], N = N[include], missrate = (Nmiss/(N+Nmiss))[include], altfreq = AF[include])
      out <- group.info[match(itmp.idx,group.info$idx),]
      n.p <- length(variant.indices)
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
      n.batchi<-length(itmp.idx)
      tmp.weight <- MAF.weights.beta.fun(AF[include], MAF.weights.beta[1], MAF.weights.beta[2])
      if(use.minor.allele) tmp.weight[AF[include] > 0.5] <- -tmp.weight[AF[include] > 0.5]
      if(!is.null(cohort.group.idx)) tmp.weight <- rep(tmp.weight, n.cohort.groups)
      U <- U*tmp.weight*(gamma(MAF.weights.beta[1])*gamma(MAF.weights.beta[2])/gamma(MAF.weights.beta[1]+MAF.weights.beta[2])/sqrt(2))
      V <- V*tmp.weight*(gamma(MAF.weights.beta[1])*gamma(MAF.weights.beta[2])/gamma(MAF.weights.beta[1]+MAF.weights.beta[2])/sqrt(2))
      LDscore<-rep(NA,n.batchi)
      # ss<-rep(NA,n.batchi)
      # vv<-rep(NA,n.batchi)
      # LDscore<-sapply(as.list(seq(1,n.batchi)),function(x) L2(x))
      for (j in 1:n.batchi){
        jidx<-itmp.idx[j]
        if (jidx %in% variant.indices) {
          jtmp.idx<-group.info[group.info$pos>=group.info$pos.start[(i-1)*nperbatch+j] & group.info$pos<=group.info$pos.end[(i-1)*nperbatch+j],]$idx
          jV<-V[variant.indices==jidx,]
          jV0<-match(jtmp.idx,variant.indices)
          jV0<-jV0[!is.na(jV0)]
          nvariants<-length(jV0)
          V0<-V[jV0,]
          # jN<-N[variant.indices==jidx]
          jN<-N[jidx]
          LDscore[j]<- (1+1/(N.randomvec-2)+1/(jN-2))*rowSums(tcrossprod((tcrossprod(jV,V0))/(jN-1)))-nvariants/(N.randomvec-2)-nvariants/(jN-2)
          # ss[j]<-sum(jV)
          # vv[j]<-sum(rowSums(V0))
          rm(jtmp.idx,jV0,V0)
        }
      }
      out$LDscore<-LDscore
      # out$ss<-ss
      # out$vv<-vv
      all.out <- rbind(all.out, out)
    }
  }
  for(i in 1:n.cohort) close(cons[[i]])
  return(all.out)
  # return(V)
}

MAF.weights.beta.fun <- function(freq, beta1, beta2) {
  freq[freq > 0.5] <- 1 - freq[freq > 0.5]
  ifelse(freq <= 0, 0, dbeta(freq, beta1, beta2))
}
