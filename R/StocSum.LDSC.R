#' Create random vectors for a glmmkin object for LD Score regresion
#' @description Generate random vectors from multivariate normal distribution with mean 0 and covariance matrix P, in which \eqn P = I-1inv((t(1)1))t(1).
#' @param obj The glmmkin object.
#' @param Z NULL. %% A list of design matrices for the random effects. The length must match the number of variance components.
#' @param N.randomvec The number of random vectors to generate (default = 1000).
#' @param group.idx NULL. %%A length N index vector showing which observation belongs to which variance group, for heteroscedastic linear mixed models (default = NULL for homoscedastic linear mixed models).
#' @param cluster.idx NULL. %% A length N index vector showing which observation belongs to which cluster (default = NULL for no clusters).
#' @param robust FALSE. %% A logical switch: whether robust variance should be used (default = FALSE).
#' @return A list of class glmmkin.randomvec
#' \item{theta}{ set to be 1. %% Variance estimates, inherited from the glmmkin object.}
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
#' obj <- LDSC.glmmkin2randomvec(nullmod)
#' }
#' @keywords random vector
#' @export

LDSC.glmmkin2randomvec <- function(obj, Z = NULL, N.randomvec = 1000, group.idx = NULL, cluster.idx = NULL, robust = FALSE) {
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


####LDSC.stat is same as G.stat


#' Calculate LDscore
#' @description Use summary statistics and stochastic statistics from G.stat to estimate LDscore.
#' @param meta.file.prefix a character vector for prefix of intermediate files (*.sample.* and *.resample.*). 
#' @param n.files an integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix))).
#' @param N.randomvec the number of replicates to simulate the random vectors in \code{StocSum.LDSC.glmmkin2randomvec} (default = 1000).
#' @param MAF.range a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
#' @param MAF.weight.beta MAF.weight.beta = c(0.5, 0.5). %% a numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies. This internal minor allele frequency weight is multiplied by the external weight given by the group.file. To turn off internal minor allele frequency weight and only use the external weight given by the group.file, use c(1, 1) to assign flat weights (default = c(1, 25)). Applied to the combined samples.
#' @param missing.cutoff the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
#' @param wind.b a positive integer define the window size in bases (default = 1000000).
#' @param use.minor.allele a logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SKAT results, but Burden (as well as SKAT-O and hybrid test to combine the burden test and SKAT) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction. Use with caution, as major/minor alleles may flip in different cohorts. In that case, minor allele will be determined based on the allele frequency in the combined samples.
#' @param auto.flip a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
#' @param nperbatch an integer for how many SNPs to be included in a batch (default = 10000). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user’s system.
#' @param ncores a positive integer indicating the number of cores to be used in parallel computing (default = 1).
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
#' GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
#' GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
#' nullmod <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"))
#' obj <- glmmkin2randomvec(nullmod)
#' out.prefix <- "test"
#' gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
#' G.stat(obj, geno.file = gdsfile, meta.file.prefix = out.prefix, MAF.range=c(0,0.5), miss.cutoff = 1)
#' xxxxxxxxxxxxxxxxxxxxxxxxxxx
#' }
#' @keywords LD Score
#' @export

LDSC.win <- function(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), N.randomvec=1000, MAF.range = c(1e-3, 0.5), MAF.weights.beta = c(0.5, 0.5), miss.cutoff = 1, wind.b=1000000, use.minor.allele = FALSE, auto.flip = FALSE, nperbatch = 100000, ncores = 1)
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
