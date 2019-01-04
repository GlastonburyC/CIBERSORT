# CIBERSORT R script v1.02 (last updated 02-02-2015)
# modified to speed up dataloading.
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#
#       Options:
#       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii) QN = Quantile normalization of input mixture (default = TRUE)
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


#dependencies
library(e1071)
library(parallel)
library(preprocessCore)
library(data.table) # CG added to speed up file loading time (x100s)
#Core algorithm
CoreAlg <- function(X, y){

    #try different values of nu
    svn_itor <- 1:3

    res <- function(i){
        if(i==1){nus <- 0.25}
        if(i==2){nus <- 0.5}
        if(i==3){nus <- 0.75}
        model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
        model
    }

    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)

    nusvm <- rep(0,svn_itor)
    corrv <- rep(0,svn_itor)

    #do cibersort
    t <- 1
    while(t <= svn_itor) {
        weights = t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights<0)]<-0
        w<-weights/sum(weights)
        u <- sweep(X,MARGIN=2,w,'*')
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- cor(k, y)
        t <- t + 1
    }

    #pick best model
    rmses <- nusvm
    mn <- which.min(rmses)
    model <- out[[mn]]

    #get and normalize coefficients
    q <- t(model$coefs) %*% model$SV
    q[which(q<0)]<-0
    w <- (q/sum(q))

    mix_rmse <- rmses[mn]
    mix_r <- corrv[mn]

    newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#do permutations
doPerm <- function(perm, X, Y){
    itor <- 1
    Ylist <- as.list(data.matrix(Y))
    dist <- matrix()

    while(itor <= perm){
        #print(itor)

        #random mixture
        yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

        #standardize mixture
        yr <- (yr - mean(yr)) / sd(yr)

        #run CIBERSORT core algorithm
        result <- CoreAlg(X, yr)

        mix_r <- result$mix_r

        #store correlation
        if(itor == 1) {dist <- mix_r}
        else {dist <- rbind(dist, mix_r)}

        itor <- itor + 1
    }
    newList <- list("dist" = dist)
}

#main function
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){

    #read in data
    X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
    Y <- fread(mixture_file, header=F, sep="\t")
    Y <- as.data.frame(Y)
    namesY = as.matrix(Y[,2])
    X <- data.matrix(X)
    Y <- data.matrix(Y)
    rownames(Y) = make.names(namesY, unique=TRUE)
    P <- perm #number of permutations

    #anti-log if max < 50 in mixture file
    if(max(Y) < 50) {Y <- 2^Y}

    #quantile normalization of mixture file
    if(QN == TRUE){
        tmpc <- colnames(Y)
        tmpr <- rownames(Y)
        Y <- normalize.quantiles(Y)
        colnames(Y) <- tmpc
        rownames(Y) <- tmpr
    }

    #intersect genes
    Xgns <- row.names(X)
    Ygns <- rownames(Y)
    YintX <- Ygns %in% Xgns
    Y <- Y[YintX,]
    XintY <- Xgns %in% row.names(Y)
    X <- X[XintY,]

    #standardize sig matrix
    X <- (X - mean(X)) / sd(as.vector(X))

    #empirical null distribution of correlation coefficients
    if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

    #print(nulldist)

    header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
    #print(header)

    output <- matrix()
    itor <- 1
    mixtures <- dim(Y)[2]
    pval <- 9999

    #iterate through mixtures
    while(itor <= mixtures){

        y <- Y[,itor]

        #standardize mixture
        y <- (y - mean(y)) / sd(y)

        #run SVR core algorithm
        result <- CoreAlg(X, y)

        #get results
        w <- result$w
        mix_r <- result$mix_r
        mix_rmse <- result$mix_rmse

        #calculate p-value
        if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

        #print output
        out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
        if(itor == 1) {output <- out}
        else {output <- rbind(output, out)}

        itor <- itor + 1

    }

    #save results
    write.table(rbind(header,output), file="CIBERSORT-Results.array_antiLOG_comparison.txt", sep="\t", row.names=F, col.names=F, quote=F)

    #return matrix object containing all results
    obj <- rbind(header,output)
    obj <- obj[,-1]
    obj <- obj[-1,]
    obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
    rownames(obj) <- colnames(Y)
    colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
    head(obj)
}
