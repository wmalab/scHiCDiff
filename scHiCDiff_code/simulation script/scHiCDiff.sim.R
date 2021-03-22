
#' generate a differential chromatin simulation in scHi-C
#'
#' 
#' This function generate the counts of the simulated samples in 3 steps:
#' 1. merge the data from single cells to get a pseudo-bulk data
#' 2. generate a pair of pseudo-bulk HiC matrices with predefined differential 
#' interections (one rep for each condition,2 conditions)
#' 3. downsample the pseudo-bulk datasets to a series of simulated single-cell-like 
#' data for each condition

#' Input:
#' @param file.path  The pathway of single cell files. All scHi-C data used in simulation
#' should be stored in this pathway. Each scHi-C file is performed as three-column format
#' containing the first interacting region of the bin pair, the second interacting region 
#' of the bin pair and the interaction frequency of the bin pair.
#' @param fold.change double. The amount of fold change.
#' @param pDiff double. The probability that an interaction will be diffrential.
#' @param noise.prop double. Remove interactions that show
#'  a dispersion value larger than this value.
#'
#' Output: a list that contains the simulated replicates and the
#' matrix of the true DCI regions.
#' The list contains the following elements:
#' Hic1.sim : A list containing the simulated scHi-C matrices of the first condition.
#' Hic2.sim : A list containing the simulated scHi-C matrices of the second condition.
#' diff.sim : A sparceMatric containing the position of the differential interactions.
#'


scHiCDiff.sim <- function(file.path,fold.change,pDiff=0.01,noise.prop=0.9){
  
library(Matrix)
library(mvtnorm)
library(rasterVis)
library(gridExtra)
library(HiTC)
library(edgeR)
library(ggsci)
  library(HiCcompare)

#step1 merge the data from single cells to get a pseudo-bulk data
merge.data <- function(rawHiC.path,rawHiC.path2,noise.prop){
  count <- c()
  for(j in 1:length(rawHiC.path)){
    if(j==1){
      dat1 <- read.table(rawHiC.path[1])
      dat1 <- as.matrix(dat1)
      dat1 <- cbind(dat1[,1],dat1[,2],dat1[,3])
      chr.length1 <- max(dat1[,1:2])
      hic.table1 <- matrix(0,nrow=chr.length1,ncol=chr.length1)
      hic.table <- sparse2full(dat1)
      start <- as.numeric(colnames(hic.table)[1])
      end <- chr.length1
      hic.table1[start:end,start:end] <- hic.table
      count[j] <- sum(dat1[,3])
    }else{
      dat2 <- read.table(rawHiC.path[j])
      dat2 <- as.matrix(dat2)
      dat2 <- cbind(dat2[,1],dat2[,2],dat2[,3])
      chr.length2 <- max(dat2[,1:2])
      hic.table2 <- matrix(0,nrow=chr.length2,ncol=chr.length2)
      hic.table <- sparse2full(dat2)
      start <- as.numeric(colnames(hic.table)[1])
      end <- chr.length2
      hic.table2[start:end,start:end] <- hic.table
      count[j] <- sum(dat2[,3])
      if(chr.length1>chr.length2){
        hic.table1[1:chr.length2,1:chr.length2] <- hic.table1[1:chr.length2,1:chr.length2]+hic.table2
      }else{
        hic.table2[1:chr.length1,1:chr.length1] <- hic.table2[1:chr.length1,1:chr.length1]+hic.table1
        hic.table1 <- hic.table2
        chr.length1 <- chr.length2
      }
    }
  }
  hic.table1[lower.tri(hic.table1)] <- 0
  hic.merge.data <- Matrix(hic.table1,sparse=TRUE)
  hic.merge.data <- as(hic.merge.data,"dgCMatrix")
  
  
  #calculate p.prop
  p.prop1 <- count/sum(count)
  for(j in 1:length(rawHiC.path2)){
    if(j==1){
      dat1 <- read.table(rawHiC.path2[1])
      dat1 <- as.matrix(dat1)
      aa1 <- dat1[,1]*10000+dat1[,2]
    }else{
      dat2 <- read.table(rawHiC.path2[j])
      dat2 <- as.matrix(dat2)
      aa2 <- dat2[,1]*10000+dat2[,2]
      posi <- which(!aa2 %in% aa1)
      if(length(posi)!=0){
        dat0 <- matrix(0,ncol=2+j,nrow=nrow(dat1)+length(posi))
        dat0[1:nrow(dat1),1:(1+j)] <- dat1
        dat0[(nrow(dat1)+1):(nrow(dat1)+length(posi)),1:2] <- dat2[posi,1:2]
        aa1 <- dat0[,1]*10000+dat0[,2]
        posi1 <- match(aa2,aa1)
        dat0[posi1,(2+j)] <- dat2[,3]
        dat1 <- dat0
      }else{
        dat0 <- matrix(0,ncol=2+j,nrow=nrow(dat1))
        dat0[1:nrow(dat1),1:(1+j)] <- dat1
        posi1 <- match(aa2,aa1)
        dat0[posi1,(2+j)] <- dat2[,3]
        dat1<- dat0
      }
    }
  }
  
  
  for(j in 1:nrow(dat1)){
    dat1[j,-c(1:2)] <- noise.prop*dat1[j,-c(1,2)]/sum(dat1[j,-c(1,2)])+(1-noise.prop)*p.prop1
  }  
  
  colnames(dat1) <- c("bin1","bin2",1:length(rawHiC.path2))
  return(list(hic.merge.data=hic.merge.data,p.prop=dat1))
}


#step2 generate a pair of pseudo-bulk HiC matrices with predefined differential interections (one rep for each condition,2 conditions)
generate.simulation <- function(Hicmat,foldDiff,pDiff) {
    hic.smry <- summary(Hicmat)
    hic.smry <- subset(hic.smry, !is.na(x))
    if (nrow(hic.smry) == 0) {
      stop("It seems that the provided merged-matrix in empty")
    }
    
    npos <- nrow(hic.smry)
    counts <- matrix(0,nrow=npos,ncol=2)
    for (i in 1:npos) {
      counts[i, 1] <- rnbinom(1,mu=round(hic.smry$x[i]),size=1/0.001)
      counts[i, 2] <- rnbinom(1,mu=round(hic.smry$x[i]),size=1/0.001)
    }
    rownames(counts) <- paste0(hic.smry$i[1:npos], "_", hic.smry$j[1:npos])
    
    
    group <- factor(c(1, 2))
    id <-  sort(sample(npos,floor(pDiff*npos)))
    counts[id,1] <- foldDiff*counts[id,1]
    pos.i <- as.numeric(gsub("_\\d+", "",rownames(counts)))
    pos.j <- as.numeric(gsub("\\d+_", "",rownames(counts)))
    Hic1 <- list()
    pos <- which(group==1)
    for (i in 1:length(pos)){
      p <- pos[i]
      Hic1[[i]] <- sparseMatrix(i=pos.i,j=pos.j,x=counts[, p],dims=dim(Hicmat))
    }
    Hic2 <- list()
    pos <- which(group==2)
    for (i in 1:length(pos)){
      p <- pos[i]
      Hic2[[i]] <- sparseMatrix(i=pos.i,j=pos.j,x=counts[, p],dims=dim(Hicmat))
    }
    change.i <- as.numeric(gsub("_\\d+", "",rownames(counts[id,])))
    change.j <- as.numeric(gsub("\\d+_", "",rownames(counts[id,])))
    change.mat <- sparseMatrix(i= change.i,j= change.j,x=rep(1,length(change.i)),dims=dim(Hicmat))
    
  
    return(list(Hic1=Hic1,Hic2=Hic2,change.mat=change.mat))
  }


#step3 downsample the pseudo-bulk datasets to a series of simulated single-cell-like data for each condition
downsampling <- function(object,p.prop){
  n <- ncol(p.prop)-2
  dat1 <- object$Hic1[[1]]
  dat1 <- Matrix(dat1,sparse=FALSE)
  dat1 <- as(dat1,"matrix")
  colnames(dat1) <- 1:nrow(dat1)
  dat1 <- full2sparse(dat1)
  dat1.sim <- array(0,dim=c(n,nrow(dat1),3))
  
  for(i in 1:nrow(dat1)){
    dat1.sim[,i,1] <- rep(dat1$region1[i],n)
    dat1.sim[,i,2] <- rep(dat1$region2[i],n)
    prob.i <- p.prop[intersect(which(p.prop[,1]==dat1$region1[i]),which(p.prop[,2]==dat1$region2[i])),-c(1:2)]
    dat1.sim[,i,3] <- c(rmultinom(1,dat1$IF[i],prob=prob.i))
  }
  
  dat2 <- object$Hic2[[1]]
  dat2 <- Matrix(dat2,sparse=FALSE)
  dat2 <- as(dat2,"matrix")
  colnames(dat2) <- 1:nrow(dat2)
  dat2 <- full2sparse(dat2)
  dat2.sim <- array(0,dim=c(n,nrow(dat2),3))
  
  for(i in 1:nrow(dat2)){
    dat2.sim[,i,1] <- rep(dat2$region1[i],n)
    dat2.sim[,i,2] <- rep(dat2$region2[i],n)
    prob.i <- p.prop[intersect(which(p.prop[,1]==dat2$region1[i]),which(p.prop[,2]==dat2$region2[i])),-c(1:2)]
    dat2.sim[,i,3] <- c(rmultinom(1,dat2$IF[i],prob=prob.i))
  }
  return(list(Hic1.sim=dat1.sim,Hic2.sim=dat2.sim))
}


#simuluation process
filename <- dir(file.path)
filename1 <- c()
for(i in 1:length(filename)){
  filename1[i] <- paste(file.path,filename[i],sep="/")
}
step.11 <- merge.data(rawHiC.path=filename1,rawHiC.path2=filename1,noise.prop=noise.prop)
step.22 <- generate.simulation(step.11$hic.merge.data,pDiff=pDiff,foldDiff=fold.change)
step.33 <- downsampling(step.22,p.prop=step.11$p.prop)
Hic1.sim <- step.33$Hic1.sim
Hic2.sim <- step.33$Hic2.sim
diff.sim <- step.22$change.mat
diff.sim <- as(diff.sim,"matrix")
colnames(diff.sim) <- 1:nrow(diff.sim)
diff.sim <- full2sparse(diff.sim)

return(list(Hic1.sim=Hic1.sim,Hic2.sim=Hic2.sim,diff.sim=diff.sim))
}