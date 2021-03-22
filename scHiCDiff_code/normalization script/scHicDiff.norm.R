#' scHiCNorm.adjust
#' 
#'
#' This function is firstly processed by scHiCNorm with NBH option to eliminate
#' the systematic biases in single-cell Hi-C data and then is divided by size factors
#' to account for the cell-specific genomic distance effect. The first part code copied from
#' scHiCNorm paper directly and is available at http://dna.cs.miami.edu/scHiCNorm.
#' 
#' 
#' 
#' Input:
#' @param bias.info.path  The pathway of the three local features (effective length, 
#' GC content and mappability of fragment ends) of all bins. The generation of these
#' items is available at http://dna.cs.miami.edu/scHiCNorm.
#' @param dat_HiC A N*N scHi-C matrix.
#'
#'
#' Output:
#' A data frame containing the normalized data. Rows are bin pairs and columns contain the 
#' following items:
#' 
#' 
#' region_1, region_2:  The interacting region of the bin pair.
#' IF: The final normalized interaction frequency of the bin pair.
#' genom.distance: The genomic distance between the two bins.
#' size.factor: The genomic-distance size factor of the bin pair.
#'




scHiCNorm.adjust <- function(bias.info.path,dat_HiC){
library(HiCcompare)

scHicNorm <- function(bias.info.path,dat_HiC,method.num=6,diag.log=FALSE){
  # load library
  library(pscl)
  library(MASS)
  library(HiCcompare)
  # read data 
  dat_feature <- read.table(bias.info.path)
  names(dat_feature) <- c("binStart", "binEnd", "density", "GC", "map")
  if (nrow(dat_feature) != nrow(dat_HiC)) {
    stop("The number of rows from the first two input files is not the same.\n", call. = FALSE)
  }
  chr.length <- nrow(dat_feature)
  # determine the method
  methods <- c("Poisson", "NB", "ZIP", "ZINB", "PH", "NBH")
  method <- methods[method.num]
  
  
  # tidy data frame: remove rows
  rows_del_feature <- which(dat_feature$density == 0 | dat_feature$GC == 0 | dat_feature$map == 0)
  rows_del_HiC <- which(apply(dat_HiC, 1, sum) == 0)
  rows_del <- c(rows_del_feature, rows_del_HiC)
  rows_del <- rows_del[!duplicated(rows_del)]
  dat_feature <- dat_feature[-rows_del,]
  dat_HiC <- dat_HiC[-rows_del,-rows_del]
  
  if (length(rows_del) > 0) {
    cat("The following rows are removed: ", rows_del, ".\n")
  }
  
  mat_HiC <- as.matrix(dat_HiC)
  vec_HiC <- mat_HiC[upper.tri(mat_HiC, diag=diag.log)]
  
  # get the feature matrix and Z-score
  mat_density <- as.matrix(log(dat_feature[, 3] %o% dat_feature[, 3]))
  mat_GC      <- as.matrix(log(dat_feature[, 4] %o% dat_feature[, 4]))
  mat_map     <- as.matrix(log(dat_feature[, 5] %o% dat_feature[, 5]))
  mat_density <- (mat_density - mean(c(mat_density))) / sd(c(mat_density))
  mat_GC      <- (mat_GC - mean(c(mat_GC))) / sd(c(mat_GC))
  
  vec_density <- mat_density[upper.tri(mat_density, diag=diag.log)]
  vec_GC      <- mat_GC[upper.tri(mat_GC, diag=diag.log)]
  vec_map     <- mat_map[upper.tri(mat_map, diag=diag.log)]
  
  fit_model <- NULL
  
  if (method == "Poisson") {
    cat("Using Poisson. \n")
    fit_model <- glm(vec_HiC ~ vec_density + vec_GC + offset(vec_map), family="poisson")  
  } else if (method == "NB") { 
    cat("Using Negative Binomial. \n")
    fit_model <- glm.nb(vec_HiC ~ vec_density + vec_GC + offset(vec_map))  
  } else if (method == "ZIP") {
    cat("Using Zero-inflated Poisson. \n")
    fit_model <- zeroinfl(vec_HiC ~ vec_density + vec_GC + offset(vec_map) | vec_density + vec_GC + vec_map, dist="poisson")
  } else if (method == "ZINB") {
    cat("Using Zero-inflated Negative Binomial. \n")
    fit_model <- zeroinfl(vec_HiC ~ vec_density + vec_GC + offset(vec_map) | vec_density + vec_GC + vec_map, dist="negbin")
  } else if (method == "PH") {
    cat("Using Poisson Hurdle. \n")
    fit_model <- hurdle(vec_HiC ~ vec_density + vec_GC + offset(vec_map) | vec_density + vec_GC + vec_map)
  } else if (method == "NBH") {
    cat("Using Negative Ninomial Hurdle. \n")
    fit_model <- hurdle(vec_HiC ~ vec_density + vec_GC + offset(vec_map) | vec_density + vec_GC + vec_map, dist="negbin")
  } else {
    stop("Please select a valid method.\n", call. = FALSE)
  }
  
  # output normalized matrix
  fitted_vals <- fit_model$fitted.values
  fitted_vals_mat <- matrix(0, nrow(mat_HiC), ncol(mat_HiC))
  fitted_vals_mat[upper.tri(fitted_vals_mat, diag=diag.log)] = fitted_vals
  fitted_vals_mat <- t(fitted_vals_mat)
  fitted_vals_mat[upper.tri(fitted_vals_mat, diag=diag.log)] = fitted_vals
  res <- round(mat_HiC / fitted_vals_mat, 4)
  if(diag.log==FALSE){
    diag(res) <- 0}
  rows.remain <- 1:chr.length
  rows.remain <- rows.remain[-rows_del]
  res1 <- matrix(0, ncol=chr.length,nrow=chr.length)
  res1[rows.remain,rows.remain] <- res
  
  return(res1)
}

 res1 <- scHicNorm(bias.info.path,dat_HiC)
 if(sum(res1)==0){
   stop("The normalization failed.\n", call. = FALSE)
 }else{
   colnames(res1) <- 1:nrow(res1)
   res2 <- full2sparse(res1)
   res2 <- cbind(res2,genom.distance=abs(res2$region1-res2$region2),size.factor=rep(0,nrow(res2)))
   gd.uniq <- sort(unique(res2$genom.distance))
   size.factor <- c()
   for(j in 1:length(gd.uniq)){
     size.factor[j] <- median(res2$IF[which(res2$genom.distance==gd.uniq[j])])
     res2$IF[which(res2$genom.distance==gd.uniq[j])] <- res2$IF[which(res2$genom.distance==gd.uniq[j])]/size.factor[j]
     res2$size.factor[which(res2$genom.distance==gd.uniq[j])] <- rep(size.factor[j],length(which(res2$genom.distance==gd.uniq[j])))
   }
 }


 return(res2)
 

}
