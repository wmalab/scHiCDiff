#' scHiCDiff.NBH
#'
#' 
#' This function is used to detect differentially chromatin interactions(DCIs) between
#' two specified groups of cells in a normalized read counts matrix of single-cell Hi-C 
#' (scHi-C) data by assuming assuming the read counts follows negative binomial hurdle
#' (NBH) distribution. 
#' 

#' Input:
#' @param count.table  A non-negative  matrix of scHi-C normalized read counts.The rows of the 
#' matrix are bin pair and columns are samples/cells.
#' @param group  A vector of factor which mentions the two condition to be compared, corresponding 
#' to the columns in the count table.

#'
#' Output:
#' 
#' A data frame containing the differential chromatin interaction (DCI) analysis results, rows are 
#' bin pairs and columns contain the following items:
#' 
#' 
#' bin_1,bin_2:  The interacting region of the bin pair
#' mu_1, mu_2, theta_1, theta_2, pi_1, pi_2: MLE of the negative binomial hurdle distribution's
#' parameters of group 1 and group 2, where mu and theta represent the mean and dispersion estimate of
#' truncated negative binomial part, and pi denotes the estimate of percentange of zero count in NBH.
#' norm_total_mean_1, norm_total_mean_2: Mean of normalized read counts of group 1 and group 2.
#' norm_foldChange: norm_total_mean_1/norm_total_mean_2.
#' chi2LR1: Chi-square statistic for hypothesis testing of H0.
#' pvalue: P value of hypothesis testing of H0 (underlying whether a bin pair is a DCI).
#' pvalue.adj.FDR: Adjusted P value of H0's pvalue using Benjamini & Hochberg's method.
#' Remark: Record of abnormal program information.
#' 
#'







scHiCDiff.NBH <- function(count.table, group){
  library(pscl)
  library(VGAM)
  library(maxLik)
  library(countreg)
  binpair <- count.table[,c(1,2)]
  counts <- count.table[,-c(1,2)]
  row.names(counts) <- 1:nrow(counts)
  
  
  binpairNum <- nrow(counts)
  sampleNum <- ncol(counts)
  gc()
  
  counts_norm <- as.matrix(counts)
  counts_norm <- ceiling(counts_norm)
  
  totalMean_1 <- rowMeans(counts[, group == levels(group)[1]])
  totalMean_2 <- rowMeans(counts[, group == levels(group)[2]])
  foldChange <- totalMean_1/totalMean_2
  All_Mean_FC <- cbind(totalMean_1, totalMean_2, foldChange)
  
  # Memory management
  remove(counts, totalMean_1, totalMean_2, foldChange)
  counts_norm <- Matrix(counts_norm, sparse = TRUE)
  gc()
  
  
  # Function of testing homogeneity of two NBH populations
  CallDCI.NBH <- function(i){
    
    # Memory management
    if(i %% 100 == 0)
      gc()
    
    # Function input and output
    counts_1 <- counts_norm[i, group == levels(group)[1]]
    counts_2 <- counts_norm[i, group == levels(group)[2]]
    results_binpair <- data.frame(row.names = row.names(counts_norm)[i], mu_1 = NA, mu_2 = NA, theta_1 = NA, theta_2 = NA, pi_1 = NA, pi_2 = NA, norm_total_mean_1 = NA, norm_total_mean_2 = NA, norm_foldChange = NA, chi2LR1 = NA, pvalue = NA, pvalue.adj.FDR = NA, Remark = NA)
    
  
  
    # MLE of parameters of NBH counts_1
    if(sum(counts_1 == 0) > 0){
      if(sum(counts_1 == 0) == length(counts_1)){
        pi_1 <- 0
        mu_1 <- 0
        theta_1 <- 1
        loglike_1 <- 0
      }else{
        options(show.error.messages = FALSE)
        nbh_try <- try(hurdle(formula = counts_1 ~ 1 | 1, dist = "negbin"), silent=TRUE)
        if('try-error' %in% class(nbh_try)){
          print("MLE of NBH failed!");
          results_binpair[1,"Remark"] <- "NBH failed!"
          return(results_binpair)
        }else{
          pi_1 <- plogis(nbh_try$coefficients$zero);names(pi_1) <- NULL
          mu_1 <- exp(nbh_try$coefficients$count);names(mu_1) <- NULL
          theta_1 <- nbh_try$theta;names(theta_1) <- NULL
          loglike_1 <- nbh_try$loglik;names(loglike_1) <- NULL
        }
      }
    }else{
      op <- options(warn=2)
      tnb_try <- try(zerotrunc(formula = counts_1 ~ 1,dist="negbin"), silent=TRUE)
      options(op)
      if('try-error' %in% class(tnb_try)){
        print("MLE of TNB failed!");
        results_binpair[1,"Remark"] <- "TNB failed!"
        return(results_binpair)
      }else{
        tnb.coef <- coef(tnb_try)
        pi_1 <- 1
        mu_1 <- exp(tnb_try$coefficients);names(mu_1) <- NULL
        theta_1 <- tnb_try$theta;names(theta_1) <- NULL
        loglike_1 <- tnb_try$loglik;names(loglike_1) <- NULL
      }
    }
    
    # MLE of parameters of NBH counts_2
    if(sum(counts_2 == 0) > 0){
      if(sum(counts_2 == 0) == length(counts_1)){
        pi_2 <- 0
        mu_2 <- 0
        theta_2 <- 1
        loglike_2 <- 0
      }else{
        options(show.error.messages = FALSE)
        nbh_try <- try(hurdle(formula = counts_2 ~ 1 | 1, dist = "negbin"), silent=TRUE)
        if('try-error' %in% class(nbh_try)){
          print("MLE of NBH failed!");
          results_binpair[1,"Remark"] <- "NBH failed!"
          return(results_binpair)
        }else{
          pi_2 <- plogis(nbh_try$coefficients$zero);names(pi_2) <- NULL
          mu_2 <- exp(nbh_try$coefficients$count);names(mu_2) <- NULL
          theta_2 <- nbh_try$theta;names(theta_2) <- NULL
          loglike_2 <- nbh_try$loglik;names(loglike_2) <- NULL
          #prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
        }
      }
    }else{
      op <- options(warn=2)
      tnb_try <- try(zerotrunc(formula = counts_2 ~ 1,dist="negbin"), silent=TRUE)
      options(op)
      if('try-error' %in% class(tnb_try)){
        print("MLE of TNB failed!");
        results_binpair[1,"Remark"] <- "TNB failed!"
        return(results_binpair)
      }else{
        tnb.coef <- coef(tnb_try)
        pi_2 <- 1
        mu_2 <- exp(tnb_try$coefficients);names(mu_2) <- NULL
        theta_2 <- tnb_try$theta;names(theta_2) <- NULL
        loglike_2 <- tnb_try$loglik;names(loglike_2) <- NULL
      }
    }
    
    # Restricted MLE under H0 (MLE of c(counts_1, counts_2))
    if(sum(c(counts_1, counts_2) == 0) > 0){
      options(show.error.messages = FALSE)
      nbh_try <- try(hurdle(formula = c(counts_1, counts_2) ~ 1 | 1, dist = "negbin"), silent=TRUE)
      if('try-error' %in% class(nbh_try)){
        print("MLE of NBH failed!");
        results_binpair[1,"Remark"] <- "NBH failed!"
        return(results_binpair)
      }else{
        pi_res <- plogis(nbh_try$coefficients$zero);names(pi_res) <- NULL
        mu_res <- exp(nbh_try$coefficients$count);names(mu_res) <- NULL
        theta_res <- nbh_try$theta;names(theta_res) <- NULL
        loglike_res <- nbh_try$loglik;names(loglike_res) <- NULL
      }
      
    }else{
      op <- options(warn=2)
      tnb_try <- try(zerotrunc(formula = c(counts_1,counts_2) ~ 1,dist="negbin"), silent=TRUE)
      options(op)
      if('try-error' %in% class(tnb_try)){
        print("MLE of TNB failed!");
        results_binpair[1,"Remark"] <- "TNB failed!"
        return(results_binpair)
      }else{
        tnb.coef <- coef(tnb_try)
        pi_res <- 1
        mu_res <- exp(tnb_try$coefficients);names(mu_res) <- NULL
        theta_res <- tnb_try$theta;names(theta_res) <- NULL
        loglike_res <- tnb_try$loglik;names(loglike_res) <- NULL
      }
    }
    
    # Log likelihood functions
    logL <- function(counts_1, mu_1, theta_1, pi_1, counts_2, mu_2, theta_2, pi_2){
      logL_1 <- sum(dhnbinom(counts_1,mu=mu_1, theta=theta_1, pi=pi_1, log = TRUE))
      logL_2 <- sum(dhnbinom(counts_2,mu=mu_2, theta=theta_2, pi=pi_2, log = TRUE))
      logL <- logL_1 + logL_2
      logL
    }
    
    
    # LRT test of H0
    chi2LR1 <- 2 *(logL(counts_1, mu_1, theta_1, pi_1, counts_2, mu_2, theta_2, pi_2) - logL(counts_1, mu_res, theta_res, pi_res, counts_2, mu_res, theta_res, pi_res))
    pvalue <- 1 - pchisq(chi2LR1, df = 3)
    
    #output
    results_binpair[1,"mu_1"] <- mu_1
    results_binpair[1,"mu_2"] <- mu_2
    results_binpair[1,"theta_1"] <- theta_1
    results_binpair[1,"theta_2"] <- theta_2
    results_binpair[1,"pi_1"] <- 1-pi_1
    results_binpair[1,"pi_2"] <- 1-pi_2
    results_binpair[1,"norm_total_mean_1"] <- mean(counts_1)
    results_binpair[1,"norm_total_mean_2"] <- mean(counts_2)
    results_binpair[1,"norm_foldChange"] <- results_binpair[1,"norm_total_mean_1"] / results_binpair[1,"norm_total_mean_2"]
    results_binpair[1,"chi2LR1"] <- chi2LR1
    results_binpair[1,"pvalue"] <- pvalue
    
    # Return results_binpair
    return(results_binpair)
  }
  
  
  # Call DCI by bin pair
  results <- matrix(data=NA, nrow = binpairNum, ncol = 15, dimnames = list(row.names(counts_norm), c("bin1","bin2","mu_1", "mu_2","theta_1", "theta_2", "pi_1", "pi_2", "norm_total_mean_1", "norm_total_mean_2", "norm_foldChange", "chi2LR1", "pvalue", "pvalue.adj.FDR","Remark")))
  results <- as.data.frame(results)
  for(i in 1:binpairNum){
    cat("\r",paste0("scHiCDiff.NBH is analyzing ", i," of ",binpairNum," bin pairs"))
    results[i,1:2] <- binpair[i,]
    results[i,3:15] <- CallDCI.NBH(i)
  }
  
  
  # Format output results
  results[,"pvalue.adj.FDR"] <- p.adjust(results[,"pvalue"], method="fdr")
  results <- results[order(results[,"chi2LR1"], decreasing = TRUE),]
  
  return(results)
  
}




