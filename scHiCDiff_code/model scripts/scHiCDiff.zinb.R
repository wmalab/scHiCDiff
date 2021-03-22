#' scHiCDiff.ZINB
#'
#' DEsingle(Miao et al., 2018) utilized a zero-inflated Negative Binomial regression model to 
#' analyze differential expression in scRNA-seq data. Inspired by this paper, we also tried to 
#' apply zero-inflated Negative Binomial model in scHi-C data to eliminate the extreme sparisity
#' problem of scHi-C data.
#'
#'
#' This function is used to detect differentially chromatin interactions(DCIs) between
#' two specified groups of cells in a normalized read counts matrix of single-cell Hi-C 
#' (scHi-C) data by assuming assuming the read counts follows zero-inflated negative binomial
#' (ZINB) distribution. 
#' DEsingle function source code(https://github.com/miaozhun/DEsingle/blob/master/R/DEsingle.R) 
#' provide useful reference to this function writing.
#' 
#' 
#' 
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
#' mu_1, mu_2, theta_1, theta_2, pi_1, pi_2: MLE of the zero-inflated negative binomial distribution's
#' parameters of group 1 and group 2, where mu and theta represent the mean and dispersion estimate of
#' negative binomial part, and pi denotes the estiamte of percentange of exact zero count in ZINB.
#' norm_total_mean_1, norm_total_mean_2: Mean of normalized read counts of group 1 and group 2.
#' norm_foldChange: norm_total_mean_1/norm_total_mean_2.
#' chi2LR1: Chi-square statistic for hypothesis testing of H0.
#' pvalue: P value of hypothesis testing of H0 (underlying whether a bin pair is a DCI).
#' pvalue.adj.FDR: Adjusted P value of H0's pvalue using Benjamini & Hochberg's method.
#' Remark: Record of abnormal program information.
#' 
#'

scHiCDiff.ZINB <- function(count.table, group){
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
  remove(counts)
  counts_norm <- Matrix(counts_norm, sparse = TRUE)
  gc()
  
  
  # Function of testing homogeneity of two ZINB populations
  CallDCI.ZINB <- function(i){
    if(i %% 100 == 0)
      gc()
    
    # Function input and output
    counts_1 <- counts_norm[i, group == levels(group)[1]]
    counts_2 <- counts_norm[i, group == levels(group)[2]]
    results_binpair <- data.frame(row.names = row.names(counts_norm)[i], mu_1 = NA, mu_2 = NA, theta_1 = NA, theta_2 = NA, pi_1 = NA, pi_2 = NA, norm_total_mean_1 = NA, norm_total_mean_2 = NA, norm_foldChange = NA, chi2LR1 = NA, pvalue = NA, pvalue.adj.FDR = NA, Remark = NA)
    
    
    # MLE of parameters of ZINB counts_1
    if(sum(counts_1 == 0) > 0){
      if(sum(counts_1 == 0) == length(counts_1)){
        pi_1 <- 1
        mu_1 <- 0
        theta_1 <- 1
      }else{
        options(show.error.messages = FALSE)
        zinb_try <- try(gamlssML(counts_1, family="ZINBI"), silent=TRUE)
        options(show.error.messages = TRUE)
        if('try-error' %in% class(zinb_try)){
          zinb_try_twice <- try(zeroinfl(formula = counts_1 ~ 1 | 1, dist = "negbin"), silent=TRUE)
          if('try-error' %in% class(zinb_try_twice)){
            print("MLE of ZINB failed!");
            results_binpair[1,"Remark"] <- "ZINB failed!"
            return(results_binpair)
          }else{
            zinb_1 <- zinb_try_twice
            pi_1 <- plogis(zinb_1$coefficients$zero);names(pi_1) <- NULL
            mu_1 <- exp(zinb_1$coefficients$count);names(mu_1) <- NULL
            theta_1 <- zinb_1$theta;names(theta_1) <- NULL
          }
        }else{
          zinb_1 <- zinb_try
          pi_1 <- zinb_1$nu;names(pi_1) <- NULL
          mu_1 <- zinb_1$mu;names(mu_1) <- NULL
          theta_1 <- 1/zinb_1$sigma;names(theta_1) <- NULL
        }
      }
    }else{
      op <- options(warn=2)
      nb_try <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
      options(op)
      if('try-error' %in% class(nb_try)){
        nb_try_twice <- try(fitdistr(counts_1, "Negative Binomial"), silent=TRUE)
        if('try-error' %in% class(nb_try_twice)){
          nb_try_again <- try(mle2(counts_1~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_1), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
          if('try-error' %in% class(nb_try_again)){
            nb_try_fourth <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
            if('try-error' %in% class(nb_try_fourth)){
              print("MLE of NB failed!");
              results_binpair[1,"Remark"] <- "NB failed!"
              return(results_binpair)
            }else{
              nb_1 <- nb_try_fourth
              pi_1 <- 0
              mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
              theta_1 <- nb_1$theta;names(theta_1) <- NULL
            }
          }else{
            nb_1 <- nb_try_again
            pi_1 <- 0
            mu_1 <- exp(nb_1@coef["logmu"]);names(mu_1) <- NULL
            theta_1 <- 1/nb_1@coef["invk"];names(theta_1) <- NULL
          }
        }else{
          nb_1 <- nb_try_twice
          pi_1 <- 0
          mu_1 <- nb_1$estimate["mu"];names(mu_1) <- NULL
          theta_1 <- nb_1$estimate["size"];names(theta_1) <- NULL
        }
      }else{
        nb_1 <- nb_try
        pi_1 <- 0
        mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
        theta_1 <- nb_1$theta;names(theta_1) <- NULL
      }
    }
    
    # MLE of parameters of ZINB counts_2
    if(sum(counts_2 == 0) > 0){
      if(sum(counts_2 == 0) == length(counts_2)){
        pi_2 <- 1
        mu_2 <- 0
        theta_2 <- 1
        
      }else{
        options(show.error.messages = FALSE)
        zinb_try <- try(gamlssML(counts_2, family="ZINBI"), silent=TRUE)
        options(show.error.messages = TRUE)
        if('try-error' %in% class(zinb_try)){
          zinb_try_twice <- try(zeroinfl(formula = counts_2 ~ 1 | 1, dist = "negbin"), silent=TRUE)
          if('try-error' %in% class(zinb_try_twice)){
            print("MLE of ZINB failed!");
            results_binpair[1,"Remark"] <- "ZINB failed!"
            return(results_binpair)
          }else{
            zinb_2 <- zinb_try_twice
            pi_2 <- plogis(zinb_2$coefficients$zero);names(pi_2) <- NULL
            mu_2 <- exp(zinb_2$coefficients$count);names(mu_2) <- NULL
            theta_2 <- zinb_2$theta;names(theta_2) <- NULL
            
          }
        }else{
          zinb_2 <- zinb_try
          pi_2 <- zinb_2$nu;names(pi_2) <- NULL
          mu_2 <- zinb_2$mu;names(mu_2) <- NULL
          theta_2 <- 1/zinb_2$sigma;names(theta_2) <- NULL
         
        }
      }
    }else{
      op <- options(warn=2)
      nb_try <- try(glm.nb(formula = counts_2 ~ 1), silent=TRUE)
      options(op)
      if('try-error' %in% class(nb_try)){
        nb_try_twice <- try(fitdistr(counts_2, "Negative Binomial"), silent=TRUE)
        if('try-error' %in% class(nb_try_twice)){
          nb_try_again <- try(mle2(counts_2~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_2), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
          if('try-error' %in% class(nb_try_again)){
            nb_try_fourth <- try(glm.nb(formula = counts_2 ~ 1), silent=TRUE)
            if('try-error' %in% class(nb_try_fourth)){
              print("MLE of NB failed!");
              results_binpair[1,"Remark"] <- "NB failed!"
              return(results_binpair)
            }else{
              nb_2 <- nb_try_fourth
              pi_2 <- 0
              mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
              theta_2 <- nb_2$theta;names(theta_2) <- NULL
            }
          }else{
            nb_2 <- nb_try_again
            pi_2 <- 0
            mu_2 <- exp(nb_2@coef["logmu"]);names(mu_2) <- NULL
            theta_2 <- 1/nb_2@coef["invk"];names(theta_2) <- NULL
          }
        }else{
          nb_2 <- nb_try_twice
          pi_2 <- 0
          mu_2 <- nb_2$estimate["mu"];names(mu_2) <- NULL
          theta_2 <- nb_2$estimate["size"];names(theta_2) <- NULL
        }
      }else{
        nb_2 <- nb_try
        pi_2 <- 0
        mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
        theta_2 <- nb_2$theta;names(theta_2) <- NULL
      }
    }
    
    # Restricted MLE under H0 (MLE of c(counts_1, counts_2))
    if(sum(c(counts_1, counts_2) == 0) > 0){
      options(show.error.messages = FALSE)
      zinb_try <- try(gamlssML(c(counts_1, counts_2), family="ZINBI"), silent=TRUE)
      options(show.error.messages = TRUE)
      if('try-error' %in% class(zinb_try)){
        zinb_try_twice <- try(zeroinfl(formula = c(counts_1, counts_2) ~ 1 | 1, dist = "negbin"), silent=TRUE)
        if('try-error' %in% class(zinb_try_twice)){
          print("MLE of ZINB failed!");
          results_binpair[1,"Remark"] <- "ZINB failed!"
          return(results_binpair)
        }else{
          zinb_res <- zinb_try_twice
          pi_res <- plogis(zinb_res$coefficients$zero);names(pi_res) <- NULL
          mu_res <- exp(zinb_res$coefficients$count);names(mu_res) <- NULL
          theta_res <- zinb_res$theta;names(theta_res) <- NULL
        }
      }else{
        zinb_res <- zinb_try
        pi_res <- zinb_res$nu;names(pi_res) <- NULL
        mu_res <- zinb_res$mu;names(mu_res) <- NULL
        theta_res <- 1/zinb_res$sigma;names(theta_res) <- NULL
      }
    }else{
      op <- options(warn=2)
      nb_try <- try(glm.nb(formula = c(counts_1, counts_2) ~ 1), silent=TRUE)
      options(op)
      if('try-error' %in% class(nb_try)){
        nb_try_twice <- try(fitdistr(c(counts_1, counts_2), "Negative Binomial"), silent=TRUE)
        if('try-error' %in% class(nb_try_twice)){
          nb_try_again <- try(mle2(c(counts_1, counts_2)~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(c(counts_1, counts_2)), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
          if('try-error' %in% class(nb_try_again)){
            nb_try_fourth <- try(glm.nb(formula = c(counts_1, counts_2) ~ 1), silent=TRUE)
            if('try-error' %in% class(nb_try_fourth)){
              print("MLE of NB failed!");
              results_binpair[1,"Remark"] <- "NB failed!"
              return(results_binpair)
            }else{
              nb_res <- nb_try_fourth
              pi_res <- 0
              mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
              theta_res <- nb_res$theta;names(theta_res) <- NULL
            }
          }else{
            nb_res <- nb_try_again
            pi_res <- 0
            mu_res <- exp(nb_res@coef["logmu"]);names(mu_res) <- NULL
            theta_res <- 1/nb_res@coef["invk"];names(theta_res) <- NULL
          }
        }else{
          nb_res <- nb_try_twice
          pi_res <- 0
          mu_res <- nb_res$estimate["mu"];names(mu_res) <- NULL
          theta_res <- nb_res$estimate["size"];names(theta_res) <- NULL
        }
      }else{
        nb_res <- nb_try
        pi_res <- 0
        mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
        theta_res <- nb_res$theta;names(theta_res) <- NULL
      }
    }
    
    # Log likelihood functions
    logL <- function(counts_1, mu_1, theta_1 ,pi_1, counts_2, mu_2, theta_2, pi_2){
      logL_1 <- sum(dzinegbin(counts_1, size = theta_1, prob = theta_1/(theta_1 + mu_1), pstr0 = pi_1, log = TRUE))
      logL_2 <- sum(dzinegbin(counts_2, size = theta_2, prob = theta_2/(theta_2 + mu_2), pstr0 = pi_2, log = TRUE))
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
    results_binpair[1,"pi_1"] <- pi_1
    results_binpair[1,"pi_2"] <- pi_2
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
    cat("\r",paste0("scHiCDiff.ZINB is analyzing ", i," of ",binpairNum," bin pairs"))
    results[i,1:2] <- binpair[i,]
    results[i,3:15] <- CallDCI.ZINB(i)
  }
  
  
  # Format output results
  results[,"pvalue.adj.FDR"] <- p.adjust(results[,"pvalue"], method="fdr")
  results <- results[order(results[,"chi2LR1"], decreasing = TRUE),]
  
  return(results)
  
}







