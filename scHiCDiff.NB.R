scHiCDiff.NB <-
function(count.table, group){
    library(pscl)
    library(gamlss)
    library(VGAM)
    library(maxLik)
    library(countreg)
    library(Matrix)
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
    
    
    # Function of testing homogeneity of two NB populations
    CallDCI.NB <- function(i){
        if(i %% 100 == 0)
            gc()
        
        # Function input and output
        counts_1 <- counts_norm[i, group == levels(group)[1]]
        counts_2 <- counts_norm[i, group == levels(group)[2]]
        results_binpair <- data.frame(row.names = row.names(counts_norm)[i], mu_1 = NA, mu_2 = NA, size_1 = NA, size_2 = NA, norm_total_mean_1 = NA, norm_total_mean_2 = NA, norm_foldChange = NA, chi2LR1 = NA, pvalue = NA, pvalue.adj.FDR = NA, Remark = NA)
        
        # Log likelihood functions
        logL <- function(counts_1, mu_1, size_1, counts_2, mu_2, size_2){
            logL_1 <- sum(dnbinom(counts_1,mu=mu_1, size=size_1, log = TRUE))
            logL_2 <- sum(dnbinom(counts_2,mu=mu_2, size=size_2, log = TRUE))
            logL <- logL_1 + logL_2
            logL
        }
        
        # MLE of parameters of NB counts_1
        if(sum(counts_1 == 0) == length(counts_1)){
            mu_1 <- 0
            size_1 <- 1
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
                            mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
                            size_1 <- nb_1$theta;names(size_1) <- NULL
                        }
                    }else{
                        nb_1 <- nb_try_again
                        mu_1 <- exp(nb_1@coef["logmu"]);names(mu_1) <- NULL
                        size_1 <- 1/nb_1@coef["invk"];names(size_1) <- NULL
                    }
                }else{
                    nb_1 <- nb_try_twice
                    mu_1 <- nb_1$estimate["mu"];names(mu_1) <- NULL
                    size_1 <- nb_1$estimate["size"];names(size_1) <- NULL
                }
            }else{
                nb_1 <- nb_try
                mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
                size_1 <- nb_1$theta;names(size_1) <- NULL
            }
        }
        
        
        # MLE of parameters of NB counts_2
        if(sum(counts_2 == 0) == length(counts_2)){
            mu_2 <- 0
            size_2 <- 1
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
                            mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
                            size_2 <- nb_2$theta;names(size_2) <- NULL
                        }
                    }else{
                        nb_2 <- nb_try_again
                        mu_2 <- exp(nb_2@coef["logmu"]);names(mu_2) <- NULL
                        size_2 <- 1/nb_2@coef["invk"];names(size_2) <- NULL
                    }
                }else{
                    nb_2 <- nb_try_twice
                    mu_2 <- nb_2$estimate["mu"];names(mu_2) <- NULL
                    size_2 <- nb_2$estimate["size"];names(size_2) <- NULL
                }
            }else{
                nb_2 <- nb_try
                mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
                size_2 <- nb_2$theta;names(size_2) <- NULL
            }
        }
        
        
        
        # Restricted MLE under H0 (MLE of c(counts_1, counts_2))
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
                        mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
                        size_res <- nb_res$theta;names(size_res) <- NULL
                        
                    }
                }else{
                    nb_res <- nb_try_again
                    mu_res <- exp(nb_res@coef["logmu"]);names(mu_res) <- NULL
                    size_res <- 1/nb_res@coef["invk"];names(size_res) <- NULL
                    
                }
            }else{
                nb_res <- nb_try_twice
                mu_res <- nb_res$estimate["mu"];names(mu_res) <- NULL
                size_res <- nb_res$estimate["size"];names(size_res) <- NULL
                
            }
        }else{
            nb_res <- nb_try
            mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
            size_res <- nb_res$theta;names(size_res) <- NULL
        }
        
        
        
        # # LRT test of H0
        chi2LR1 <- 2 *(logL(counts_1, mu_1, size_1, counts_2, mu_2, size_2) - logL(counts_1, mu_res, size_res, counts_2, mu_res, size_res))
        pvalue <- 1 - pchisq(chi2LR1, df = 2)
        
        # Format output
        
        results_binpair[1,"mu_1"] <- mu_1
        results_binpair[1,"mu_2"] <- mu_2
        results_binpair[1,"size_1"] <- size_1
        results_binpair[1,"size_2"] <- size_2
        results_binpair[1,"norm_total_mean_1"] <- mean(counts_1)
        results_binpair[1,"norm_total_mean_2"] <- mean(counts_2)
        results_binpair[1,"norm_foldChange"] <- results_binpair[1,"norm_total_mean_1"] / results_binpair[1,"norm_total_mean_2"]
        results_binpair[1,"chi2LR1"] <- chi2LR1
        results_binpair[1,"pvalue"] <- pvalue
        
        # Return results_binpair
        return(results_binpair)
    }
    
    # Call DCI by bin pair
    results <- matrix(data=NA, nrow = binpairNum, ncol = 13, dimnames = list(row.names(counts_norm), c("bin1","bin2","mu_1", "mu_2","size_1", "size_2", "norm_total_mean_1", "norm_total_mean_2", "norm_foldChange", "chi2LR1", "pvalue", "pvalue.adj.FDR","Remark")))
    results <- as.data.frame(results)
    for(i in 1:binpairNum){
        cat("\r",paste0("scHiCDiff.NB is analyzing ", i," of ",binpairNum," bin pairs"))
        results[i,1:2] <- binpair[i,]
        results[i,3:13] <- CallDCI.NB(i)
    }
    
    
    # Format output results
    results[,"pvalue.adj.FDR"] <- p.adjust(results[,"pvalue"], method="fdr")
    results <- results[order(results[,"chi2LR1"], decreasing = TRUE),]
    
    return(results)
}
