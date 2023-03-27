scHiCDiff.KS <-
function(count.table,group){
    library(twosamples)
    binpair <- count.table[,c(1,2)]
    counts <- count.table[,-c(1,2)]
    row.names(counts) <- 1:nrow(counts)
    binpairNum <- nrow(counts)
    
    callDCI.KS <- function(i){
        
        counts_1 <- count.table[i,group == levels(group)[1]]
        counts_2 <- count.table[i,group == levels(group)[2]]
        results_binpair <- data.frame(row.names = row.names(counts)[i], stat = NA, pvalue = NA, pvalue.adj = NA)
        res.binpair <- ks_test(counts_1,counts_2)
        results_binpair[1,"stat"] <- res.binpair[1]
        results_binpair[1,"pvalue"] <- res.binpair[2]
        
        return(results_binpair)
    }
    
    results <- matrix(data=NA, nrow=binpairNum,ncol=5, dimnames=list(row.names(counts),c("bin1","bin2","stat","pvalue","pvalue.adj")))
    results <- as.data.frame(results)
    for(i in 1:binpairNum){
        cat("\r",paste0("scHiCDiff.KS is analyzing",i,"of",binpairNum,"bin paris"))
        results[i,1:2] <- binpair[i,]
        results[i,3:5] <- callDCI.KS(i)
    }
    
    results[,"pvalue.adj"] <- p.adjust(results[,"pvalue"],method="fdr")
    results <- results[order(results[,"stat"],decreasing = TRUE),]
    
    return(results)
}
