scHiCDiff.common.DCI <- function(result.ks,result.cvm,result.nb,result.zinb,pvalue){
  result.cvm <- result.cvm[result.cvm$pvalue.adj<pvalue,]
  result.ks <- result.ks[result.ks$pvalue.adj<pvalue,]
  result.nb <- result.nb[result.nb$pvalue.adj.FDR<pvalue,]
  result.zinb <- result.zinb[result.zinb$pvalue.adj.FDR<pvalue,]
  posi.cvm <- result.cvm$bin1*100000+result.cvm$bin2
  posi.ks <- result.ks$bin1*100000+result.ks$bin2
  posi.nb <- result.nb$bin1*100000+result.nb$bin2
  posi.zinb <- result.zinb$bin1*100000+result.zinb$bin2
  posi <- intersect(intersect(intersect(posi.cvm,posi.ks),posi.nb),posi.zinb)
  bin1 <- floor(posi/100000)
  bin2 <- posi-bin1*100000
  bin <- data.frame(bin1,bin2)
  return(bin1)
}