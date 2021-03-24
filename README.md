# scHiCDiff

We presented a novel statistical method scHiCDiff to detect differential chromatin interactions (DCIs) between two Hi-C experiments at single-cell level. Here, we introduced 3 models to capture the DCIs: Negative Binomial, zero-inflated Negative Binomial(ZINB) and Negative Binomial Hurdle (NBH) regression model. The first NB model is commonly applied in bulk HiC data; and the last two are models specifically designed to eliminate the effects of extreme sparsity and heterogeneity in sciHi-C data. Once the models were fitted, they were performed a rigorous likelihood ratio test to capture the bin pairs showing significant changes in contact counts.


# Installation

To accelerate data processing and use as less memory as possible, scHiCDiff requires the Matrix pacakges. For specific Hi-C data processing, we tend to use the HiTC, HiCcompare pacakges. In addition, to fit the regression models, we also need the R packages ggsci, VGAM etc.

Thus, with the installation of packages Matrix, mvtnorm, HiTC, HiCcompare, edgeR, ggsci, pscl, VGAM, maxLik, countreg and gamlss, the source code can be performed under R language version 4.0.2.

```
require(Matrix)
require(mvtnorm)
require(HiTC)
require(HiCcompare)
require(edgeR)
require(ggsci)
require(pscl)
require(VGAM)
require(maxLik)
require(countreg)
require(gamlss)
```


# Sample examples

## generate simulated scHi-Cs


```
data.file <- "path/sim.test.data"
simRes <- scHiCDiff.sim(data.file,fold.change=5)

simRes$Hic1.sim
SimRes$Hic2.sim
SimRes$diff.sim
```


## find differential chromatin interactions (DCIs) in real data

The data getting from chr11 of oocyte and zygote cells with resolution=200kb (Flyamer et.al.) were untilized as sample data. In the sample data file, it lists all bin pairs with at least one non-zero counts in one of cell types. The first two columns represent the interacting region of each listed bin pair, then followed 86 columns denote the normalized read counts for oocyte cells and the last 34 columns denote the normalized read counts for zygote cells. 

```
count.table <- read.table(paste("path/oocyte.zygote.filtered.chr11.txt")
count.table <- as.matrix(count.table)
group <- factor(c(rep(1,86), rep(2,34)))
result.nb <- scHiCDiff.NB(count.table,group)
result.zinb <- scHiCDiff.ZINB(count.table,group)
result.nbh <- scHiCDiff.NBH(count.table,group)

result.nb[1:5,]
result.zinb[1:5,]
result.nbh[1:5,]
```
