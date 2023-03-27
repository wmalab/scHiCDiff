# scHiCDiff

scHiCDiff is a novel statistical algorithm to detect differential chromatin interactions (DCIs) between two Hi-C experiments at single-cell level. Here, we introduced 4 ways to capture the DCIs: two non-parametric tests (Kolmogorov–Smirnov test/ Cramér-von Mises test) and parametric likelihood ratio test with two regression models (Negative Binomial/ Zero-inflated Negative Binomial). Non-parametric tests are advantageous by allowing us detecting DCIs without any assumption on data distribution; negative binomial(NB) is the most common assumption for interaction counts in bulk Hi-C parametric approaches, while zero-inflated Negative Binomial(ZINB) regression models is specially designated for the interaction comparison at single-cell level by taking the excessive zeros feature into consideration.


# Installation

To install and load the developmental version of scHiCDiff in R:

```

install.packages("path/scHiCDiff_1.0.tar.gz", repos = NULL, type ="source")
library(scHiCDiff)

```


# Usage

The functions in scHiCDiff can be classified as two types: The first type is the normalization function (scHiCDiff.sim) and the other type is the detection function (scHiCDiff.KS, scHiCDiff.CVM, scHiCDiff.NB and scHiCDiff.ZINB). 

## Normalization Function

The inputs of the normalization function scHiCDiff.norm are illustrated below:

``` 
bias.info.path      The pathway of the three local features (effective length,GC content 
                    and mappability of fragment ends) of all bins. The generation of these 
                    items is available at http://dna.cs.miami.edu/scHiCNorm.
dat_HiC             A N*N scHi-C matrix.
```

The function returns the normalized Hi-C matrix.




## Detection Functions

The inputs for all detection functions are illustrated below:

```
count.table      A non-negative  matrix of scHi-C normalized read counts.The rows of the 
                 matrix are bin pair and columns are samples/cells.
group            A vector of factor which mentions the two condition to be compared, 
                 corresponding to the columns in the count table.
```

The detection function will return a data frame containing the differential chromatin interaction (DCI) analysis results, rows are bin pairs and columns lists the related statistics.

The outputs for the two parametric models are listed below:

```
bin_1,bin_2          The interacting region of the bin pair.
mu_1,mu_2,           MLE of the parameters of NB/ZINB/NBH of group 1 and group 2,
theta_1,theta_2      where mu and theta represent the mean and dispersion estimate of
(pi_1,pi_2)          negative binomial, pi denotes the estimate of zero percentange
norm_total_mean_1,   Mean of normalized read counts of group 1 and group 2.
norm_total_mean_2
norm_foldChange      norm_total_mean_1/norm_total_mean_2.
chi2LR1              Chi-square statistic for hypothesis testing of H0.
pvalue               P value of hypothesis testing of H0 (underlying whether a bin pair 
                     is a DCI).
pvalue.adj.FDR       Adjusted P value of H0's pvalue using Benjamini & Hochberg's method.
Remark               Record of abnormal program information.
```
The outputs for the non-parametric tests are shown below:

```
bin_1,bin_2          The interacting region of the bin pair.
test.statistic       The statistic given by KS/CVM test.
pvalue               P value of hypothesis testing of H0 (underlying whether a bin pair 
                     is a DCI).
pvalue.adj           Adjusted P value of H0's pvalue using Benjamini & Hochberg's method.
```

Example: The data getting from chr19 of oocyte and zygote cells with resolution=200kb (Flyamer et.al.) were untilized as sample data. In the sample data file, it lists all bin pairs with at least one non-zero counts in one of cell types. The first two columns represent the interacting region of each listed bin pair, then followed 89 columns denote the normalized read counts for oocyte cells and the last 54 columns denote the normalized read counts for zygote cells.


```
count.table <- read.table(paste("path/sampledata/oocyte.zygote.chr19.txt")
count.table <- as.matrix(count.table)
group <- factor(c(rep(1,89), rep(2,54)))
result.ks <- scHiCDiff.KS(count.table,group)
result.cvm <- scHiCDiff.CVM(count.table,group)
result.nb <- scHiCDiff.NB(count.table,group)
result.zinb <- scHiCDiff.ZINB(count.table,group)
#common DCIs identified by four methods
common.DCIs <- scHiCDiff.common.DCI(result.ks,result.cvm,result.nb,result.zinb)
```







