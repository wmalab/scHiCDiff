# scHiCDiff

scHiCDiff is a novel statistical algorithm to detect differential chromatin interactions (DCIs) between two Hi-C experiments at single-cell level. Here, we introduced 5 ways to capture the DCIs: two non-parametric tests (Kolmogorov–Smirnov test/ Cramér-von Mises test) and parametric likelihood ratio test with three regression models (Negative Binomial/ Zero-inflated Negative Binomial/ Negative Binomial Hurdle). Non-parametric tests are advantageous by allowing us detecting DCIs without any assumption on data distribution; negative binomial(NB) is the most common assumption for interaction counts in bulk Hi-C parametric approaches, while zero-inflated Negative Binomial(ZINB) and Negative Binomial Hurdle (NBH) regression models are specially designated for the interaction comparison at single-cell level by taking the excessive zeros feature into consideration.


# Installation

To accelerate data processing and use as less memory as possible, scHiCDiff requires the Matrix pacakges. For specific Hi-C data processing, we tend to use the HiTC, HiCcompare pacakges. For non-parametric tests, we utilize the R package twosamples to perform. In addition, to fit the regression models, we also need the R packages ggsci, VGAM etc.

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
require(twosamples)

```


# Usage

The functions in scHiCDiff can be classified as two types: The first type is the simulation function (scHiCDiff.sim) and the other type is the detection function (scHiCDiff.KS, scHiCDiff.CVM, scHiCDiff.NB, scHiCDiff.ZINB and scHiCDiff.NBH). 

## Simulation Function

The inputs of the simulation function scHiCDiff.sim are illustrated below:

```
file.path       The pathway of single cell files. All scHi-C data used in simulation
                should be stored in this pathway. Each scHi-C file is performed as 
                three-column format containing the first interacting region of the bin 
                pair, the second interacting region of the bin pair and the interaction 
                frequency of the bin pair.
fold.change     The amount of fold change.
resolution      The resolution of singel-cell HiC data, eg:200kb will input 200,000
sample.num      The number of single cells tending to generate in each condition.
                (<= the number of inputted singel cells)
pDiff           The probability that an interaction will be differential.
```

The function returns a list that contains the simulated replicates and the matrix of the true DCI regions. The list contains the following elements:

```
Hic1.sim        A list containing the simulated scHi-C matrices of the first condition.
Hic2.sim        A list containing the simulated scHi-C matrices of the second condition.
diff.sim        A sparceMatric containing the position of the differential interactions.
```


Simulation Example: The simulation test data is a dataset with 8 single-cells getting from chr1 of Diploid ESC cultured with 2i in Nagano et al. with resolution=200kb.

```
data.file <- "path/sampledata/sim.test.data"
simRes <- scHiCDiff.sim(data.file,fold.change=5)
```



## Detection Functions

The inputs for all detection functions are illustrated below:

```
count.table      A non-negative  matrix of scHi-C normalized read counts.The rows of the 
                 matrix are bin pair and columns are samples/cells.
group            A vector of factor which mentions the two condition to be compared, 
                 corresponding to the columns in the count table.
```

The detection function will return a data frame containing the differential chromatin interaction (DCI) analysis results, rows are bin pairs and columns lists the related statistics.

The outputs for the three parametric models are listed below:

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
pvalue.adj.FDR       Adjusted P value of H0's pvalue using Benjamini & Hochberg's method.
```

Example: The data getting from chr11 of oocyte and zygote cells with resolution=200kb (Flyamer et.al.) were untilized as sample data. In the sample data file, it lists all bin pairs with at least one non-zero counts in one of cell types. The first two columns represent the interacting region of each listed bin pair, then followed 86 columns denote the normalized read counts for oocyte cells and the last 34 columns denote the normalized read counts for zygote cells.


```
count.table <- read.table(paste("path/sampledata/oocyte.zygote.filtered.chr11.txt")
count.table <- as.matrix(count.table)
group <- factor(c(rep(1,86), rep(2,34)))
result.ks <- scHiCDiff.KS(count.table,group)
result.cvm <- scHiCDiff.CVM(count.table,group)
result.nb <- scHiCDiff.NB(count.table,group)
result.zinb <- scHiCDiff.ZINB(count.table,group)
result.nbh <- scHiCDiff.NBH(count.table,group)
```







