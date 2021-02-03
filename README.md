# geneticOffsetR
Custom R functions from:  Fitzpatrick MC, Chhatre VE, Soolanayakanahally RY, Keller SR (in review) Experimental support for genomic prediction of climate maladaptation using the machine learning approach Gradient Forests. *Molecular Ecology Resources*



#### Here's a simple example for estimating p-values for outlier detection using GF modeling. We also provide a function for calculating genetic offsets ([Fitzpatrick & Keller 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.12376)), but also see Andy Gougherty's [repository](https://github.com/agougher/poplarAdaptiveOffset) that provides some additional functionality for offset analyses described in [Gougherty et al. (2021) Nature Climate Change](https://www.nature.com/articles/s41558-020-00968-6).


File path to allele freq data. Contains data for ~6500 SNPs, of which ~1500 are intergenic.
```{r}
filePath <- "~/lociTable.csv"
```

Load environmental table (row for each population, column for each variable).
```{r}
load("~/envTable.RData")
```

Names of variables used in fitting GF:
```{r}
vars <- c("y", "elev", "bio2", "bio10", "bio11", "bio18", "bio19")
```

Names of SNPs to be used for outlier detection (~1,500 intergenic SNPs in this example).
```{r}
load("~/igenSNPs.RData")
```

Use the `runGF` function to fit GF models to ~6,500 SNPs individually, in parallel. Function returns only the importance (R2) values from the GF models, not the entire model object. Note that the allele frequency data provided contain the intergenic SNPs, so these are fitted using GF along with the other SNPs.

```{r}
gfSNPs <- runGF(freqFilePath=filePath,
                 envTab=envPop, 
                 vars=vars, 
                 cores=12, 
                 ntree=500)
```

Next, extract the importance (R2) values and format as a table.
```{r}
gfSNPs.R2 <- gfR2tab(gfSNPs)
```

Lastly, calculate p-values using the ~1,500 intergenic SNPs.
```{r}
pVals <- findOL(runGFObj=gfSNPs, 
                     summaryTab=gfSNPs.R2, 
                     nmSNPs=iGen)
```
