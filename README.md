# geneticOffsetR
Custom R functions and data sets from:  Fitzpatrick MC, Chhatre VE, Soolanayakanahally RY, Keller SR (2021), Experimental support for genomic prediction of climate maladaptation using the machine learning approach Gradient Forests. *Mol Ecol Resour.* https://doi.org/10.1111/1755-0998.13374

***

## DATA
## List of files in the /data folder and brief explanation of each:

- **balsam_core336inds_42pops.vcf.gz:** *bgzipped vcf file containing SNP genotypes for 336 Populus balsamifera from the core ancesty deme*
- **GF_GBS_predict_gardens.R:** *R script for analysis of height grwoth as a function of gnetic and climate offsets*
- **core336.bayenv.gardenOffset.csv:** *Genetic offset predictions based on outliers identified by Bayenv2*
- **core336.Bayenv-lfmm.gardenOffset.csv:** *Genetic offset predictions based on outliers identified by Bayenv2 and LFMM*
- **core336.bf.pc1.emp.txt:** *List of outlier SNPs with empirical p-values from Bayenv2 analysis using climate PC1*
- **core336.bf.pc2.emp.txt:** *List of outlier SNPs with empirical p-values from Bayenv2 analysis using climate PC2*
- **core336.gf.allOls.gardenOffset.csv:** *Genetic offset predictions based on outliers identified by all outlier methods (Bayenv2, LFMM, GF-X, GF-Raw)*
- **core336.gf.rawFreqs.gardenOffset.csv:** *Genetic offset predictions based on outliers identified by GF-Raw*	
- **core336.gf.stdFreqs.gardenOffset.csv:** *Genetic offset predictions based on outliers identified by GF-X*
- **core336.lfmm.gardenOffset.csv:** *Genetic offset predictions based on outliers identified by LFMM*
- **core336.500_RANDOM.gardenOffset.csv:** *Genetic offset predictions based on 999 bootstrap replicates of 500 randomly sampled SNPs*
- **core336.mahal.gardenOffset.csv:** *Mahalanobis multivariate climate distances between populations and common garden sites*
- **gardenQuery.RData:** *Vermont common garden data*
- **IH_phenology_2015.csv:** *Indian Head common garden data*
- **manhattan_data_core_pc1.txt** *List of outlier SNPs with empirical p-values from LFMM analysis using climate PC1*
- **manhattan_data_core_pc2.txt** *List of outlier SNPs with empirical p-values from LFMM analysis using climate PC2*
- **snpR2.csv:** *GF R^2 values for modeled SNPs*

***

## R SCRIPTS

The R functions allow fitting GF to large numbers of SNPs and help with outlier detection. We also provide a function for calculating genetic offsets ([Fitzpatrick & Keller 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.12376)), but also see Andy Gougherty's [repository](https://github.com/agougher/poplarAdaptiveOffset) that provides additional functionality for genetic offset analyses as described in [Gougherty et al. (2021) Nature Climate Change](https://www.nature.com/articles/s41558-020-00968-6).  

***

**Here's a simple example for estimating p-values for outlier detection using GF modeling on 6000+ SNPs. See the comments within the R functions for more details.**

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
