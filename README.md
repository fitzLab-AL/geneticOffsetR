# geneticOffsetR
Custom R functions from:  Fitzpatrick MC, Chhatre VE, Soolanayakanahally RY, Keller SR (in review) Experimental support for genomic prediction of climate maladaptation using the machine learning approach Gradient Forests.

```{r}

# path to allele freq data
filePath <- "/Users/mfitzpatrick/code/geneticOffsetR/lociTable.csv"

# env data for each population
load("/Users/mfitzpatrick/code/geneticOffsetR/envTable.RData")

# variables to use in fitting GF
vars <- c("y", "elev", "bio2", "bio10", "bio11", "bio18", "bio19")

# names of intergenic SNPs for outlier detection
load("/Users/mfitzpatrick/code/geneticOffsetR/igenSNPs.RData")

# fit GF models to 6000+ SNPs individually, in parallel and return
# importance (R2) values
gfSNPs <- runGF(freqFilePath=filePath,
                 envTab=envPop, 
                 vars=vars, 
                 cores=12, 
                 ntree=500)

# extract importance (R2) values and format table
gfSNPs.R2 <- gfR2tab(gfSNPs)

# perform outlier detection using ~1500 intergenic SNPs
pVals <- findOL(runGFObj=gfSNPs, 
                     summaryTab=gfSNPs.R2, 
                     nmSNPs=iGen)
```
