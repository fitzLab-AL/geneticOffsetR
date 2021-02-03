################################################################################
# Supporting R functions used in:
# Fitzpatrick et al. (in review) Experimental support for genomic prediction of 
# climate maladaptation using the machine learning approach Gradient Forests.
# Molecular Ecology Resources.
#
# written by MC Fitzpatrick at the Appalachian Lab, Frostburg, MD
#
# Code is provided as is, without support 
################################################################################


# gradientForest Batch Function ------------------------------------------------
# 
# This function will fit GF models to a large number of SNPs,
# optionally in parallel. Function also can be used to fit a random
# subset of SNPs. To save memory, it returns ONLY importance (R2) values,
# not the GF model object. The output from this function can be fed into
# the outlier detection functions below.


runGF <- function(freqFilePath, envTab, vars, ntree=500, rmRare=F, SNPs= "all", 
                  cores=NA, indLoci=T, nRand=0){
  
  ### CHECK THE FORMATTING OF THE EXAMPLE DATA & MAKE SURE YOUR DATA
  ### TABLES ARE THE SAME (OR EDIT CODE TO MATCH YOUR FORMATTING)
  
  # freqFilePath = path to allele freq. file
  # envTab = table of env data, with a row for each population
  # vars = vector of variable names to be used in gf modeling
  # ntree = the number of trees to fit using GF
  # rmRare = remove rare alleles?
  # SNPs = run all SNPs or subset (for subset, provide vector of SNP IDs)?
  # cores = # number of processors to use
  # indLoci = run loci individually or all in the same model?
  # nRand = number of random SNPs to model if nRand > 0
  
  require(data.table)
  require(gradientForest)
  require(parallel)
  require(foreach)
  require(doParallel)
  
  # create custom object to hold output 
  gfOutObj <- setClass("gfOutObj", slots = c(alFreq="data.frame", imp="list"))
  
  # read in allele frequency table
  alFreq <- fread(freqFilePath, header=T, data.table=F)
  
  # names of SNPs
  # see formatting of example data
  loci <- alFreq$chr_snp
  
  # transpose data to pops x alleles
  alFreq <- t(alFreq[,-1])
  colnames(alFreq) <- loci
  
  # set NA values to 0
  #if(sum(is.na(alFreq))>0){alFreq[is.na(alFreq)] <- 0}
  
  # rearrange populations and check to see 
  # if they are in same order in the env table and the SNP table
  # Check formatting of both tables for column names, etc.
  envTab <- envTab[match(rownames(alFreq), envTab$pop_code),]
  if(sum(rownames(alFreq)!=envTab$pop_code)!=0){
    stop("allele and env tables do not align!")
  }
  
  # run all SNPs or a user defined subset?
  if(length(SNPs)>1){alFreq <- alFreq[, match(SNPs, colnames(alFreq))]}
  
  # run random subset of SNPs?
  if(nRand>0){alFreq <- alFreq[, match(sample(loci, nRand), colnames(alFreq))]}
  
  # remove SNPs that are polymorphic in fewer than five populations?
  if(rmRare==T){
    alFreq  <- alFreq[,which(unlist(lapply(apply(alFreq, 2, unique), length))>=6)]
  }
  
  # check to make sure user has set the number of cores (min = 1)
  if(indLoci==T & is.na(cores)){
    stop(paste0("If fitting GF models individually to each SNP (indLoci=T), you need to provide the number of processors to use (i.e., cores = 4). You have a maximum of ", detectCores(), " cores available."))
  }
  
  # run in parallel if fitting SNPs individually
  if(indLoci==T & !is.na(cores)){
    # fit gf model to each SNP individually
    cl <- makeCluster(cores, setup_strategy = "sequential")
    registerDoParallel(cl)
    
    gfMods <- foreach(k=1:ncol(alFreq), .verbose=F, .packages=c("gradientForest")) %dopar%{
      locus <- data.frame(alFreq[,k])
      names(locus) <- colnames(alFreq)[k]
      gf.mod <- gradientForest(data.frame(envTab[, vars], locus), 
                               predictor.vars=vars, response.vars=colnames(alFreq)[k], 
                               corr.threshold=0.5, ntree=ntree, trace=F)
      # extract importance values
      if(!is.null(gf.mod)){
        imps <- importance(gf.mod)
        imps <- imps[order(names(imps))]
        data.frame(imps, SNP = colnames(alFreq)[k])
      }
    }
    
    stopCluster(cl)
    return(gfOutObj(alFreq = data.frame(alFreq), imp = gfMods))
  } else {
    # run all SNPs at once if not fitting individually
    gf.mod <- gradientForest(data.frame(envTab[, vars], alFreq), 
                             predictor.vars=vars, response.vars=colnames(alFreq), 
                             corr.threshold=0.5, ntree=ntree, trace=T)
    return(gf.mod)
  }
}
##////////////////////////////////////////////////////////////////////////////##


# Prepare importance values ----------------------------------------------------
# Function to prepare the output from the runGF batch function for
# outlier detection. The only input to the function is the object
# returned by the runGF function. It returns a data.frame with a row for each
# SNP with R2>0 and a column with the R2 for each predictor variable.
########### build output GF data frames
gfR2tab <- function(gfMods.list){
  gfMods.list <- gfMods.list@imp
  i=1
  while(is.null(gfMods.list[[i]])){i=i+1}
  tab <- do.call(rbind, gfMods.list)
  vrNm <- rep(row.names(tab)[1:nrow(gfMods.list[[i]])], 
              nrow(tab)/nrow(gfMods.list[[i]]))
  tab <- data.frame(variable=vrNm, tab)
  tab <- reshape2::dcast(tab, SNP~variable, value.var="imps")
  totalR2 <- rowSums(tab[,-1])
  return(data.frame(tab, totalR2=totalR2))}
##////////////////////////////////////////////////////////////////////////////##


# Estimate p-values & find outliers --------------------------------------------
findOL <- function(runGFObj, summaryTab, alFreqs, nmSNPs){
  
  # runGFObj = object from runGF function
  # summaryTab = object from gfR2tab function
  # nmSNPs = names of SNPs from which to develop empirical null model

  # get name of SNP if it has a positive R2
  models <- runGFObj@imp
  posR2 <- unlist(lapply(models, function(x){
    return(as.character(unique(x[,2])))}))
  
  # Find which loci have R2 < 0 (no GF model for those) & assign R2=0
  alFreqs <- runGFObj@alFreq
  negR2 <- !(colnames(alFreqs) %in% posR2)
  negR2 <- colnames(alFreqs)[negR2] 
  noGF <- data.frame(matrix(0, nrow=length(negR2), ncol=ncol(summaryTab)))
  colnames(noGF) <- colnames(summaryTab)           
  noGF$SNP <- negR2
  summaryTab <- rbind(summaryTab, noGF)
  summaryTab <- summaryTab[match(colnames(alFreqs), summaryTab$SNP),]
  
  #populate results dataframe with quantiles of R2 values
  # using R2's from intergenic SNPs
  iGen.R2 <- summaryTab[which(summaryTab[,1] %in% nmSNPs),-1]
  pV <- sapply(1:nrow(summaryTab), function(x, summaryTab, nmSNPs, iGen.R2){
    snps2Rank <- rbind(summaryTab[x,-1], iGen.R2)
    P <- apply(snps2Rank, 2, function(y){
      rankSNP <- frank(y)
      return(1-rankSNP[1]/length(rankSNP))
    })}, summaryTab, nmSNPs, iGen.R2)
  
  # format output as data.frame
  pV <- t(pV)
  colnames(pV) <- paste("pValue_", colnames(pV), sep="")
  pV <- data.frame(summaryTab, pV)
  return(pV)
}
##////////////////////////////////////////////////////////////////////////////##


# calculate genetic offset -----------------------------------------------------


########## calculate adaptive offset for populations in space or time
offsetDist <- function(gfMod, vars, env2, combined=F,
                       pops = envPop$pop_code, weighted=FALSE){
  #gfMod = gf model for prediction
  #vars = names of env variables
  #env2 = new environment (new place / time)
  
  # calculates the proportion of SNPs with an R2=0
  # this is used to adjust the transformed values to reflect
  # these GF models with R2=0.
  noR2prop <- length(gfMod$result)/length(gfMod$Y)
  
  if(weighted){
    transEnv2 <- predict(gfMod, env2[,vars])*noR2prop #new env
    transEnv1 <- predict(gfMod)*noR2prop #current env
  } else {
    transEnv2 <- predict(gfMod, env2[,vars]) #new env
    transEnv1 <- predict(gfMod) #current env
  }
  
  if(combined){
    transEnv1 <- transEnv1[1:42,names(transEnv2)]
    # transEnv1 <- data.frame(pop_code=rep(pops, length(unique(gfMod$X$gf.name))),
    #                                      transEnv1)
    # transEnv1 <- aggregate(.~pop_code, data=transEnv1, mean)[,-1]
  }
  
  #calculate Euclidean distance in transformed env space
  num <- nrow(transEnv1)
  dOut <- lapply(1:num, function(x, tEnv1, tEnv2){
    as.numeric(pdist(tEnv1[x,],  tEnv2)@dist)}, tEnv2=transEnv2, tEnv1=transEnv1)
  return(dOut)
}
##////////////////////////////////////////////////////////////////////////////##