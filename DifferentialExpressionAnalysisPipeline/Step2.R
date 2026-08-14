# Step 2

tic()

##--------------------
## Only looking at kids that have blood samples at birth and have measurements of the 
## covariates we will use in differential expression analysis:
## (RTI stands for RNA Timepoint of Interest, which is birth)
MetadataRTI = Metadata[BarcodesOfBirthKids,]   
MetadataColumnsForCovariates = MetadataRTI[, CovariateArray]
LogicalsOfNonNA = complete.cases(MetadataColumnsForCovariates)
IndicesOfNonNA = which(LogicalsOfNonNA == TRUE)

BarcodesOfNonNA = BarcodesOfBirthKids[IndicesOfNonNA]
MetadataNonNA = Metadata[BarcodesOfNonNA,]

##---------------------
## Making the model for DESeq2:
AllCovariatesForDE = c(TechVars, CovariateArray)
NumberOfCovariates = length(AllCovariatesForDE)

CurrentString = paste("~ ", AllCovariatesForDE[1], sep = '')
if (NumberOfCovariates > 1) {
  for (h in 2:NumberOfCovariates) {
    CurrentString = paste(CurrentString, ' + ', AllCovariatesForDE[h], sep = '')
  }
}
CurrentString = paste(CurrentString, ' + CaseOrControl', sep = '')
print(CurrentString)


##---------------------
## Preparing to run DESeq2:

if (PrefilterOption == 1) {
  ddsForNorm = DESeqDataSetFromMatrix(countData = OrderedCountsMatrix[,AllNewbornBarcodes],
                             colData = Metadata[AllNewbornBarcodes,],
                             design = ~ CaseOrControl )
  ddsForNorm <- estimateSizeFactors(ddsForNorm)
  NormCountsMatrix <- counts(ddsForNorm, normalized = TRUE)
  keep <- rowSums(NormCountsMatrix >= 10) >= 5
}

dds = DESeqDataSetFromMatrix(countData = OrderedCountsMatrix[,BarcodesOfNonNA],
                             colData = MetadataNonNA,
                             design = as.formula(CurrentString))

if (PrefilterOption == 1) {
  dds <- dds[keep,]  
}

##-------------------

if (OutlierWeightingOption == 1) {
  
  EnsemblsFromNoWeighting = rownames(CooksMatrix)
  EnsemblsFromWeighting = rownames(dds)
  WeightsMatrix = matrix(1, nrow(dds), ncol(dds))
  rownames(WeightsMatrix) = EnsemblsFromWeighting
  colnames(WeightsMatrix) = colnames(dds)
  
  for (q in 1:length(EnsemblsFromNoWeighting)) {
    CurrentEnsembl = EnsemblsFromNoWeighting[q]
    if( (CurrentEnsembl %in%  EnsemblsFromWeighting) == TRUE) {
      IndsToZero = which(CooksMatrix[CurrentEnsembl,] > quantile_value)
      WeightsMatrix[CurrentEnsembl,IndsToZero] = 10^(-9)
    }
  }
  assay(dds, "weights") = WeightsMatrix
  
}

##----------------------
## Running DESeq2:

dds = DESeq(dds) 
ddsOriginal = dds

IndicesOfConvergence = which(mcols(dds)$betaConv)
dds = dds[IndicesOfConvergence, ]

ResultsDF = data.frame(results(dds, contrast = ContrastVariables, pAdjustMethod = "fdr", cooksCutoff=TRUE))

CoefficientMatrix = coef(dds) 
StandardErrorMatrix = coef(dds, SE=TRUE)

##------------------------------
## Post-filtering (removing genes for which there are not at least 44 samples containing at least 10 (normalized) counts):
## RE stands for Reasonable Expression:  
IndicesOfRE <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 44  
ddsRE <- dds[IndicesOfRE,]  
ResultsDFRE =  data.frame(results(ddsRE, contrast = ContrastVariables, pAdjustMethod = "fdr", cooksCutoff=TRUE)) #FALSE "bonferroni"

toc()

##------------------------------
## Saving:

save(ResultsDF, ResultsDFRE, CoefficientMatrix, StandardErrorMatrix, dds, file = FileNameStep3)

