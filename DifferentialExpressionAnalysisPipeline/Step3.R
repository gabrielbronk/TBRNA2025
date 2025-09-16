## Script for Step 3

## This step determines which covariates are correlated with which RNAs. 
## I do so using the N covariates shown in step (1) to be correlated with 
## the outcome. I run DESeq2 N times, each time including a different one 
## of these covariates in the DESeq2 model (the model contains only that 
## covariate, not the TB outcome). Thus we obtain a list of RNAs that are 
## correlated with each covariate. 

##------------------------------------------------------------------------------

CovariateDEList = list()

DESeqResultsFromBackwardsSelection$Name = rownames(DESeqResultsFromBackwardsSelection)
if (PValueOption == 1) {
  DEGIndices = which(DESeqResultsFromBackwardsSelection$padj < 0.05)
} else if (PValueOption == 2) {
  DEGIndices = which(DESeqResultsFromBackwardsSelection$pvalue < 0.001)
} 
DEGNames = DESeqResultsFromBackwardsSelection[DEGIndices,'Name']
NumberOfDEGs = length(DEGNames)

##------------------------------------------------------------------------------
## Now I round the corrected counts so that they can be used in DESeq2 again:
##------------------------------------------------------------------------------

RoundedCountsMatrixRTI = round(CorrectedCountsMatrix)
NAInds = is.na(RoundedCountsMatrixRTI)
CountsMatrixRTI = RoundedCountsMatrixRTI
CountsMatrixRTI[NAInds] = 0  

CurrentCountsMatrixRTI = CountsMatrixRTI[,BarcodesOfInterest]
MetadataRTI = Metadata[1:NumberOfKidsOfInterest,]

##==============================================================================
## Running DESeq2 to evaluate whether the covariate is correlated with the RNAs:
##==============================================================================

RNACovariateMatrix = matrix(1000,NumberOfDEGs,NumberOfSignificantCovariates)
CaseVsControlResultsList = list()
CaseVsControlResultsListFewerDEGs = list()

for (CovariateNum in 1:NumberOfSignificantCovariates) {
  CurrentCovariateName = CovariateArray[CovariateNum]
  
  CurrentString = paste('~ ', CurrentCovariateName, sep = '')
  IndicesNonNA = which(is.na(MetadataRTI[,CurrentCovariateName]) == FALSE)
  
  CurrentCountsMatrixRTINonNA = CurrentCountsMatrixRTI[,IndicesNonNA]
  MetadataRTINonNA = MetadataRTI[IndicesNonNA,]
  
  dds = DESeqDataSetFromMatrix(countData = CurrentCountsMatrixRTINonNA, 
                               colData = MetadataRTINonNA, 
                               design = as.formula(CurrentString))   
  dds = DESeq(dds)
  ddsOriginal = dds
  dds = dds[which(mcols(dds)$betaConv),]
  ## It's possible some RNAs's GLMs did not converge so it will give NAs, but in the later steps of the pipeline
  ## I'll just play it safe and assume those nonconverged RNAs correlate with the covariate and will include that
  ## covariate in the DESeq2 model for those RNAs.

  if (CovariateClassArray[CovariateNum] == 0) {
    CurrentContrastVariables = ContrastMatrix[CovariateNum,] 
    ResultsDF = data.frame(results(dds, contrast = CurrentContrastVariables, pAdjustMethod = "fdr"))
  } else if (CovariateClassArray[CovariateNum] == 1) {
    ResultsDF = data.frame(results(dds, name = CurrentCovariateName, pAdjustMethod = "fdr"))
  }
  RNACovariateMatrix[,CovariateNum] = ResultsDF[DEGNames,'padj']
  CovariateDEList[[CovariateNum]] = ResultsDF
}

RNACovariateDF = as.data.frame(RNACovariateMatrix)
colnames(RNACovariateDF) = CovariateArray
rownames(RNACovariateDF) = DEGNames

RNACovariateDFLogicals = RNACovariateDF
RNACovariateDFLogicals[RNACovariateDF < 0.05] = 1
RNACovariateDFLogicals[RNACovariateDF > 0.05] = 0


##------------------------------------------------------------------------------
## Saving the Output, and Checking There are No Issues:
##------------------------------------------------------------------------------

## Checking to see if there are any NAs - for instance, if there is a DEG we found in step 2
## for which a covariate's GLM did not converge, then we would get an NA:
NARowsAndColumns = which(is.na(RNACovariateMatrix)==TRUE, arr.ind = TRUE)

save(RNACovariateDF, RNACovariateDFLogicals, DEGNames, 
     CovariateDEList, NARowsAndColumns, NAInds,
     file = FileNameStep3)

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Checks:
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("Any NAs in the DESeq2 results for covariates correlations with RNAs?")
print(NARowsAndColumns)

