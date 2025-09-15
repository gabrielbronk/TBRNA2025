## Script for Step 7:

## This script assembles the final DEG list, using the results from step (6), 
## where the p-value for each DEG is from the run of DESeq2 in step (6) 
## that used that DEG's confounding covariate set in the DESeq2 model.


##==============================================================================
## Making the MatrixOfConfounding:
##==============================================================================

DEGNames = rownames(CausaMatrix)
NumberOfDEGs = length(DEGNames)
  
OriginalDEGResultsDF = DEGList[[IndexOfBestTrial]]


MatrixOfConfounding = matrix(0,dim(CausaMatrix)[1],dim(CausaMatrix)[2])
IndicesOfConfounding = which(CausaMatrix > NominalPValueThreshold & CausaMatrix < 1000)
MatrixOfConfounding[IndicesOfConfounding] = 1


##==============================================================================
## Making Final DEG List:
##==============================================================================

IDArray = ConfoundingCovariateSetsDF[,1]
AssembledDEGMatrix = matrix(0,NumberOfDEGs,6)

for (i in 1:NumberOfDEGs) {
  
  ID = c()
  for (j in 1:NumberOfSignificantCovariates) {
    ID = paste(ID, MatrixOfConfounding[i,j], sep = "")
  }
  ## Now we find which index of IDArray this ID is:
  IndexOfID = which(IDArray == ID)
  
  CurrentDEResultsDF = CaseVsControlResultsList[[IndexOfID]]
  CurrentDEG = DEGNames[i]
  
  ## If there were confounders, then we use the DE analysis results from DESeq2 using confounders:
  if (sum(MatrixOfConfounding[i,]) > 0) {
    CurrentRow = as.numeric(CurrentDEResultsDF[CurrentDEG,])
    ## However, if there were not confounders, then we use the original DE analysis results:
  } else {
    CurrentRow = as.numeric(OriginalDEGResultsDF[CurrentDEG,])
  }
  AssembledDEGMatrix[i,] = CurrentRow
  
}

AssembledDEGDF = as.data.frame(AssembledDEGMatrix)
colnames(AssembledDEGDF) = colnames(OriginalDEGResultsDF)
rownames(AssembledDEGDF) = DEGNames

AssembledDEGsAndConfoundersMatrix = cbind(AssembledDEGMatrix, MatrixOfConfounding)
rownames(AssembledDEGMatrix) = DEGNames

##==============================================================================
## Now we see which of the original DEGs are still significant now that we have 
## included confounding covariates in DESeq2:
##==============================================================================

if (PValueOption == 1) {
SurvivingDEGIndices = which(AssembledDEGDF$padj < 0.05)
} else if (PValueOption == 2) {
  SurvivingDEGIndices = which(AssembledDEGDF$pvalue < 0.001)
}
NumberOfSurvivingDEGs = length(SurvivingDEGIndices)

FinalDEGDF = AssembledDEGDF[SurvivingDEGIndices,]
SurvivingDEGNames = rownames(FinalDEGDF)

MatrixOfConfoundingFinal = MatrixOfConfounding[SurvivingDEGIndices,]
NumberOfConfoundersFinal = rowSums(MatrixOfConfoundingFinal)


save(FinalDEGDF, MatrixOfConfoundingFinal, AssembledDEGsAndConfoundersMatrix,
     NumberOfSurvivingDEGs, SurvivingDEGNames, file = FileNameStep7)

