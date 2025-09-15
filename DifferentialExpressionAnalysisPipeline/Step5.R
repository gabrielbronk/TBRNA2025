## Script for Step 5:

## This script takes the results of Step 3 and determines which sets of 
## confounding covariates were found in Step 3. The output is ConfoundingCovariateSetsDF.

NumberOfDEGs = dim(RNACovariateDFLogicals)[1]
NumberOfCovariates = length(CovariateArray)

##==============================================================================
## Determining the Combinations of Confounding Covariates:
##==============================================================================

MatrixOfConfounding = matrix(0,dim(RNACovariateDFLogicals)[1],dim(RNACovariateDFLogicals)[2])
IndicesOfConfounding = which(RNACovariateDFLogicals == 1 | is.na(RNACovariateDFLogicals) == TRUE)
MatrixOfConfounding[IndicesOfConfounding] = 1

IDMatrix = matrix(0,NumberOfDEGs,1)
TallyMatrix = matrix(0,NumberOfDEGs,1)
Row = 1
for (i in 1:NumberOfDEGs) {
  ID = c()
  for (j in 1:NumberOfCovariates) {
    ID = paste(ID, MatrixOfConfounding[i,j], sep = "")
  }
  FoundRowOfIDMatrix = which(IDMatrix == ID)
  if (length(FoundRowOfIDMatrix) == 0) {
    IDMatrix[Row] = ID
    TallyMatrix[Row] = 1
    Row = Row + 1
  } else {
    TallyMatrix[FoundRowOfIDMatrix] = TallyMatrix[FoundRowOfIDMatrix] + 1
  }
}

IDTallyDF = data.frame('ID' = IDMatrix, 'Tally' = TallyMatrix)
IndicesToKeep = which(IDTallyDF[,1] != "0")
ConfoundingCovariateSetsDF = IDTallyDF[IndicesToKeep,]


save(ConfoundingCovariateSetsDF, MatrixOfConfounding, file=FileNameStep5)
