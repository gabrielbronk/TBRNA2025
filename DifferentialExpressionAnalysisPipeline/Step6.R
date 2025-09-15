## Script for Step 6:

## In this script we run DESeq2 again, this time running it M times, 
## each time with a different confounding covariate set from step (5).

ListOfCovariateDEResults = list()

##==============================================================================
## First we round the corrected counts so that they can be used in DESeq2 again:
##==============================================================================

RoundedCountsMatrixRTI = round(CorrectedCountsMatrix)
NAInds = is.na(RoundedCountsMatrixRTI)
CountsMatrixRTI = RoundedCountsMatrixRTI
CountsMatrixRTI[NAInds] = 0

##==============================================================================
## Running DESeq2 to determine differentially expressed RNAs for case vs. control:
##==============================================================================

CaseVsControlResultsList = list()
NumberOfConfounderSets = dim(ConfoundingCovariateSetsDF)[1]

for (ConfounderSetNum in 1:NumberOfConfounderSets) {
  CurrentConfounderSet = ConfoundingCovariateSetsDF[ConfounderSetNum, 1]
  
  ## Converting the Current Covariate Set from a binary string (e.g. '010001')
  ## to an array (e.g. c(0,1,0,0,0,1)):
  CurrentConfounderSetArray = c()
  for (i in 1:NumberOfSignificantCovariates) {
    CurrentConfounderSetArray = c(CurrentConfounderSetArray,
                                  substr(CurrentConfounderSet, i, i))
  }
  
  IndicesOfCurrentConfounders = which(CurrentConfounderSetArray == '1')
  CurrentConfounders = CovariateArray[IndicesOfCurrentConfounders]
  NumberOfCurrentConfounders = length(CurrentConfounders)
  
  if (NumberOfCurrentConfounders == 0) {
    CurrentString = '~ CaseOrControl'
    
  } else if (NumberOfCurrentConfounders == 1) {
    CurrentString = paste("~ ", CurrentConfounders, ' + CaseOrControl', sep = '')
    
  } else if (NumberOfCurrentConfounders > 1) {
    CurrentString = paste("~ ", CurrentConfounders[1], sep = '')
    
    for (h in 2:NumberOfCurrentConfounders) {
      CurrentString = paste(CurrentString, ' + ', CurrentConfounders[h], sep = '')
    }
    CurrentString = paste(CurrentString, ' + CaseOrControl', sep = '')
    
  }
  
  ##----------------------------------------------------------------------------
  ## Now I run DE Analysis:
  
  MetadataRTI = Metadata[BarcodesOfInterest,]
  
  if (NumberOfCurrentConfounders == 0) {
    BarcodesOfNonNA = BarcodesOfInterest
  } else {
    MetadataColumnsForConfounders = MetadataRTI[, CurrentConfounders]
    LogicalsOfNonNA = complete.cases(MetadataColumnsForConfounders)
    IndicesOfNonNA = which(LogicalsOfNonNA == TRUE)
    BarcodesOfNonNA = BarcodesOfInterest[IndicesOfNonNA]
  }
  
  
  if (NumberOfCurrentConfounders == 0) {
    CaseVsControlResultsList[[ConfounderSetNum]] = 0 ## I put this because we don't need to do DESeq2 
    ## since we already have the DE results for no confounders.
  } else {
  dds = DESeqDataSetFromMatrix( countData = CountsMatrixRTI[, BarcodesOfNonNA],
        colData = MetadataRTI[BarcodesOfNonNA, ], design = as.formula(CurrentString) )
  
  dds = DESeq(dds)
  ddsOriginal = dds
  NewIndicesOfConvergence = which(mcols(dds)$betaConv)
  dds = dds[NewIndicesOfConvergence, ]
  
  # NewIndicesOfNonconvergence = which(mcols(dds)$betaConv == FALSE)
  # AllConsideredRNANames = rownames(CountsMatrixRTI)
  # NewNonconvergedRNANames = AllConsideredRNANames[NewIndicesOfNonconvergence]
  ## There's a small chance this might get rid of some DEGs from the previous analysis,
  ## but it is what it is because if DESeq2 can't figure out the coefficients,
  ## then I can't know whether it's a DEG or not (except maybe by running DESeq for more trials).
  ##########################################
  
  ResultsDF = data.frame(results(dds, contrast = ContrastVariables, pAdjustMethod = "fdr"))
  CaseVsControlResultsList[[ConfounderSetNum]] = ResultsDF
  
  if (SectionToRun == 'DAlternate') {
  ## I just added this in for Step6Alternate (Don't pay  attention to this if just doing regular step 6):
    for (q in 1:NumberOfSignificantCovariates){
    ListOfCovariateDEResults[[q]] = data.frame(results(dds, contrast = ContrastMatrix[q,], pAdjustMethod = "fdr"))
  }
  }
    
  }
  
}

save(CaseVsControlResultsList,  ListOfCovariateDEResults, file = FileNameStep6)
