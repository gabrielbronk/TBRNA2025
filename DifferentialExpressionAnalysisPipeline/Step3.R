## Script for Step 3:

## In this script we run DESeq2 with the technical variables selected in Step 2,
## In addition to technical variables, this script also can use demographic and clinical covariates.

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

  CurrentConfounders = CovariateArray
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
  
  
  if (DCCovariatesOption == 0) {
    CurrentString = '~ CaseOrControl'
  }
  
  ##----------------------------------------------------------------------------
  ## Now I run DE Analysis:
  
  MetadataRTI = Metadata[BarcodesOfBirthKids,]   
  
  if (NumberOfCurrentConfounders == 0) {
    BarcodesOfNonNA = BarcodesOfBirthKids
  } else if (DCCovariatesOption == 0) {
    BarcodesOfNonNA = BarcodesOfBirthKids
  } else {
    MetadataColumnsForConfounders = MetadataRTI[, CurrentConfounders]
    LogicalsOfNonNA = complete.cases(MetadataColumnsForConfounders)
    IndicesOfNonNA = which(LogicalsOfNonNA == TRUE)
    BarcodesOfNonNA = BarcodesOfBirthKids[IndicesOfNonNA]
  }
  
  
  dds = DESeqDataSetFromMatrix( countData = CountsMatrixRTI[, BarcodesOfNonNA],
        colData = MetadataRTI[BarcodesOfNonNA, ], design = as.formula(CurrentString) )
  
  dds = DESeq(dds) 
  ddsOriginal = dds
  
  NewIndicesOfConvergence = which(mcols(dds)$betaConv)
  dds = dds[NewIndicesOfConvergence, ]
  
  ## Just for inspecting non-convergence:
  NonconvergedIndices = which(mcols(ddsOriginal)$betaConv == FALSE)
  NumberOfNonConvergedIndices = length(NonconvergedIndices)
  
  ## Cooks cut off argument is set to FALSE so that the p-values will be produced (rather than
  ## just getting NA flags wherever there were outliers. The outliers were removed and replaced,
  ## so as stated in the DESeq2 documentation, there will not be false positives).
  ResultsDF = data.frame(results(dds, contrast = ContrastVariables, pAdjustMethod = "fdr", cooksCutoff=FALSE))
 
  IndsOfDEGs = which(ResultsDF$padj<0.05)
  NumberOfDEGs = length(IndsOfDEGs)
  CoefficientMatrix = coef(dds) 
  StandardErrorMatrix = coef(dds, SE=TRUE)
  
  if (DCCovariatesOption == 0) {
    FileNameStep3 = paste(Folder, '/', StringForOutputFileName, 'OutputStep3NoDCCovariates.RData', sep = "")
  }
  
  save(ResultsDF, NumberOfDEGs, CoefficientMatrix, StandardErrorMatrix, file = FileNameStep3) 

  
  

