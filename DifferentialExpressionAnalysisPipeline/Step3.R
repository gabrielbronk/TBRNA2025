## Script for Step 6 One Only:

## In this script we run DESeq2, this time running it with all demographic and clinical
## covariates. 

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
  
  
  if (NoCovariatesOption == 1) {
    CurrentString = '~ CaseOrControl'
  }
  
  ##----------------------------------------------------------------------------
  ## Now I run DE Analysis:
  
  MetadataRTI = Metadata[BarcodesOfBirthKids,]   ##BarcodesOfInterest
  
  if (NumberOfCurrentConfounders == 0) {
    BarcodesOfNonNA = BarcodesOfBirthKids
  } else if (NoCovariatesOption == 1) {
    BarcodesOfNonNA = BarcodesOfBirthKids
  } else {
    MetadataColumnsForConfounders = MetadataRTI[, CurrentConfounders]
    LogicalsOfNonNA = complete.cases(MetadataColumnsForConfounders)
    IndicesOfNonNA = which(LogicalsOfNonNA == TRUE)
    BarcodesOfNonNA = BarcodesOfBirthKids[IndicesOfNonNA]
  }
  
  
  dds = DESeqDataSetFromMatrix( countData = CountsMatrixRTI[, BarcodesOfNonNA],
        colData = MetadataRTI[BarcodesOfNonNA, ], design = as.formula(CurrentString) )
  
  
  # dds <- estimateSizeFactors(dds)
  # dds <- estimateDispersions(dds)
  # dds <- nbinomWaldTest(dds, maxit=500)
  # 
  
  
  
  
  
  
  dds = DESeq(dds) #
  ddsOriginal = dds
  
  
  NewIndicesOfConvergence = which(mcols(dds)$betaConv)
  NumberOfBadIndices = dim(CountsMatrixRTI)[1]-length(NewIndicesOfConvergence)
  dds = dds[NewIndicesOfConvergence, ]
  
  
  ## Just for inspecting non-convergence:
  NonconvergedIndices = which(mcols(ddsOriginal)$betaConv == FALSE)
  #CountsMatrixUsed = CountsMatrixRTI[, BarcodesOfNonNA]
  #CountsMatrixUsed[NonconvergedIndices,]
  NumberOfNonConvergedIndices = length(NonconvergedIndices)
  
  
  # NewIndicesOfNonconvergence = which(mcols(dds)$betaConv == FALSE)
  # AllConsideredRNANames = rownames(CountsMatrixRTI)
  # NewNonconvergedRNANames = AllConsideredRNANames[NewIndicesOfNonconvergence]
  ## There's a small chance this might get rid of some DEGs from the previous analysis,
  ## but it is what it is because if DESeq2 can't figure out the coefficients,
  ## then I can't know whether it's a DEG or not (except maybe by running DESeq for more trials).

  
  ## Cooks cut off argument is set to FALSE so that the p-values will be produced (rather than
  ## just getting NA flags wherever there were outliers. The outliers were removed and replaced,
  ## so there's no need to worry about getting a false positive - this was stated by the DESeq2 documentation).
  ResultsDF = data.frame(results(dds, contrast = ContrastVariables, pAdjustMethod = "fdr", cooksCutoff=FALSE))
 
  IndsOfDEGs = which(ResultsDF$padj<0.05)
  NumberOfDEGs = length(IndsOfDEGs)
  CoefficientMatrix = coef(dds) 
  StandardErrorMatrix = coef(dds, SE=TRUE)
  
  
  
  ## Ensuring no p-values came up as NA due to some fluke in DESeq2:
  NumberOfFlukeNAs = length(which(is.na(ResultsDF$pvalue) == TRUE))
  print('NumberOfFlukeNAs')
  print(NumberOfFlukeNAs)
  
  
  ## Inspecting the non-converged RNAs (to see if fold change was 0/0):
  ResultsDFOriginal = data.frame(results(ddsOriginal, contrast = ContrastVariables, pAdjustMethod = "fdr", cooksCutoff=FALSE))
  NumberOfFoldChangeNAs = length(which(is.na(ResultsDFOriginal$log2FoldChange) == TRUE))
  print('NumberOfFoldChangeNAs')
  print(NumberOfFoldChangeNAs)
  print('NumberOfNonConvergedIndices')
  print(NumberOfNonConvergedIndices)
  print('NumberOfBadIndices')
  print(NumberOfBadIndices)
  
  #Just for inspecting NA fold changes:
  #NAFoldChangeIndices = which(is.na(ResultsDFOriginal$log2FoldChange) == TRUE)
  #NAFoldChangeCountsMatrixUsed = CountsMatrixUsed[NAFoldChangeIndices,]
  
  
  if (NoCovariatesOption == 1) {
    FileNameStep3 = paste(Folder, '/', StringForOutputFileName, 'OutputStep3NoDCCovariates.RData', sep = "")
  }
  
  if (SwappedCovariatesOption > 0) {
    FileNameStep3 =  paste(Folder, '/', StringForOutputFileName, 'OutputStep3CovariatesFrom', CovariateNameString,
    '.RData', sep = "")
    if (NoCovariatesOption == 1) {
      FileNameStep3 =  paste(Folder, '/', StringForOutputFileName, 'OutputStep3CovariatesFrom', CovariateNameString,
                             'NoDCCovariates.RData', sep = "")
    }  
  }
  
  save(ResultsDF, NumberOfDEGs, CoefficientMatrix, StandardErrorMatrix, file = FileNameStep3) #

  
  

