## Script for Step 2

## Backward selection of technical variables with what I call "sameness score"
## and
## Performing Differential Expression Analysis (using DESeq2) after correcting the counts matrix for 
## sequencing batch and technical variables; 

## This takes about 24 hours to complete. First it corrects the counts matrix for sequencing batch and
## technical variables. Next it uses the corrected counts to perform differential expression analysis 
## (using DESeq2) - it simply uses case vs control in the DESeq2 model and does not include any 
## covariates in the model in this step (the technical variables have been used in counts correction, as 
## I just mentioned, but demographic/clinical covariates are not used in this step). The script performs 
## backwards selection to find the set of technical variables that, when used to correct the counts 
## matrix, results in the most differentially expressed genes (DEGs) for the case vs. control. We then 
## select that set of technical variables (the one that yielded the most DEGs), and we obtain a list of 
## putative DEGs (putative because we do not know yet if they are a result of confounding, which we 
## will control for in steps below.)

  
##==============================================================================
## Initializing arrays:
##==============================================================================

  ## These are the covariates you will be using in backwards selection:
  CovariateArray = c(TechnicalVariables, 'RIN')
  NumberOfCovariates = length(CovariateArray)
  CovariateClassArray = rep(1,1,NumberOfCovariates)
  ## 0 indicates a 2-level factor. 1 indicates a continuous variable. 2 Indicates a multi-level factor.
  
  Outcome = 'CaseOrControl'
  
  if (TypeOfBatch == '') {
    NumberOfEntries = NumberOfCovariates
  } else {
    NumberOfEntries = NumberOfCovariates+1
  }

  ## NNSP stands for Number of Nominally Significant P-values:
  NNSPMatrix = matrix(10^9,NumberOfEntries,NumberOfEntries)
  
  ## NSQ stands for Number of Significant q-values:
  NSQMatrix = matrix(10^9,NumberOfEntries,1)
  
  SamenessScoreMatrix = matrix(-1,NumberOfEntries,NumberOfCovariates)
  
  ## By "Non Zero RNAs", I just mean this is all RNAs that are expressed 
  ## (since I previously filtered out RNAs that have entirely 0 counts):
  NumberOfNonZeroRNAs = dim(CurrentCountsMatrix)[1]
  
  ## This records the number of covariates that have that RNA as nominally significant, and it does that
  ## for every trial of Backwards Selection (i.e. for each row):
  RecordOfSameness = matrix(0,NumberOfEntries,NumberOfNonZeroRNAs)
  ReferenceRNAs = 1:NumberOfNonZeroRNAs
  
  DEGList = list()
  
  RecordOfRemovedCovariates = c() 
  CurrentCovariateArray = CovariateArray 
  CurrentCovariateClassArray = CovariateClassArray

  
  ##==============================================================================
  ## Starting Backwards Selection:
  ##==============================================================================

    ##--------------------------------------------------------------------------
    ## A Loop to determine the Sameness Scores of each covariate:
    ##--------------------------------------------------------------------------
    for (TrialNum in 1:(NumberOfCovariates+1)) { 
      
      print(TrialNum)
      
      ##------------------------------------------------------------------------
      ## Creating Design Formula:
      ##------------------------------------------------------------------------
      CurrentNumberOfCovariates = length(CurrentCovariateArray)
      
      if (CurrentNumberOfCovariates == 0 & TypeOfBatch == '') {
        ## This stops the loop.
      } else {
        if (CurrentNumberOfCovariates > 0) {
          CurrentString = paste("~", TypeOfBatch)
          for (w in 1:length(CurrentCovariateArray)) {
            CurrentString = paste(CurrentString, "+", CurrentCovariateArray[w])
          }
        } else if (CurrentNumberOfCovariates == 0) {
          CurrentString = paste("~ ", TypeOfBatch, sep = '')
        }
        ##------------------------------------------------------------------------------
        
        
        ##==============================================================================
        ## Correcting the Counts Matrix for Technical Variables:
        ##==============================================================================
        
        ## 0 means I'm using no batch but I am using technical variables.
        ## 1 means I'm using Sequencing Batch and technical variables.
        ## 2 means I'm using Sequencing Batch but no technical variables.
        if (CurrentNumberOfCovariates == 0) {
          BatchOption = 2
        } else if (CurrentNumberOfCovariates > 0) {
          if (TypeOfBatch == '') {
            BatchOption = 0 
          } else {
            BatchOption = 1
          }
        }
        
        ## VC stands for Variables for Correction:
        MetadataVC = Metadata[,c(TypeOfBatch,CurrentCovariateArray)]
        if (CurrentNumberOfCovariates == 0) {
          MetadataVC = data.frame(TypeOfBatch = Metadata[,c(TypeOfBatch)])
          colnames(MetadataVC) = TypeOfBatch
          rownames(MetadataVC) = rownames(Metadata)
        }  
        
        
        
      if (CorrectionOption == 2) {
       
        SelectionString = 'Sameness'
  
        FileNameCurrentCountsMatrix = 
          paste('CorrectedCountsMatrices/', 'SequencingBatch', SelectionString, 'Trial', as.character(TrialNum), '.RData', sep = "")
        
        load(FileNameCurrentCountsMatrix)  
        
      } else {
        ## Correcting counts matrix if you have not created it previously:
        
        
        ##----------------------------------------------------------------------
        ## Making the DESeq2 Object:
        ##----------------------------------------------------------------------
        
        print("DE Analysis of Technical Variables Section Has Started")
        dds = DESeqDataSetFromMatrix(countData = CurrentCountsMatrix, colData = MetadataVC,
                                     design = as.formula(CurrentString) )   
        
        
        ##----------------------------------------------------------------------
        ## Filtering Out RNAs with Low Counts
        ##----------------------------------------------------------------------
        
        ## Now I filter out RNAs that have low counts. This is because if I don't filter,
        ## it gives a warning that some of the RNAs' GLMs did not converge. Therefore, I presume 
        ## those RNAs' coefficients would not be reliable and hence we could not use those 
        ## coefficients to correct the counts matrix. Therefore, we should just exclude those 
        ## non-converged RNAs from our analysis. From https://support.bioconductor.org/p/65091/ Michael Love (the DESeq2 author) suggests to
        ## exclude any RNAs that for which there are not at least 4 samples with normalized counts of 
        ## 10 or greater. This will get rid of most of the RNAs that have this convergence issue.
        
        dds <- estimateSizeFactors(dds)
        nc <- counts(dds, normalized = TRUE)
        ReasonableExpressionFilter <- rowSums(nc >= 10) >= 4
        dds <- dds[ReasonableExpressionFilter, ]
        
        ## Since I am using this filter, I also need to apply it to the counts matrix:
        CurrentCountsMatrixFiltered = CurrentCountsMatrix[ReasonableExpressionFilter,]
        
        ##----------------------------------------------------------------------
        ## Running DESeq2:
        ##----------------------------------------------------------------------
        
        dds = DESeq(dds)
        ddsOriginal = dds
        
        ##----------------------------------------------------------------------
        ## Filtering Out RNAs whose GLMs did not Converge
        ##----------------------------------------------------------------------
        
        ## Here I narrow down RNAs based on which ones converged because if an RNA didn't converge, then 
        ## I think the coefficients may be wrong and the corrected counts for that RNA would therefore be wrong.
        IndicesOfConvergence = which(mcols(dds)$betaConv)
        dds = ddsOriginal[IndicesOfConvergence, ]
        
        NumberOfSamplesForCorrection =  NumberOfSamples
        
        ## Since I am narrowing down RNAs based on which ones converged, I 
        ## do so on the counts matrix also:
        NewCountsMatrix = CurrentCountsMatrixFiltered[IndicesOfConvergence,1:NumberOfSamplesForCorrection]
        
        ##----------------------------------------------------------------------
        ## Preparing for Counts Matrix Correction:
        ##----------------------------------------------------------------------
        
        MySizeFactors = sizeFactors(dds)
        CoefficientMatrixOriginal = coef(dds) 
        CoefficientMatrix = CoefficientMatrixOriginal
        
        NewCountsMatrixOriginal = NewCountsMatrix
        NumberOfRNAs = dim(NewCountsMatrix)[1]
        ## Note these say RTI (which stands for RNA Timepoint of Interest), but we're correcting the counts 
        ## based on all 1815 samples, so RTI is a misnomer.
        MetadataRTI = Metadata[1:NumberOfSamplesForCorrection,]
        MetadataVCRTI = MetadataVC[1:NumberOfSamplesForCorrection,]
        if (CurrentNumberOfCovariates == 0) {
          MetadataVCRTI = data.frame(TypeOfBatch = MetadataVC[1:NumberOfSamplesForCorrection,c(TypeOfBatch)])
          colnames(MetadataVCRTI) = TypeOfBatch
          rownames(MetadataVCRTI) = rownames(Metadata[1:NumberOfSamplesForCorrection,])
        }
        
        ##----------------------------------------------------------------------
        ## Correcting the Counts Matrix:
        ##----------------------------------------------------------------------
        
        print("Counts Correction Section Begun")
        
        ## FCC stands for For Counts Correction:
        if (BatchOption == 0) {
          MetadataFCC = MetadataVCRTI
          CoefficientMatrixFCC = CoefficientMatrix[, 2:(1 + CurrentNumberOfCovariates)]
          ## The "1+" is there because the first column in the CoefficientMatrix is the intercept.
        } else  if (BatchOption == 1) {
          MetadataFCC = MetadataVCRTI[, -1] ## The -1 removes the Sequencing Batch info
          CoefficientMatrixFCC = CoefficientMatrix[, 7:(6 + CurrentNumberOfCovariates)]
          ## 7 is here because there is 1 column for the intercept and 5 columns for Sequencing batches.
          BatchCoefficientMatrixFCC = CoefficientMatrix[, 2:6]
        } else  if (BatchOption == 2) {
          BatchCoefficientMatrixFCC = CoefficientMatrix[, 2:6] 
        }
        
        if (BatchOption == 0 | BatchOption == 1) {
          ## Correcting for Technical Variables:
          for (j in 1:NumberOfSamplesForCorrection) {
            if (class(MetadataFCC)=="data.frame") {
              CurrentCovariateValues = as.numeric(MetadataFCC[j, ])
            } else {
              CurrentCovariateValues = as.numeric(MetadataFCC[j])
            }
            
            for (i in 1:NumberOfRNAs) {
              if (class(CoefficientMatrixFCC)=="numeric") {
                CurrentCoefficentArray = as.numeric(CoefficientMatrixFCC[i])
              } else {
                CurrentCoefficentArray = as.numeric(CoefficientMatrixFCC[i, ])
              }
              
              if (BatchOption == 1) {
                CorrectionFactor = 2 ^ (CurrentCoefficentArray %*% CurrentCovariateValues) 
              }
              if (BatchOption == 0) {
                CorrectionFactor = (2 ^ (CurrentCoefficentArray %*% CurrentCovariateValues))*MySizeFactors[j]
              }
              NewCountsMatrix[i, j] = NewCountsMatrixOriginal[i, j]/CorrectionFactor
            }
          }
        }
        
        ##-------------------------------------------
        CorrectedCountsMatrix = NewCountsMatrix
        
        if (BatchOption == 1 | BatchOption == 2) {
          if (TypeOfBatch == 'SequencingBatch') {
            MinimumBatch = 0
          } else if (TypeOfBatch == 'LibraryBatch') {
              MinimumBatch = 1
          }
          
          
          ## Correcting for Batch:
          for (j in 1:NumberOfSamplesForCorrection) {
            
            if (TypeOfBatch == 'SequencingBatch') {
              CurrentBatch = as.numeric(MetadataVCRTI[j,TypeOfBatch]) - 1 ## This is the sequencing batch
            ## The -1 is there because as.numeric adds 1 because there's a batch 0.
            } else if (TypeOfBatch == 'LibraryBatch') {
              CurrentBatch = as.numeric(MetadataVCRTI[j,TypeOfBatch]) 
            }
            
            for (i in 1:NumberOfRNAs) {
              ## The first entry is batch 1vs0. The 2nd entry is batch 2vs0. etc.
              BatchCoefficient = BatchCoefficientMatrixFCC[i, CurrentBatch]  
              
              if (CurrentBatch == MinimumBatch) {
                BatchCorrectionFactor = 1
              } else {
                BatchCorrectionFactor = 2 ^ (BatchCoefficient)
              }
              CorrectedCountsMatrix[i, j] = NewCountsMatrix[i, j] / (BatchCorrectionFactor*MySizeFactors[j])
            }
          }
        }
         
      if (CorrectionOption == 1) {
      ## Saving the corrected counts matrix:
      FileNameCurrentCountsMatrix = 
        paste('CorrectedCountsMatrices/', TypeOfBatch, 'SamenessTrial', as.character(TrialNum), '.RData', sep = "")
      
      save(CorrectedCountsMatrix, file = FileNameCurrentCountsMatrix)
      }  
      }
        ##======================================================================
        ## A For Loop to get the nominally significant RNAs for each covariate:
        ##======================================================================
        
        if (CorrectionOption < 2) {
        NominalRNAsList = list()
        for (CovariateNum in 1:NumberOfCovariates) {
          if ((CovariateNum %in% RecordOfRemovedCovariates) == FALSE){
            if (CovariateClassArray[CovariateNum] == 1) {
              CurrentCovariateName = CovariateArray[CovariateNum]
              CovariateResultsDF = data.frame(results(dds, name = CurrentCovariateName, pAdjustMethod = "fdr"))
              NominalRNAs = which(CovariateResultsDF$pvalue < 0.05)
              NNSP = length(NominalRNAs)
              NominalRNAsList[[CovariateNum]] = NominalRNAs
              NNSPMatrix[TrialNum, CovariateNum] = NNSP
            }
          } else {
            NominalRNAsList[[CovariateNum]] = 1 ## This is just so that the code doesn't give an error later.
          }
        }
        
        ##======================================================================
        ## A For Loop to get the nominally significant RNAs for the batch:
        ##======================================================================
       
        AllNominalIndices = c()
        
        if (TypeOfBatch == 'SequencingBatch') {
          ## For Sequencing Batch:
          for (SB in 1:5) {
          CurrentContrastVariables = c("SequencingBatch", as.character(SB), '0')
          CovariateResultsDF = data.frame(results(dds, contrast = CurrentContrastVariables,
                                                  pAdjustMethod = "fdr"))
          CurrentNominalIndices = which(CovariateResultsDF$pvalue < 0.05)
          AllNominalIndices = c(AllNominalIndices, CurrentNominalIndices)
          UniqueNominalIndices = unique(AllNominalIndices)
          }
          } else if (TypeOfBatch == 'LibraryBatch') {
          ## For Library Batch:
          for (SB in 2:21) {
          CurrentContrastVariables = c("LibraryBatch", as.character(SB), '1')
          CovariateResultsDF = data.frame(results(dds, contrast = CurrentContrastVariables,
                                                    pAdjustMethod = "fdr"))
          CurrentNominalIndices = which(CovariateResultsDF$pvalue < 0.05)
          AllNominalIndices = c(AllNominalIndices, CurrentNominalIndices)
          UniqueNominalIndices = unique(AllNominalIndices)
          }
        }
        
        }
        
        ##======================================================================
        ## Rounding the corrected counts so that they can be used in DESeq2 again:
        ##======================================================================
        
        ## Here we narrow down the corrected counts matrix from 
        ## all 1815 samples to just the samples at the birth timepoint 
        ## before we do the case vs. control DE analysis.
        
        BarcodesForDESeq = c(BarcodesOfTestKids,BarcodesOfControlKids,BarcodesOfOtherBirthKids)
        RoundedCountsMatrixRTI = round(CorrectedCountsMatrix[,BarcodesForDESeq])
        
        NAInds = is.na(RoundedCountsMatrixRTI)
        CountsMatrixRTI = RoundedCountsMatrixRTI
        CountsMatrixRTI[NAInds] = 0
        
        MetadataRTIReduced = Metadata[BarcodesForDESeq,]
        ##======================================================================
        ## Executing the DESeq2 function for Case Vs Control:
        ##======================================================================
        
        print("DE Analysis of Case vs Control Has Started")
        
        dds = DESeqDataSetFromMatrix(countData = CountsMatrixRTI, colData = MetadataRTIReduced, 
                                     design = ~ CaseOrControl)  
        
        
        dds = DESeq(dds)
        ddsOriginal = dds
        ## Then this line further reduces narrows dds down to only the rows that have betaConv = TRUE
        ## (i.e. the RNAs whose GLMs have converged):
        dds <- dds[which(mcols(dds)$betaConv),]
        
        ###======================================================================
        ## Getting p-values for Case vs Control:
        ##======================================================================
        print("Results")
        
        CovariateResultsDF = data.frame(results(dds, contrast = ContrastVariables, pAdjustMethod = "fdr"))
        
        if (CorrectionOption < 2) {
        NNSP = length(which(CovariateResultsDF$pvalue < 0.05)) 
        NNSPMatrix[TrialNum,(NumberOfCovariates+1)] = NNSP
        }
        
        ## Getting q-values for the outcome of interest:
        NSQ = length(which(CovariateResultsDF$padj < 0.05)) 
        NSQMatrix[TrialNum] = NSQ
        

        
        DEGList[[TrialNum]] = CovariateResultsDF
        
  
        ##======================================================================
        ## Determining the Sameness score for each covariate:
        ##======================================================================
       
        if (CorrectionOption < 2) {
        CovariateNumbersArray = 1:NumberOfCovariates
        
        for (CovariateNum in CovariateNumbersArray) {
          ReferenceNominalRNAs = NominalRNAsList[[CovariateNum]]
          ## This sees how much overlap there is of the nominally significant DE RNAs of the 
          ## Reference Covariate with each of the other covariates and adds up all the overlap 
          ## to make the SamenessScore:
          SamenessScore = 0
          
          ## First we look at the overlap with the batch if using a batch:
          if (nchar(TypeOfBatch) > 0) {
            NominalRNAsLogicals = ReferenceNominalRNAs %in% UniqueNominalIndices   
            SamenessScore = SamenessScore + length(which(NominalRNAsLogicals == TRUE))
          }
          
          for (g in CovariateNumbersArray[-CovariateNum]) {
            ComparisonNominalRNAs = NominalRNAsList[[g]]
            NominalRNAsLogicals = ReferenceNominalRNAs %in% ComparisonNominalRNAs 
            SamenessScore = SamenessScore + length(which(NominalRNAsLogicals == TRUE))
          }
          SamenessScoreMatrix[TrialNum,CovariateNum] = SamenessScore
        }
        
        ##------------------------------------------------------------------------------
        ## Determining which covariate has the highest sameness score, which is the covariate we'll remove:
        if (length(TVNumbersToNeverExclude) > 0){
          HighestSamenessScore = max(SamenessScoreMatrix[TrialNum, -TVNumbersToNeverExclude])
        } else {
          HighestSamenessScore = max(SamenessScoreMatrix[TrialNum,])
        }
        IndexToRemove = which(SamenessScoreMatrix[TrialNum,] == HighestSamenessScore)
        RecordOfRemovedCovariates = c(RecordOfRemovedCovariates, IndexToRemove)
        CurrentCovariateArray = CovariateArray[-RecordOfRemovedCovariates]
        
        ##------------------------------------------------------------------------------
        ## This isn't necessary, but here I make a histogram to visualize how this algorithm works:
        
        RecordOfSamenessRow = matrix(0,1,NumberOfNonZeroRNAs)  
        for (g in CovariateNumbersArray) {
          ComparisonNominalRNAs = NominalRNAsList[[g]]
          NominalRNAsLogicals = ReferenceRNAs %in% ComparisonNominalRNAs 
          NominalRNAsIndices = which(NominalRNAsLogicals == TRUE)
          RecordOfSamenessRow[NominalRNAsIndices] = RecordOfSamenessRow[NominalRNAsIndices]+1
        }
        ## Also including overlap with batch if using batch:
        if (nchar(TypeOfBatch) > 0) {
          NominalRNAsLogicals = ReferenceRNAs %in% UniqueNominalIndices 
          NominalRNAsIndices = which(NominalRNAsLogicals == TRUE)
          RecordOfSamenessRow[NominalRNAsIndices] = RecordOfSamenessRow[NominalRNAsIndices]+1
        }
        RecordOfSameness[TrialNum,] = RecordOfSamenessRow
        
        
        print('Trial')
        print(TrialNum)
        
        print('Index to remove')
        print(IndexToRemove)
        }
        print('Trial')
        print(TrialNum)
        
        print('NSQ')
        print(NSQ)
        
      }
    } 
 
  ##============================================================================  
  ## The Backwards Selection trials are finished.
  ##============================================================================
  
    

  
  ##------------------------------------------------------------------------------
  ## Reformatting the NNSP Matrix:
  if (CorrectionOption < 2) {
  NotApplicableIndices = which(NNSPMatrix == 10^9)
  FinalNNSPMatrix = NNSPMatrix
  FinalNNSPMatrix[NotApplicableIndices] == 'NA'
  NNSPDF = as.data.frame(FinalNNSPMatrix)
  colnames(NNSPDF) = c(CovariateArray, Outcome)
  print(NNSPDF)
  
  ##------------------------------------------------------------------------------
  
  FinalRecordOfRemovedCovariates = RecordOfRemovedCovariates[1:NumberOfCovariates]
  }
  IndexOfBestTrial = which(NSQMatrix == max(NSQMatrix))
  
  
  if (CorrectionOption < 2) {
  save(IndexOfBestTrial, DEGList, CovariateResultsDF, NNSPDF, NSQMatrix,
       RecordOfSameness, SamenessScoreMatrix,
       FinalRecordOfRemovedCovariates, NominalRNAsList, file = FileNameStep2)
  } else {
    save(IndexOfBestTrial, DEGList, CovariateResultsDF, NSQMatrix,
           file = FileNameStep2)
  }
  