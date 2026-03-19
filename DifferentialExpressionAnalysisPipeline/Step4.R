## Script for Step 4

## This step takes a minute to complete. This step compares the DEGs found when 
## including demographic/clinical covariates or not including demographic/clinical 
## covariates. As discussed in our publication, only RNAs found to be differentially 
## expressed in both cases can be confidently identified as differentially expressed. 
## This step also finds RNAs that do not have an adjusted p-value < 0.05 but are still 
## likely to be differentially expressed (see the Methods section of our publication).

  
  ##------------------------------------------------------------------------------  
  ## Seeing which DEGs (when using demographic/clinical covariates) were also significant when not using demographic/clinical covariates:
  
  load(paste(Folder, '/', StringForOutputFileName, 'OutputStep3NoDCCovariates.RData', sep = ""))
  NoCovariatesResultsDF = ResultsDF
  NCCoefficientMatrix = CoefficientMatrix
  NCStandardErrorMatrix = StandardErrorMatrix
  
  load(paste(Folder, '/', StringForOutputFileName, 'OutputStep3.RData', sep = ""))
  CovariatesResultsDF = ResultsDF
  CCoefficientMatrix = CoefficientMatrix
  CStandardErrorMatrix = StandardErrorMatrix
  
  ## NC stands for No Covariates:
  NCEnsemblIDs = rownames(NoCovariatesResultsDF)
  NCSignificantIndices = which(NoCovariatesResultsDF$padj < 0.05)
  NCSignificantEnsemblIDs = NCEnsemblIDs[NCSignificantIndices]
  
  ## C stands for Covariates:
  CEnsemblIDs = rownames(CovariatesResultsDF)
  CSignificantIndices = which(CovariatesResultsDF$padj < 0.05)
  CSignificantEnsemblIDs = CEnsemblIDs[CSignificantIndices]
  
  ## This creates the column in Supplementary Tables 1A-3A called "Consistently Differentially Expressed?":
  ConsistentDELogicals = CSignificantEnsemblIDs %in% NCSignificantEnsemblIDs
  ConsistentEnsemblIDs = CSignificantEnsemblIDs[ConsistentDELogicals]
  InconsistentEnsemblIDs = CSignificantEnsemblIDs[!ConsistentDELogicals]
  
  NumberOfCRNAs = length(CEnsemblIDs)
  NAArray = rep('NA', NumberOfCRNAs) 
  
  ##------------------------------------------------------------------------------
  
  ## Here, the code inspects the DEGs when not using demographic/clinical covariates to see if  
  ## coefficients change only minimally when demographic/clinical covariates are added.
  ## This section creates the column in Supplementary Tables 1B-3B called
  ## "DE with Demographic/Clinical Covariates?":
  
  ## Beta refers to the coefficient estimated by DESeq2 for the outcome variable (where the outcome variable is TB disease 
  ## vs. no disease, or TB progressor vs. non-progressor, or infected vs. uninfected). Each element in 
  ## these arrays is a beta corresponding to each of the RNAs:
  NCBetas = NCCoefficientMatrix[NCSignificantEnsemblIDs,'CaseOrControl_B_vs_A']
  CBetas = CCoefficientMatrix[NCSignificantEnsemblIDs,'CaseOrControl_B_vs_A']
  
  ## These are the standard errors (SE) of the betas:
  NCSE = NCStandardErrorMatrix[NCSignificantEnsemblIDs,'SE_CaseOrControl_B_vs_A']
  CSE = CStandardErrorMatrix[NCSignificantEnsemblIDs,'SE_CaseOrControl_B_vs_A']
  
  ## Q is the percent changes in betas when demographic/clinical covariates are added to the DE model:
  Q = 100*(CBetas - NCBetas)/NCBetas
  
  ## This is the standard errors in the percent changes:
  dQ = abs(Q)*sqrt((CSE/CBetas)^2 + (NCSE/NCBetas)^2 )
  
  ExtremeQ = Q-dQ
  ## We consider the RNA to be differentially Expressed (DE) or likely DE if the ExtremeQ is less negative than -10:
  IndicesOfDE = which(ExtremeQ > -10)
  IndicesOfNotDE = which(ExtremeQ <= -10)
  DEEnsemblIDs = NCSignificantEnsemblIDs[IndicesOfDE]
  NotDEEnsemblIDs = NCSignificantEnsemblIDs[IndicesOfNotDE]
  
  ##--------------
  ## If the coefficient does decrease but the RNA is still found to be differentially expressed (p-adj < 0.05),
  ## then we include the RNA as DE (and remove the RNA from the "Not DE" list):
  ExtraDEIndices = which((NotDEEnsemblIDs %in% ConsistentEnsemblIDs)==TRUE)
  ExtraDEEnsemblIDs = NotDEEnsemblIDs[ExtraDEIndices]
  
  Indices = which((NotDEEnsemblIDs %in% ExtraDEEnsemblIDs) == FALSE)
  NotDEEnsemblIDs = NotDEEnsemblIDs[Indices]
  
  ##---------------
  ## Here we segregate the DEGs from the likely DEGs:
  LikelyDELogicals = DEEnsemblIDs %in% ConsistentEnsemblIDs
  LikelyDEEnsemblIDs = DEEnsemblIDs[which(LikelyDELogicals==FALSE)]
  print(LikelyDEEnsemblIDs)

  NumberOfLikelyDERNAs = length(LikelyDEEnsemblIDs)
  print('NumberOfLikelyDERNAs')
  print(NumberOfLikelyDERNAs)
  

##------------------------------------------------------------------------------
## Now we assemble the final DEG lists:

for (DCCovariatesOption in c(0, 1)) {

  if (DCCovariatesOption == 1) {
    FileNameStep3 = paste(Folder,
                          '/',
                          StringForOutputFileName,
                          'OutputStep3.RData',
                          sep = "")
  } else if (DCCovariatesOption == 0) {
    FileNameStep3 = paste(Folder,
                          '/',
                          StringForOutputFileName,
                          'OutputStep3NoDCCovariates.RData',
                          sep = "")
  }
  ## This loads the DEG list:
  load(FileNameStep3)
  
  
  MyEnsIDs = rownames(ResultsDF)
  DETable = cbind(MyEnsIDs, ResultsDF)
  colnames(DETable)[1] = "ENSEMBL"
  
  MyFullNames = bitr(
    MyEnsIDs,
    fromType = "ENSEMBL",
    toType = c("SYMBOL", "GENENAME"),
    OrgDb = "org.Hs.eg.db"
  )
  TableForPaper = merge(DETable, MyFullNames, by = "ENSEMBL", all.x = TRUE)
  
  
  
  
  
  
  ##-------------------------
  ## Looking at the likely differentially expressed RNAs:
  FullNamesLikelyDEEnsemblIDs = bitr(
    LikelyDEEnsemblIDs,
    fromType = "ENSEMBL",
    toType = c("SYMBOL", "GENENAME"),
    OrgDb = "org.Hs.eg.db"
  )
  
  
  
  
  
  ##-------------------------
  ## Removing duplicated rows:
  Duplicates = which(duplicated(TableForPaper[, 1]) == TRUE)
  print("These indices are duplicates:")
  print(Duplicates)
  TableForPaperNoDuplicates = TableForPaper[!duplicated(TableForPaper[, 1]), ]
  
  ##-------------------------
  ## Ordering by adjusted p-value:
  OrderedIndices = order(TableForPaperNoDuplicates[, 'padj'])
  OrderedTableForPaper = TableForPaperNoDuplicates[OrderedIndices,]
  
  ##-------------------------
  ## Making the right number of significant figures:
  OrderedTableForPaper[, 2:7] = signif(OrderedTableForPaper[, 2:7], 2)
  
  for (q in 2:7) {
    OrderedTableForPaper[, q] =  formatC(
      OrderedTableForPaper[, q],
      digits = 2,
      format = "G",
      flag = " "
    )
  }
  
  
  FinalTableForPaper = OrderedTableForPaper[, c(1, 8, 9, 6, 7, 3, 4, 5, 2)]
  rownames(FinalTableForPaper) = FinalTableForPaper[, 1]
  
  ##-------------------------
  if (DCCovariatesOption == 1) {
    FinalTableForPaper[, 10] = NAArray
    colnames(FinalTableForPaper)[10] = 'Consistently Differentially Expressed?'
    FinalTableForPaper[ConsistentEnsemblIDs, 'Consistently Differentially Expressed?'] = 'Yes'
    FinalTableForPaper[InconsistentEnsemblIDs, 'Consistently Differentially Expressed?'] = 'No'
    
  }
  ##-------------------------
  if (DCCovariatesOption == 0) {
    NumberOfNCRNAs = length(NCEnsemblIDs)
    NAArray = rep('NA', NumberOfNCRNAs)
    
    FinalTableForPaper[, 10] = NAArray
    colnames(FinalTableForPaper)[10] = 'DE with Demographic/Clinical Covariates?'
    FinalTableForPaper[ConsistentEnsemblIDs, 'DE with Demographic/Clinical Covariates?'] = 'Yes'
    FinalTableForPaper[LikelyDEEnsemblIDs, 'DE with Demographic/Clinical Covariates?'] = 'Likely'
    FinalTableForPaper[NotDEEnsemblIDs, 'DE with Demographic/Clinical Covariates?'] = 'No'
    
    FinalTableForPaper[, 11] = NAArray
    colnames(FinalTableForPaper)[11] = 'Percent Change in Coefficient'
    FinalTableForPaper[, 12] = NAArray
    colnames(FinalTableForPaper)[12] = 'Standard Error of Percent Change in Coefficient'
    
    for (a in 1:length(NCSignificantEnsemblIDs)) {
      CurrentEnsID = NCSignificantEnsemblIDs[a]
      FinalTableForPaper[CurrentEnsID, 11] = Q[CurrentEnsID]
      FinalTableForPaper[CurrentEnsID, 12] = dQ[CurrentEnsID]
    }
    
    
  }
  ##----------------------
  
  
  
  if (DCCovariatesOption == 0) {
    FileNameStep4 = paste(StringForFileNameStep4, 'NoDCCovariates.csv',
                          sep = "")
  } else if (DCCovariatesOption == 1) {
    FileNameStep4 = paste(StringForFileNameStep4, '.csv',
                          sep = "")
  }

  write.csv(FinalTableForPaper, FileNameStep4, row.names = FALSE)
  

}
