# Step 3

  
  NormCountsMatrix = counts(dds, normalized = TRUE, replaced = TRUE)
  
  IndsToLookAt = which(ResultsDFRENoWeighting$padj < PadjThreshold)
  NumberOfIndsToLookAt = length(IndsToLookAt)
  NumberOfRNAs = dim(ResultsDFRENoWeighting)[1]
  
  ##==============================================================================
  ## This sees how many outlier samples there are for each DE RNA:
  ##==============================================================================
  
  CooksIndexArrayMatrix = matrix(0,NumberOfRNAs,2)
  for (i in 1:NumberOfRNAs) {
    CurrentEnsemblID = rownames(ResultsDFRENoWeighting[i,])
    CurrentCooksArray = CooksMatrix[CurrentEnsemblID, ]
    
    CooksIndexArrayMatrix[i,1] = CurrentEnsemblID
    CooksIndexArrayMatrix[i,2] = length(which(CurrentCooksArray > quantile_value))
  }
  rownames(CooksIndexArrayMatrix) = CooksIndexArrayMatrix[,1]
  
  ##==============================================================================
  ## This makes a table to compare the DE results with or without outlier removal:
  ##==============================================================================
  if  (OutlierWeightingOption == 1) {
    DEGNoWeightingInds = which(ResultsDFRENoWeighting$padj < PadjThreshold)
    DEGEnsembls = rownames(ResultsDFRENoWeighting)[DEGNoWeightingInds]
    ComparisonTable = cbind(ResultsDFRENoWeighting[DEGEnsembls,c(2,5,6)],
                            ResultsDFRE[DEGEnsembls,c(2,5,6)],
                            CooksIndexArrayMatrix[DEGEnsembls,2])
    FlaggingColumn = matrix(0,length(DEGEnsembls),1)
    
    for (i in 1:length(DEGEnsembls)) {
      MyNumber = ComparisonTable[i,1]*ComparisonTable[i,4]
      if (is.na(MyNumber) == FALSE){
        if (MyNumber < 0) {
          FlaggingColumn[i] = 1
        }
      }
    }
    ComparisonTable$Flagging = FlaggingColumn
    
    MyEnsIDs = rownames(ComparisonTable)
    DETable = cbind(MyEnsIDs,ComparisonTable)
    colnames(DETable)[1] = "ENSEMBL"
    
    MyFullNames = bitr(MyEnsIDs, fromType="ENSEMBL", toType=c("SYMBOL", "GENENAME"), OrgDb="org.Hs.eg.db")
    ComparisonTableWithNames = merge(DETable,MyFullNames,by="ENSEMBL",all.x=TRUE)
    ## Removing duplicated rows:
    Duplicates = which(duplicated(ComparisonTableWithNames[,1])==TRUE)
    CTable = ComparisonTableWithNames[!duplicated(ComparisonTableWithNames[,1]),]
      ##CTable is an abbreviation for ComparisonTableWithNamesNoDuplicates
    
    
    
    
    
    
    ##============================================================================
    ## This section makes the DE table for the paper (i.e. Table 2):
    ##============================================================================
    
    MainTable = as.data.frame(matrix(0,1,6))
    NumOfAssembledDEGs = 0
    for (i in 1:NumberOfIndsToLookAt){
      PadjReg = CTable[i,4]
      PadjOR = CTable[i,7] ## OR stands for Outliers Removed
      
      MyNumber = ComparisonTable[i,1]*ComparisonTable[i,4]
      if (is.na(MyNumber) == FALSE){ ## This ensures the log fold changes are not NA.
      if (MyNumber > 0) { ## This ensures the fold changes are in the same direction regardless of whether outliers are removed.
      if ((is.na(PadjReg) == FALSE) & (is.na(PadjOR) == FALSE)  ){ ## This ensures the adjusted p-values are not NA. 
        
        if ((PadjReg < PadjThreshold) & (PadjOR < PadjThreshold)) {
          NumOfAssembledDEGs = NumOfAssembledDEGs + 1
          MainTable[NumOfAssembledDEGs,1] = CTable[i,1]
          MainTable[NumOfAssembledDEGs,2:3] = CTable[i,10:11]
          if (PadjReg > PadjOR){
            MainTable[NumOfAssembledDEGs,4:6] = CTable[i,2:4]
          } else {
            MainTable[NumOfAssembledDEGs,4:6] = CTable[i,5:7]
          }
        }
        
      }
      }
      }
    }
      
    IndsForOrdering = order(MainTable[,6])
    OrderedMainTable = MainTable[IndsForOrdering,]
    ## Making the right number of significant figures:
    OrderedMainTable[,4:6] = signif(OrderedMainTable[,4:6],2)
    for (q in 4:6) {
      OrderedMainTable[,q] =  formatC(OrderedMainTable[,q], digits = 2, format = "G", flag = " ")
    }   
    colnames(OrderedMainTable) = c("Ensembl ID",	"Gene Symbol",	"Gene Name",	"log2(Fold Change)", "P-Value",	"Adjusted P-Value")	
    
    ## Saving the DE Results:
    write.csv(OrderedMainTable, FileNameStep3Main, row.names = FALSE)
    
    
    
    
    
    ##============================================================================
    ## This section makes a table for the Supplementary Information, showing
    ## DE analysis results when outliers are not excluded:
    ##============================================================================
    MyEnsIDs = rownames(ResultsDFRENoWeighting)
    DETable = cbind(MyEnsIDs, ResultsDFRENoWeighting, CooksIndexArrayMatrix[,2])
    colnames(DETable)[1] = "ENSEMBL"
    
    MyFullNames = bitr(
      MyEnsIDs,
      fromType = "ENSEMBL",
      toType = c("SYMBOL", "GENENAME"),
      OrgDb = "org.Hs.eg.db"
    )
    TableForPaper = merge(DETable, MyFullNames, by = "ENSEMBL", all.x = TRUE)
    
    ##-------------------------
    ## Removing duplicated rows:
    Duplicates = which(duplicated(TableForPaper[, 1]) == TRUE)
    print("These indices are duplicates:")
    print(Duplicates)
    TableForPaperNoDuplicates = TableForPaper[!duplicated(TableForPaper[, 1]), ]
    
    ##-------------------------
    ## Ordering by adjusted p-value:
    OrderedIndices = order(TableForPaperNoDuplicates[, 'pvalue'])
    OrderedTableForPaper = TableForPaperNoDuplicates[OrderedIndices,]
    
    ##-------------------------
    ## Making the right number of significant figures:
    OrderedTableForPaper[, 2:7] = signif(OrderedTableForPaper[, c(2:7)], 2)
    
    for (q in 2:7) {
      OrderedTableForPaper[, q] =  formatC(
        OrderedTableForPaper[, q],
        digits = 2,
        format = "G",
        flag = " "
      )
    }
    
    FinalTableForPaper = OrderedTableForPaper[, c(1, 9, 10, 6, 7, 3, 4, 5, 2, 8)]
    rownames(FinalTableForPaper) = FinalTableForPaper[, 1]
    
    
    ## Saving the DE Results:
    write.csv(FinalTableForPaper, FileNameStep3Supplement, row.names = FALSE)
    
    
    
    
    
    ##============================================================================
    ## This section makes a table for the Supplementary Information, showing
    ## DE analysis results when outliers are excluded. OR stands for Outliers Removed:
    ##============================================================================
    MyEnsIDsOR = rownames(ResultsDFRE)
    DETableOR = cbind(MyEnsIDsOR, ResultsDFRE)
    colnames(DETableOR)[1] = "ENSEMBL"
    
    MyFullNamesOR = bitr(
      MyEnsIDsOR,
      fromType = "ENSEMBL",
      toType = c("SYMBOL", "GENENAME"),
      OrgDb = "org.Hs.eg.db"
    )
    TableForPaperOR = merge(DETableOR, MyFullNamesOR, by = "ENSEMBL", all.x = TRUE)
    
    ##-------------------------
    ## Removing duplicated rows:
    DuplicatesOR = which(duplicated(TableForPaperOR[, 1]) == TRUE)
    print("These indices are duplicates:")
    print(DuplicatesOR)
    TableForPaperNoDuplicatesOR = TableForPaperOR[!duplicated(TableForPaperOR[, 1]), ]
    
    ##-------------------------
    ## Ordering by adjusted p-value:
    OrderedIndicesOR = order(TableForPaperNoDuplicatesOR[, 'padj'])
    OrderedTableForPaperOR = TableForPaperNoDuplicatesOR[OrderedIndicesOR,]
    
    ##-------------------------
    ## Making the right number of significant figures:
    OrderedTableForPaperOR[, 2:7] = signif(OrderedTableForPaperOR[, c(2:7)], 2)
    
    for (q in 2:7) {
      OrderedTableForPaperOR[, q] =  formatC(
        OrderedTableForPaperOR[, q],
        digits = 2,
        format = "G",
        flag = " "
      )
    }
    
    FinalTableForPaperOR = OrderedTableForPaperOR[, c(1, 8, 9, 6, 7, 3, 4, 5, 2)]
    rownames(FinalTableForPaperOR) = FinalTableForPaperOR[, 1]
    
    
    ## Saving the DE Results:
    write.csv(FinalTableForPaperOR, FileNameStep3SupplementOutliersRemoved, row.names = FALSE)
    
  }
  
  
  
