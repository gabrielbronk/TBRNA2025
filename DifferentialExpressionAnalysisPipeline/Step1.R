## Script for Step 1

## Determining if any covariates are correlated with the outcome 
## (where the outcome is TB progressor vs non progressor, or TB infection vs uninfected, 
## or TB disease vs no disease, or TB infection or disease vs no infection).

##==============================================================================
## Seeing if Sequencing Batch correlates with the outcome:

  SequencingBatchTable = matrix(0,2,6)  
  for (z in 1:NumberOfKidsOfInterest) {
    
    if (Metadata[z,"CaseOrControl"]=='A') {
      CurrentRow = 1
    } else {
      CurrentRow = 2
    }
    CurrentColumn = as.numeric(Metadata[z,"SequencingBatch"]) 
    SequencingBatchTable[CurrentRow, CurrentColumn] = SequencingBatchTable[CurrentRow, CurrentColumn] + 1
  }
  
  ## Performing Chi-Square Test:
  ChiSquareResultSequencing = chisq.test(SequencingBatchTable,simulate.p.value = TRUE)
  FractionsSequencing = SequencingBatchTable[1,]/colSums(SequencingBatchTable)
  FractionControl = length(which(Metadata[,"CaseOrControl"]=='A'))/NumberOfKidsOfInterest
  
  
  ##==============================================================================
  ## Seeing if Library Batch correlates with the outcome:
  
  LibraryBatchTable = matrix(0,2,21)  
  for (z in 1:NumberOfKidsOfInterest) {
    
    if (Metadata[z,"CaseOrControl"]=='A') {
      CurrentRow = 1
    } else {
      CurrentRow = 2
    }
    CurrentColumn = as.numeric(Metadata[z,"LibraryBatch"]) 
    LibraryBatchTable[CurrentRow, CurrentColumn] = LibraryBatchTable[CurrentRow, CurrentColumn] + 1
  }
  
  ZeroIndices = which(colSums(LibraryBatchTable)==0)
  
  ## Performing Chi-Square Test:
  ChiSquareResultLibrary = chisq.test(LibraryBatchTable[,-ZeroIndices],simulate.p.value = TRUE)
  FractionsLibrary = SequencingBatchTable[1,]/colSums(LibraryBatchTable)
  FractionControl = length(which(Metadata[,"CaseOrControl"]=='A'))/NumberOfKidsOfInterest
  
  
##==============================================================================
## Seeing if any covariates correlate with the outcome, 
## Looking at Binary Variables:
##==============================================================================

  IndicesOfTest = which(ArrayOfCaseOrControl == 'B')
  IndicesOfControl = which(ArrayOfCaseOrControl == 'A')
  
  CovariatesArray = c('sex', 'child_hiv_birth', 'ethnicity', 'MaternalDepression', 'MaternalPTSD', 
                      'work', 'Prenatal_alcohol_composite', 'preterm', 'Food_insecurity_ANC',
                      'education_binary', 'household_income_binary', 'Prenatal_smoking_composite')
  
  ## This is one level of the covariates:
  CovariatesArray0 = c('Male', 'HIV unexposed', 'Black african', 0, 0, 
                       'working', 'No Exposure', 'Full term', 'Percieved food secure',
                       'completed secondary', 'middle', 'No')
  
  ## This is the other level of the covariates:
  CovariatesArray1 = c('Female', 'HIV exposed uninfected', 'Mixed ancestory', 1, 1, 
                       'not working', 'Exposure', 'Early term', 'Percieved food insecure',
                       'some secondary', 'poor', 'Yes')
  
  NumberOfCovsForTable = length(CovariatesArray)
  
  FractionTable = matrix(0,4,NumberOfCovsForTable)
  FractionDF = as.data.frame(FractionTable)
  colnames(FractionDF) = CovariatesArray
  rownames(FractionDF) = c('Fraction in Asymptomatic', 'Fraction in TB', 'Chi Square P Value', 'Mean of Fractions')
  
  for (d in 1:NumberOfCovsForTable) {
    NumberOf0inTest = length(which(Metadata[IndicesOfTest,CovariatesArray[d]] == CovariatesArray0[d]))
    NumberOf1inTest = length(which(Metadata[IndicesOfTest,CovariatesArray[d]] == CovariatesArray1[d]))
    FractionOf1InTest = NumberOf1inTest/(NumberOf0inTest + NumberOf1inTest)
    
    NumberOf0inControl = length(which(Metadata[IndicesOfControl,CovariatesArray[d]] == CovariatesArray0[d]))
    NumberOf1inControl = length(which(Metadata[IndicesOfControl,CovariatesArray[d]] == CovariatesArray1[d]))  
    FractionOf1InControl = NumberOf1inControl/(NumberOf0inControl + NumberOf1inControl)
    
    print(NumberOf0inTest)
    print(NumberOf1inTest)
    print(NumberOf0inControl)
    print(NumberOf1inControl)
    print(NumberOf0inTest + NumberOf1inTest + NumberOf0inControl + NumberOf1inControl)
    
    ##----------------------------------------------------------------------------
    ## Performing Chi-Square Test:
    
    CurrentTable = matrix(0,2,2)
    CurrentTable[1,1] = NumberOf0inControl
    CurrentTable[2,1] = NumberOf1inControl
    CurrentTable[1,2] = NumberOf0inTest
    CurrentTable[2,2] = NumberOf1inTest
    ChiSquareResult = chisq.test(CurrentTable,simulate.p.value = TRUE)
    
    MeanFraction = (FractionOf1InTest + FractionOf1InControl)/2
    
    FractionDF[1,d] = FractionOf1InControl
    FractionDF[2,d] = FractionOf1InTest
    FractionDF[3,d] = ChiSquareResult$p.value
    FractionDF[4,d] = MeanFraction
  }                     
  

##==============================================================================
## Seeing if any covariates correlate with the outcome, 
## Looking at Ordinal Variables (i.e. more than two levels):
##==============================================================================

  ## RTI stands for RNA Timepoint of Interest (since these are the kids' blood samples at the timepoint of interest):
  MetadataRTI = Metadata[1:(NumberOfTestKids + NumberOfControlKids),]
  ArrayOfCaseOrControlNarrowed = c( rep(1,NumberOfTestKids), 
                                    rep(0,NumberOfControlKids) )
  
  CovariateArray = c('household_income', 'education')
  
  KendallPValueTable = matrix(1000,length(CovariateArray),1)
  for (i in 1:length(CovariateArray)){
    
    KendallObject = cor.test(ArrayOfCaseOrControlNarrowed,  MetadataRTI[,CovariateArray[i]], method=c("kendall"))
    KendallPValueTable[i,1] = KendallObject$p.value
  }
  rownames(KendallPValueTable) = CovariateArray
  colnames(KendallPValueTable) = 'p Value'

##==============================================================================
## Seeing if any covariates correlate with the outcome,
## this time looking at continuous covariates:
##==============================================================================

  MetadataRTI[,'CaseOrControl'] = ArrayOfCaseOrControlNarrowed
  
  CovariateArray = c('birth_weight', 'gestation_delivery','number_in_household')
  #CovariateArray = c(TechnicalVariables, 'RIN')
  ## This initializes the AnovaPValueTable by filling it with 1000's:
  AnovaPValueTable = matrix(1000,length(CovariateArray),1)
  for (w in 1:length(CovariateArray)) {
    MyFormula = as.formula( paste('CaseOrControl ~ ', as.character(CovariateArray[w]), sep='') )
    AnovaObject = aov(MyFormula, data = MetadataRTI)
    Output = summary(AnovaObject)
    AnovaPValueTable[w] = Output[[1]][1,5]
  }
  rownames(AnovaPValueTable) = CovariateArray
  colnames(AnovaPValueTable) = 'p Value'
  
##==============================================================================
## Seeing if any covariates correlate with the outcome,
## this time looking at the technical variables (which are continuous):
##==============================================================================
  
  CovariateArray = c(TechnicalVariables, 'RIN')
  ## This initializes the AnovaPValueTable by filling it with 1000's:
  ## TV stands for Technical Variables:
  TTestTableTV = matrix(1000,length(CovariateArray),1)
  
  ## I use t-test now but I didn't bother to change the name of the object from ANOVAPValueTable
  for (w in 1:length(CovariateArray)) {
    MyFormula = as.formula( paste(as.character(CovariateArray[w]), ' ~ CaseOrControl', sep='') )
    TTestObject = t.test(MyFormula, data = MetadataRTI, alternative = c("two.sided"))
    TTestTableTV[w] = TTestObject$p.value
  }
  rownames(TTestTableTV) = CovariateArray
  colnames(TTestTableTV) = 'p Value'
  
  
