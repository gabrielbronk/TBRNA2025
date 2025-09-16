## Script for Step 8:

## Using the RNAs found in the final DEG list from step (7), 
## for each of these RNAs, I perform Bayesian logistic regression 
## to determine the relative risk of the outcome (e.g. progression to TB disease) 
## for individuals with the highest expression levels vs individuals with the 
## lowest expression levels.

## Note: THIS SCRIPT ASSUMES THERE IS ONLY ONE COVARIATE IN CovariateNames. This 
## step 8 script is not written for more than one confounding covariate.

FittingOption = 2
## Use Option 1 to just calculate the odds ratios. This takes just a minute.
## Use Option 2 to calculate the odds ratios and the relative risks (this takes about 24 hours).

ParallelOption = 1
## Use 0 to not run parallel jobs.
## Use 1 to run parallel jobs.

SavingOption = 1
## Use Option 0 to do nothing.
## Use Option 1 to save the output.

library('data.table')

## WARNING: You may need to adjust the range and spacing of the parameters I and B 
## (or if you make a sequence that is over the log space, then make sure 
## to multiply the posterior by the appropriate volume of parameter space). 
## The bottom of this script will output information that tells you whether
## the range and intervals need adjustment:

NumberOfIsForZero = 5000
NumberOfBsForZero = 15000

NumberOfIsForOne = 250
NumberOfBsForOne = 1000      
NumberOfSs = 300

NumOfBins = 300

## This loads up the CorrectedCountsMatrix, which is the counts matrix for all blood samples. 
## The counts matrix was corrected for technical variables and batch (see Step2Fewer.Rmd or Step2Sameness.Rmd):
NewCountsMatrixOriginal = CorrectedCountsMatrix[,c(BarcodesOfTestKids, BarcodesOfControlKids)]
MetadataRTIOriginal = Metadata[c(BarcodesOfTestKids, BarcodesOfControlKids),]
YOriginal = c(rep(1,length(BarcodesOfTestKids)), rep(0,length(BarcodesOfControlKids)))

AllDEGsResultsDF = as.matrix(matrix(10^-6,NumberOfSurvivingDEGs,10)) #4#24#18
AllDEGsResultsDFGLM = as.matrix(matrix(10^-6,NumberOfSurvivingDEGs,12)) 
## 10^-6 is just a place holder.



##==============================================================================
## Fitting the Data with Logistic Regression and a Jeffreys Prior, and then 
## Obtaining the Confidence Intervals on the Relative Risk:
##==============================================================================

if (ParallelOption == 0) {
  LoopNumbers = 1:NumberOfSurvivingDEGs
} else if (ParallelOption == 1) {
  TaskID = Sys.getenv("SLURM_ARRAY_TASK_ID") 
  TaskIDNumber = as.integer(TaskID)
  LoopNumbers = TaskIDNumber
}


tic()
for (DEGNum in LoopNumbers){ 
  
  ##----------------------------------------------------------------------------
  ## Preparation for model fitting:
  ##----------------------------------------------------------------------------
  
  CovariateIndices = which(MatrixOfConfoundingFinal[DEGNum,]==1)
  NumberOfCovariatesTried = length(CovariateArray)
  NumberOfCovariatesForRNA = length(CovariateIndices)
  
  CovariateNames = c()
  if (NumberOfCovariatesForRNA > 0) {
  for (u in 1:NumberOfCovariatesForRNA) {
    if (CovariateClassArray[CovariateIndices[u]] == 0) {
      CovariateNames[u] = paste(CovariateArray[CovariateIndices[u]],'Numeric',sep="") 
    } else if (CovariateClassArray[CovariateIndices[u]] == 1) {
      CovariateNames[u] = CovariateArray[CovariateIndices[u]]
    }
  }
  }
  
  
##------------------------------------------------------------------------------
    ## This subsection removes samples if they don't have measurements of the covariates 
    ## that I am controling for in this regression:
    d = as.numeric(MetadataRTIOriginal[,CovariateNames])
    NACovariateInds = which(is.na(d)==TRUE)
    if (length(NACovariateInds) > 0) {
      d = d[-NACovariateInds]
      NewCountsMatrix = NewCountsMatrixOriginal[,-NACovariateInds]
      MetadataRTI = MetadataRTIOriginal[-NACovariateInds,]
      Y = YOriginal[-NACovariateInds]
      NumberOfTestKids = sum(Y)
      NumberOfControlKids = length(Y) - sum(Y)
    } else {
      NewCountsMatrix = NewCountsMatrixOriginal
      MetadataRTI = MetadataRTIOriginal
      Y = YOriginal
      NumberOfTestKids = sum(Y)
      NumberOfControlKids = length(Y) - sum(Y)
    }
    NumOfKids = length(Y)
    ##-----------------------
    
    d1 = d[1:NumberOfTestKids]
    d0 = d[(NumberOfTestKids+1):(NumberOfTestKids + NumberOfControlKids)]
  
  
    XData = t(NewCountsMatrix[SurvivingDEGNames[DEGNum],])
    x = as.numeric(XData[,1])
    x1 = x[1:NumberOfTestKids]
    x0 = x[(NumberOfTestKids+1):(NumberOfTestKids + NumberOfControlKids)]
    xmax = max(x)  
    xmin = min(x) 
    xmean = mean(x)
    Xtest = mean(x1)
    Xcontrol = mean(x0)
    
    OrderedX = x[order(x)]
    x2 = OrderedX[ceiling(0.02*NumOfKids)]
    x3 = OrderedX[ceiling(0.03*NumOfKids)]
    x97 = OrderedX[ceiling(0.97*NumOfKids)]
    x98 = OrderedX[ceiling(0.98*NumOfKids)]
##------------------------------------------------------------------------------
  
  
  
##==============================================================================  
## Using the BRGLM R package to perform logistic regression:
##==============================================================================
  Y = c(rep(1,NumberOfTestKids), rep(0,NumberOfControlKids))
  Data = cbind(data.frame('x' = x, 'y' = Y), MetadataRTI[,CovariateNames] )
  colnames(Data) = c('x','y',CovariateNames)
  
  ## Creating Regression Formula:
  CurrentString = "y ~ x"
  if (NumberOfCovariatesForRNA > 0) {
    for (u in 1:NumberOfCovariatesForRNA) {
      CurrentString = paste(CurrentString, "+", CovariateNames[u])
    }
  }
  
  BayesianRegressionModel = brglm(as.formula(CurrentString), family = binomial('logit'), data = Data, 
                                  method = "brglm.fit")
  
  ConfIntObject = confint(BayesianRegressionModel)
  
  Int = as.numeric(BayesianRegressionModel$coefficients['(Intercept)'])
  IntBottom = ConfIntObject['(Intercept)', 1]
  IntTop = ConfIntObject['(Intercept)', 2]
  IntTopDistance = IntTop - Int
  IntBottomDistance = Int - IntBottom 
  
  Beta = as.numeric(BayesianRegressionModel$coefficients['x'])
  BetaBottom = ConfIntObject['x', 1]
  BetaTop = ConfIntObject['x', 2]
  BetaTopDistance = BetaTop - Beta
  BetaBottomDistance = Beta - BetaBottom
  
  
  if (NumberOfCovariatesForRNA > 0) {
  BetaCoef = as.numeric(BayesianRegressionModel$coefficients[3])
  BetaCoefBottom = ConfIntObject[3, 1]
  BetaCoefTop = ConfIntObject[3, 2]
  BetaCoefTopDistance = BetaCoefTop - BetaCoef
  BetaCoefBottomDistance = BetaCoef - BetaCoefBottom
  }
  
  LogORMax = log10(exp(Beta*(xmax-xmin)))
  LogORMaxBottom = log10(exp(BetaBottom*(xmax-xmin)))
  LogORMaxTop = log10(exp(BetaTop*(xmax-xmin)))
  LogORTopDistance = LogORMaxTop - LogORMax
  LogORBottomDistance = LogORMax - LogORMaxBottom 
  
  AllDEGsResultsDFGLM[DEGNum,1] = LogORMax
  AllDEGsResultsDFGLM[DEGNum,2] = LogORMaxBottom
  AllDEGsResultsDFGLM[DEGNum,3] = LogORMaxTop
  
  AllDEGsResultsDFGLM[DEGNum,4] = Int
  AllDEGsResultsDFGLM[DEGNum,5] = IntBottom
  AllDEGsResultsDFGLM[DEGNum,6] = IntTop
   
  AllDEGsResultsDFGLM[DEGNum,7] = Beta
  AllDEGsResultsDFGLM[DEGNum,8] = BetaBottom
  AllDEGsResultsDFGLM[DEGNum,9] = BetaTop
  
  if (NumberOfCovariatesForRNA > 0) {
  AllDEGsResultsDFGLM[DEGNum,10] = BetaCoef
  AllDEGsResultsDFGLM[DEGNum,11] = BetaCoefBottom
  AllDEGsResultsDFGLM[DEGNum,12] = BetaCoefTop  
  }


  
  
  
  
##==============================================================================
## My custom fitting script starts here:
##==============================================================================
  
if (FittingOption == 2) {

##==============================================================================
## For 0 Covariates:
##==============================================================================
if (NumberOfCovariatesForRNA == 0) {  
  
  NumberOfIs = NumberOfIsForZero
  NumberOfBs = NumberOfBsForZero
  
  set.seed(DEGNum)
  IArray = runif(NumberOfIs, (IntBottom-IntBottomDistance), (IntTop + IntTopDistance) )
  BArray = runif(NumberOfBs, (BetaBottom-BetaBottomDistance), (BetaTop + BetaTopDistance) )
 
  
  ##----------------------------------------------------------------------------
  ## This section defines functions that will be part of the Jeffreys Prior:
  ##----------------------------------------------------------------------------
  
  f0 = function(x, I, B) {
    sum( ( 2 + exp(I + B*x) + exp(-(I + B*x)) )^(-1) )
  }
  
  f1 = function(x, I, B) {
    sum(x*( ( 2 + exp(I + B*x) + exp(-(I + B*x)) )^(-1) ))
  }
  
  f2 = function(x, I, B) {
    sum((x^2)*( ( 2 + exp(I + B*x) + exp(-(I + B*x)) )^(-1) ))
  }
  
  ##----------------------------------------------------------------------------
  ## This section evaluates the values of the posterior probability distribution 
  ## for relative risk for gene expression xmax vs gene expression xmin (normalized, 
  ## batch and technical variable-corrected counts):
  ##----------------------------------------------------------------------------
  
  ## ECO stands for "Exponentiated numbers Cut Off" (this cut off improves the
  ## script's ability to computer numbers numerically):
  ECO = 50
  NumberOfParameterSets = NumberOfIs*NumberOfBs
  
  ParametersAndLogPosteriorList = list()
  ## There will be a separate list entry for each I.
  
  tic()
  
  ##----------------------------------------------------------------------------
  ## Calculating the values of loglikelihoods and priors:
  ##----------------------------------------------------------------------------
  for (INum in 1:NumberOfIs) {  
    ParametersAndLogPosteriorDF = as.data.frame(matrix(0,NumberOfBs,4))#
    I = IArray[INum]
    
    for (BNum in 1:NumberOfBs) {  
      B = BArray[BNum] 
      
      x1NumbersArray = -(I + B*x1)
      HighIndices1 = which(x1NumbersArray > ECO)
      Highx1 = x1[HighIndices1]
      LowIndices1 = which(x1NumbersArray <= ECO)
      Lowx1 = x1[LowIndices1]
      
      x0NumbersArray = I + B*x0
      HighIndices0 = which(x0NumbersArray > ECO)
      Highx0 = x0[HighIndices0]
      LowIndices0 = which(x0NumbersArray <= ECO)
      Lowx0 = x0[LowIndices0]
      
      LogLikelihood = -sum( log(1 + exp(-(I + B*Lowx1))) ) -sum(-(I + B*Highx1))  -sum( log(1 + exp(I + B*Lowx0)) ) -sum(I + B*Highx0)
      Prior = sqrt( f0(x, I, B)*f2(x, I, B) - (f1(x, I, B))^2 )

      LogPosterior = LogLikelihood + log(Prior) 
      RRMax = (exp(I) + exp(-B*xmin))/(exp(I) + exp(-B*xmax))
      RR97 = (exp(I) + exp(-B*x3))/(exp(I) + exp(-B*x97))
      RR98 = (exp(I) + exp(-B*x2))/(exp(I) + exp(-B*x98))
      
      
      
      ParametersAndLogPosteriorDF[BNum,1] = LogPosterior
      ParametersAndLogPosteriorDF[BNum,2] = RRMax 
      ParametersAndLogPosteriorDF[BNum,3] = RR97 
      ParametersAndLogPosteriorDF[BNum,4] = RR98
    }
    ParametersAndLogPosteriorList[[INum]] = ParametersAndLogPosteriorDF
    
  }
  TimeForFitting = toc()
  print("Time For Fitting:")
  print(TimeForFitting)
  
  ##----------------------------------------------------------------------------
  ## Determining the value of the BestLogPosterior so that it can be subtracted 
  ## from every logposterior value:
  ##----------------------------------------------------------------------------
  tic()
  
  BestLogPosteriorArray = matrix(0,1,NumberOfIs)
  for (INum in 1:NumberOfIs) { 
    CurrentDF = ParametersAndLogPosteriorList[[INum]]
    BestLogPosteriorArray[INum] = max(CurrentDF[,1], na.rm = TRUE) 
  }
  BestLogPosterior = max(BestLogPosteriorArray)
  
  ## Subtracting off the BestLogPosterior from each log posterior value, and then replacing the logposteriors with
  ## posteriors, and also writing the log10(RR)s instead of the RRs:
  NormalizationFactor = 0
  for (INum in 1:NumberOfIs) { 
    CurrentDF = ParametersAndLogPosteriorList[[INum]]
    CurrentPosteriorArray = exp(CurrentDF[,1] - BestLogPosterior)
    NormalizationFactor = NormalizationFactor + sum(CurrentPosteriorArray)
    
    CurrentDF[,1] = CurrentPosteriorArray
    CurrentDF[,2] = log10(CurrentDF[,2])
    CurrentDF[,3] = log10(CurrentDF[,3])
    CurrentDF[,4] = log10(CurrentDF[,4])
    colnames(CurrentDF) = c('Posterior', 'LogRRMax','LogRR97','LogRR98')
    ParametersAndLogPosteriorList[[INum]] = CurrentDF
  }
 
  TimeForAssembling = toc()
  print('TimeForAssembling')
  print(TimeForAssembling)
  
  ##============================================================================
  ## This section finds the 95% CIs for RRMax by integrating over all of 
  ## parameter space to get the probability density function for RRmax:
  ##============================================================================
  tic()
  
  LowerLimit = LogORMaxBottom - 3*LogORBottomDistance
  UpperLimit = LogORMaxTop + 3*LogORTopDistance
  BinSize = (UpperLimit - LowerLimit)/NumOfBins
  
  MyProbDistRRMax = matrix(0,1,NumOfBins)
  MyProbDistRR97 = matrix(0,1,NumOfBins)
  MyProbDistRR98 = matrix(0,1,NumOfBins)
  for (INum in 1:NumberOfIs) { 
    CurrentDF = ParametersAndLogPosteriorList[[INum]]
    for (BNum in 1:NumberOfBs) { 
      CurrentPosterior = CurrentDF[BNum,'Posterior']/NormalizationFactor 
      
      CurrentLogRRMax = CurrentDF[BNum,'LogRRMax']  
      CurrentBin = ceiling((CurrentLogRRMax-LowerLimit)/BinSize)
      MyProbDistRRMax[CurrentBin] = MyProbDistRRMax[CurrentBin] + CurrentPosterior
      
      CurrentLogRR97 = CurrentDF[BNum,'LogRR97']  
      CurrentBin = ceiling((CurrentLogRR97-LowerLimit)/BinSize)
      MyProbDistRR97[CurrentBin] = MyProbDistRR97[CurrentBin] + CurrentPosterior
      
      CurrentLogRR98 = CurrentDF[BNum,'LogRR98']  
      CurrentBin = ceiling((CurrentLogRR98-LowerLimit)/BinSize)
      MyProbDistRR98[CurrentBin] = MyProbDistRR98[CurrentBin] + CurrentPosterior
    }
  }

  Bins = seq(LowerLimit,(UpperLimit-BinSize/2),BinSize) + BinSize/2

  
  ##---------
  IndicesOfDecreasingPosterior = order(MyProbDistRRMax, decreasing = TRUE)
  ## RPDF stands for RRMax and Posterior Data Frame:
  OrderedProbDist = MyProbDistRRMax[IndicesOfDecreasingPosterior]
  OrderedRRs = Bins[IndicesOfDecreasingPosterior]
  
  CumulativeSumProbDist = cumsum(OrderedProbDist)
  IndicesToKeep = which(CumulativeSumProbDist < 0.95)
  ## FI stands for From Integration:
  RRMaxNFCITopFI = max(OrderedRRs[IndicesToKeep])
  RRMaxNFCIBottomFI = min(OrderedRRs[IndicesToKeep])
  RRMaxBestFI = OrderedRRs[1]
  
  
  ##---------
  IndicesOfDecreasingPosterior = order(MyProbDistRR97, decreasing = TRUE)
  ## RPDF stands for RRMax and Posterior Data Frame:
  OrderedProbDist = MyProbDistRR97[IndicesOfDecreasingPosterior]
  OrderedRRs = Bins[IndicesOfDecreasingPosterior]
  
  CumulativeSumProbDist = cumsum(OrderedProbDist)
  IndicesToKeep = which(CumulativeSumProbDist < 0.95)
  ## FI stands for From Integration:
  RR97NFCITopFI = max(OrderedRRs[IndicesToKeep])
  RR97NFCIBottomFI = min(OrderedRRs[IndicesToKeep])
  RR97BestFI = OrderedRRs[1]
  
  
  ##---------
  IndicesOfDecreasingPosterior = order(MyProbDistRR98, decreasing = TRUE)
  ## RPDF stands for RRMax and Posterior Data Frame:
  OrderedProbDist = MyProbDistRR98[IndicesOfDecreasingPosterior]
  OrderedRRs = Bins[IndicesOfDecreasingPosterior]
  
  CumulativeSumProbDist = cumsum(OrderedProbDist)
  IndicesToKeep = which(CumulativeSumProbDist < 0.95)
  ## FI stands for From Integration:
  RR98NFCITopFI = max(OrderedRRs[IndicesToKeep])
  RR98NFCIBottomFI = min(OrderedRRs[IndicesToKeep])
  RR98BestFI = OrderedRRs[1]
  
  
  
  print(RRMaxBestFI)
  print(RRMaxNFCIBottomFI)
  print(RRMaxNFCITopFI)
  
  TimeForFittingPartThree = toc()
  print("Time For Integration:")
  print(TimeForFittingPartThree)
  ##--------------------------------------------------------------------------
  
  
  ResultsDF = as.matrix(matrix(0,3,1))
  ResultsDF[1,] = RRMaxBestFI
  ResultsDF[2,] = c(RRMaxNFCITopFI)
  ResultsDF[3,] = c(RRMaxNFCIBottomFI)

  colnames(ResultsDF) = c('Relative Risk Max')
  rownames(ResultsDF) = c('MAP', 'Top of 95% CI', 'Bottom of 95% CI')
  
  AllDEGsResultsDF[DEGNum,] = c(RRMaxBestFI, RRMaxNFCIBottomFI, RRMaxNFCITopFI, 
                                RR97BestFI, RR97NFCIBottomFI, RR97NFCITopFI,
                                RR98BestFI, RR98NFCIBottomFI, RR98NFCITopFI,
                                BinSize)

}
  
 
  
  ##============================================================================
  ## For 1 Covariate:
  ##============================================================================
  if (NumberOfCovariatesForRNA == 1) {  
    
    NumberOfIs = NumberOfIsForOne
    NumberOfBs = NumberOfBsForOne
    
    
    set.seed(DEGNum)
    IArray = runif(NumberOfIs, (IntBottom-IntBottomDistance), (IntTop + IntTopDistance) )
    BArray = runif(NumberOfBs, (BetaBottom-BetaBottomDistance), (BetaTop + BetaTopDistance) )
    SArray = runif(NumberOfSs, (BetaCoefBottom-BetaCoefBottomDistance), (BetaCoefTop + BetaCoefTopDistance) )
    
    NumberOfIs = length(IArray)
    NumberOfBs = length(BArray)
    NumberOfSs = length(SArray)
    
   
    ##----------------------------------------------------------------------------
    ## This section defines functions that will be part of the Jeffreys Prior:
    ##----------------------------------------------------------------------------
    
    f0 = function(x, d, I, B, S) {
      sum( ( 2 + exp(I + B*x + S*d) + exp(-(I + B*x + S*d)) )^(-1) )
    }
    
    f1 = function(x, d, I, B, S) {
      sum(x*( ( 2 + exp(I + B*x + S*d) + exp(-(I + B*x + S*d)) )^(-1) ))
    }
    
    f2 = function(x, d, I, B, S) {
      sum((x^2)*( ( 2 + exp(I + B*x + S*d) + exp(-(I + B*x + S*d)) )^(-1) ))
    }
    
    ##----------------------------------------------------------------------------
    ## This section evaluates the values of the posterior probability distribution 
    ## for relative risk for gene expression xmax vs gene expression xmin (normalized, 
    ## batch and technical variable-corrected counts):
    ##----------------------------------------------------------------------------
    
    ## ECO stands for "Exponentiated numbers Cut Off" (this cut off improves the
    ## script's ability to computer numbers numerically):
    ECO = 50
  
    NumberOfParameterSets = NumberOfIs*NumberOfBs*NumberOfSs
    ## RR stands for Relative Risk:
    ParametersAndLogPosteriorList = list()
    ## There will be a separate list entry for each I.
    
    tic()
    
    ##----------------------------------------------------------------------------
    ## Calculating the values of loglikelihoods and priors:
    ##----------------------------------------------------------------------------
    
    j = 1
    for (INum in 1:NumberOfIs) {  
      I = IArray[INum]
      
      for (SNum in 1:NumberOfSs) {  
        S = SArray[SNum] 
        
        ParametersAndLogPosteriorDF = as.data.frame(matrix(0,NumberOfBs,4))
        for (BNum in 1:NumberOfBs) { 
            
            B = BArray[BNum] 
       
        
        ##----------------------------------------
        x1NumbersArray = -(I + B*x1 + S*d1)
        HighIndices1 = which(x1NumbersArray > ECO)
        Highx1 = x1[HighIndices1]
        Highd1 = d1[HighIndices1]
        LowIndices1 = which(x1NumbersArray <= ECO)
        Lowx1 = x1[LowIndices1]
        Lowd1 = d1[LowIndices1]
        
        x0NumbersArray = I + B*x0 + S*d0
        HighIndices0 = which(x0NumbersArray > ECO)
        Highx0 = x0[HighIndices0]
        Highd0 = d0[HighIndices0]
        LowIndices0 = which(x0NumbersArray <= ECO)
        Lowx0 = x0[LowIndices0]
        Lowd0 = d0[LowIndices0]
        
        LogLikelihood = (-sum( log(1 + exp(-(I + B*Lowx1 + S*Lowd1 ))) ) -sum(-(I + B*Highx1 + S*Highd1)) 
                        -sum( log(1 + exp(I + B*Lowx0 + S*Lowd0)) ) -sum(I + B*Highx0 + S*Highd0) )

        Prior = sqrt( f0(x, d, I, B, S)*f2(x, d, I, B, S) - (f1(x, d, I, B, S))^2 )

        
        LogPrior = log(Prior)
        LogPosterior = LogLikelihood + LogPrior 
        
        


        RRMax = sum( (1 + exp( -(I + B*xmax + S*d ) ))^(-1) )/sum( (1 + exp( -(I + B*xmin + S*d ) ))^(-1) )
        RR97 = sum( (1 + exp( -(I + B*x97 + S*d ) ))^(-1) )/sum( (1 + exp( -(I + B*x3 + S*d ) ))^(-1) )
        RR98 = sum( (1 + exp( -(I + B*x98 + S*d ) ))^(-1) )/sum( (1 + exp( -(I + B*x2 + S*d ) ))^(-1) )


        
        
        ParametersAndLogPosteriorDF[BNum,1] = LogPosterior     
        ParametersAndLogPosteriorDF[BNum,2] = RRMax            
        ParametersAndLogPosteriorDF[BNum,3] = RR97            
        ParametersAndLogPosteriorDF[BNum,4] = RR98           
        
        }
        
        ParametersAndLogPosteriorList[[j]] = ParametersAndLogPosteriorDF
        j = j+1
      }
      
  }
    
    TimeForFittingPartOne = toc()
    print("Time For Fitting One:")
    print(TimeForFittingPartOne)


    ##--------------------------------------------------------------------------
    ##--------------------------------------------------------------------------
    ## Determining the value of the BestLogPosterior so that it can be subtracted 
    ## from every logposterior value:
    ##--------------------------------------------------------------------------
    tic()
    
    BestLogPosteriorArray = matrix(0,1,NumberOfIs*NumberOfSs)
    for (j in 1:(NumberOfIs*NumberOfSs)) { 
      CurrentDF = ParametersAndLogPosteriorList[[j]]
      BestLogPosteriorArray[j] = max(CurrentDF[,1], na.rm = TRUE) 
    }
    BestLogPosterior = max(BestLogPosteriorArray, na.rm = TRUE) #
    
    ## Subtracting off the BestLogPosterior from each log posterior value, and then replacing the logposteriors with
    ## posteriors, and also writing the log10(RR)s instead of the RRs:
    NormalizationFactor = 0
    for (j in 1:(NumberOfIs*NumberOfSs)) { 
      CurrentDF = ParametersAndLogPosteriorList[[j]]
      CurrentPosteriorArray = exp(CurrentDF[,1] - BestLogPosterior)
      NormalizationFactor = NormalizationFactor + sum(CurrentPosteriorArray)
      
      CurrentDF[,1] = CurrentPosteriorArray
      CurrentDF[,2] = log10(CurrentDF[,2])
      CurrentDF[,3] = log10(CurrentDF[,3])
      CurrentDF[,4] = log10(CurrentDF[,4])
      colnames(CurrentDF) = c('Posterior', 'LogRRMax','LogRR97','LogRR98')
      ParametersAndLogPosteriorList[[j]] = CurrentDF
    }
    
    TimeForAssembling = toc()
    print('TimeForAssembling')
    print(TimeForAssembling)
    
    ##============================================================================
    ## This section finds the 95% CIs for RRMax by integrating over all of 
    ## parameter space to get the probability density function for RRmax (this might not actually
    ## change the result found above):
    ##============================================================================
    tic()
    
    LowerLimit = LogORMaxBottom - 3*LogORBottomDistance
    UpperLimit = LogORMaxTop + 3*LogORTopDistance
    BinSize = (UpperLimit - LowerLimit)/NumOfBins
    
    MyProbDistRRMax = matrix(0,1,NumOfBins)
    MyProbDistRR97 = matrix(0,1,NumOfBins)
    MyProbDistRR98 = matrix(0,1,NumOfBins)
    for (j in 1:(NumberOfIs*NumberOfSs)) {
      CurrentDF = ParametersAndLogPosteriorList[[j]]
      for (BNum in 1:NumberOfBs) { 
        CurrentPosterior = CurrentDF[BNum,'Posterior']/NormalizationFactor 
        
        CurrentLogRRMax = CurrentDF[BNum,'LogRRMax']  
        CurrentBin = ceiling((CurrentLogRRMax-LowerLimit)/BinSize)
        MyProbDistRRMax[CurrentBin] = MyProbDistRRMax[CurrentBin] + CurrentPosterior
        
        CurrentLogRR97 = CurrentDF[BNum,'LogRR97']  
        CurrentBin = ceiling((CurrentLogRR97-LowerLimit)/BinSize)
        MyProbDistRR97[CurrentBin] = MyProbDistRR97[CurrentBin] + CurrentPosterior
        
        CurrentLogRR98 = CurrentDF[BNum,'LogRR98']  
        CurrentBin = ceiling((CurrentLogRR98-LowerLimit)/BinSize)
        MyProbDistRR98[CurrentBin] = MyProbDistRR98[CurrentBin] + CurrentPosterior
      }
    }
    
    Bins = seq(LowerLimit,(UpperLimit-BinSize/2),BinSize) + BinSize/2

    
    ##---------
    IndicesOfDecreasingPosterior = order(MyProbDistRRMax, decreasing = TRUE)
    ## RPDF stands for RRMax and Posterior Data Frame:
    OrderedProbDist = MyProbDistRRMax[IndicesOfDecreasingPosterior]
    OrderedRRs = Bins[IndicesOfDecreasingPosterior]
    
    CumulativeSumProbDist = cumsum(OrderedProbDist)
    IndicesToKeep = which(CumulativeSumProbDist < 0.95)
    ## FI stands for From Integration:
    RRMaxNFCITopFI = max(OrderedRRs[IndicesToKeep])
    RRMaxNFCIBottomFI = min(OrderedRRs[IndicesToKeep])
    RRMaxBestFI = OrderedRRs[1]
    
    
    ##---------
    IndicesOfDecreasingPosterior = order(MyProbDistRR97, decreasing = TRUE)
    ## RPDF stands for RRMax and Posterior Data Frame:
    OrderedProbDist = MyProbDistRR97[IndicesOfDecreasingPosterior]
    OrderedRRs = Bins[IndicesOfDecreasingPosterior]
    
    CumulativeSumProbDist = cumsum(OrderedProbDist)
    IndicesToKeep = which(CumulativeSumProbDist < 0.95)
    ## FI stands for From Integration:
    RR97NFCITopFI = max(OrderedRRs[IndicesToKeep])
    RR97NFCIBottomFI = min(OrderedRRs[IndicesToKeep])
    RR97BestFI = OrderedRRs[1]
    
    
    ##---------
    IndicesOfDecreasingPosterior = order(MyProbDistRR98, decreasing = TRUE)
    ## RPDF stands for RRMax and Posterior Data Frame:
    OrderedProbDist = MyProbDistRR98[IndicesOfDecreasingPosterior]
    OrderedRRs = Bins[IndicesOfDecreasingPosterior]
    
    CumulativeSumProbDist = cumsum(OrderedProbDist)
    IndicesToKeep = which(CumulativeSumProbDist < 0.95)
    ## FI stands for From Integration:
    RR98NFCITopFI = max(OrderedRRs[IndicesToKeep])
    RR98NFCIBottomFI = min(OrderedRRs[IndicesToKeep])
    RR98BestFI = OrderedRRs[1]
    
    
    
    print(RRMaxBestFI)
    print(RRMaxNFCIBottomFI)
    print(RRMaxNFCITopFI)
    
    TimeForFittingPartThree = toc()
    print("Time For Integration:")
    print(TimeForFittingPartThree)
    ##--------------------------------------------------------------------------
    
    
   
    
    AllDEGsResultsDF[DEGNum,] = c(RRMaxBestFI, RRMaxNFCIBottomFI, RRMaxNFCITopFI, 
                                  RR97BestFI, RR97NFCIBottomFI, RR97NFCITopFI,
                                  RR98BestFI, RR98NFCIBottomFI, RR98NFCITopFI,
                                  BinSize)
    
   
    
  }  
  
  
  
  
  
   
  
  
  
  
  
  }

}

colnames(AllDEGsResultsDF) = c('Relative Risk Max', 'Bottom of 95% CI', 'Top of 95% CI', 
                               'Relative Risk 97', 'Bottom of 95% CI', 'Top of 95% CI', 
                               'Relative Risk 98', 'Bottom of 95% CI', 'Top of 95% CI', 
                               'Bin Size')

rownames(AllDEGsResultsDF) = SurvivingDEGNames


colnames(AllDEGsResultsDFGLM) = c('Log10 OR Max', 'Bottom of 95% CI', 'Top of 95% CI', 
                               'Intercept', 'Bottom of 95% CI', 'Top of 95% CI',
                               'Beta1', 'Bottom of 95% CI', 'Top of 95% CI',
                               'Beta2', 'Bottom of 95% CI', 'Top of 95% CI')

rownames(AllDEGsResultsDFGLM) = SurvivingDEGNames

## You look at this after just running FittingOption = 1 (i.e. just the BRGLM R package fit) and
## it allows you to get a sense of what the parameter value ranges are that should be 
# tried in my custom fitting section (FittingOption = 2):
## The entries are the most extreme values of the following:
## I Bottom 95% CI, I Top 95% CI, B Bottom 95% CI, B Top 95% CI, S Bottom 95% CI, S Top 95% CI
RangeMaxAndMins = c( min(AllDEGsResultsDFGLM[,5]),
                    max(AllDEGsResultsDFGLM[,6]),
                    min(AllDEGsResultsDFGLM[,8]),
                    max(AllDEGsResultsDFGLM[,9]),
                    min(AllDEGsResultsDFGLM[,11]),
                    max(AllDEGsResultsDFGLM[,12]) )
## WARNING: If you see 10^-6 in either of the S columns, then inspect that column of                   
## AllDEGsResultsDFGLM to get the real value because 10^-6 was just a place holder.

if (ParallelOption == 0){
  FileNameForSaving = FileNameStep8
} else {
 
  FileNameForSaving = paste(Folder, '/', StringForOutputFileName, 
                            'OutputStep8ForProbTesting', as.character(TaskIDNumber), '.RData', sep = "")
                            
}


if (SavingOption == 1) {
  save(AllDEGsResultsDF, AllDEGsResultsDFGLM, RangeMaxAndMins,
       MyProbDistRRMax, MyProbDistRR97, MyProbDistRR98, Bins,
       file = FileNameForSaving)
}
