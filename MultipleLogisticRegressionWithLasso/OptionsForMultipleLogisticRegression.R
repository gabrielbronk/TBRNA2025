
##============================================
## Choose which options you want to use here:
##============================================

JobArrayOption = 0
## Set JobArrayOption = 0 if not running this script in parallel (this is just for a single fold of outer cross validation). This is just if you want to see how this script runs on RStudio on your laptop, for instance.
## Set JobArrayOption = 1 if running this script with parallel jobs using a Job Array on a slurm job scheduler on a computer cluster. Parallelization is necessary because each of the 20 folds of outer cross validation takes about 10 hours to run. Therefore, you will have to launch this script as parallel jobs in order to run a full nested cross validation in a reasonable amount of time.

SaveFileOption = 1
## Use SaveFileOption = 0 to not save the output.
## Use SaveFileOption = 1 to save the output in an RData file.

NameOfOutputFile = 'TBNestedCrossValidation'
## Every output file will be named this, but at the end of the name of each file it will also include "Job1", "Job2", "Job3" etc. depending on which parallel job it is.

AnalysisOption = 1
## Use AnalysisOption = 1 for predicting progressors vs. Mtb-infected non-progressors.
## Use AnalysisOption = 2 for predicting progressors vs everyone else (Mtb-infected non-progressors and uninfected individuals).
## Use AnalysisOption = 3 for predicting Mtb-infected (both progressors and non-progressors) vs uninfected individuals.

CrossValidationOption = 1
## Use CrossValidationOption = 1 to do nested cross validation, which assesses the performance of the model.
## Use CrossValidationOption = 2 to do just inner cross validation. Here we do inner cross validation on the whole dataset to do hyperparameter tuning and model fitting, with no data held out for outer folds.

## Number of folds for outer cross validation (OCV):
kOCV = 20 

## Number of folds for inner cross validation:
kin = 10  

## OCV stands for Outer Cross Validation. This is the number of OCVs to run in each parallel job:
NumberOfOCVsPerJob = 1   



##------------------------------------------------------------------------------
## Hyperparameters for Lasso Regularization:
##------------------------------------------------------------------------------

## Here are the hyperparameters of the multiple logistic regression model with lasso regularization.
## The jhyperparameter values below are the ones tried in our study. These values do not need to be changed
## (but they can be changed if you wish).

## CutOff is the probability value at which to consider a model prediction to be a positive prediction. 
## In the manuscript, we refer to this as the classification threshold T.
## For example, if CutOff is 0.5, then any prediction greater than 0.5 is considered positive
## (i.e. TB progressor), and any prediction less than 0.5 is considered negative (non-progressor). 
## We Vary the CutOff rather than setting it to just 0.5 because we may want to use a different CutOff
## if we are maximizing parformance metrics such as positive predictive value (PPV) or negative predictive
## value (NPV). If you are only interested in AUC, the script will skip using the CutOffArray.

## Lambda is a parameter that controls the how many variables (in our case RNAs) are included in the
## model. See the glmnet R package manual for more information. NumberOfLambdas and MinimumLambdaRatio
## determine the lambda values to try. The values below were chosen so that the model would use no more
## than about 6 RNAs (because having a model with a small number of RNAs would be easier to use in a
## clinical setting). 

## U is the number of RNAs to start with, which then gets narrowed down by the lasso regularization. For example, we take the top 40 (or top 60 or top 80 etc) most significantly differentially expressed RNAs and feed them into glmnet.

CutOffArray = seq(0.02, 0.98, by = 0.02) 
NumberOfLambdas = 13  
MinimumLambdaRatio = 0.717
UArray = c(40, 60, 80, 100, 120, 130, 150, 175, 200, 250, 300, 400, 500, 700, 1000)


##------------------------------------------------------------------------------

NumberOfCutOffs = length(CutOffArray)
NumberOfUs = length(UArray)
MaxNumberOfInputRNAs = max(UArray)

TotalNumberOfJobs = kOCV / NumberOfOCVsPerJob