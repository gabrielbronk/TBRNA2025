# TBRNA2025
Scripts for "Host RNA biomarkers in cord blood predict Mycobacterium tuberculosis infection or TB disease during childhood"

Prior to using these scripts, please first refer to the Supplementary Methods for our paper. 

To run the same analyses as performed in our paper, refer to the Read Me files in each subfolder. The subfolders are:

1)	MultipleLogisticRegressionWithLasso

This folder contains the scripts for the model that predicts whether children will progress to TB disease if they become infected with Mycobacterium tuberculosis. This model is a multiple logistic regression model with lasso regularization. 

2)	DifferentialExpressionAnalysisPipeline

This folder contains the pipeline to perform our differential expression analysis methodology as well as the logistic regression that computes relative risk estimates.

3)	Data

This folder contains all data needed to run the above analyses, except RNA expression count matrices.

