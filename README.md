# Higgs -> ZZ using a Boosted Decision Tree, with Atlas Open Data.

## Introduction

This framework makes use of [ROOT](https://root.cern.ch/) in C language, and the [Atlas Open Data](http://opendata.atlas.cern). The purpouse of this framework is to train a [TMVA](https://root.cern/manual/tmva/) tool, called Boosted Decision tree, to identify and discriminate the signal and background of the [H->ZZ](http://opendata.atlas.cern/release/2020/documentation/physics/FL2.html) process from simulated and real samples. 

## Analysis

Files SignalToRoot and BackgroundToRoot store the signal and background data from the [AOD](http://opendata.atlas.cern/release/2020/documentation/datasets/intro.html) into root files.

The created root files are reorganised and separated with the mergeDataTest, mergeDataTrain and mergeTestData, in order to have different samples for training the BDT, and to test it. 

Once is done, TMVAClassification uses the train samples to train the BDT, and TMVAClassificationApplication uses the datasets from the training, and the test samples, to discriminate the data with the criteria.

## About files in this Framework

+ The SignalToRoot and BackgroundToRoot files take data from AOD and, filter and store it in root files in the RootFiles folder.
+ The mergeDataTest, mergeDataTrain and mergeTestData files, prepare the root files for being processed by TMVA.
+ TMVAClassification will create the datasets for the testing the samples with TMVACLassificationApplication (and  TMVAClassificationApplicationData for data) This returns the multivariate discriminant, for reject background and keep signal.
+ BDTAnalysis is the notebook for plotting and changing the parameters of TMVA from the resulting samples.

## To run:

(Only necessary if changes are made to the analysis and filtering section, elsewhere, root files are ready upload to be read, including the results from the TMVA training.)

+ 1) run the  Samples Processing
```
source run_processing_samples.sh 
```
+ 2) Train the TMVA and then test it
```
source run_tmva.sh
```

## To visualize the results

Go to the BDTAnalysis notebook, change the bdt parameters and plot the different variables, with the already processed root files, or new loaded files.



