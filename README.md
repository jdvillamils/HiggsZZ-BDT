# Higgs -> ZZ using a Boosted Decision Tree, with Atlas Open Data.

## Introduction

This framework makes use of [ROOT](https://root.cern.ch/) in C language, and the [Atlas Open Data](http://opendata.atlas.cern). The purpouse of this framework is to train a [TMVA](https://root.cern/manual/tmva/) tool, called Boosted Decision tree, to identify and discriminate the signal and background of the [H->ZZ](http://opendata.atlas.cern/release/2020/documentation/physics/FL2.html) process from simulated and real samples. 

## Analysis

Files SignalToRoot and BackgroundToRoot store the signal and background data from the [AOD](http://opendata.atlas.cern/release/2020/documentation/datasets/intro.html) into root files.

The created root files are reorganised and separated with the mergeDataTest, mergeDataTrain and mergeTestData, in order to have different samples for training the BDT, and to test it. 

Once is done, TMVAClassification uses the train samples to train the BDT, and TMVAClassificationApplication uses the datasets from the training, and the test samples, to discriminate the data with the criteria.

## To run:

Only neccesary if changes are made to the analysis and filtering section, elsewhere, root files are ready to be read, including the results from the TMVA training.

+ 1) Run the SignalToRoot and BackgroundToRoot files, this will create the root files in the RootFiles folder.
+ 2) Run the mergeDataTest, mergeDataTrain and mergeTestData files, in order to prepare the root files for being processed by TMVA.
+ 3) Run TMVAClassification, this will create the datasets for the testing.
+ 4) Run TMVACLassificationApplication. This results in the finel multivariate discriminant, for reject background and keep signal.
+ 5) The plots and information of the process can be seen in the BDTAnalysis notebook.

To add to the previous data (simulated), the real data meassured by ATLAS:

+ 1) Run DataToRoot.C
+ 2) In TMVAClassificationApplication, comment lines 269, 313 and 427 and uncomment lines 270, 314 and 428. Then run the script.
+ 3) Add the TMVAppData.root into BDTAnalysis notebook

## To visualize the results

Go to the BDTAnalysis notebook, change the bdt parameters and plot the different variables, with the already processed root files.

P.S.
All files are run from the terminal as
```
root -l File.C 
```

