echo Starting Signal and Background processing 

root -q SignalToRoot1.C

root -q SignalToRoot2.C

root -q SignalToRoot3.C

root -q BackgroundToRoot1.C

root -q BackgroundToRoot2.C

root -q BackgroundToRoot3.C

root -q BackgroundToRoot4.C

root -q BackgroundToRoot5.C

echo Starting Data processing

root -q DataToRoot.C

echo Preparing Training and testing samples

root -q mergeDataTest.C

root -q mergeDataTrain.C

root -q mergeTestData.C

echo Done! Samples are ready to be processed by TMVA 


