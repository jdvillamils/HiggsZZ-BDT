echo Training TMVA

root -q TMVAClassification.C

echo Testing TMVA over simulated samples

root -q TMVAClassificationApplication.C

echo Testing TMVA over data samples

root -q TMVAClassificationApplicationData.C

echo All data has been processed succesfully