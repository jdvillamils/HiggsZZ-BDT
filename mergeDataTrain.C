#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <time.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <iostream>
#include <stdio.h>

int mergeDataTrain(void){
    
    
    
    
    //----------------------
        //Train
    //---------------------

    
    TChain* sigchain = new TChain("SignalTrain");
    TChain* bacchain = new TChain("BackgroundTrain");
    bacchain->AddFile("RootFiles/DataBackground.root");
    sigchain->AddFile("RootFiles/DataSignal.root");
    
    Float_t foursignaltrain;
    Float_t fourbacktrain;
    Float_t ptsignaltrain;
    Float_t ptbacktrain;
    Float_t invmassz1;
    Float_t invmassz2;
    Float_t invmassz1b;
    Float_t invmassz2b;
    Float_t fourrap;
    Float_t fourrapb;
    Float_t foure;
    Float_t foureb;
    Float_t wstrain;
    Float_t wbtrain;
    sigchain->SetBranchAddress("FourLepSystemMTrain", &foursignaltrain);
    sigchain->SetBranchAddress("FourLepSystemptTrain", &ptsignaltrain);
    sigchain->SetBranchAddress("WeightSignalTrain", &wstrain);
    bacchain->SetBranchAddress("FourLepSystemMTrain", &fourbacktrain);
    bacchain->SetBranchAddress("FourLepSystemptTrain", &ptbacktrain);
    bacchain->SetBranchAddress("WeightBackTrain", &wbtrain);
    
    sigchain->SetBranchAddress("InvMassZ1Train", &invmassz1);
    bacchain->SetBranchAddress("InvMassZ1Train", &invmassz1b);
    sigchain->SetBranchAddress("InvMassZ2Train", &invmassz2);
    bacchain->SetBranchAddress("InvMassZ2Train", &invmassz2b);
    sigchain->SetBranchAddress("FourLepRapidityTrain", &fourrap);
    bacchain->SetBranchAddress("FourLepRapidityTrain", &fourrapb);
    sigchain->SetBranchAddress("FourLepSystemETrain", &foure);
    bacchain->SetBranchAddress("FourLepSystemETrain", &foureb);
    
    
    TFile *target = new TFile("RootFiles/DataTrain.root","RECREATE");
    
    TTree *signal = new TTree("Signal","Signal from samples");
    TTree *background = new TTree("Background","Background from samples");
    
    Float_t fourtrainsignal;
    Float_t fourtrainback;
    Float_t pttrainsignal;
    Float_t pttrainback;
    Float_t invmassz1sig;
    Float_t invmassz1bac;
    Float_t invmassz2sig;
    Float_t invmassz2bac;
    Float_t fourrapsig;
    Float_t fourrapbac;
    Float_t fouresig;
    Float_t fourebac;
    Float_t ws;
    Float_t wb;
    
    signal->Branch("FourLepSystemM", &fourtrainsignal);
    signal->Branch("FourLepSystempt", &pttrainsignal);
    signal->Branch("InvMassZ1", &invmassz1sig);
    signal->Branch("InvMassZ2", &invmassz2sig);
    signal->Branch("FourLepRapidity", &fourrapsig);
    signal->Branch("FourLepSystemE", &fouresig);
    signal->Branch("Weight", &ws);
    background->Branch("FourLepSystemM", &fourtrainback);
    background->Branch("FourLepSystempt", &pttrainback);
    background->Branch("InvMassZ1", &invmassz1bac);
    background->Branch("InvMassZ2", &invmassz2bac);
    background->Branch("FourLepRapidity", &fourrapbac);
    background->Branch("FourLepSystemE", &fourebac);
    background->Branch("Weight", &wb);
    
    int nentries, nbytes, k;

    nentries = (Int_t)sigchain->GetEntries();
    for (k = 0; k < nentries; k++)
    {
        nbytes = sigchain->GetEntry(k);
        fourtrainsignal=foursignaltrain;
        pttrainsignal=ptsignaltrain;
        invmassz1sig=invmassz1;
        invmassz2sig=invmassz2;
        fourrapsig=fourrap;
        fouresig=foure;
        ws=wstrain;
        signal->Fill();
    }
    
    int nentriesb, nbytesb, j;

    nentriesb = (Int_t)bacchain->GetEntries();
    for (j = 0; j < nentriesb; j++)
    {
        nbytesb = bacchain->GetEntry(j);
        fourtrainback=fourbacktrain;
        pttrainback=ptbacktrain;
        invmassz1bac=invmassz1b;
        invmassz2bac=invmassz2b;
        fourrapbac=fourrapb;
        fourebac=foureb;
        wb=wbtrain;
        background->Fill();
    }
    signal->Write();
    background->Write();
    target->Close();
    return 0;
}
    