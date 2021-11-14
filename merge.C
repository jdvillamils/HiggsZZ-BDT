//General combination of test and train samples
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

int merge(void){
    
    
    
    
    //----------------------
        //Train
    //---------------------
    /*{
    
    TChain* sigchain = new TChain("SignalTrain");
    TChain* bacchain = new TChain("BackgroundTrain");
    bacchain->AddFile("DataBackground.root");
    sigchain->AddFile("DataSignal.root");
    
    Float_t foursignaltrain;
    Float_t fourbacktrain;
    Float_t ptsignaltrain;
    Float_t ptbacktrain;
    sigchain->SetBranchAddress("FourLepSystemMTrain", &foursignaltrain);
    sigchain->SetBranchAddress("FourLepSystemptTrain", &ptsignaltrain);
    bacchain->SetBranchAddress("FourLepSystemMTrain", &fourbacktrain);
    bacchain->SetBranchAddress("FourLepSystemptTrain", &ptbacktrain);
    
    
    TFile *target = new TFile("DataTrain.root","RECREATE");
    
    TTree *signal = new TTree("Signal","Signal from samples");
    TTree *background = new TTree("Background","Background from samples");
    
    Float_t fourtrainsignal;
    Float_t fourtrainback;
    Float_t pttrainsignal;
    Float_t pttrainback;
    
    signal->Branch("FourLepSystemM", &fourtrainsignal);
    signal->Branch("FourLepSystempt", &pttrainsignal);
    background->Branch("FourLepSystemM", &fourtrainback);
    background->Branch("FourLepSystempt", &pttrainback);
    
    int nentries, nbytes, k;

    nentries = (Int_t)sigchain->GetEntries();
    for (k = 0; k < nentries; k++)
    {
        nbytes = sigchain->GetEntry(k);
        fourtrainsignal=foursignaltrain;
        pttrainsignal=ptsignaltrain;
        signal->Fill();
    }
    
    int nentriesb, nbytesb, j;

    nentriesb = (Int_t)bacchain->GetEntries();
    for (j = 0; j < nentriesb; j++)
    {
        nbytesb = bacchain->GetEntry(j);
        fourtrainback=fourbacktrain;
        pttrainback=ptbacktrain;
        background->Fill();
    }
    signal->Write();
    background->Write();
    target->Close();
    return 0;
//  }*/
    
    
    //----------------
        //Test
    //-----------------
    
 //  /*{    
    TChain* sigchain = new TChain("SignalTest");
    TChain* bacchain = new TChain("BackgroundTest");
    bacchain->AddFile("DataBackground.root");
    sigchain->AddFile("DataSignal.root");
    
    Float_t foursignaltest;
    Float_t fourbacktest;
    Float_t ptsignaltest;
    Float_t ptbacktest;
    sigchain->SetBranchAddress("FourLepSystemMTest", &foursignaltest);
    sigchain->SetBranchAddress("FourLepSystemptTest", &ptsignaltest);
    bacchain->SetBranchAddress("FourLepSystemMTest", &fourbacktest);
    bacchain->SetBranchAddress("FourLepSystemptTest", &ptbacktest);
    
    
    TFile *target = new TFile("DataTest.root","RECREATE");
    
    TTree *signal = new TTree("Signal","Signal from samples");
    TTree *background = new TTree("Background","Background from samples");
    
    Float_t fourtestsignal;
    Float_t fourtestback;
    Float_t pttestsignal;
    Float_t pttestback;
    
    signal->Branch("FourLepSystemM", &fourtestsignal);
    signal->Branch("FourLepSystempt", &pttestsignal);
    background->Branch("FourLepSystemM", &fourtestback);
    background->Branch("FourLepSystempt", &pttestback);
    
    int nentries, nbytes, k;

    nentries = (Int_t)sigchain->GetEntries();
    for (k = 0; k < nentries; k++)
    {
        nbytes = sigchain->GetEntry(k);
        fourtestsignal=foursignaltest;
        pttestsignal=ptsignaltest;
        signal->Fill();
    }
    
    int nentriesb, nbytesb, j;

    nentriesb = (Int_t)bacchain->GetEntries();
    for (j = 0; j < nentriesb; j++)
    {
        nbytesb = bacchain->GetEntry(j);
        fourtestback=fourbacktest;
        pttestback=ptbacktest;
        background->Fill();
    }
    signal->Write();
    background->Write();
    target->Close();
    return 0;
//   }*/
}
    