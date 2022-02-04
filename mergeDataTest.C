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

int mergeDataTest(void){
  
    //----------------
        //Test
    //-----------------
 
    TChain* sigchain = new TChain("SignalTest");
    TChain* bacchain = new TChain("BackgroundTest");
    sigchain->AddFile("RootFiles/DataSignal1.root");
    sigchain->AddFile("RootFiles/DataSignal2.root");
    sigchain->AddFile("RootFiles/DataSignal3.root");
    bacchain->AddFile("RootFiles/DataBackground1.root");
    bacchain->AddFile("RootFiles/DataBackground2.root");
    bacchain->AddFile("RootFiles/DataBackground3.root");
    bacchain->AddFile("RootFiles/DataBackground4.root");
    bacchain->AddFile("RootFiles/DataBackground5.root");
    
    Float_t foursignaltest;
    Float_t fourbacktest;
    Float_t ptsignaltest;
    Float_t ptbacktest;
    Float_t invmassz1;
    Float_t invmassz2;
    Float_t invmassz1b;
    Float_t invmassz2b;
    Float_t fourrap;
    Float_t fourrapb;
    Float_t foure;
    Float_t foureb;
    Float_t wstest;
    Float_t wbtest;
    sigchain->SetBranchAddress("FourLepSystemMTest", &foursignaltest);
    sigchain->SetBranchAddress("FourLepSystemptTest", &ptsignaltest);
    sigchain->SetBranchAddress("WeightSignalTest", &wstest);
    bacchain->SetBranchAddress("FourLepSystemMTest", &fourbacktest);
    bacchain->SetBranchAddress("FourLepSystemptTest", &ptbacktest);
    bacchain->SetBranchAddress("WeightBackTest", &wbtest);
    
    sigchain->SetBranchAddress("InvMassZ1Test", &invmassz1);
    bacchain->SetBranchAddress("InvMassZ1Test", &invmassz1b);
    sigchain->SetBranchAddress("InvMassZ2Test", &invmassz2);
    bacchain->SetBranchAddress("InvMassZ2Test", &invmassz2b);
    sigchain->SetBranchAddress("FourLepRapidityTest", &fourrap);
    bacchain->SetBranchAddress("FourLepRapidityTest", &fourrapb);
    sigchain->SetBranchAddress("FourLepSystemETest", &foure);
    bacchain->SetBranchAddress("FourLepSystemETest", &foureb);
    
    
    TFile *target = new TFile("RootFiles/DataTest.root","RECREATE");
    
    TTree *signal = new TTree("Signal","Signal from samples");
    TTree *background = new TTree("Background","Background from samples");
    
    Float_t fourtestsignal;
    Float_t fourtestback;
    Float_t pttestsignal;
    Float_t pttestback;
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
    
    signal->Branch("FourLepSystemM", &fourtestsignal);
    signal->Branch("FourLepSystempt", &pttestsignal);
    signal->Branch("InvMassZ1", &invmassz1sig);
    signal->Branch("InvMassZ2", &invmassz2sig);
    signal->Branch("FourLepRapidity", &fourrapsig);
    signal->Branch("FourLepSystemE", &fouresig);
    signal->Branch("Weight", &ws);
    background->Branch("FourLepSystemM", &fourtestback);
    background->Branch("FourLepSystempt", &pttestback);
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
        fourtestsignal=foursignaltest;
        pttestsignal=ptsignaltest;
        invmassz1sig=invmassz1;
        invmassz2sig=invmassz2;
        fourrapsig=fourrap;
        fouresig=foure;
        ws=wstest;
        signal->Fill();
    }
    
    int nentriesb, nbytesb, j;

    nentriesb = (Int_t)bacchain->GetEntries();
    for (j = 0; j < nentriesb; j++)
    {
        nbytesb = bacchain->GetEntry(j);
        fourtestback=fourbacktest;
        pttestback=ptbacktest;
        invmassz1bac=invmassz1b;
        invmassz2bac=invmassz2b;
        fourrapbac=fourrapb;
        fourebac=foureb;
        wb=wbtest;
        background->Fill();
    }
    signal->Write();
    background->Write();
    target->Close();
    return 0;

}
    