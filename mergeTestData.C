//Combines the signal and background trees for future testing 

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

int mergeTestData(void){
    
    TChain* sigchain = new TChain("Signal");
    TChain* bacchain = new TChain("Background");
    sigchain->AddFile("RootFiles/DataTest.root");
    bacchain->AddFile("RootFiles/DataTest.root");
    
    Float_t fours;
    Float_t fourb;
    Float_t pts;
    Float_t ptb;
    Float_t invz1;
    Float_t invz1b;
    Float_t invz2;
    Float_t invz2b;
    Float_t fourrap;
    Float_t fourrapb;
    Float_t foure;
    Float_t foureb;
    
    sigchain->SetBranchAddress("FourLepSystemM", &fours);
    bacchain->SetBranchAddress("FourLepSystemM", &fourb);
    sigchain->SetBranchAddress("FourLepSystempt", &pts);
    bacchain->SetBranchAddress("FourLepSystempt", &ptb);
    sigchain->SetBranchAddress("InvMassZ1", &invz1);
    bacchain->SetBranchAddress("InvMassZ1", &invz1b);
    sigchain->SetBranchAddress("InvMassZ2", &invz2);
    bacchain->SetBranchAddress("InvMassZ2", &invz2b);
    sigchain->SetBranchAddress("FourLepRapidity", &fourrap);
    bacchain->SetBranchAddress("FourLepRapidity", &fourrapb);
    sigchain->SetBranchAddress("FourLepSystemE", &foure);
    bacchain->SetBranchAddress("FourLepSystemE", &foureb);
    
    
    TFile *target = new TFile("RootFiles/DataTestmerge.root","RECREATE");
    
    TTree *data = new TTree("Data","Data from samples, excluding Train data");
    
    Float_t fourdata;
    Float_t ptdata;
    Float_t invmassz1data;
    Float_t invmassz2data;
    Float_t fourrapdata;
    Float_t fouredata;
    
    data->Branch("FourLepSystemM", &fourdata);
    data->Branch("FourLepSystempt", &ptdata);
    data->Branch("InvMassZ1", &invmassz1data);
    data->Branch("InvMassZ2", &invmassz2data);
    data->Branch("FourLepRapidity", &fourrapdata);
    data->Branch("FourLepSystemE", &fouredata);
    
    int nentries, nbytes, k;

    nentries = (Int_t)sigchain->GetEntries();
    for (k = 0; k < nentries; k++)
    {
        nbytes = sigchain->GetEntry(k);
        fourdata=fours;
        ptdata=pts;
        invmassz1data=invz1;
        invmassz2data=invz2;
        fourrapdata=fourrap;
        fouredata=foure;
        data->Fill();
    }
    
    int nentriesb, nbytesb, j;

    nentriesb = (Int_t)bacchain->GetEntries();
    for (j = 0; j < nentriesb; j++)
    {
        nbytesb = bacchain->GetEntry(j);
        fourdata=fourb;
        ptdata=ptb;
        invmassz1data=invz1b;
        invmassz2data=invz2b;
        fourrapdata=fourrapb;
        fouredata=foureb;
        data->Fill();
    }
    data->Write();
    target->Close();
    return 0;

    
    
}