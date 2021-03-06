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



int SignalToRoot3(void){


    //Get Data from repository    
    TString path = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/";
    
    
    //Creates a Chain for store SIGNAL info
    TChain* schain = new TChain("mini");

    //schain->AddFile(path+"MC/mc_345060.ggH125_ZZ4lep.4lep.root");
    //schain->AddFile(path+"MC/mc_341947.ZH125_ZZ4lep.4lep.root");
    schain->AddFile(path+"MC/mc_341964.WH125_ZZ4lep.4lep.root");
    //schain->AddFile(path+"MC/mc_344235.VBFH125_ZZ4lep.4lep.root");//Problematic

    //Define variables to use in SIGNAL
    
    vector<float> *Lep_pts;
    Bool_t TrigEs;
    Bool_t TrigMs;
    UInt_t Lep_ns;
    vector<float> *Lep_phis;
    vector<float> *Lep_etas;
    vector<float> *Lep_Es;
    vector<float> *Lep_ptcone30s;
    vector<float> *Lep_etcone20s;
    vector<unsigned int> *Lep_types;
    vector<int> *Lep_charges;
    vector<float>   *lep_trackd0pvunbiaseds;
    vector<float>   *lep_tracksigd0pvunbiaseds;
    vector<float>   *Lep_z0s;
    UInt_t          jet_ns;
    vector<float>   *jet_pts;
    vector<float>   *jet_etas;
    Float_t         mcWeight;
    Float_t         scaleFactor_PILEUP;
    Float_t         scaleFactor_ELE;
    Float_t         scaleFactor_MUON;
    Float_t         scaleFactor_LepTRIGGER;
   
    
 //------------------------------------------------   
   TBranch        *b_TrigEs;   //!
   TBranch        *b_TrigMs;   //!   //!
   TBranch        *b_Lep_ns;   //!
   TBranch        *b_Lep_pts;   //!
   TBranch        *b_Lep_etas;   //!
   TBranch        *b_Lep_phis;   //!
   TBranch        *b_Lep_Es;   //!
   TBranch        *b_Lep_z0s;   //!
   TBranch        *b_Lep_charges;   //!
   TBranch        *b_Lep_types; 
   TBranch        *b_Lep_ptcone30s;   //!
   TBranch        *b_Lep_etcone20s;   //!
   TBranch        *b_lep_trackd0pvunbiaseds;   //!
   TBranch        *b_lep_tracksigd0pvunbiaseds;   //!
   TBranch        *b_jet_ns;   //!
   TBranch        *b_jet_pts;   //!
   TBranch        *b_jet_etas;
   TBranch        *b_mcWeight;   //!
   TBranch        *b_scaleFactor_PILEUP;   //!
   TBranch        *b_scaleFactor_ELE;   //!
   TBranch        *b_scaleFactor_MUON;   //!
   TBranch        *b_scaleFactor_LepTRIGGER;     //!
//---------------------------------------------------
    
    //Setting branch address for signal and background
    schain->SetBranchAddress("lep_pt", &Lep_pts, &b_Lep_pts);
    schain->SetBranchAddress("trigE", &TrigEs, &b_TrigEs);
    schain->SetBranchAddress("trigM", &TrigMs, &b_TrigMs);
    schain->SetBranchAddress("lep_n", &Lep_ns, &b_Lep_ns);
    schain->SetBranchAddress("lep_phi", &Lep_phis, &b_Lep_phis);
    schain->SetBranchAddress("lep_eta", &Lep_etas, &b_Lep_etas);
    schain->SetBranchAddress("lep_E", &Lep_Es, &b_Lep_Es);
    schain->SetBranchAddress("lep_ptcone30", &Lep_ptcone30s, &b_Lep_ptcone30s);
    schain->SetBranchAddress("lep_etcone20", &Lep_etcone20s, &b_Lep_etcone20s);
    schain->SetBranchAddress("lep_type", &Lep_types, &b_Lep_types);
    schain->SetBranchAddress("lep_charge", &Lep_charges, &b_Lep_charges);
    schain->SetBranchAddress("lep_trackd0pvunbiased", &lep_trackd0pvunbiaseds, &b_lep_trackd0pvunbiaseds);
    schain->SetBranchAddress("lep_tracksigd0pvunbiased", &lep_tracksigd0pvunbiaseds, &b_lep_tracksigd0pvunbiaseds);
    schain->SetBranchAddress("lep_z0", &Lep_z0s, &b_Lep_z0s);
    schain->SetBranchAddress("jet_n", &jet_ns, &b_jet_ns);
    schain->SetBranchAddress("jet_pt", &jet_pts, &b_jet_pts);
    schain->SetBranchAddress("jet_eta", &jet_etas, &b_jet_etas);
    schain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
    schain->SetBranchAddress("scaleFactor_PILEUP", &scaleFactor_PILEUP, &b_scaleFactor_PILEUP);
    schain->SetBranchAddress("scaleFactor_ELE", &scaleFactor_ELE, &b_scaleFactor_ELE);
    schain->SetBranchAddress("scaleFactor_MUON", &scaleFactor_MUON, &b_scaleFactor_MUON);
    schain->SetBranchAddress("scaleFactor_LepTRIGGER", &scaleFactor_LepTRIGGER, &b_scaleFactor_LepTRIGGER);

    
    //output root file for TRAINING
    TFile *target = new TFile("RootFiles/DataSignal3.root","RECREATE");
    
    TTree *signaltrain = new TTree("SignalTrain","Signal for training");
    TTree *signaltest  = new TTree("SignalTest" , "Signa for testing");
    
    
    //Output/interest Variables
    Float_t fourlepsystems; //Train
    Float_t fourlepsystems0; //Test -----------0->Test
    Float_t fourlepsystempts; //Train
    Float_t fourlepsystempts0; //Test -----------0->Test
    Float_t invmassz1;
    Float_t invmassz2;
    Float_t invmassz10;
    Float_t invmassz20;
    Float_t fourrap;
    Float_t fourrap0;
    Float_t foure;
    Float_t foure0;
    Float_t wstrain;
    Float_t wstest;
    
    signaltrain->Branch("FourLepSystemMTrain", &fourlepsystems);
    signaltest->Branch("FourLepSystemMTest", &fourlepsystems0);
    signaltrain->Branch("FourLepSystemptTrain", &fourlepsystempts);
    signaltest->Branch("FourLepSystemptTest", &fourlepsystempts0);
    signaltrain->Branch("InvMassZ1Train", &invmassz1);
    signaltest->Branch("InvMassZ1Test", &invmassz10);
    signaltrain->Branch("InvMassZ2Train", &invmassz2);
    signaltest->Branch("InvMassZ2Test", &invmassz20);
    signaltrain->Branch("FourLepRapidityTrain", &fourrap);
    signaltest->Branch("FourLepRapidityTest", &fourrap0);
    signaltrain->Branch("FourLepSystemETrain", &foure);
    signaltest->Branch("FourLepSystemETest", &foure0);
    signaltrain->Branch("WeightSignalTrain", &wstrain);
    signaltest->Branch("WeightSignalTest", &wstest);
    
    int counter=0;


/*{
----------------------------------------
           SIGNAL SECTION
----------------------------------------
}*/
        
   int nentries, nbytes, k;

    nentries = (Int_t)schain->GetEntries();
    for (k = 0; k < nentries; k++)
    {
        nbytes = schain->GetEntry(k);
        if(TrigEs || TrigMs)
        {
            
        int goodlep_index[Lep_ns];
        int goodlep_n = 0;
        int lep_index =0;
            
        

            
            
        for(int i=0; i<Lep_ns; i++)
                {
        TLorentzVector leptemp;  leptemp.SetPtEtaPhiE(Lep_pts->at(i)/1000., Lep_etas->at(i), Lep_phis->at(i), Lep_Es->at(i)/1000.);
            
        if( Lep_pts->at(i) > 5000. && TMath::Abs(Lep_etas->at(i)) < 2.5 && ( (Lep_ptcone30s->at(i)/Lep_pts->at(i)) < 0.3) && ( (Lep_etcone20s->at(i) / Lep_pts->at(i)) < 0.3 ) ) 
                    {
		// electron
		if ( Lep_types->at(i) == 11 && Lep_pts->at(i) > 7000. && TMath::Abs(Lep_etas->at(i)) <2.47 ) 
                        {
		  if( TMath::Abs(lep_trackd0pvunbiaseds->at(i))/lep_tracksigd0pvunbiaseds->at(i) < 5 && TMath::Abs(Lep_z0s->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) 
                            {
		    goodlep_n = goodlep_n + 1;
		    goodlep_index[lep_index] = i;
		    lep_index++;
                            }
                        }
		//muon
		if ( Lep_types->at(i) == 13) 
                        {
		  if( TMath::Abs(lep_trackd0pvunbiaseds->at(i))/lep_tracksigd0pvunbiaseds->at(i) < 3 && TMath::Abs(Lep_z0s->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) 
                            {
		    goodlep_n = goodlep_n + 1;
		    goodlep_index[lep_index] = i;
		    lep_index++;
                            }
                        }
                    } // if of lepton classification
                }   //lepton for loop
                      
            if(goodlep_n == 4 )
	    {
	      
	      int goodlep1_index = goodlep_index[0];
	      int goodlep2_index = goodlep_index[1];
	      int goodlep3_index = goodlep_index[2];
	      int goodlep4_index = goodlep_index[3];
	      
	      //first lepton pT > 25 GeV, second > 15 GeV and third > 10 GeV		      
	      if (Lep_pts->at(goodlep1_index) > 25000. && Lep_pts->at(goodlep2_index) > 15000. && Lep_pts->at(goodlep3_index) > 10000. ) 
		{ 		 
              
        // TLorentzVector definitions
		  TLorentzVector Lepton_1  = TLorentzVector();
		  TLorentzVector Lepton_2  = TLorentzVector();
		  TLorentzVector Lepton_3  = TLorentzVector();
		  TLorentzVector Lepton_4  = TLorentzVector();
		  
		  Lepton_1.SetPtEtaPhiE(Lep_pts->at(goodlep1_index), Lep_etas->at(goodlep1_index), Lep_phis->at(goodlep1_index),Lep_Es->at(goodlep1_index));
		  Lepton_2.SetPtEtaPhiE(Lep_pts->at(goodlep2_index), Lep_etas->at(goodlep2_index), Lep_phis->at(goodlep2_index),Lep_Es->at(goodlep2_index));
		  Lepton_3.SetPtEtaPhiE(Lep_pts->at(goodlep3_index), Lep_etas->at(goodlep3_index), Lep_phis->at(goodlep3_index),Lep_Es->at(goodlep3_index));
		  Lepton_4.SetPtEtaPhiE(Lep_pts->at(goodlep4_index), Lep_etas->at(goodlep4_index), Lep_phis->at(goodlep4_index),Lep_Es->at(goodlep4_index));
		  
		  
		  // minimisation of difference from the Z mass
		  float delta_Z1=0; 
		  float delta_Z2=0; 
		  float InvMassZ1=0; 
		  float InvMassZ2=0;
		  float delta_Z1_1=0; float delta_Z1_2=0; float delta_Z1_3=0;
		  float delta_Z2_1=0; float delta_Z2_2=0; float delta_Z2_3=0;
		  float InvMassZ1_1=0; float InvMassZ1_2=0; float InvMassZ1_3=0;
		  float InvMassZ2_1=0; float InvMassZ2_2=0; float InvMassZ2_3=0;
		  float sum_ZZ1=0; float sum_ZZ2=0; float sum_ZZ3=0;
		  
		  // final values
		  float InvMassZ1_min=0; float InvMassZ2_min=0; float sum_ZZ_fin=0;
		  
		  
		  float sum_charges = Lep_charges->at(goodlep1_index) + Lep_charges->at(goodlep2_index) + Lep_charges->at(goodlep3_index) + Lep_charges->at(goodlep4_index);			 
		  
		  // step-by-step
		  // opposite charge leptons
		  if ( sum_charges == 0  ) 
		    {
         
          
		      int sum_types  = Lep_types->at(goodlep1_index) + Lep_types->at(goodlep2_index) + Lep_types->at(goodlep3_index) + Lep_types->at(goodlep4_index) ;
              
              // type e=11, mu=13
		      // begin case e+e-e+e- or mu+mu-mu+mu-
		      if ( sum_types == 44 || sum_types == 52  )
			{
			  if ( Lep_types->at(goodlep1_index) == Lep_types->at(goodlep2_index) && ( (Lep_charges->at(goodlep1_index) * Lep_charges->at(goodlep2_index)) < 0 )  )
			    {
		              InvMassZ1_1=(Lepton_1+Lepton_2).Mag()/1000.;
			      InvMassZ2_1=(Lepton_3+Lepton_4).Mag()/1000.;
			      delta_Z1_1 =  TMath::Abs(InvMassZ1_1 - 91.18); 
			      delta_Z2_1 =  TMath::Abs(InvMassZ2_1 - 91.18);
			    }
			  if ( Lep_types->at(goodlep1_index) == Lep_types->at(goodlep3_index)  && ( (Lep_charges->at(goodlep1_index) * Lep_charges->at(goodlep3_index)) < 0 ) )
			    {
			      InvMassZ1_2=(Lepton_1+Lepton_3).Mag()/1000.;
			      InvMassZ2_2=(Lepton_2+Lepton_4).Mag()/1000.;
			      delta_Z1_2 =  TMath::Abs(InvMassZ1_2 - 91.18); 
			      delta_Z2_2 =  TMath::Abs(InvMassZ2_2 - 91.18);
			    }
			  if ( Lep_types->at(goodlep1_index) == Lep_types->at(goodlep4_index)  && ( (Lep_charges->at(goodlep1_index) * Lep_charges->at(goodlep4_index)) < 0 ) )
			    {
			      InvMassZ1_3=(Lepton_1+Lepton_4).Mag()/1000.;
			      InvMassZ2_3=(Lepton_2+Lepton_3).Mag()/1000.;
			      delta_Z1_3 =  TMath::Abs(InvMassZ1_3 - 91.18); 
			      delta_Z2_3 =  TMath::Abs(InvMassZ2_3 - 91.18);
			    }

			  if(delta_Z1_1 < delta_Z2_1) { InvMassZ1_min = InvMassZ1_1; InvMassZ2_min = InvMassZ2_1;}
                          if(delta_Z2_1 < delta_Z1_1) { InvMassZ1_min = InvMassZ2_1; InvMassZ2_min = InvMassZ1_1;}

			  if(delta_Z1_2 < delta_Z2_2) { InvMassZ1_min = InvMassZ1_2; InvMassZ2_min = InvMassZ2_2;}
                          if(delta_Z2_2 < delta_Z1_2) { InvMassZ1_min = InvMassZ2_2; InvMassZ2_min = InvMassZ1_2;}

			  if(delta_Z1_3 < delta_Z2_3) { InvMassZ1_min = InvMassZ1_3; InvMassZ2_min = InvMassZ2_3;}
                          if(delta_Z2_3 < delta_Z1_3) { InvMassZ1_min = InvMassZ2_3; InvMassZ2_min = InvMassZ1_3;}

			} // cases of eeee or mumumumu
              
               if ( sum_types == 48 )
			{
			  
			  if ( Lep_types->at(goodlep1_index) == Lep_types->at(goodlep2_index)  && ( (Lep_charges->at(goodlep1_index) * Lep_charges->at(goodlep2_index)) < 0 ) )
			    {
			      InvMassZ1=(Lepton_1+Lepton_2).Mag()/1000.;
			      InvMassZ2=(Lepton_3+Lepton_4).Mag()/1000.;
			      delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
			      delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
			    }
			  if ( Lep_types->at(goodlep1_index) == Lep_types->at(goodlep3_index)  && ( (Lep_charges->at(goodlep1_index) * Lep_charges->at(goodlep3_index)) < 0 ) )
			    {
			      InvMassZ1=(Lepton_1+Lepton_3).Mag()/1000.;
			      InvMassZ2=(Lepton_2+Lepton_4).Mag()/1000.;
			      delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
			      delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
			    }
			  if ( Lep_types->at(goodlep1_index) == Lep_types->at(goodlep4_index)  && ( (Lep_charges->at(goodlep1_index) * Lep_charges->at(goodlep4_index)) < 0 ) )
			    {
			      InvMassZ1=(Lepton_1+Lepton_4).Mag()/1000.;
			      InvMassZ2=(Lepton_2+Lepton_3).Mag()/1000.;
			      delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
			      delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
			    }
			  
			  if(delta_Z1 < delta_Z2) { InvMassZ1_min = InvMassZ1; InvMassZ2_min = InvMassZ2;}
			  if(delta_Z2 < delta_Z1) { InvMassZ1_min = InvMassZ2; InvMassZ2_min = InvMassZ1;}
			} // eemumu overe
              
              
              if ( (sum_types == 44 || sum_types == 52 || sum_types == 48) )
			{
			  
			  TLorentzVector FourLepSystem = TLorentzVector();
			  FourLepSystem = Lepton_1 + Lepton_2 + Lepton_3 + Lepton_4;
			  float FourLepSystem_M = FourLepSystem.Mag()/1000.;
			  float FourLepSystem_pt = FourLepSystem.Pt()/1000.;
			  float FourLepSystem_y = FourLepSystem.Rapidity();
              float FourLepSystem_E = FourLepSystem.E()/1000.;
			 

                          //Preselection of good jets
			  int goodjet_n = 0;
			  int goodjet_index = 0;
				  /*{
			  if (jet_ns > 0)
			    {
			      for(unsigned int p=0; p<jet_ns; p++)
				{
				  if(jet_pts->at(p) > 30000. && TMath::Abs(jet_etas->at(p)) < 4.4)
				    {
				      goodjet_n++;
				      goodjet_index = p;
				    }
				}
			    }
                  }*/
                  
                          Float_t xsec_weight = (10000.*0.0003769)/(149400.);
                          Float_t scaleFactor = xsec_weight*scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER*scaleFactor_PILEUP;
      
                            Float_t m_mcWeight = mcWeight;
            
                            Float_t weight = scaleFactor*m_mcWeight;

                      
                  counter++;
                  if(counter==1) //Fill Train Data 33%
                  {
                      fourlepsystems=FourLepSystem_M;
                      fourlepsystempts=FourLepSystem_pt;
                      invmassz1=InvMassZ1_min;
                      invmassz2=InvMassZ2_min;
                      fourrap=FourLepSystem_y;
                      foure=FourLepSystem_E;
                      wstrain=weight;
                      signaltrain->Fill();
                  }
                  if(counter==2)
                  {
                      fourlepsystems0=FourLepSystem_M;
                      fourlepsystempts0=FourLepSystem_pt;
                      invmassz10=InvMassZ1_min;
                      invmassz20=InvMassZ2_min;
                      fourrap0=FourLepSystem_y;
                      foure0=FourLepSystem_E;
                      wstest=weight;
                      signaltest->Fill();
                  }
                   if(counter==3)
                  {
                      counter=0;
                      fourlepsystems0=FourLepSystem_M;
                      fourlepsystempts0=FourLepSystem_pt;
                      invmassz10=InvMassZ1_min;
                      invmassz20=InvMassZ2_min;
                      fourrap0=FourLepSystem_y;
                      foure0=FourLepSystem_E;
                      wstest=weight;
                      signaltest->Fill();
                  }
              }
          }
        }
        }
        }
    }
   
    signaltrain->Write();
    signaltest->Write();
    target->Close();
    

    
    return 0;  
}
