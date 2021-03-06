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

int HZZSamplesToRoot(void){
    
    //Get Data from repository    
    TString path = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/";
    
    
    //Creates a Chain for store SIGNAL info
    //TChain* schain = new TChain("mini");
    //Creates a chain for store BACKGROUND info
    TChain* bchain = new TChain("mini");
    
    bchain->AddFile(path+"MC/mc_361106.Zee.4lep.root");
    bchain->AddFile(path+"MC/mc_361107.Zmumu.4lep.root");
    bchain->AddFile(path+"MC/mc_410000.ttbar_lep.4lep.root");
    //bchain->AddFile(path+"MC/mc_363490.llll.4lep.root");
    //bchain->AddFile(path+"MC/mc_361108.Ztautau.4lep.root"); //not at HZZ JN
    /*{
    schain->AddFile(path+"MC/mc_345060.ggH125_ZZ4lep.4lep.root");
    schain->AddFile(path+"MC/mc_341947.ZH125_ZZ4lep.4lep.root");
    schain->AddFile(path+"MC/mc_341964.WH125_ZZ4lep.4lep.root");
    schain->AddFile(path+"MC/mc_344235.VBFH125_ZZ4lep.4lep.root");
}*/
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
    
    //Define variables to use in BACKGROUND
    vector<float> *Lep_ptb;
    Bool_t TrigEb;
    Bool_t TrigMb;
    UInt_t Lep_nb;
    vector<float> *Lep_phib;
    vector<float> *Lep_etab;
    vector<float> *Lep_Eb;
    vector<float> *Lep_ptcone30b;
    vector<float> *Lep_etcone20b;
    vector<unsigned int> *Lep_typeb;
    vector<int> *Lep_chargeb;
    vector<float>   *lep_trackd0pvunbiasedb;
    vector<float>   *lep_tracksigd0pvunbiasedb;
    vector<float>   *Lep_z0b;
    UInt_t          jet_nb;
    vector<float>   *jet_ptb;
    vector<float>   *jet_etab;
    
 //------------------------------------------------   
  /*{ TBranch        *b_TrigEs;   //!
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
   }*/
   TBranch        *b_TrigEb;   //!
   TBranch        *b_TrigMb;   //!   //!
   TBranch        *b_Lep_nb;   //!
   TBranch        *b_Lep_ptb;   //!
   TBranch        *b_Lep_etab;   //!
   TBranch        *b_Lep_phib;   //!
   TBranch        *b_Lep_Eb;   //!
   TBranch        *b_Lep_z0b;   //!
   TBranch        *b_Lep_chargeb;   //!
   TBranch        *b_Lep_typeb; 
   TBranch        *b_Lep_ptcone30b;   //!
   TBranch        *b_Lep_etcone20b;   //!
   TBranch        *b_lep_trackd0pvunbiasedb;   //!
   TBranch        *b_lep_tracksigd0pvunbiasedb;   //!
   TBranch        *b_jet_nb;   //!
   TBranch        *b_jet_ptb;   //!
   TBranch        *b_jet_etab;
//---------------------------------------------------
    
    //Setting branch address for signal and background
 /*{   schain->SetBranchAddress("lep_pt", &Lep_pts, &b_Lep_pts);
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
    }*/
    bchain->SetBranchAddress("lep_pt", &Lep_ptb, &b_Lep_ptb);
    bchain->SetBranchAddress("trigE", &TrigEb, &b_TrigEb);
    bchain->SetBranchAddress("trigM", &TrigMb, &b_TrigMb);
    bchain->SetBranchAddress("lep_n", &Lep_nb, &b_Lep_nb);
    bchain->SetBranchAddress("lep_phi", &Lep_phib, &b_Lep_phib);
    bchain->SetBranchAddress("lep_eta", &Lep_etab, &b_Lep_etab);
    bchain->SetBranchAddress("lep_E", &Lep_Eb, &b_Lep_Eb);
    bchain->SetBranchAddress("lep_ptcone30", &Lep_ptcone30b, &b_Lep_ptcone30b);
    bchain->SetBranchAddress("lep_etcone20", &Lep_etcone20b, &b_Lep_etcone20b);
    bchain->SetBranchAddress("lep_type", &Lep_typeb, &b_Lep_typeb);
    bchain->SetBranchAddress("lep_charge", &Lep_chargeb, &b_Lep_chargeb);
    bchain->SetBranchAddress("lep_trackd0pvunbiased", &lep_trackd0pvunbiasedb, &b_lep_trackd0pvunbiasedb);
    bchain->SetBranchAddress("lep_tracksigd0pvunbiased", &lep_tracksigd0pvunbiasedb, &b_lep_tracksigd0pvunbiasedb);
    bchain->SetBranchAddress("lep_z0", &Lep_z0b, &b_Lep_z0b);
    bchain->SetBranchAddress("jet_n", &jet_nb, &b_jet_nb);
    bchain->SetBranchAddress("jet_pt", &jet_ptb, &b_jet_ptb);
    bchain->SetBranchAddress("jet_eta", &jet_etab, &b_jet_etab);
    


    
    //output root file
    TFile *target = new TFile("Data11.root","RECREATE");
    
    TTree *signal = new TTree("Signal","Signal from samples");
    TTree *background = new TTree("Background","Background from samples");
    
    
    Float_t invmassz1mins;
    Float_t invmassz2mins;
    Float_t invmassz1minb;
    Float_t invmassz2minb;
    /*{
    vector<float> *Lepton_pts;
    vector<float> *Lepton_ptb;
    }*/
    //signal->Branch("met"   , &mets   , "met/F"   );
    //signal->Branch("InvMassZ1_min", &invmassz1mins);
    //signal->Branch("InvMassZ2_min", &invmassz2mins);
    background->Branch("InvMassZ1_min", &invmassz1minb);
    background->Branch("InvMassZ2_min", &invmassz2minb);


/*{
----------------------------------------
           SIGNAL SECTION
----------------------------------------
}*/
        
 /*{   int nentries, nbytes, k;

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
			 

                          //Preselection of good jets
			  int goodjet_n = 0;
			  int goodjet_index = 0;
				  
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
                  
                  //here goes the stuff
                  invmassz1mins=InvMassZ1_min;
                  invmassz2mins=InvMassZ2_min;
                  signal->Fill();
                  
              }
          }
        }
        }
        
        }
        //mets=missing_et;
        //Lepton_pts=Lep_pts;
        //signal->Fill();
    }
    }*/

/*{
----------------------------------------
        BACKGROUND SECTION
----------------------------------------
}*/

    int nentriesb, nbytesb, b;

    nentriesb = (Int_t)bchain->GetEntries();
    for (b = 0; b < nentriesb; b++)
    {
        nbytesb = bchain->GetEntry(b);
        if(TrigEb || TrigMb)
        {
            
        int goodlep_indexb[Lep_nb];
        int goodlep_nb = 0;
        int lep_indexb =0;
            
        for(int i=0; i<Lep_nb; i++)
                {
        TLorentzVector leptempb;  leptempb.SetPtEtaPhiE(Lep_ptb->at(i)/1000., Lep_etab->at(i), Lep_phib->at(i), Lep_Eb->at(i)/1000.);
            
        if( Lep_ptb->at(i) > 5000. && TMath::Abs(Lep_etab->at(i)) < 2.5 && ( (Lep_ptcone30b->at(i)/Lep_ptb->at(i)) < 0.3) && ( (Lep_etcone20b->at(i) / Lep_ptb->at(i)) < 0.3 ) ) 
                    {
		// electron
		if ( Lep_typeb->at(i) == 11 && Lep_ptb->at(i) > 7000. && TMath::Abs(Lep_etab->at(i)) <2.47 ) 
                        {
		  if( TMath::Abs(lep_trackd0pvunbiasedb->at(i))/lep_tracksigd0pvunbiasedb->at(i) < 5 && TMath::Abs(Lep_z0b->at(i)*TMath::Sin(leptempb.Theta())) < 0.5) 
                            {
		    goodlep_nb = goodlep_nb + 1;
		    goodlep_indexb[lep_indexb] = i;
		    lep_indexb++;
                            }
                        }
		//muon
		if ( Lep_typeb->at(i) == 13) 
                        {
		  if( TMath::Abs(lep_trackd0pvunbiasedb->at(i))/lep_tracksigd0pvunbiasedb->at(i) < 3 && TMath::Abs(Lep_z0b->at(i)*TMath::Sin(leptempb.Theta())) < 0.5) 
                            {
		    goodlep_nb = goodlep_nb + 1;
		    goodlep_indexb[lep_indexb] = i;
		    lep_indexb++;
                            }
                        }
                    } // if of lepton classification
                }   //lepton for loop
                      
            if(goodlep_nb == 4 )
	    {
	      
	      int goodlep1_indexb = goodlep_indexb[0];
	      int goodlep2_indexb = goodlep_indexb[1];
	      int goodlep3_indexb = goodlep_indexb[2];
	      int goodlep4_indexb = goodlep_indexb[3];
	      
	      //first lepton pT > 25 GeV, second > 15 GeV and third > 10 GeV		      
	      if (Lep_ptb->at(goodlep1_indexb) > 25000. && Lep_ptb->at(goodlep2_indexb) > 15000. && Lep_ptb->at(goodlep3_indexb) > 10000. ) 
		{ 		 
              
        // TLorentzVector definitions
		  TLorentzVector Lepton_1b  = TLorentzVector();
		  TLorentzVector Lepton_2b  = TLorentzVector();
		  TLorentzVector Lepton_3b  = TLorentzVector();
		  TLorentzVector Lepton_4b  = TLorentzVector();
		  
		  Lepton_1b.SetPtEtaPhiE(Lep_ptb->at(goodlep1_indexb), Lep_etab->at(goodlep1_indexb), Lep_phib->at(goodlep1_indexb),Lep_Eb->at(goodlep1_indexb));
		  Lepton_2b.SetPtEtaPhiE(Lep_ptb->at(goodlep2_indexb), Lep_etab->at(goodlep2_indexb), Lep_phib->at(goodlep2_indexb),Lep_Eb->at(goodlep2_indexb));
		  Lepton_3b.SetPtEtaPhiE(Lep_ptb->at(goodlep3_indexb), Lep_etab->at(goodlep3_indexb), Lep_phib->at(goodlep3_indexb),Lep_Eb->at(goodlep3_indexb));
		  Lepton_4b.SetPtEtaPhiE(Lep_ptb->at(goodlep4_indexb), Lep_etab->at(goodlep4_indexb), Lep_phib->at(goodlep4_indexb),Lep_Eb->at(goodlep4_indexb));
		  
		  
		  // minimisation of difference from the Z mass
		  float delta_Z1b=0; 
		  float delta_Z2b=0; 
		  float InvMassZ1b=0; 
		  float InvMassZ2b=0;
		  float delta_Z1_1b=0; float delta_Z1_2b=0; float delta_Z1_3b=0;
		  float delta_Z2_1b=0; float delta_Z2_2b=0; float delta_Z2_3b=0;
		  float InvMassZ1_1b=0; float InvMassZ1_2b=0; float InvMassZ1_3b=0;
		  float InvMassZ2_1b=0; float InvMassZ2_2b=0; float InvMassZ2_3b=0;
		  float sum_ZZ1b=0; float sum_ZZ2b=0; float sum_ZZ3b=0;
		  
		  // final values
		  float InvMassZ1_minb=0; float InvMassZ2_minb=0; float sum_ZZ_finb=0;
		  
		  
		  float sum_chargesb = Lep_chargeb->at(goodlep1_indexb) + Lep_chargeb->at(goodlep2_indexb) + Lep_chargeb->at(goodlep3_indexb) + Lep_chargeb->at(goodlep4_indexb);			 
		  
		  // step-by-step
		  // opposite charge leptons
		  if ( sum_chargesb == 0  ) 
		    {
         
          
		      int sum_typesb  = Lep_typeb->at(goodlep1_indexb) + Lep_typeb->at(goodlep2_indexb) + Lep_typeb->at(goodlep3_indexb) + Lep_typeb->at(goodlep4_indexb) ;
              
              // type e=11, mu=13
		      // begin case e+e-e+e- or mu+mu-mu+mu-
		      if ( sum_typesb == 44 || sum_typesb == 52  )
			{
			  if ( Lep_typeb->at(goodlep1_indexb) == Lep_typeb->at(goodlep2_indexb) && ( (Lep_chargeb->at(goodlep1_indexb) * Lep_chargeb->at(goodlep2_indexb)) < 0 )  )
			    {
		          InvMassZ1_1b=(Lepton_1b+Lepton_2b).Mag()/1000.;
			      InvMassZ2_1b=(Lepton_3b+Lepton_4b).Mag()/1000.;
			      delta_Z1_1b =  TMath::Abs(InvMassZ1_1b - 91.18); 
			      delta_Z2_1b =  TMath::Abs(InvMassZ2_1b - 91.18);
			    }
			  if ( Lep_typeb->at(goodlep1_indexb) == Lep_typeb->at(goodlep3_indexb)  && ( (Lep_chargeb->at(goodlep1_indexb) * Lep_chargeb->at(goodlep3_indexb)) < 0 ) )
			    {
			      InvMassZ1_2b=(Lepton_1b+Lepton_3b).Mag()/1000.;
			      InvMassZ2_2b=(Lepton_2b+Lepton_4b).Mag()/1000.;
			      delta_Z1_2b =  TMath::Abs(InvMassZ1_2b - 91.18); 
			      delta_Z2_2b =  TMath::Abs(InvMassZ2_2b - 91.18);
			    }
			  if ( Lep_typeb->at(goodlep1_indexb) == Lep_typeb->at(goodlep4_indexb)  && ( (Lep_chargeb->at(goodlep1_indexb) * Lep_chargeb->at(goodlep4_indexb)) < 0 ) )
			    {
			      InvMassZ1_3b=(Lepton_1b+Lepton_4b).Mag()/1000.;
			      InvMassZ2_3b=(Lepton_2b+Lepton_3b).Mag()/1000.;
			      delta_Z1_3b =  TMath::Abs(InvMassZ1_3b - 91.18); 
			      delta_Z2_3b =  TMath::Abs(InvMassZ2_3b - 91.18);
			    }

			  if(delta_Z1_1b < delta_Z2_1b) { InvMassZ1_minb = InvMassZ1_1b; InvMassZ2_minb = InvMassZ2_1b;}
                          if(delta_Z2_1b < delta_Z1_1b) { InvMassZ1_minb = InvMassZ2_1b; InvMassZ2_minb = InvMassZ1_1b;}

			  if(delta_Z1_2b < delta_Z2_2b) { InvMassZ1_minb = InvMassZ1_2b; InvMassZ2_minb = InvMassZ2_2b;}
                          if(delta_Z2_2b < delta_Z1_2b) { InvMassZ1_minb = InvMassZ2_2b; InvMassZ2_minb = InvMassZ1_2b;}

			  if(delta_Z1_3b < delta_Z2_3b) { InvMassZ1_minb = InvMassZ1_3b; InvMassZ2_minb = InvMassZ2_3b;}
                          if(delta_Z2_3b < delta_Z1_3b) { InvMassZ1_minb = InvMassZ2_3b; InvMassZ2_minb = InvMassZ1_3b;}

			} // cases of eeee or mumumumu
              
               if ( sum_typesb == 48 )
			{
			  
			  if ( Lep_typeb->at(goodlep1_indexb) == Lep_typeb->at(goodlep2_indexb)  && ( (Lep_chargeb->at(goodlep1_indexb) * Lep_chargeb->at(goodlep2_indexb)) < 0 ) )
			    {
			      InvMassZ1b=(Lepton_1b+Lepton_2b).Mag()/1000.;
			      InvMassZ2b=(Lepton_3b+Lepton_4b).Mag()/1000.;
			      delta_Z1b =  TMath::Abs(InvMassZ1b - 91.18); 
			      delta_Z2b =  TMath::Abs(InvMassZ2b - 91.18);
			    }
			  if ( Lep_typeb->at(goodlep1_indexb) == Lep_typeb->at(goodlep3_indexb)  && ( (Lep_chargeb->at(goodlep1_indexb) * Lep_chargeb->at(goodlep3_indexb)) < 0 ) )
			    {
			      InvMassZ1b=(Lepton_1b+Lepton_3b).Mag()/1000.;
			      InvMassZ2b=(Lepton_2b+Lepton_4b).Mag()/1000.;
			      delta_Z1b =  TMath::Abs(InvMassZ1b - 91.18); 
			      delta_Z2b =  TMath::Abs(InvMassZ2b - 91.18);
			    }
			  if ( Lep_typeb->at(goodlep1_indexb) == Lep_typeb->at(goodlep4_indexb)  && ( (Lep_chargeb->at(goodlep1_indexb) * Lep_chargeb->at(goodlep4_indexb)) < 0 ) )
			    {
			      InvMassZ1b=(Lepton_1b+Lepton_4b).Mag()/1000.;
			      InvMassZ2b=(Lepton_2b+Lepton_3b).Mag()/1000.;
			      delta_Z1b =  TMath::Abs(InvMassZ1b - 91.18); 
			      delta_Z2b =  TMath::Abs(InvMassZ2b - 91.18);
			    }
			  
			  if(delta_Z1b < delta_Z2b) { InvMassZ1_minb = InvMassZ1b; InvMassZ2_minb = InvMassZ2b;}
			  if(delta_Z2b < delta_Z1b) { InvMassZ1_minb = InvMassZ2b; InvMassZ2_minb = InvMassZ1b;}
			} // eemumu overe
              
              
              if ( (sum_typesb == 44 || sum_typesb == 52 || sum_typesb == 48) )
			{
			  
			  TLorentzVector FourLepSystemb = TLorentzVector();
			  FourLepSystemb = Lepton_1b + Lepton_2b + Lepton_3b + Lepton_4b;
			  float FourLepSystem_Mb = FourLepSystemb.Mag()/1000.;
			  float FourLepSystem_ptb = FourLepSystemb.Pt()/1000.;
			  float FourLepSystem_yb = FourLepSystemb.Rapidity();
			 
/*{
                          //Preselection of good jets
			  int goodjet_nb = 0;
			  int goodjet_indexb = 0;
				  
			  if (jet_nb > 0)
			    {
			      for(unsigned int h=0; h<jet_nb; h++)
				{
				  if(jet_ptb->at(h) > 30000. && TMath::Abs(jet_etab->at(h)) < 4.4)
				    {
				      goodjet_nb++;
				      goodjet_indexb = h;
				    }
				}
			    }
                 }*/ 
                  //here goes the stuff
                  invmassz1minb=InvMassZ1_minb;
                  invmassz2minb=InvMassZ2_minb;
                  background->Fill();
                  
              }
          }
        }
        }
        
        }
        //mets=missing_et;
        //Lepton_pts=Lep_pts;
        //signal->Fill();
    }
    

    signal->Write();
    background->Write();
    target->Close();
    return 0;  
}
