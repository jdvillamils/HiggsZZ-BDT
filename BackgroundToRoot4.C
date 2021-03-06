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


int BackgroundToRoot4(void){
    
    TString path = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/4lep/";
    
    
    
    //Creates a Chain for store BACKGROUND info
    TChain* bchain = new TChain("mini");
    
    //Backgrounds from ZZ
    //bchain->AddFile(path+"MC/mc_363490.llll.4lep.root");
    //bchain->AddFile(path+"MC/mc_363356.ZqqZll.4lep.root");
    //bchain->AddFile(path+"MC/mc_363492.llvv.4lep.root");
    
    //minor bkgs
    bchain->AddFile(path+"MC/mc_361106.Zee.4lep.root");
    //bchain->AddFile(path+"MC/mc_361107.Zmumu.4lep.root");
    //bchain->AddFile(path+"MC/mc_410000.ttbar_lep.4lep.root"); 
    /*bchain->AddFile(path+"/mc_361108.Ztautau.4lep.root"); //not at HZZ JN
    bchain->AddFile(path+"/mc_361101.Wplusmunu.4lep.root");
    bchain->AddFile(path+"/mc_361100.Wplusenu.4lep.root");
    bchain->AddFile(path+"/mc_361102.Wplustaunu.4lep.root");
    bchain->AddFile(path+"/mc_361103.Wminusenu.4lep.root");
    bchain->AddFile(path+"/mc_361104.Wminusmunu.4lep.root");
    bchain->AddFile(path+"/mc_361105.Wminustaunu.4lep.root");
    bchain->AddFile(path+"/mc_410011.single_top_tchan.4lep.root");
    bchain->AddFile(path+"/mc_410012.single_antitop_tchan.4lep.root");
    bchain->AddFile(path+"/mc_410013.single_top_wtchan.4lep.root");
    bchain->AddFile(path+"/mc_410014.single_antitop_wtchan.4lep.root");
    bchain->AddFile(path+"/mc_410025.single_top_schan.4lep.root");
    bchain->AddFile(path+"/mc_410026.single_antitop_schan.4lep.root");
    bchain->AddFile(path+"/mc_363491.lllv.4lep.root");
    bchain->AddFile(path+"/mc_363493.lvvv.4lep.root");
    bchain->AddFile(path+"/mc_363358.WqqZll.4lep.root");
    bchain->AddFile(path+"/mc_363489.WlvZqq.4lep.root");
    */
    
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
    Float_t         mcWeight;
    Float_t         scaleFactor_PILEUP;
    Float_t         scaleFactor_ELE;
    Float_t         scaleFactor_MUON;
    Float_t         scaleFactor_LepTRIGGER;
    


    //Setting branch address for signal and background
    bchain->SetBranchAddress("lep_pt", &Lep_ptb);//, &b_Lep_ptb);
    bchain->SetBranchAddress("trigE", &TrigEb);//, &b_TrigEb);
    bchain->SetBranchAddress("trigM", &TrigMb);//, &b_TrigMb);
    bchain->SetBranchAddress("lep_n", &Lep_nb);//, &b_Lep_nb);
    bchain->SetBranchAddress("lep_phi", &Lep_phib);//, &b_Lep_phib);
    bchain->SetBranchAddress("lep_eta", &Lep_etab);//, &b_Lep_etab);
    bchain->SetBranchAddress("lep_E", &Lep_Eb);//, &b_Lep_Eb);
    bchain->SetBranchAddress("lep_ptcone30", &Lep_ptcone30b);//, &b_Lep_ptcone30b);
    bchain->SetBranchAddress("lep_etcone20", &Lep_etcone20b);//, &b_Lep_etcone20b);
    bchain->SetBranchAddress("lep_type", &Lep_typeb);//, &b_Lep_typeb);
    bchain->SetBranchAddress("lep_charge", &Lep_chargeb);//, &b_Lep_chargeb);
    bchain->SetBranchAddress("lep_trackd0pvunbiased", &lep_trackd0pvunbiasedb);//, &b_lep_trackd0pvunbiasedb);
    bchain->SetBranchAddress("lep_tracksigd0pvunbiased", &lep_tracksigd0pvunbiasedb);//, &b_lep_tracksigd0pvunbiasedb);
    bchain->SetBranchAddress("lep_z0", &Lep_z0b);//, &b_Lep_z0b);
    bchain->SetBranchAddress("jet_n", &jet_nb);//, &b_jet_nb);
    bchain->SetBranchAddress("jet_pt", &jet_ptb);//, &b_jet_ptb);
    bchain->SetBranchAddress("jet_eta", &jet_etab);//, &b_jet_etab);
    bchain->SetBranchAddress("mcWeight", &mcWeight);//, &b_mcWeight);
    bchain->SetBranchAddress("scaleFactor_PILEUP", &scaleFactor_PILEUP);//, &b_scaleFactor_PILEUP);
    bchain->SetBranchAddress("scaleFactor_ELE", &scaleFactor_ELE);//, &b_scaleFactor_ELE);
    bchain->SetBranchAddress("scaleFactor_MUON", &scaleFactor_MUON);//, &b_scaleFactor_MUON);
    bchain->SetBranchAddress("scaleFactor_LepTRIGGER", &scaleFactor_LepTRIGGER);//, &b_scaleFactor_LepTRIGGER);
    
    
    
    //output root file
    
    TFile *target = new TFile("RootFiles/DataBackground4.root","RECREATE");
    
    TTree *backtrain = new TTree("BackgroundTrain","Background for training");
    TTree *backtest  = new TTree("BackgroundTest" , "Background for testing");

    //Output/interest Variables
    Float_t fourlepsystemb; //Train
    Float_t fourlepsystemb0; //Test -----------0->Test
    Float_t fourlepsystemptb; //Train
    Float_t fourlepsystemptb0; //Test -----------0->Test
    Float_t invmassz1;
    Float_t invmassz2;
    Float_t invmassz10;
    Float_t invmassz20;
    Float_t fourrap;
    Float_t fourrap0;
    Float_t foure;
    Float_t foure0;
    Float_t wbtrain;
    Float_t wbtest;


    backtrain->Branch("FourLepSystemMTrain", &fourlepsystemb);
    backtest->Branch("FourLepSystemMTest", &fourlepsystemb0);
    backtrain->Branch("FourLepSystemptTrain", &fourlepsystemptb);
    backtest->Branch("FourLepSystemptTest", &fourlepsystemptb0);
    backtrain->Branch("InvMassZ1Train", &invmassz1);
    backtest->Branch("InvMassZ1Test", &invmassz10);
    backtrain->Branch("InvMassZ2Train", &invmassz2);
    backtest->Branch("InvMassZ2Test", &invmassz20);
    backtrain->Branch("FourLepRapidityTrain", &fourrap);
    backtest->Branch("FourLepRapidityTest", &fourrap0);
    backtrain->Branch("FourLepSystemETrain", &foure);
    backtest->Branch("FourLepSystemETest", &foure0);
    backtrain->Branch("WeightBackTrain", &wbtrain);
    backtest->Branch("WeightBackTest", &wbtest);
    
    int counter=0;



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
              float FourLepSystem_Eb = FourLepSystemb.E()/1000.;
			 

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
                   Float_t xsec_weight = (10000.*1950.5295)/(150277594200.);
                   Float_t scaleFactor = xsec_weight*scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER*scaleFactor_PILEUP;
      
        Float_t m_mcWeight = mcWeight;
            
        Float_t weight = scaleFactor*m_mcWeight;
                  
                  //here goes the stuff
                  counter++;
                  if(counter==1) //Fill Train Data 33%
                  {
                      fourlepsystemb=FourLepSystem_Mb;
                      fourlepsystemptb=FourLepSystem_ptb;
                      invmassz1=InvMassZ1_minb;
                      invmassz2=InvMassZ2_minb;
                      fourrap=FourLepSystem_yb;
                      foure=FourLepSystem_Eb;
                      wbtrain=weight;
                      backtrain->Fill();
                  }
                  if(counter==2)
                  {
                      fourlepsystemb0=FourLepSystem_Mb;
                      fourlepsystemptb0=FourLepSystem_ptb;
                      invmassz10=InvMassZ1_minb;
                      invmassz20=InvMassZ2_minb;
                      fourrap0=FourLepSystem_yb;
                      foure0=FourLepSystem_Eb;
                      wbtest=weight;
                      backtest->Fill();
                      
                  }
                  if(counter==3)
                  { 
                      fourlepsystemb0=FourLepSystem_Mb;
                      fourlepsystemptb0=FourLepSystem_ptb;
                      invmassz10=InvMassZ1_minb;
                      invmassz20=InvMassZ2_minb;
                      fourrap0=FourLepSystem_yb;
                      foure0=FourLepSystem_Eb;
                      wbtest=weight;
                      backtest->Fill();
                      counter=0;
                  }

                  
              }
          }
        }
        }
        
        }
    }
    
    backtrain->Write();
    backtest->Write();
    target->Close();
    return 0;  
}