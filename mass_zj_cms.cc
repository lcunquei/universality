// main01.cc is a phart of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TList.h"
#include "TVector3.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"  
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh" 
#include "fastjet/contrib/Njettiness.hh"
#include <ctime>
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <valarray>


using namespace Pythia8;
//___________________________________________________________________
//               FUNCTIONS

bool EtaCut(fastjet::PseudoJet fjJet, double etaMin, double etaMax);
Double_t RelativePhi(Double_t mphi,Double_t vphi);


//___________________________________________________________________


int main(int argc, char* argv[]) {


   Int_t cislo = -1;                 //unique number for each file
   Int_t tune  = -1;                 //pythia tune
   Int_t  had=1;                     //hadronization on/off 
   Int_t charged=0;                  //full or track-based jets
   Int_t unev=1;                     //underlying event (ISR+MPI) 


   if(argc!=8){  
     std::cout<<"Usage:"<<std::endl<<"./pygen <PythiaTune> <Number> <nEvts> <jetR> <charged> <unev> <had>"<<endl;
      return 0;
   }
   tune  = atoi(argv[1]);
   cislo = atoi(argv[2]);
   
   had= atoi(argv[7]);  
  
   charged = atoi(argv[5]);
   unev    = atoi(argv[6]);
 
   Int_t nEvent= atoi(argv[3]);   //(Int_t) 1e3 + 1.0;
   TString name;
   //__________________________________________________________________________
   //                        ANALYSIS SETTINGS

   double jetParameterR   = (double) atof(argv[4]); //jet R
   double trackLowPtCut   = 0.; 
   double trackEtaCut     = 3;
   double deltajet_muon = 0.4;
   double etacut_muon =2.4;
   double ptcut_muon=10;    //2 muons with pT>10 GeV
   double zmass_cut_min=70;
   double zmass_cut_max=110;
   double zptcut=40;
   //analyis cuts from https://arxiv.org/pdf/1702.01060.pdf  

   Float_t ptHatMin=30;
   Float_t ptHatMax=1000;
   
   TRandom3 *r3=new TRandom3();
   //__________________________________________________________________________
   //                        PYTHIA SETTINGS

   // Generator. Process selection. LHC initialization. Histogram.
   Pythia pythia;
   pythia.readString("Beams:idA = 2212"); //beam 1 proton
   pythia.readString("Beams:idB = 2212"); //beam 2 proton
   pythia.readString("Beams:eCM = 13000");
   pythia.readString(Form("Tune:pp = %d",tune));  //tune 1-13    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d",cislo));
  
  pythia.readString("WeakBosonAndParton:qg2gmZq = on");
  pythia.readString("WeakBosonAndParton:qqbar2gmZg  = on");
  pythia.readString("WeakZ0:gmZmode = 2");

  if(ptHatMin<0 || ptHatMax <0){     
      pythia.readString("PhaseSpace:pTHatMin = 0."); // <<<<<<<<<<<<<<<<<<<<<<<
   }else{
      name = Form("PhaseSpace:pTHatMin = %f", (Float_t) ptHatMin);
      pythia.readString(name.Data()); 
      name = Form("PhaseSpace:pTHatMax = %f", (Float_t) ptHatMax);
      pythia.readString(name.Data()); 
   }

   if(unev==0){
     pythia.readString("PartonLevel:MPI = off");
     pythia.readString("PartonLevel:ISR = off");
   }
   if(had==0){
     pythia.readString("HadronLevel:Hadronize=off");}

 
   pythia.readString("111:mayDecay  = off"); //pi0s  
   pythia.readString("310:mayDecay  = off"); //K0s
   pythia.readString("3122:mayDecay = off"); //labda0
   pythia.readString("3112:mayDecay = off"); //sigma-
   pythia.readString("3212:mayDecay = off"); //sigma0
   pythia.readString("3222:mayDecay = off"); //sigma+
   pythia.readString("3312:mayDecay = off"); //xi-
   pythia.readString("3322:mayDecay = off"); //xi+
   pythia.readString("3334:mayDecay = off"); //omega-

  
   pythia.init();
 
   //___________________________________________________ 
   //                      FASTJET  SETTINGS

   double etamin_Sig = - trackEtaCut + jetParameterR; //signal jet eta range
   double etamax_Sig = - etamin_Sig;


   fastjet::Strategy strategy = fastjet::Best;
   fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
   fastjet::JetDefinition *jetDefAKT_Sig = NULL;

   jetDefAKT_Sig = new fastjet::JetDefinition(fastjet::antikt_algorithm, 
                                              jetParameterR, 
                                              recombScheme, 
                                              strategy);
 
   // Fastjet input
    std::vector<fastjet::PseudoJet> fjInputs;
    std::vector<fastjet::PseudoJet> fjPartMax;
    std::vector<fastjet::PseudoJet> fjMuonCandidates;
    std::vector<fastjet::PseudoJet> fjMuonCandidates_ordered;
  
   TH1D *histotriggergamma = new TH1D("histotriggergamma","histotriggergamma", 200,0.0,400);
   histotriggergamma->Sumw2();


   
   TProfile* fHistXsection = new TProfile("fHistXsection", "fHistXsection", 1, 0, 1);
   fHistXsection->GetYaxis()->SetTitle("xsection");

   TH1F* fHistTrials = new TH1F("fHistTrials", "fHistTrials", 1, 0, 1);
   fHistTrials->GetYaxis()->SetTitle("trials");


   const Int_t nVar = 5;
   TTree *fTreeObservables = new TTree("variables", "variables");
   TString *fShapesVarNames = new TString [nVar];
   
   float fShapesVar[5];
   
   fShapesVarNames[0] = "ptZ";
   fShapesVarNames[1] = "ptJet";
   fShapesVarNames[2] = "mass";
   fShapesVarNames[3] = "angle";
   fShapesVarNames[4] = "Zmass";
  
   
     for(Int_t ivar=0; ivar < nVar; ivar++){
      fTreeObservables->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));}

  
  
   for(int iEvent = 0; iEvent < nEvent; iEvent++){
      if(!pythia.next()) continue;
      
      fjInputs.resize(0);
      fjMuonCandidates.resize(0);
      fjMuonCandidates_ordered.resize(0);
      int in=-1;
      double ptmax=0;
      Double_t index=0; 
      Double_t fourvec[4];
      Double_t muon[4];
   
     for(Int_t i = 0; i < pythia.event.size(); ++i){
             if(pythia.event[i].isFinal()){
             if(pythia.event[i].pT() < trackLowPtCut) continue;                
             if(TMath::Abs(pythia.event[i].eta()) > trackEtaCut) continue;  
             if(TMath::Abs(pythia.event[i].id())==13){
             muon[0]=pythia.event[i].px();
             muon[1]=pythia.event[i].py();
             muon[2]=pythia.event[i].pz();
             muon[3]=pythia.event[i].e();
	     
             
             fastjet::PseudoJet partmuon(muon);
             if(partmuon.perp()<ptcut_muon) continue;
             if(TMath::Abs(partmuon.eta())<etacut_muon) continue;
	      partmuon.set_user_index(i);
             fjMuonCandidates.push_back(partmuon);}}}

	     

      fjMuonCandidates_ordered = sorted_by_pt(fjMuonCandidates);
      if(fjMuonCandidates_ordered.size()<2) continue;
      int index1=fjMuonCandidates_ordered[0].user_index();
      int index2=fjMuonCandidates_ordered[1].user_index();
       if(-1*pythia.event[index1].id()!=pythia.event[index2].id()) continue; 
      fastjet::PseudoJet zet;
      zet.reset(0,0,0,0);
      zet=fjMuonCandidates_ordered[0]+fjMuonCandidates_ordered[1];
      if(zet.m()>zmass_cut_max || zet.m()<zmass_cut_min)  continue;
      if(zet.perp()<zptcut) continue;  
 
     for(Int_t i = 0; i < pythia.event.size(); ++i){
             if(pythia.event[i].isFinal()){
             if(charged==1) if(!pythia.event[i].isCharged()) continue;         
             if(pythia.event[i].pT() < trackLowPtCut) continue;                
	     if(TMath::Abs(pythia.event[i].eta()) > trackEtaCut) continue;     
             if(i==fjMuonCandidates_ordered[0].user_index() || i==fjMuonCandidates_ordered[1].user_index()) continue;
             fourvec[0]=pythia.event[i].px();
             fourvec[1]=pythia.event[i].py();
             fourvec[2]=pythia.event[i].pz();
             fourvec[3]=pythia.event[i].e(); 
             fastjet::PseudoJet particle(fourvec);
	     fjInputs.push_back(particle);
            
	     }}
             

      
     
      
    
      histotriggergamma->Fill(zet.perp());     
    
     
      vector<fastjet::PseudoJet> inclusiveJets_Sig;
      fastjet::ClusterSequence clustSeq_Sig(fjInputs, *jetDefAKT_Sig);
      inclusiveJets_Sig = clustSeq_Sig.inclusive_jets(1.); 


          for(unsigned int ijet = 0; ijet < inclusiveJets_Sig.size(); ijet++){
          fastjet::PseudoJet fjJet = inclusiveJets_Sig.at(ijet);
	  Double_t dphi=RelativePhi(zet.phi(),fjJet.phi()); 
          if(fjJet.delta_R(fjMuonCandidates[0])<=deltajet_muon) continue;
          if(fjJet.delta_R(fjMuonCandidates[1])<=deltajet_muon) continue; 
	  if(TMath::Abs(fjJet.eta())>etamax_Sig) continue;
          double jetmass=TMath::Sqrt(fjJet.e()*fjJet.e()-fjJet.perp()*fjJet.perp()-fjJet.pz()*fjJet.pz());
    cout<<"jetmass"<< jetmass<<endl;	  
   
   fShapesVar[0]=zet.perp();
   fShapesVar[1]=fjJet.perp();
   fShapesVar[2]=jetmass;
   fShapesVar[3]=dphi;
   fShapesVar[4]=zet.m();
   fTreeObservables->Fill();
   
	  }
      	           
  

   }// End of event loop.

   //____________________________________________________
   //          SAVE OUTPUT

   TString tag = Form("RecoilJetMassANTIKT%02d",TMath::Nint(jetParameterR*10));

   TFile* outFile = new TFile(Form("%s_tune%d_c%d_charged%d_unev%d_hadv%d.root",tag.Data(), tune,cislo,charged,unev,had), "RECREATE");
   
    outFile->cd();
    histotriggergamma->Write();
    fTreeObservables->Write(); 
    outFile->Close();


   pythia.stat();
   return 0;
}

//_________________________________________________________________________ 
//_________________________________________________________________________ 
//_________________________________________________________________________ 
//_________________________________________________________________________ 
bool EtaCut(fastjet::PseudoJet fjJet, double etaMin, double etaMax) {
   if(fjJet.eta() > etaMax || fjJet.eta() < etaMin){
      return false;
   }else{
      return true;
   }
}
//_________________________________________________________________________ 
Double_t RelativePhi(Double_t mphi,Double_t vphi){
   //Get relative azimuthal angle of two particles -pi to pi
   if      (vphi < -TMath::Pi()) vphi += TMath::TwoPi();
   else if (vphi > TMath::Pi())  vphi -= TMath::TwoPi();

   if      (mphi < -TMath::Pi()) mphi += TMath::TwoPi();
   else if (mphi > TMath::Pi())  mphi -= TMath::TwoPi();

   Double_t dphi = mphi - vphi;
   if      (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
   else if (dphi > TMath::Pi())  dphi -= TMath::TwoPi();

   return dphi;//dphi in [-Pi, Pi]
}

