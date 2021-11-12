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
   double trackLowPtCut   = 0.; //GeV
   double trackEtaCut     = 2.;
   double zetacut =2.4;
   double zpt=54;
   Float_t ptHatMin=50;
   Float_t ptHatMax=1000;
   
   TRandom3 *r3=new TRandom3();
   //__________________________________________________________________________
   //                        PYTHIA SETTINGS

   // Generator. Process selection. LHC initialization. Histogram.
   Pythia pythia;
   pythia.readString("Beams:idA = 2212"); //beam 1 proton
   pythia.readString("Beams:idB = 2212"); //beam 2 proton
   pythia.readString("Beams:eCM = 2760.");
   pythia.readString(Form("Tune:pp = %d",tune));  //tune 1-13    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1

   pythia.readString("Random:setSeed = on");
   pythia.readString(Form("Random:seed = %d",cislo));
  //pythia.readString("HardQCD:all = on");
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
            pythia.readString("PhaseSpace:bias2Selection = on");
         pythia.readString("PhaseSpace:bias2SelectionPow = 5.");
         pythia.readString("PhaseSpace:bias2SelectionRef = 8.");
        pythia.readString(name.Data());




   }

   if(unev==0){
     pythia.readString("PartonLevel:MPI = off");
     pythia.readString("PartonLevel:ISR = off");
     }
     if(had==0){
     pythia.readString("HadronLevel:Hadronize=off");}

   pythia.readString("23:mayDecay  = off"); //Z0s  
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

  
   TH1D *histotriggerZ = new TH1D("histotriggerZ","histotriggerZ", 200,0.0,400);
   histotriggerZ->Sumw2();


   
   TProfile* fHistXsection = new TProfile("fHistXsection", "fHistXsection", 1, 0, 1);
   fHistXsection->GetYaxis()->SetTitle("xsection");

   TH1F* fHistTrials = new TH1F("fHistTrials", "fHistTrials", 1, 0, 1);
   fHistTrials->GetYaxis()->SetTitle("trials");


   const Int_t nVar = 8;
   TTree *fTreeObservables = new TTree("variables", "variables");
   TString *fShapesVarNames = new TString [nVar];
   
   float fShapesVar[8];
   
   fShapesVarNames[0] = "ptGamma";
   fShapesVarNames[1] = "ptJet";
   fShapesVarNames[2] = "mass";
   fShapesVarNames[3] = "angle";
   fShapesVarNames[4] = "asym";
   fShapesVarNames[5] = "veto";
   fShapesVarNames[6] = "flavour";
   fShapesVarNames[7] = "weight";
     for(Int_t ivar=0; ivar < nVar; ivar++){
      fTreeObservables->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/F", fShapesVarNames[ivar].Data()));}

  
  
   for(int iEvent = 0; iEvent < nEvent; iEvent++){
      if(!pythia.next()) continue;
      
      fjInputs.resize(0);
      int in=-1;
      double ptmax=0;
      double emax=0;
      Double_t index=0; 
      Double_t fourvec[4];
      int indexveto=-1;
      double ptveto=0;
       double weight=pythia.info.weight();
      for(Int_t i = 0; i < pythia.event.size(); ++i){
             if(pythia.event[i].isFinal()){
	     if(charged==1 && TMath::Abs(pythia.event[i].id())!=23) if(!pythia.event[i].isCharged()) continue;         
             if(pythia.event[i].pT() < trackLowPtCut) continue;                
	     if(TMath::Abs(pythia.event[i].eta()) > trackEtaCut) continue;     
             fourvec[0]=pythia.event[i].px();
             fourvec[1]=pythia.event[i].py();
             fourvec[2]=pythia.event[i].pz();
             fourvec[3]=pythia.event[i].e(); 
             fastjet::PseudoJet particle(fourvec);
	     fjInputs.push_back(particle);
	     if(TMath::Abs(pythia.event[i].id())==23 && particle.perp()>ptmax){
	       ptmax=particle.perp();
	       emax=particle.e();
	       in=i;
	     }
            
	     }
             

      }
      if(in==-1) continue;
     
      if(TMath::Abs(pythia.event[in].eta())>trackEtaCut) continue;
      
      if(ptmax<zpt) continue;
    
      histotriggerZ->Fill(ptmax);     
      vector<fastjet::PseudoJet> inclusiveJets;
      
      fastjet::ClusterSequence clustSeq_Sig(fjInputs, *jetDefAKT_Sig);
      inclusiveJets = sorted_by_pt(clustSeq_Sig.inclusive_jets(1.));
      double dphijz=0;
      int indexjet=-1;
      double jetptmax=0;
      //get the hardest recoil jet
      for(Int_t n=0;n<inclusiveJets.size();n++){
       
        fastjet::PseudoJet fjJet = inclusiveJets.at(n);
        dphijz=RelativePhi(pythia.event[in].phi(),fjJet.phi());
         if(TMath::Abs(fjJet.eta())>etamax_Sig) continue;
	 if(TMath::Abs(dphijz)<=2) continue;
 	 if(fjJet.perp()>jetptmax){jetptmax=fjJet.perp();
	  indexjet=n;}

	

      }
      
      if(indexjet==-1) continue;
      //get the second hardest recoil jet
      for(int m=indexjet+1;m<inclusiveJets.size();m++){
	fastjet::PseudoJet subleadJet = inclusiveJets.at(m);
	if(TMath::Abs(subleadJet.eta()>etamax_Sig)) continue;
        dphijz=RelativePhi(pythia.event[in].phi(),subleadJet.phi());
       if(TMath::Abs(dphijz)<=2) continue;
	indexveto=m;
	break;
      }
	 
      if(indexveto==-1) ptveto=-1;
      if(indexveto!=-1) ptveto=inclusiveJets.at(indexveto).perp();
            
 
      double jetmass=0;
      double asym=inclusiveJets[indexjet].perp()/ptmax;
      double jetmass=TMath::Sqrt(fjJet.e()*fjJet.e()-fjJet.perp()*fjJet.perp()-fjJet.pz()*fjJet.pz());     
         
   fShapesVar[0]=ptmax;
   fShapesVar[1]=inclusiveJets[indexjet].perp();
   fShapesVar[2]=jetmass;
   fShapesVar[3]=dphijz;
   fShapesVar[4]=asym;
   fShapesVar[5]=ptveto; 
   fShapesVar[6]=pythia.event[6].id();
   fShapesVar[7]=weight;
   fTreeObservables->Fill();
   
	  
      	           
  

   }// End of event loop.

   //____________________________________________________
   //          SAVE OUTPUT

   TString tag = Form("RecoilJetMassANTIKT%02d",TMath::Nint(jetParameterR*10));

   TFile* outFile = new TFile(Form("%s_tune%d_c%d_charged%d_unev%d_hadv%d.root",tag.Data(), tune,cislo,charged,unev,had), "RECREATE");
   
    outFile->cd();
    histotriggerZ->Write();
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

