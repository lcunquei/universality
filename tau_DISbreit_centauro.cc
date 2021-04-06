// main01.cc is a part of the PYTHIA event generator.
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
#include "fastjet/JetDefinition.hh"

#include "fastjet/Selector.hh"
#include "fastjet/contrib/Centauro.hh"
#include "fastjet/EECambridgePlugin.hh"

#include <ctime>
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <valarray>


using namespace Pythia8;

   
   
    
//___________________________________________________________________
//-------------------------------------------------------------------
// Lorentz transformations (Use For Breit Frame)
//-------------------------------------------------------------------

fastjet::PseudoJet operator-(const fastjet::PseudoJet &p)
{
  return fastjet::PseudoJet(-p.px(), -p.py(), -p.pz(), p.E());
}

// boost the vector p with the boost given by (bx,by,bz)
fastjet::PseudoJet boost(const fastjet::PseudoJet p, 
		const double &bx, const double &by, const double &bz){
  double b2 = bx*bx + by*by + bz*bz;
  //assert(b2 < 1.0);
  if(b2 >= 1.0)
    {
      b2 = 0.999;
      //cout << "Bad in boost" << endl;
      //cout << bx << " " << by << " " << bz << endl;
      //cout << " " << endl;
    }
  //cout << "Matt's b2 = " << b2 << endl;
  double gamma = 1.0/sqrt(1.0 - b2);
  double bp = bx*p.px() + by*p.py() + bz*p.pz();
  double gamma2 = (b2 > 0.0 ? (gamma - 1.0)/b2 : 0.0);

  return fastjet::PseudoJet(p.px() + gamma2*bp*bx + gamma*bx*p.E(),
		   p.py() + gamma2*bp*by + gamma*by*p.E(),
		   p.pz() + gamma2*bp*bz + gamma*bz*p.E(),
		   gamma*(p.E() + bp));
}

// boost the vector p with the boost given by b (i.e. (bx/bE,by/bE,bz/bE))
fastjet::PseudoJet boost(const fastjet::PseudoJet p, const fastjet::PseudoJet b){
  return boost(p, b.px()/b.E(), b.py()/b.E(), b.pz()/b.E());
}

// rotation around the x axis
fastjet::PseudoJet rotateX(const fastjet::PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return fastjet::PseudoJet(p.px(),
		   cp*p.py()-sp*p.pz(),
		   sp*p.py()+cp*p.pz(),
		   p.E());
}

// rotation around the y axis
fastjet::PseudoJet rotateY(const fastjet::PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return fastjet::PseudoJet(cp*p.px()+sp*p.pz(),
		   p.py(),
		   cp*p.pz()-sp*p.px(),
		   p.E());
}

// rotation around the z axis
fastjet::PseudoJet rotateZ(const fastjet::PseudoJet p, const double& psi){
  double sp = sin(psi);
  double cp = cos(psi);
  return fastjet::PseudoJet(cp*p.px()-sp*p.py(),
		   sp*p.px()+cp*p.py(),
		   p.pz(),
		   p.E());
}


int main(int argc, char* argv[]) {


   Int_t cislo = -1;                 //unique number for each file
   Int_t tune  = 6;                 //pythia tune to be used for DIS
   
   Int_t had = 1;               //hadronisation on/off
   Int_t charged=0;                  //full or track-based jets
   Int_t unev=1;                     //underlying event on/off(ISR+MPI) 


   if(argc!=8){  
   cout<<"Usage:"<<endl<<"./pygen <PythiaTune> <Number> <nEvts> <jetR> <charged> <unev> <had>"<<endl;
      return 0;
   }
   tune  = atoi(argv[1]);
   cislo = atoi(argv[2]);
  
   had  = atof(argv[7]);
   charged = atoi(argv[5]);
   unev    = atoi(argv[6]);
 
   
   Int_t nEvent= atoi(argv[3]);   //(Int_t) 1e3 + 1.0;
   TString name;
   //__________________________________________________________________________
   //                        ANALYSIS SETTINGS

   double jetParameterR   = (double) atof(argv[4]); //jet R
   double trackLowPtCut   = 0.; //GeV
   double trackEtaCut     = 5;
   

   //__________________________________________________________________________
   //                        PYTHIA SETTINGS

   // Generator. 
   Pythia pythia;
   pythia.readString("Beams:frameType = 2");
   pythia.readString("Beams:idB = 11"); //beam 1 electron
   pythia.readString("Beams:idA = 2212"); //beam 2 proton
   pythia.readString("Beams:eB =18");
   pythia.readString("Beams:eA =275");
   pythia.readString(Form("Tune:pp = %d",tune));  //tune 1-13    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1

   
  // Beam energies, minimal Q2, number of events to generate.
  double eProton   = 275.;
  double eElectron = 18;
  double Q2min     = 900.;
 

   
   pythia.readString("Random:setSeed = on");
   pythia.readString(Form("Random:seed = %d",cislo));


   // Set up DIS process within some phase space.
  // Neutral current (with gamma/Z interference).
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  // Uncomment to allow charged current.
  //pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  // Phase-space cut: minimal Q2 of process.
  pythia.settings.parm("PhaseSpace:Q2Min", Q2min);

  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");


   
    pythia.readString("11:mayDecay  = off"); //pi0s
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
   //                      FASTJET  SETTINGS and centauro algorithm definition

   double etamin_Sig = - trackEtaCut + jetParameterR; //signal jet eta range
   double etamax_Sig = - etamin_Sig;



  fastjet::contrib::CentauroPlugin * centauro_plugin = new fastjet::contrib::CentauroPlugin(1.0);
  fastjet::JetDefinition jet_def(centauro_plugin);
  std::vector<fastjet::PseudoJet> fjInputs;


  
   

   //___________________________________________________ 
   //                HISTOGRAMS

   
  

    TH1D *histoJet = new TH1D("histoJet","histoJet", 50,0,50);
   histoJet->Sumw2();
   
    TH1D *histoJetEta = new TH1D("histoJetEta","histoJetEta", 200,-10,10);
   histoJetEta->Sumw2();




   TH1F* fHistXsection = new TH1F("fHistXsection", "fHistXsection", 1, 0, 1);
   fHistXsection->GetYaxis()->SetTitle("xsection");

   TH1F* fHistTrials = new TH1F("fHistTrials", "fHistTrials", 1, 0, 1);
   fHistTrials->GetYaxis()->SetTitle("trials");


    double Wmax = sqrt(4.* eProton * eElectron);
    TH1D *Qhist=new TH1D("Q [GeV]", "q",100, 0., 50.);
    TH1D *Whist=new TH1D("W [GeV]", "w",100, 0., Wmax);
    TH1D *xhist=new TH1D("x","x",100, 0., 1.);
    TH1D *yhist=new TH1D("y", "y",100, 0., 1.);
    TH1D *pTehist=new TH1D("pT of scattered electron [GeV]","pt", 100, 0., 50.);
    TH1D *pTrhist=new TH1D("particle pt in breit frame", "ptr",100, 0., 50.);

  Qhist->Sumw2();
  Whist->Sumw2();
  xhist->Sumw2();
  yhist->Sumw2();
  pTehist->Sumw2();
  pTrhist->Sumw2();

  Int_t count=0;

   // Begin event loop. 
   for(int iEvent = 0; iEvent < nEvent; iEvent++){
      if(!pythia.next()) continue;
       count=count+1;

      Double_t weight=pythia.info.sigmaGen();
       fHistXsection->Fill(weight);
       fHistTrials->Fill(count);
      fjInputs.resize(0);
      Double_t index=0; 
      Double_t fourvec[4];
      
       // Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
    Vec4 pProton = pythia.event[1].p();
    Vec4 peIn    = pythia.event[4].p();
    Vec4 peOut   = pythia.event[6].p();
    Vec4 pPhoton = peIn - peOut;

    // Q2, W2, Bjorken x, y.
    double Q2    = - pPhoton.m2Calc();
    double W2    = (pProton + pPhoton).m2Calc();
    double x     = Q2 / (2. * pProton * pPhoton);
    double y     = (pProton * pPhoton) / (pProton * peIn);


    // Set Up Boost to Breit Frame
    //const Particle* part2 = pythia.event[1]; // Incoming Proton
    //const Particle* part4 = pythia.event[3]; // Virtual Photon

      fastjet::PseudoJet proton(pythia.event[1].px(),pythia.event[1].py(),pythia.event[1].pz(),pythia.event[1].e());
      fastjet::PseudoJet gamma(pythia.event[3].px(),pythia.event[3].py(),pythia.event[3].pz(),pythia.event[3].e());
      

      fastjet::PseudoJet boost_vector = -(gamma + 2.0*x*proton);

      double bxL = boost_vector.px()/boost_vector.E(); // Check for valid boost
      double byL = boost_vector.py()/boost_vector.E();
      double bzL = boost_vector.pz()/boost_vector.E();

      if(bxL*bxL + byL*byL + bzL*bzL >= 1.0) 
	{
         continue;
	}



    
    // Fill kinematics histograms.
    Qhist->Fill( sqrt(Q2) );
    Whist->Fill( sqrt(W2) );
    xhist->Fill( x );
    yhist->Fill( y );
    pTehist->Fill( pythia.event[6].pT() );

    fastjet::PseudoJet boosted_proton = boost(proton,boost_vector);
      Double_t phi_p = boosted_proton.phi();
      Double_t theta_p = TMath::ATan2(boosted_proton.perp(),boosted_proton.pz());

      

    
    for(Int_t i = 0; i < pythia.event.size(); ++i){
      if(pythia.event[i].statusAbs()==43){
	     if(charged==1) if(!pythia.event[i].isCharged()) continue;  
             fourvec[0]=pythia.event[i].px();
             fourvec[1]=pythia.event[i].py();
             fourvec[2]=pythia.event[i].pz();
             fourvec[3]=pythia.event[i].e(); 
             fastjet::PseudoJet particle(fourvec);
             fastjet::PseudoJet boosted_particle = boost(particle,boost_vector);
	     boosted_particle = rotateZ(boosted_particle, -phi_p);
	     boosted_particle = rotateY(boosted_particle, -theta_p);
	     
	     fjInputs.push_back(boosted_particle);
             pTrhist->Fill( boosted_particle.perp());
            

      }



    }
      

      vector<fastjet::PseudoJet> inclusiveJets;
     
      fastjet::ClusterSequence clustSeq_Sig(fjInputs, jet_def);
          
      inclusiveJets = sorted_by_E(clustSeq_Sig.inclusive_jets(0));
      if(inclusiveJets.size()==0) continue;
      histoJet->Fill(inclusiveJets[0].perp());
      histoJetEta->Fill(inclusiveJets[0].eta());
     
          
          
	 	  
      	          }// End of event loop.

   //____________________________________________________
   //          SAVE OUTPUT

   TString tag = Form("DISs%02d",TMath::Nint(jetParameterR*10));
   TFile* outFile = new TFile(Form("%s_tune%d_c%d_charged%d_unev%d_had%d.root",tag.Data(),tune,cislo,charged,unev,had), "RECREATE");
   
    outFile->cd();
    
    

     pTrhist->Write();
     Qhist->Write();
     Whist->Write();
     xhist->Write();
     fHistXsection->Write();
     fHistTrials->Write();
     yhist->Write();
     pTehist->Write();
     histoJet->Write();
     histoJetEta->Write();
     outFile->Close();
   

   pythia.stat();
   return 0;
}


