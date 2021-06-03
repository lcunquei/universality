#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TString.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TTree.h>
#include <TMath.h>
#endif
void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text, Float_t textSize = 0.06)
{
  TLatex* latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextSize(textSize);
  latex->SetTextColor(color);
  latex->SetTextFont(42);
  latex->Draw();
}



void plot(){
  Int_t option=0; 
  TCanvas *can[10];
   TCanvas *can2[6];
   TCanvas *can3[6];
    TCanvas *can4[6];
    TCanvas *can5[6];
  TPad *pad[11];
  TPad *pad2[11];
  TPad *pad3[11];
  TPad *pad4[11];
  TPad *pad5[11];
  TLegend *lego;
  TPad *pad0;
   TLegend *legon;
  TH1D *recoilmass[11];
  TH1D *recoilmpt[11];
  TH1D *trig[11];
    TH1D *tau1[11];
  TH1D *tau2[11];
   TH1D *tau3[11];
    TH1D *tau4[11];


   TH1D *taub1_noh[6];
  TH1D *taub2_noh[6];
   TH1D *taub3_noh[6];
    TH1D *taub4_noh[6];
    
  for(Int_t j=0;j<11;j++){

         recoilmass[j]=new TH1D(Form("recoilmass[%d]",j),Form("recoilmass[%d]",j),100,0,100);
       recoilmass[j]->Sumw2();
recoilmpt[j]=new TH1D(Form("recoilmpt[%d]",j),Form("recoilmpt[%d]",j),100,-6,0);
 recoilmpt[j]->Sumw2();

tau1[j]=new TH1D(Form("tau1[%d]",j),Form("tau1[%d]",j),100,-10,0);
 tau1[j]->Sumw2();

 tau2[j]=new TH1D(Form("tau2[%d]",j),Form("tau2[%d]",j),100,-10,0);
 tau2[j]->Sumw2();

 tau3[j]=new TH1D(Form("tau3[%d]",j),Form("tau3[%d]",j),100,-10,0);
 tau3[j]->Sumw2();

 tau4[j]=new TH1D(Form("tau4[%d]",j),Form("tau4[%d]",j),100,-10,0);
 tau4[j]->Sumw2();





 
  }



  gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
   
    TString fname[20];
    TFile *input;
   TFile *input2;
    TFile *input3;
    TTree *mc;
     TTree *mc2;
     Float_t ptGamma, ptJet,mass,angle,asym;
     
    
    fname[0]="MassR02_had0.root";
    fname[1]="MassR04_had0.root";
    fname[2]="MassR06_had0.root";
    fname[3]="MassR08_had0.root";
    fname[4]="MassR1_had0.root";
    fname[5]="MassR02_had1.root";
    fname[6]="MassR04_had1.root";
    fname[7]="MassR06_had1.root";
    fname[8]="MassR08_had1.root";
    fname[9]="MassR1_had1.root";




    
    
     fname[10]="TauR02_had0.root";
      fname[11]="TauR04_had0.root";
       fname[12]="TauR06_had0.root";
        fname[13]="TauR08_had0.root";
	fname[14]="TauR1_had0.root";

     fname[15]="TauR02_had1.root";
      fname[16]="TauR04_had1.root";
       fname[17]="TauR06_had1.root";
       	fname[18]="TauR08_had1.root";
         fname[19]="TauR1_had1.root";
     
    Int_t nEv=0;



		  for(Int_t i=0;i<10;i++){

                    input=TFile::Open(fname[i]);
		    mc=(TTree*)input->Get("variables");
		    mc->SetBranchAddress("ptGamma", &ptGamma);
      mc->SetBranchAddress("ptJet", &ptJet);
      mc->SetBranchAddress("mass", &mass);
      mc->SetBranchAddress("angle", &angle);
      mc->SetBranchAddress("asym", &asym);
      trig[i]=(TH1D*)input->Get("histotriggergamma");
		    
		    
        nEv=mc->GetEntries();
        for(int iEntry=0; iEntry< nEv; iEntry++){
        mc->GetEntry(iEntry);
	if(ptJet<30) continue;
        if(TMath::Abs(angle)<=2) continue;	
      	if(ptGamma>300 || ptGamma<100) continue;
        recoilmass[i]->Fill(mass); 
        recoilmpt[i]->Fill(log(mass/ptJet));

	}
        int norm = trig[i]->Integral(51,150);
        recoilmass[i]->Scale(1./norm);
        recoilmpt[i]->Scale(1./norm);
        recoilmass[i]->SetLineColor(i+1);
        recoilmpt[i]->SetLineColor(i+1);
		  }
       
		
		  Float_t tau,zcut,Q2,ptjet,size;

       for(Int_t i=10;i<20;i++){

                    input=TFile::Open(fname[i]);
		    mc2=(TTree*)input->Get("variables");

      mc2->SetBranchAddress("tau", &tau);
      mc2->SetBranchAddress("ptjet", &ptjet);
      mc2->SetBranchAddress("zcut", &zcut);
      mc2->SetBranchAddress("Q2", &Q2);
       mc2->SetBranchAddress("size", &size);
      
		    
		    
        nEv=mc2->GetEntries();
        for(int iEntry=0; iEntry< nEv; iEntry++){
        mc2->GetEntry(iEntry);
	if(size==1) continue;
	
		

        if(zcut<0.2) tau1[i]->Fill(log(tau));
	if(zcut<0.4) tau2[i]->Fill(log(tau)); 
        if(zcut<0.6) tau3[i]->Fill(log(tau));
	if(zcut<0.8) tau4[i]->Fill(log(tau));

         

	
	}
	tau1[i]->Scale(1./tau1[i]->Integral(1,-1));
        tau2[i]->Scale(1./tau2[i]->Integral(1,-1));
	tau3[i]->Scale(1./tau3[i]->Integral(1,-1));
	tau4[i]->Scale(1./tau4[i]->Integral(1,-1));

 tau1[i]->SetLineWidth(3);
 tau2[i]->SetLineWidth(3);
 tau3[i]->SetLineWidth(3);
 tau4[i]->SetLineWidth(3);  
       }
       





		  

 	  	    Int_t i=3;

     can[i]= new TCanvas(Form("canvas%d",i),Form("canvas%d",i) ,1100,1100);
     can[i]->SetTicks();
     can[i]->cd();
     pad[i] = new TPad("pad0","This is pad0",0.,0.,1,1.);
     pad[i]->SetFillColor(0);
     pad[i]->SetMargin(0.15,0.9,0.25,0.9);
     pad[i]->Draw();
     pad[i]->SetTicks(1,1);
     pad[i]->cd();
     recoilmpt[0]->GetYaxis()->SetTitleOffset(1.15);
     recoilmpt[0]->GetXaxis()->SetTitleOffset(0.9);
   
   recoilmpt[0]->GetXaxis()->SetLabelFont(42);
   recoilmpt[0]->GetYaxis()->SetLabelFont(42);
   recoilmpt[0]->GetXaxis()->SetLabelSize(0.04);
   recoilmpt[0]->GetYaxis()->SetLabelSize(0.04);
    
  recoilmpt[0]->GetXaxis()->SetTitleFont(42);
  recoilmpt[0]->GetYaxis()->SetTitleFont(42);
  
  recoilmpt[0]->GetXaxis()->SetTitleSize(0.065);
  recoilmpt[0]->GetYaxis()->SetTitleSize(0.055);
  recoilmpt[0]->SetMarkerSize(1);
  recoilmpt[0]->SetMarkerSize(1);

   recoilmpt[0]->SetMarkerSize(1);
  recoilmpt[0]->SetMarkerSize(1);
  recoilmpt[0]->GetYaxis()->SetRangeUser(0,0.04);

              
    recoilmpt[0]->Draw("");
    recoilmpt[1]->Draw("same");
    recoilmpt[2]->Draw("same");

    recoilmpt[3]->Draw("same");
    recoilmpt[4]->Draw("same");
    recoilmpt[5]->Draw("same");

    recoilmpt[6]->Draw("same");
    recoilmpt[7]->Draw("same");
    recoilmpt[8]->Draw("same");

    recoilmpt[9]->Draw("same");
    recoilmpt[10]->Draw("same");
    recoilmpt[11]->Draw("same");


    
    
    DrawLatex(0.2, 0.8, 1, "Z-jet,pp 5 TeV R=0.2", 0.03);
     DrawLatex(0.2, 0.75, 1, "Pythia8 Tune 4C", 0.03);
 	lego = new TLegend(0.7, 0.6, 0.8, 0.85);
          lego->SetBorderSize(0);
          lego->SetTextSize(0.02);
          lego->SetTextFont(42);



         
   lego->AddEntry(recoilmpt[0],"R=0.2", "PEL");
   lego->AddEntry(recoilmpt[1],"R=0.4", "PEL");
   lego->AddEntry(recoilmpt[2],"R=0.6", "PEL");
   lego->AddEntry(recoilmpt[1],"R=0.8", "PEL");
   lego->AddEntry(recoilmpt[2],"R=1", "PEL");

   lego->Draw();
   lego->SetFillColor(0);

    
    recoilmpt[0]->GetXaxis()->SetTitle("log(m_{jet}/p_{T,jet})");
    recoilmpt[0]->GetYaxis()->SetTitle("1/N_{Z} dN/d(log(m_{jet}/p_{T,jet}))");
    can[i]->SaveAs("massRecoil_R.pdf");







    
}
           


     	       
       
    
