//
//  mass_res_2018.c
//  
//
//  Created by Kalpanie Liyanage on 3/12/21.
//

#include <stdio.h>
#include "TFile.h"
#include "TFrame.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TList.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TText.h"
#include <TString.h>
#include "TMath.h"

#include "TFitResult.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"


Double_t DSCB(Double_t *x, Double_t *par){
    
    // par[0] = normalization
    // par[1] = mean
    // par[2] = sigma
    // par[3] = alphaL
    // par[4] = alphaR
    // par[5] = nL
    // par[6] = nR
    
    
    double A = (x[0] - par[1]) / par[2];
    
    if(A < -par[3]){
        //      double B = par[3] / par[5];
        //      double C = B*(1/B - par[3] - A);
        //      return par[0]*std::exp(-0.5 * par[3]*par[3]) * pow(C, -par[5]);
        double A1 = pow(par[5]/par[3],par[5])*std::exp(-0.5*par[3]*par[3]);
        double B1 = (par[5]/par[3])-par[3]-A;
        return par[0] * A1 / pow(B1, par[5]);
    }
    
    else if(A >= -par[3] && A <= par[4])
        return par[0]*std::exp(-0.5 * A*A);
    
    else{
        //      double B = par[4] / par[6];
        //      double C = B*(1/B - par[4] + A);
        //      return par[0]*std::exp(-0.5 * par[4]*par[4]) * pow(C, -par[6]);
        double A1 = pow(par[6]/par[4],par[6])*exp(-0.5*par[4]*par[4]);
        double B1 = par[6]/par[4]-par[4]+A;
        return par[0] * A1 / pow(B1, par[6]);
    }
}


void mass_res_2018(){
    
    gStyle->SetOptStat(0); //turn off stst box
    gStyle->SetOptFit(111); // fit parameters: chi2,parameter value,parameter error
    
    gStyle->SetGridColor(kGray+2); //grid color
    gStyle->SetGridWidth(1);
    
    
    TFile *file1 = new TFile("out.root", "RECREATE");
    TList *l = new TList();
    
    
    TFile *mc = TFile::Open("DY_effi_2017.root");
    
    mc->cd("Our2018OppSignResolutionVertex");
    
    const int M = 15;
    double mrange[M] = {120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500};
    
    //gROOT->LoadMacro("DSCB.c");
    //if("DSCB.c") std::cout<<"DSCB succefully open"<<std::endl;
    
    const int eta = 3;
    bool in=false;
    bool bb=false;
    bool be=true;
    
    
    TH1D* myhist[M-1]; //To create my histograms(ProjectionY)
    TH2F *h1;
    
    
    TString etas[eta] =  {
        "",
        "_BB",
        "_BE",
    };
    
        
    TString region;
    if(in) region = etas[0];
    if(bb) region = etas[1];
    if(be) region = etas[2];
        
    h1 = (TH2F*)gDirectory->Get("DileptonMassResVMass_2d"+region);
    int nbins = h1->GetNbinsX();
    TAxis *axis = h1->GetXaxis();
    std::cout<<nbins<<std::endl;
    
    
    //creating a histogram for each mass bin
    for(int i=0; i<M-1; i++){
        
        TCanvas *c = new TCanvas("DATA_MC", "DATA_MC",  1000, 1000);
        c->cd();
        
        int bmin = axis->FindBin(mrange[i]);
        int bmax = axis->FindBin(mrange[i+1]);
        std::cout<<bmin<<"  "<<bmax<<std::endl;
        TString hist = "residual"+region+"_"+mrange[i]+"_"+mrange[i+1];
        myhist[i] = h1->ProjectionY(hist, bmin, bmax);
        
        double fit_min = myhist[i]->GetMean()-2.0*myhist[i]->GetRMS();
        double fit_max = myhist[i]->GetMean()+1.7*myhist[i]->GetRMS();
        myhist[i]->Rebin(5);
        myhist[i]->GetXaxis()->SetRangeUser(fit_min,fit_max);
        myhist[i]->SetTitle("");
        myhist[i]->GetXaxis()->SetTitle("m_{RECO} / m_{GEN} - 1");
        myhist[i]->SetLineColor(kBlack);
        myhist[i]->SetMarkerStyle(20);
        myhist[i]->SetMarkerSize(0.8);
        
        myhist[i]->Draw("PE");
        c->Update();
        
        
        
        //myhist[i]->GetXaxis()->SetLimits(fit_min, fit_max);
        
        //double bin_width = 0.01; //default value
        //double Nbins = (fit_max-fit_min)/bin_width;
        //double rebin = 1000.0/Nbins;ßß
        //myhist[i]->Rebin(rebin);
        //std::cout<<hist<<"  "<<myhist[i]->GetNbinsX()<<"    "<<fit_min<<"   "<<fit_max<<"   "<<Nbins<<std::endl;
       
        //////////////////////////Fitting///////////////////////////
        
        TF1 *gauss = new TF1("gauss","gaus",fit_min,fit_max);
        gauss->SetParameters(0,myhist[i]->GetMean(),myhist[i]->GetRMS());
        myhist[i]->Fit("gauss","M0R+");
        gauss->SetLineColor(kRed);
        gauss->SetLineWidth(2);
        //gauss->Draw("SAME");
        c->Update();
        
        TF1 *funct = new TF1("funct",DSCB,fit_min,fit_max,7);
        //funct->SetParameters(myhist[i]->Integral(fit_min,fit_max),myhist[i]->GetMean(),myhist[i]->GetRMS(), 1.5, 2.0, 2.0, 20.0);
        funct->SetParameters(gauss->GetParameter(0),myhist[i]->GetMean(),myhist[i]->GetRMS(), 1.3, 1.5, 1.5, 20.0);
        //funct->FixParameter(3,1.58);
        //funct->FixParameter(4,1.8);
        //funct->FixParameter(5,2.2);
        //funct->FixParameter(6,20.0);
        funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR","nL","nR");
        myhist[i]->Fit("funct","M0R");
        funct->SetLineColor(kBlue);
        funct->SetLineWidth(2);
        funct->Draw("SAME");
        c->Update();
        
        TPaveText* tText1 = new TPaveText(0.2, 0.80, 0.4, 0.9, "brNDC");
        tText1->SetBorderSize(0);
        tText1->SetFillColor(0);
        tText1->SetFillStyle(0);
        TString name = "2017"+region+" "+mrange[i]+" < m < "+mrange[i+1];
        TText *t1 = tText1->AddText(name);
        //tText1->SetTextSize(0.005);
        tText1->Draw();
        c->Update();
        
        
        
        //myhist[i]->SaveAs("ProjectionY_2018_"+hist+".root","root");
        c->SaveAs("ProjectionY_2017_"+hist+".root","root");
        
        
        
        
    }
   
    file1->Write();
    file1->Close();
    mc->Close();
    
}
    
    
    
    

