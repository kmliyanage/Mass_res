//
//  mass_resolution.c
//  
//
//  Created by Kalpanie Liyanage on 3/26/21.
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


void mass_resolution(){
    
    
    
    
    const int M = 14;
    
    double x_co[M] = {
        160.0,
        250.0,
        350.0,
        500.0,
        700.0,
        900.0,
        1150.0,
        1450.0,
        1800.0,
        2250.0,
        2800.0,
        3450.0,
        4150.0,
        5000.0,
    };
    
    double x_er[M] = {
        40.0,
        50.0,
        50.0,
        100.0,
        100.0,
        100.0,
        150.0,
        150.0,
        200.0,
        250.0,
        300.0,
        350.0,
        350.0,
        500.0,
    };
    
    double y_co_2018_BB[M] = {
        0.01136,
        0.01411,
        0.01782,
        0.02024,
        0.02396,
        0.02691,
        0.03046,
        0.03336,
        0.03528,
        0.03854,
        0.04016,
        0.04294,
        0.0453,
        0.04715,
    };
    
    
    double y_er_2018_BB[M] = {
        0.00015,
        0.00013,
        0.0004,
        0.00019,
        0.00042,
        0.00025,
        0.00036,
        0.00031,
        0.00037,
        0.00037,
        0.00046,
        0.00037,
        0.0005,
        0.00037,
    };
    
    double y_co_2017_BB[M] = {
        0.01132,
        0.01386,
        0.01766,
        0.01966,
        0.02437,
        0.02683,
        0.03022,
        0.03395,
        0.03528,
        0.03829,
        0.03951,
        0.04414,
        0.04515,
        0.04855,
    };
    
    
    double y_er_2017_BB[M] = {
        0.00014,
        0.00013,
        0.00029,
        0.00017,
        0.00034,
        0.00026,
        0.00044,
        0.00033,
        0.00039,
        0.00034,
        0.00036,
        0.00041,
        0.00051,
        0.00036,
    };
    
    double y_co_2018_BE[M] = {
        0.01692,
        0.01865,
        0.02166,
        0.02479,
        0.03031,
        0.03209,
        0.0358,
        0.03777,
        0.04006,
        0.04389,
        0.04698,
        0.05151,
        0.05209,
        0.05573,
    };
    
    
    double y_er_2018_BE[M] = {
        0.00032,
        0.00022,
        0.00045,
        0.00032,
        0.00055,
        0.00042,
        0.0007,
        0.00047,
        0.0008,
        0.00062,
        0.00069,
        0.00065,
        0.00094,
        0.00068,
    };
    
    double y_co_2017_BE[M] = {
        0.01639,
        0.01831,
        0.02107,
        0.02414,
        0.02869,
        0.03119,
        0.034,
        0.03862,
        0.03994,
        0.0443,
        0.04764,
        0.05032,
        0.05261,
        0.05479,
    };
    
    
    double y_er_2017_BE[M] = {
        0.00034,
        0.00025,
        0.00057,
        0.00032,
        0.0006,
        0.00037,
        0.001,
        0.00047,
        0.00062,
        0.0006,
        0.00058,
        0.00069,
        0.00096,
        0.00072,
    };
    
    TGraphErrors* resolution_BB_2018 = new TGraphErrors(M, x_co, y_co_2018_BB, x_er, y_er_2018_BB);
    TGraphErrors* resolution_BB_2017 = new TGraphErrors(M, x_co, y_co_2017_BB, x_er, y_er_2017_BB);
    TGraphErrors* resolution_BE_2018 = new TGraphErrors(M, x_co, y_co_2018_BE, x_er, y_er_2018_BE);
    TGraphErrors* resolution_BE_2017 = new TGraphErrors(M, x_co, y_co_2017_BE, x_er, y_er_2017_BE);
    
    
    TCanvas *c1 = new TCanvas("BE", "BE", 1000, 800);
    c1->cd();
    
    gStyle->SetOptStat(0); //turn off stst box
    c1->SetGrid();
    gStyle->SetGridColor(kGray+2); //grid color
    gStyle->SetGridWidth(1);
    c1->Update();
    
    resolution_BE_2018->SetLineColor(kBlue);
    resolution_BE_2018->SetLineWidth(1);
    resolution_BE_2018->SetMarkerStyle(20);
    resolution_BE_2018->SetMarkerColor(kBlue);
    resolution_BE_2018->SetMarkerSize(0.1);
    resolution_BE_2018->Draw("APE");
    
    
    resolution_BE_2017->SetLineColor(kRed+1);
    resolution_BE_2017->SetLineWidth(1);
    resolution_BE_2017->SetMarkerStyle(20);
    resolution_BE_2017->SetMarkerColor(kRed+1);
    resolution_BE_2017->SetMarkerSize(0.1);
    resolution_BE_2017->Draw("SAMEP");
    c1->Update();
    
    resolution_BE_2018->GetXaxis()->SetTitle("M [GeV]");
    resolution_BE_2018->GetYaxis()->SetTitle("Mass resolution");
    //resolution_BE_2018->GetYaxis()->SetLabelFont(43);
    //resolution_BE_2018->GetYaxis()->SetLabelSize(15);
    resolution_BE_2018->SetMinimum(0.0);
    resolution_BE_2018->SetMaximum(0.08);
    resolution_BE_2018->SetTitle("");
    c1->Update();
    
    
    TPaveText* tText1 = new TPaveText(0.55,0.8,0.75,0.85, "brNDC");
    tText1->SetBorderSize(0);
    tText1->SetFillColor(0);
    tText1->SetFillStyle(0);
    TText *t1 = tText1->AddText("TunePNew BE");
    tText1->SetTextFont(42);
    tText1->SetTextSize(0.03);
    tText1->Draw();
    
    TLegend *l1 = new TLegend(0.6,0.7,0.75,0.8);
    l1->SetBorderSize(0);
    l1->SetLineWidth(0);
    l1->SetLineStyle(0);
    l1->SetLineColor(0);
    l1->SetTextFont(42);
    l1->SetTextSize(0.03);
   // l1->SetMargin(0.12);
    l1->AddEntry(resolution_BE_2018, "2018", "lp2");
    l1->AddEntry(resolution_BE_2017, "2017", "lp2");
    l1->Draw();
    
    TPaveText* tText2 = new TPaveText(0.75, 0.92, 0.93, 0.93, "brNDC");
    tText2->SetBorderSize(0);
    tText2->SetFillColor(0);
    tText2->SetFillStyle(0);
    TText *t2 = tText2->AddText("CMS (13TeV)");
    tText2->SetTextFont(42);
    tText2->SetTextSize(0.03);
    tText2->Draw();
    
    c1->Update();
    
}
