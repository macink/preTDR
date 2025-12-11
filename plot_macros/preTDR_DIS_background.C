#include "RiceStyle.h"
#include "ePIC_style.C"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>

using namespace std;

void preTDR_DIS_background(TString phiFile, TString disFile)
{
    TFile* phi_t_file = TFile::Open(phiFile, "READ");
    TFile* DIS_noVetoes = TFile::Open(disFile, "READ");

	TH1D* h_t_MC = (TH1D*) phi_t_file->Get("h_t_MC");
	TH1D* h_t_REC_wRES_cut_pi12 = (TH1D*) phi_t_file->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_phi_sartre_events = (TH1D*) phi_t_file->Get("h_Nevents");

    TH1D* h_DIS_noVetoes_proj = (TH1D*) DIS_noVetoes->Get("h_t_REC_wRES_cut_pi12");
    TH1D* h_DIS_events = (TH1D*) DIS_noVetoes->Get("h_Nevents");

    double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double sigma_DIS = 50184.8; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries();
    double nEvents_DIS = h_DIS_events->GetEntries();
    double preTDR_lumi = 50761.4; //nb^-1
    int rebin_width = 2;   

    TCanvas* c14_213 = new TCanvas("c14_213","c14_213",1,1,1200,800);
	gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    h_t_MC->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
	h_t_MC->GetXaxis()->SetTitle("|t| [GeV/c]^{2}");
    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
    h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");

    h_t_REC_wRES_cut_pi12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_wRES_cut_pi12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_wRES_cut_pi12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_wRES_cut_pi12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_wRES_cut_pi12->SetBinError(i, newError);
    }
	h_t_REC_wRES_cut_pi12->SetMarkerStyle(30);
	h_t_REC_wRES_cut_pi12->SetMarkerColor(kP8Pink);
    h_t_REC_wRES_cut_pi12->SetLineColor(kP8Pink);
	h_t_REC_wRES_cut_pi12->Draw("PEsame");

    h_DIS_noVetoes_proj->Rebin(rebin_width);
    h_DIS_noVetoes_proj->Scale(sigma_DIS*(1/nEvents_DIS)*(1/dN_dt)/rebin_width*((M_PI/2)/(M_PI/12)));
    double simu_DIS_lumi = nEvents_DIS/sigma_DIS; // nb^-1 
    double ratio_DIS_preTDR = preTDR_lumi/simu_DIS_lumi;
    for (int i = 1; i <= h_DIS_noVetoes_proj->GetNbinsX(); ++i) 
    {
        double binError = h_DIS_noVetoes_proj->GetBinError(i);
        double newError = binError / sqrt(ratio_DIS_preTDR);
        h_DIS_noVetoes_proj->SetBinError(i, newError);
    }
    h_DIS_noVetoes_proj->SetMarkerStyle(2);
    h_DIS_noVetoes_proj->SetMarkerColor(kP8Cyan);
    h_DIS_noVetoes_proj->SetLineColor(kP8Cyan);
    h_DIS_noVetoes_proj->Draw("PEsame");
        
	TLegend *w14_213 = new TLegend(0.6,0.68,0.72,0.9);
	w14_213->AddEntry(h_t_MC,"#phi: MC", "L");
    w14_213->AddEntry(h_t_REC_wRES_cut_pi12,"#phi: #theta_{max}= #pi/12", "P");    
    w14_213->AddEntry(h_DIS_noVetoes_proj,"DIS", "P");
	w14_213->SetBorderSize(0);   
    w14_213->SetFillStyle(0);
    w14_213->SetTextSize(30);
	w14_213->SetTextFont(45);
    w14_213->Draw("same");	

    c14_213->cd(1);  
    gPad->SetTopMargin(0.08); 
    TLatex* title14_213 = new TLatex();
    title14_213->SetNDC(); 
    title14_213->SetTextSize(0.05);
    title14_213->SetTextAlign(22);  
    //title14_213->DrawLatex(0.5, 0.97, "|t| Distribution with DIS No Vetoes");  

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.13, 0.85, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(.038);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.18, 0.85, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(30);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.17, 0.75, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(30);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.15, 0.8, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(30);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same"); 

    TLatex* r421111 = new TLatex(0.15, 0.15, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r421111->SetNDC();
	r421111->SetTextSize(25);
	r421111->SetTextFont(43);
	r421111->SetTextColor(kBlack);
	r421111->Draw("same");

    c14_213->Print("./plot_t_DIS_noVetoes_preTDR.pdf");
}
