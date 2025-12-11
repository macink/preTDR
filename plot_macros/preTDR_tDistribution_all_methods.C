#include "RiceStyle.h"
#include "ePIC_style.C"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>

using namespace std;

void preTDR_tDistribution_all_methods(TString phiFile)
{
    TFile* coherent_phi_file = TFile::Open(phiFile, "READ");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

	//t distribution
	TH1D* h_t_MC = (TH1D*) coherent_phi_file->Get("h_t_MC");
	TH1D* h_t_REC_L = (TH1D*) coherent_phi_file->Get("h_t_REC_EEMC_cut");
	TH1D* h_t_REC_proj12 = (TH1D*) coherent_phi_file->Get("h_t_REC_wRES_cut_pi12");
    TH1F* h_phi_sartre_events = (TH1F*) coherent_phi_file->Get("h_Nevents");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,800);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", "d#sigma/d|t| [nb/(GeV/c)^{2}] ", 100,0,0.18,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-2, 1e6);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.2);
	base1->GetYaxis()->SetTitleOffset(1.5);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	double dN_dt = 0.002;
    double sigma_coherent = 459.05; // nb
    double nEvents_coherent = h_phi_sartre_events->GetEntries(); // 6.36679 M
    double preTDR_lumi = 50761.4; // nb^-1
    int rebin_width = 2;   

    h_t_MC->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    double simu_coherent_lumi = nEvents_coherent/sigma_coherent; // nb^-1 
    double ratio_coherent_phi_preTDR = preTDR_lumi/simu_coherent_lumi;
    for (int i = 1; i <= h_t_MC->GetNbinsX(); ++i) 
    {
        double binError = h_t_MC->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_MC->SetBinError(i, newError);
    }
	h_t_MC->SetLineStyle(1);   
	h_t_MC->SetLineWidth(1);   
	h_t_MC->SetLineColor(kBlack);
	h_t_MC->Draw("HISTsame");
	h_t_MC->Draw("same");

    h_t_REC_L->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt));
    for (int i = 1; i <= h_t_REC_L->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_L->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_L->SetBinError(i, newError);
    }
	h_t_REC_L->SetMarkerStyle(20); // method L RECO
	h_t_REC_L->SetMarkerColor(kP8Blue);
	h_t_REC_L->Draw("PEsame");

    h_t_REC_proj12->GetYaxis()->SetTitle("d#sigma/d|t| [nb/(GeV/c)^{2}]");
    h_t_REC_proj12->Scale(sigma_coherent*(1/nEvents_coherent)*(1/dN_dt)*((M_PI/2)/(M_PI/12)));
    for (int i = 1; i <= h_t_REC_proj12->GetNbinsX(); ++i) 
    {
        double binError = h_t_REC_proj12->GetBinError(i);
        double newError = binError / sqrt(ratio_coherent_phi_preTDR);
        h_t_REC_proj12->SetBinError(i, newError);
    }
	h_t_REC_proj12->SetMarkerStyle(30);
	h_t_REC_proj12->SetMarkerColor(kP8Pink);
	h_t_REC_proj12->SetLineColor(kP8Pink);
	h_t_REC_proj12->Draw("PEsame");

	TLatex* r44 = new TLatex(0.2, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.85");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");
	
	TLatex* r44_0 = new TLatex(0.2, 0.79, "  |y_{"+vm_label+"}|<3.5, |M_{inv} #minus M_{"+vm_label+"}| < 0.02 GeV");
	r44_0->SetNDC();
	r44_0->SetTextSize(20);
	r44_0->SetTextFont(43);
	r44_0->SetTextColor(kBlack);
	r44_0->Draw("same");
	
    // Add labels
    TLatex* ep = new TLatex(0.18, 0.33, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.23, 0.33, " Simulation 25.10.2");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.18, 0.23, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.18, 0.28, "eAu #rightarrow e'Au'#phi, 10x100 GeV");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLatex* r142111 = new TLatex(0.18, 0.18, "L_{int} = 10 fb^{-1}/A"); // preTDR lumi
	r142111->SetNDC();
	r142111->SetTextSize(25);
	r142111->SetTextFont(43);
	r142111->SetTextColor(kBlack);
	r142111->Draw("same");
	
	TLegend *w7 = new TLegend(0.58,0.7,0.73,0.85);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(20);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
	w7->AddEntry(h_t_REC_L, "Sartre "+vm_label+" Method L RECO", "P");
    w7->AddEntry(h_t_REC_proj12, "Sartre "+vm_label+" RECO #theta_{max} = #pi/12", "P");
	w7->Draw("same");

	c1->Print("./plot_t_dist_preTDR_allMethods.pdf");
}
