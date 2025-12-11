#include "RiceStyle.h"
#include "ePIC_style.C"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>

using namespace std;

void preTDR_t_resolution(TString phiFile)
{
	TFile* coherent_phi_file = TFile::Open(phiFile, "READ");

	TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";
	TH2D* h_t_res = (TH2D*) coherent_phi_file->Get("h_t_res_proj_percent_pi12");

	TCanvas* c2 = new TCanvas("c2","c2",1,1,1000,800);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.13);
	gPad->SetRightMargin(0.1);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.13);

	TH1D* base1 = makeHist("base1", "", "|t| [GeV/c]^{2}", " #delta t/|t| (resolution) ", 100,0,0.2,kBlack);
	base1->GetYaxis()->SetRangeUser(1e-3, 1000);
	base1->GetXaxis()->SetTitleColor(kBlack);
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base1,1.2,1.6);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.8);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.7);
	base1->GetXaxis()->SetNdivisions(5,5,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();
	TH2D* h_res = (TH2D*) h_t_res;

	TH1D* h_res_1D = new TH1D("h_res_1D","",100,0,0.2);
	for(int ibin=0;ibin<h_res->GetNbinsX();ibin++)
	{
		TH1D* tmp=h_res->ProjectionY("tmp",ibin+1,ibin+1);
		double sigma = tmp->GetStdDev();
		double sigmaerror = tmp->GetStdDevError();
		h_res_1D->SetBinContent(ibin+1, sigma);
		h_res_1D->SetBinError(ibin+1, sigmaerror);
	}

	h_res_1D->SetMarkerSize(1.6);
	h_res_1D->SetMarkerColor(kBlack);
	h_res_1D->SetLineColor(kBlack);
	h_res_1D->SetMarkerStyle(20);

	h_res_1D->Fit("pol0","RMS0","",0.011,0.019);//first  dip
	h_res_1D->Fit("pol0","RMS0","",0.045,0.053);//second dip
	h_res_1D->Fit("pol0","RMS0","",0.095,0.102);//third  dip
    h_res_1D->Fit("pol0","RMS0","",0.167,0.175);//fourth  dip
	h_res_1D->Draw("Psame");

	TLatex* r42 = new TLatex(0.15, 0.91, "eAu 10x100 GeV");
	r42->SetNDC();
	r42->SetTextSize(25);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.8,0.91, "ePIC");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44_2 = new TLatex(0.17, 0.18, vm_label+" #rightarrow "+daug_label );
	r44_2->SetNDC();
	r44_2->SetTextSize(30);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	TPad* drawPad = new TPad("pad_etalab_11","pad_etalab_11",0.16,0.53,0.47,0.83);
	drawPad->SetLeftMargin(0.08);
	drawPad->SetRightMargin(0.08);
	drawPad->SetTopMargin(0.0);
	drawPad->SetBottomMargin(0.08);
	drawPad->Draw("same");
	drawPad->SetTicks();
	drawPad->SetLogz(1);
	drawPad->cd();
	TH1D* base2 = makeHist("base2", "", "MC", " resolution ", 100,0,0.2,kBlack);
	base2->GetYaxis()->SetRangeUser(-10, 1.5);
	base2->GetXaxis()->SetTitleColor(kBlack);
	base2->GetXaxis()->SetLabelColor(kBlack);
	base2->GetYaxis()->SetLabelColor(kBlack);
	base2->GetXaxis()->SetTitle("Resolution");
	base2->GetYaxis()->SetTitle("MC");
	TGaxis::SetMaxDigits(3);
	fixedFontHist1D(base2,3,3);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1);
	base2->GetXaxis()->SetNdivisions(5,5,0);
	base2->GetYaxis()->SetNdivisions(5,5,0);
	base2->Draw();

	h_res->Draw("colzsame");
	gPad->Update(); 
    TPaletteAxis *palette = (TPaletteAxis*)h_res->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9); 
    palette->SetX2NDC(0.95); 
    palette->SetY1NDC(0.08); 
    palette->SetY2NDC(0.99); 
    gPad->Modified();
    gPad->Update();

	c2->Print("./plot_t_resolution.pdf");
}
