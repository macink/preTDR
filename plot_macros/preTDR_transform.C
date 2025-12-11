#include "RiceStyle.h"
#include "ePIC_style.C"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <vector>

using namespace std;

void preTDR_transform(TString phiFile)
{
    TFile* coherent_phi_file = TFile::Open(phiFile, "READ");

	TH1D* hdsigmadt_MC = (TH1D*)coherent_phi_file->Get("h_t_MC");
	TH1D* hdsigmadt_REC = (TH1D*)coherent_phi_file->Get("h_t_REC_EEMC_cut");
	TH1D* hdsigmadt_REC_new = (TH1D*)coherent_phi_file->Get("h_t_REC_wRES_cut_pi12");

	const int nTrials = 1000; 
	const double hbarc = 0.197;
	const double t_cut = 0.2;
	const double bmin = -12;
	const double bmax = 12;
	const int noOfBins = 300;

	int nbins = hdsigmadt_MC->GetNbinsX();
	TRandom3 randGen(0); // Random seed

	TH1D* hF_b_MC_2d = new TH1D("hF_b_MC_2d", "", noOfBins, bmin, bmax);
	TH1D* hF_b_REC_2d = new TH1D("hF_b_REC_2d", "", noOfBins, bmin, bmax);
	TH1D* hF_b_REC_new_2d = new TH1D("hF_b_REC_new_2d", "", noOfBins, bmin, bmax);

	vector<vector<double>> trials_MC(noOfBins, vector<double>(nTrials));
	vector<vector<double>> trials_REC(noOfBins, vector<double>(nTrials));
	vector<vector<double>> trials_REC_new(noOfBins, vector<double>(nTrials));

	for (int trial = 0; trial < nTrials; ++trial) 
	{
    	for (int j = 1; j <= noOfBins; ++j) 
    	{
        	double b_2d = hF_b_MC_2d->GetBinCenter(j);
        	double prefactor = 1.0 / (2*TMath::Pi());

        	double F_b_MC_2d = 0, F_b_REC_2d = 0, F_b_REC_new_2d = 0;

        	for (int i = 1; i <= nbins; ++i) 
        	{
            	double tBinWidth = hdsigmadt_MC->GetBinWidth(i);
            	double t = hdsigmadt_MC->GetBinCenter(i);
            	double delta = sqrt(fabs(t));

            	// Gaussian sampling
            	double dsigmadt_MC = randGen.Gaus(hdsigmadt_MC->GetBinContent(i), sqrt(60)*hdsigmadt_MC->GetBinError(i)) / 1e7;
            	double dsigmadt_REC = randGen.Gaus(hdsigmadt_REC->GetBinContent(i), sqrt(60)*hdsigmadt_REC->GetBinError(i)) / 1e7;
            	double dsigmadt_REC_new = randGen.Gaus(hdsigmadt_REC_new->GetBinContent(i), sqrt(60)*hdsigmadt_REC_new->GetBinError(i)) / 1e7;

            	double bessel = TMath::BesselJ0(b_2d * delta / hbarc);

            	if (t > t_cut) continue;

            	double amp_MC = dsigmadt_MC > 0 ? sqrt(dsigmadt_MC) : 0;
            	double amp_REC = dsigmadt_REC > 0 ? sqrt(dsigmadt_REC) : 0;
            	double amp_REC_new = dsigmadt_REC_new > 0 ? sqrt(dsigmadt_REC_new) : 0;

            	if (t > 0.014) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
            	if (t > 0.048) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
            	if (t > 0.098) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }
                if (t > 0.17) { amp_MC *= -1; amp_REC *= -1; amp_REC_new *= -1; }

            	F_b_MC_2d += amp_MC * bessel * tBinWidth / 2;
            	F_b_REC_2d += amp_REC * bessel * tBinWidth / 2;
            	F_b_REC_new_2d += amp_REC_new * bessel * tBinWidth / 2;
        	}

        	F_b_MC_2d *= prefactor / hbarc;
        	F_b_REC_2d *= prefactor / hbarc;
        	F_b_REC_new_2d *= prefactor / hbarc;

        	trials_MC[j - 1][trial] = F_b_MC_2d;
        	trials_REC[j - 1][trial] = F_b_REC_2d;
        	trials_REC_new[j - 1][trial] = F_b_REC_new_2d;
    	}
	}

	for (int j = 1; j <= noOfBins; ++j) 
	{
    	auto& vec_MC = trials_MC[j - 1];
    	auto& vec_REC = trials_REC[j - 1];
    	auto& vec_REC_new = trials_REC_new[j - 1];

    	double mean_MC = TMath::Mean(nTrials, vec_MC.data());
    	double std_MC = TMath::RMS(nTrials, vec_MC.data());
    	double mean_REC = TMath::Mean(nTrials, vec_REC.data());
    	double std_REC = TMath::RMS(nTrials, vec_REC.data());
    	double mean_REC_new = TMath::Mean(nTrials, vec_REC_new.data());
    	double std_REC_new = TMath::RMS(nTrials, vec_REC_new.data());

    	hF_b_MC_2d->SetBinContent(j, mean_MC);
    	hF_b_MC_2d->SetBinError(j, std_MC);

    	hF_b_REC_2d->SetBinContent(j, mean_REC);
    	hF_b_REC_2d->SetBinError(j, std_REC);

    	hF_b_REC_new_2d->SetBinContent(j, mean_REC_new);
    	hF_b_REC_new_2d->SetBinError(j, std_REC_new);

    	double mc = hF_b_MC_2d->GetBinContent(j);
    	double rec = hF_b_REC_2d->GetBinContent(j);
    	double rec_new = hF_b_REC_new_2d->GetBinContent(j);

    	double emc = hF_b_MC_2d->GetBinError(j);
    	double erec = hF_b_REC_2d->GetBinError(j);
    	double erec_new = hF_b_REC_new_2d->GetBinError(j);

    	double ratio_rec_mc = (mc != 0) ? rec / mc : 0;
    	double ratio_rec_new_mc = (mc != 0) ? rec_new / mc : 0;
    	double error_ratio = 0;

		if (mc != 0) 
		{
    		error_ratio = ratio_rec_new_mc * sqrt(pow(erec_new / rec_new, 2) + pow(emc / mc, 2));
		}
	}

    TCanvas* c1 = new TCanvas("c1","c1",800,800);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.12);  
    gPad->SetLogx(0);
	gStyle->SetOptStat(0);
        
    hF_b_MC_2d->Scale(1.0 / hF_b_MC_2d->Integral("width"));
    hF_b_MC_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_MC_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_MC_2d->GetXaxis()->SetTitleOffset(1.2);
    hF_b_MC_2d->GetYaxis()->SetTitleOffset(1.5);
    hF_b_MC_2d->GetYaxis()->SetRangeUser(-0.07, 0.16);
    hF_b_MC_2d->SetLineColor(kBlack);
    hF_b_MC_2d->SetLineWidth(4);
    hF_b_MC_2d->GetXaxis()->SetLabelSize(0.04);  
    hF_b_MC_2d->GetYaxis()->SetLabelSize(0.04);

	hF_b_REC_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_REC_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_REC_2d->Scale(1.0 / hF_b_REC_2d->Integral("width"));
    hF_b_REC_2d->GetYaxis()->SetRangeUser(-0.07, 0.16);
    hF_b_REC_2d->SetMarkerStyle(20); 
    hF_b_REC_2d->SetMarkerColor(kP8Blue);
    hF_b_REC_2d->SetLineColor(kP8Blue);

	hF_b_REC_new_2d->GetXaxis()->SetTitle("b [fm]");
    hF_b_REC_new_2d->GetYaxis()->SetTitle("F(b)/#scale[0.6]{#int} F(b) db");
    hF_b_REC_new_2d->Scale(1.0 / hF_b_REC_new_2d->Integral("width"));
    hF_b_REC_new_2d->SetMarkerStyle(30); 
    hF_b_REC_new_2d->SetMarkerColor(kP8Pink);
    hF_b_REC_new_2d->SetLineColor(kP8Pink);

    hF_b_REC_2d->Draw(); 
    hF_b_REC_new_2d->Draw("PEsame"); 
    hF_b_MC_2d->SetMarkerStyle(0);
    hF_b_MC_2d->Draw("Lsame");

    TString vm_label="#phi";
	TString daug_label="K^{+}K^{-}";

    // Add labels
    TLatex* ep = new TLatex(0.22, 0.3, "ePIC");
    ep->SetNDC();
    ep->SetTextFont(62);
    ep->SetTextSize(0.030);
    ep->SetTextColor(kBlack);
    ep->Draw("same");

    TLatex* r421 = new TLatex(0.28, 0.3, " Simulation 25.10.2, 10x100 GeV");
	r421->SetNDC();
	r421->SetTextSize(25);
	r421->SetTextFont(43);
	r421->SetTextColor(kBlack);
	r421->Draw("same");

	TLatex* label_5 = new TLatex(0.22, 0.2, Form("%s #rightarrow %s", vm_label.Data(), daug_label.Data()));
	label_5->SetNDC();
	label_5->SetTextSize(25);
	label_5->SetTextFont(43);
	label_5->SetTextColor(kBlack);
	label_5->Draw("same");

    TLatex* r4211 = new TLatex(0.22, 0.25, "eAu #rightarrow e'Au'#phi");
	r4211->SetNDC();
	r4211->SetTextSize(25);
	r4211->SetTextFont(43);
	r4211->SetTextColor(kBlack);
	r4211->Draw("same");

    TLegend *leg = new TLegend(0.68,0.7,0.72,0.85);
    leg->AddEntry(hF_b_MC_2d, " MC", "L");
    leg->AddEntry(hF_b_REC_2d, " Method L", "PE"); 
    leg->AddEntry(hF_b_REC_new_2d, " Projection method", "PE"); 
    leg->SetBorderSize(0);      
    leg->SetFillStyle(0);       
    leg->SetTextFont(45);       
    leg->SetTextSize(20);    
    leg->Draw("same");

    c1->Print("./plot_transform_Fb.pdf");
}
