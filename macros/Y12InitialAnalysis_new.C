//This macro is a compendium of various plotting routines on the ppY12 dataset.
//For the most part, it reproduces plots from Raghav's internal group update on 11/30/2017 and his JetCorr update on 2/27/2018.

#include <iostream>
#include <string>
#include <math.h>
//#include "RooUnfold.h"
#include "Plots_new.h"

using namespace std;

void JetQA (TFile* inFile, const string flag1, const string flag2, const string out, const string filetype) {
  //~~~~~~~~~~~~~~CANVASES,TREES,HISTOGRAMS~~~~~~~~~~~~~~//
  TCanvas *cpt = MakeCanvas("cpt", "y",800,800);
  TCanvas *ceta = MakeCanvas("ceta","0",800,800);
  TCanvas *cphi = MakeCanvas("cphi","0",800,800);
  TCanvas *cm = MakeCanvas("cm","0",800,800);
    
  vector<TH1D*> pts = {(TH1D*) inFile->Get("pt_p_var_bin"), (TH1D*) inFile->Get("pt_g_var_bin"), (TH1D*) inFile->Get("pt_d_var_bin")};
  vector<TH1D*> etas = {(TH1D*) inFile->Get("eta_p"), (TH1D*) inFile->Get("eta_g"), (TH1D*) inFile->Get("eta_d")};
  vector<TH1D*> phis = {(TH1D*) inFile->Get("phi_p"), (TH1D*) inFile->Get("phi_g"), (TH1D*) inFile->Get("phi_d")};
  vector<TH1D*> ms = {(TH1D*) inFile->Get("m_p"), (TH1D*) inFile->Get("m_g"), (TH1D*) inFile->Get("m_d")};
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TLEGEND~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TLegend * tpt = TitleLegend(0.244, 0.15, 0.7, 0.35);
  tpt->AddEntry((TObject*)0,(flag1 + " " + flag2).c_str(), "");
  tpt->AddEntry(pts[0], "PYTHIA6", "p");
  tpt->AddEntry(pts[1], "PYTHIA6+GEANT", "p");
  tpt->AddEntry(pts[2], "STAR data", "p");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRETTIFYING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                  
  Prettify1D(pts[0], kBlue, kOpenCircle, 1, kBlue, "p^{jet}_{T} [GeV/c]", "prob.", 5, 100, -1,-1);
  Prettify1D(pts[1], kRed, kOpenSquare, 1, kRed, "p^{jet}_{T} [GeV/c]", "prob.", 5, 100, -1,-1);
  Prettify1D(pts[2], kBlack, kFullStar, 1, kBlack, "p^{jet}_{T} [GeV/c]", "prob.", 5, 100, -1,-1);
  Prettify1D(etas[0], kBlue, kOpenCircle, 1, kBlue, "#eta^{jet}", "prob.", -1,1, 0,0.05);
  Prettify1D(etas[1], kRed, kOpenSquare, 1, kRed, "#eta^{jet}", "prob.", -1,1, 0,0.05);
  Prettify1D(etas[2], kBlack, kFullStar, 1, kBlack, "#eta^{jet}", "prob.", -1,1, 0,0.05);
  Prettify1D(phis[0], kBlue, kOpenCircle, 1, kBlue, "#phi^{jet}", "prob.", 0,2*M_PI, 0,0.035);
  Prettify1D(phis[1], kRed, kOpenSquare, 1, kRed, "#phi^{jet}", "prob.", 0, 2*M_PI, 0,0.035);
  Prettify1D(phis[2], kBlack, kFullStar, 1, kBlack, "#phi^{jet}", "prob.", 0, 2*M_PI, 0,0.035);
  Prettify1D(ms[0], kBlue, kOpenCircle, 1, kBlue, "M^{jet} [GeV/c^{2}]", "prob.", 0, 10, -1,-1);
  Prettify1D(ms[1], kRed, kOpenSquare, 1, kRed, "M^{jet} [GeV/c^{2}]", "prob.", 0, 10, -1,-1);
  Prettify1D(ms[2], kBlack, kFullStar, 1, kBlack, "M^{jet} [GeV/c^{2}]", "prob.", 0, 10, -1,-1);
 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                         
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DRAWING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                   
  cpt->cd(); pts[0]->Draw("same"); pts[1]->Draw("same"); pts[2]->Draw("same"); tpt->Draw("same");
  ceta->cd(); etas[0]->Draw("same"); etas[1]->Draw("same"); etas[2]->Draw("same"); tpt->Draw("same");
  cphi->cd(); phis[0]->Draw("same"); phis[1]->Draw("same"); phis[2]->Draw("same"); tpt->Draw("same");
  cm->cd(); ms[0]->Draw("same"); ms[1]->Draw("same"); ms[2]->Draw("same"); tpt->Draw("same");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                         
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SAVING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                   
  
  cpt->SaveAs((out + "pt_" + flag1 + "_" + flag2 + filetype).c_str());
  ceta->SaveAs((out + "eta_" + flag1 + "_" + flag2 + filetype).c_str());
  cphi->SaveAs((out + "phi_" + flag1 + "_" + flag2 + filetype).c_str());
  cm->SaveAs((out + "m_" + flag1 + "_" + flag2 + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//    
  
  return;
}


//Plots jet observable resolutions in slices of pT! M, pT, zg, Rg
void Slices(TFile *inFile, const string out, const string filetype) {
  TH2D * h_dummy = new TH2D("h_dummy","",1,0,1,1,0,1);
  TH2D *ratPt = (TH2D*) inFile->Get("ratioPtvPyPt");
  TH2D *ratM = (TH2D*) inFile->Get("ratioMvPyPt");
  TH2D *ratZg = (TH2D*) inFile->Get("ratioZgvPyPt");
  TH2D *ratRg = (TH2D*) inFile->Get("ratioRgvPyPt");
  TH2D *ratPtg = (TH2D*) inFile->Get("ratioPtgvPyPt");
  TH2D *ratMg = (TH2D*) inFile->Get("ratioMgvPyPt");
  const int nHists = 5;
  double pt_bins[nHists+1] = {2,3,4,5,7,11}; //pT is usually 11 5 GeV bins from 5 to 60 GeV, so e.g. bin 0 = 5, bin 1 = 10, etc.
  
  //plots slices in pT of jet observable resolutions
  ObservablePtSlices (h_dummy, h_dummy, ratPt, out, filetype, "p_{T}^{det-jet} / p_{T}^{gen-jet}", 1, 0,pt_bins);
  ObservablePtSlices (h_dummy, h_dummy, ratM, out, filetype, "M_{jet}^{det} / M_{jet}^{gen}", 1, 0,pt_bins);
  ObservablePtSlices (h_dummy, h_dummy, ratZg, out, filetype, "z_{g}^{det} / z_{g}^{gen}", 1, 1,pt_bins);
  ObservablePtSlices (h_dummy, h_dummy, ratRg, out, filetype, "R_{g}^{det} / R_{g}^{gen}", 1, 1,pt_bins);
  ObservablePtSlices (h_dummy, h_dummy, ratPtg, out, filetype, "p_{T,g}^{det} / p_{T,g}^{gen}", 1, 1,pt_bins);
  ObservablePtSlices (h_dummy, h_dummy, ratMg, out, filetype, "M_{g}^{det} / M_{g}^{gen}", 1, 1,pt_bins);
  
  return;
}

void Resolutions(TFile *inFile, const string out, const string filetype) {
  TH2D *delPt = (TH2D*) inFile->Get("deltaPtvPyPt");
  TH2D *delM = (TH2D*) inFile->Get("deltaMvPyPt");
  TH2D *delZg = (TH2D*) inFile->Get("deltaZgvPyPt");
  TH2D *delRg = (TH2D*) inFile->Get("deltaRgvPyPt");
  TH2D *delPtg = (TH2D*) inFile->Get("deltaPtgvPyPt");
  TH2D *delMg = (TH2D*) inFile->Get("deltaMgvPyPt");
  
  //plots jet observable resolutions v. gen-level pT and associated mean & RMS for each pT bin
  Resolution(delPt, out, filetype, "#Delta p_{T} (Det - Gen) / p_{T}^{gen}", 0);
  Resolution(delM, out, filetype, "#Delta M (Det - Gen) / M^{gen}", 0);
  Resolution(delZg, out, filetype, "#Delta z_{g} (Det - Gen) / z_{g}^{gen}", 1);
  Resolution(delRg, out, filetype, "#Delta R_{g} (Det - Gen) / R_{g}^{gen}", 1);
  Resolution(delPtg, out, filetype, "#Delta p_{T,g} (Det - Gen) / p_{T,g}^{gen}", 1);
  Resolution(delMg, out, filetype, "#Delta M_{g} (Det - Gen) / M_{g}^{gen}", 1);
  
  return;
}



//Plots jet observables! M, pT, zg, Rg
void Observables(TFile *inFile, const string out, const string filetype) {
  vector<TH2D*> mpt = {(TH2D*) inFile->Get("m_v_pt_d"), (TH2D*) inFile->Get("m_v_pt_p"), (TH2D*) inFile->Get("m_v_pt_g")};
  vector<TH2D*> chfracpt = {(TH2D*) inFile->Get("ch_frac_v_pt_d"), (TH2D*) inFile->Get("ch_frac_v_pt_p"), (TH2D*) inFile->Get("ch_frac_v_pt_g")};
  vector<TH2D*> zgpt = {(TH2D*) inFile->Get("zg_v_pt_d"), (TH2D*) inFile->Get("zg_v_pt_p"), (TH2D*) inFile->Get("zg_v_pt_g")};
  vector<TH2D*> rgpt = {(TH2D*) inFile->Get("rg_v_pt_d"), (TH2D*) inFile->Get("rg_v_pt_p"), (TH2D*) inFile->Get("rg_v_pt_g")};
  vector<TH2D*> ptgpt = {(TH2D*) inFile->Get("ptg_v_pt_d"), (TH2D*) inFile->Get("ptg_v_pt_p"), (TH2D*) inFile->Get("ptg_v_pt_g")};
  vector<TH2D*> mgpt = {(TH2D*) inFile->Get("mg_v_pt_d"), (TH2D*) inFile->Get("mg_v_pt_p"), (TH2D*) inFile->Get("mg_v_pt_g")};
  const int nHists = 5;
  double det_pt_bins[nHists+1] = {0,1,2,3,5,9}; //pT is usually 11 5 GeV bins from 5 to 60 GeV, so e.g. bin 0 = 5, bin 1 = 10, etc.  
  //plots the mass for various pT ranges
  ObservablePtSlices(mpt[0], mpt[1], mpt[2], out, filetype, "M_{jet} [GeV/c^{2}]", 0, 0, det_pt_bins);
  //plots the ch fraction for various pT ranges
  ObservablePtSlices(chfracpt[0], chfracpt[1], chfracpt[2], out, filetype, "ch. frac.", 0, 0, det_pt_bins);
  //plots the zg for various pT ranges
  ObservablePtSlices(zgpt[0], zgpt[1], zgpt[2], out, filetype, "z_{g}", 0, 0, det_pt_bins);
  //plots the rg for various pT ranges
  ObservablePtSlices(rgpt[0], rgpt[1], rgpt[2], out, filetype, "R_{g}", 0, 0, det_pt_bins);
  //plots the ptg for various pT ranges
  ObservablePtSlices(ptgpt[0], ptgpt[1], ptgpt[2], out, filetype, "p_{T,g}", 0, 0, det_pt_bins);
  //plots the mg for various pT ranges
  ObservablePtSlices(mgpt[0], mgpt[1], mgpt[2], out, filetype, "M_{g}", 0, 0, det_pt_bins);  

  return;
}


void Responses(TFile* matchFile, TFile* inFile, const string out, const string filetype) {
  TH1D* pt_p = (TH1D*) inFile->Get("pt_p"); TH1D* pt_g = (TH1D*) inFile->Get("pt_g");
  TH1D* m_p = (TH1D*) inFile->Get("m_p"); TH1D* m_g = (TH1D*) inFile->Get("m_g");
  TH1D* zg_p = (TH1D*) inFile->Get("zg_p"); TH1D* zg_g = (TH1D*) inFile->Get("zg_g");
  TH1D* rg_p = (TH1D*) inFile->Get("rg_p"); TH1D* rg_g = (TH1D*) inFile->Get("rg_g");
  TH1D* ptg_p = (TH1D*) inFile->Get("ptg_p"); TH1D* ptg_g = (TH1D*) inFile->Get("ptg_g");
  TH1D* mg_p = (TH1D*) inFile->Get("mg_p"); TH1D* mg_g = (TH1D*) inFile->Get("mg_g");

  Response(matchFile, "pt_res_coarse", pt_p, pt_g, "p_{T}^{j} [GeV/c]", out, filetype);
  Response(matchFile, "m_response", m_p, m_g, "M^{j} [GeV/c^{2}]", out, filetype);
  Response(matchFile, "zg_response", zg_p, zg_g, "z_{g}", out, filetype);
  Response(matchFile, "rg_response", rg_p, rg_g, "R_{g}", out, filetype);
  Response(matchFile, "ptg_response", ptg_p, ptg_g, "p_{T,g}", out, filetype);
  Response(matchFile, "mg_response", mg_p, mg_g, "M_{g}", out, filetype);
  
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void Y12InitialAnalysis_new () {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILES&STYLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~//    
  
  string dir = "~/jetmass/";
  string gein = "out/sim/ge/";
  string pyin = "out/sim/py/";
  string datain = "out/data/";
  string matchin = "out/matching/";
  string in = "macros/hists/";
  string infile = "hists.root";
  string file = "full.root";
  string out = "~/jetmass/plots/Y12ResolutionsandResponses/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";
  
  //  TFile* pyFile = new TFile( (dir + pyin + file).c_str(), "READ");
  //TFile* geFile = new TFile( (dir + gein + file).c_str(), "READ");
  //TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  TFile* matchFile = new TFile( (dir + matchin + file).c_str(), "READ");
  TFile* inFile = new TFile( (dir + in + infile).c_str(), "READ");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //plots pT, eta, phi, and mass
  //JetQA(inFile, flag1, flag2, out, filetype);
  //plots the pT response matrix and the simulation pT spectra for fine-grained and coarse-grained binning
  //Responses(matchFile, inFile, out, filetype);
    Observables(inFile, out, filetype); //plots jet observables: M, pT, zg, Rg, Ptg, Mg
  //  Slices(inFile, out, filetype); //plots the resolution of the observables in bins of pT
  //Resolutions(inFile, out, filetype); //plots 2D observable resolution with pT and associated means & RMSs
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
