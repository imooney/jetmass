//This macro is a compendium of various plotting routines on the ppY12 dataset.
//For the most part, it reproduces plots from Raghav's internal group update on 11/30/2017 and his JetCorr update on 2/27/2018.

#include <iostream>
#include <string>
//#include "RooUnfold.h"
#include "Plots.h"

using namespace std;

void JetQA (TFile* simFile, TFile* dataFile, const string flag1, const string flag2, const string flag3, const string out, const string filetype) { 
  //~~~~~~~~~~~~~~CANVASES,TREES,HISTOGRAMS~~~~~~~~~~~~~~//
  
  TCanvas *cpt = MakeCanvas("cpt", "y",800,800);
  TCanvas *ceta = MakeCanvas("ceta","0",800,800);
  TCanvas *cphi = MakeCanvas("cphi","0",800,800);
  TCanvas *cm = MakeCanvas("cm","0",800,800);
  
  TH3D *ptetaphi_py = (TH3D*) simFile->Get(("PtEtaPhi_" + flag1 + "_" + flag2 + "_" + flag3 + "_" + "py").c_str());
  TH3D *ptetaphi_ge = (TH3D*) simFile->Get(("PtEtaPhi_" + flag1 + "_" + flag2 + "_" + flag3 + "_" + "ge").c_str());
  TH1D *m_py = (TH1D*) simFile->Get(("m_" + flag1 + "_" + flag2 + "_" + flag3 + "_" + "py").c_str());
  TH1D *m_ge = (TH1D*) simFile->Get(("m_" + flag1 + "_" + flag2 + "_" + flag3 + "_" + "ge").c_str());
  
  TH3D *ptetaphi_dat= (TH3D*) dataFile->Get(("PtEtaPhi_" + flag1 + "_" + flag2 + "_" + flag3).c_str());
  TH1D *m_dat = (TH1D*) dataFile->Get(("m_" + flag1 + "_" + flag2 + "_" + flag3).c_str());
  
  TTree *d_incl = (TTree*) dataFile->Get("incl");
  TTree *p_incl = (TTree*) simFile->Get("py_inclTree");
  TTree *g_incl = (TTree*) simFile->Get("ge_inclTree");
  
  double pT, dummy_d_weight, g_pt, g_weight, p_pt, p_weight;
  
  d_incl->SetBranchAddress("Pt", &pT); p_incl->SetBranchAddress("Pt",&p_pt); g_incl->SetBranchAddress("Pt",&g_pt);
  p_incl->SetBranchAddress("weight", &p_weight); g_incl->SetBranchAddress("weight",&g_weight);
  dummy_d_weight = 1.0;
  
  //making variable bin size histogram                                                                                                                                  
  const int nBins_pt = 8;
  double edges[nBins_pt + 1] = {5,10,15,20,25,30,40,60,100};
  
  TH1D * d_pT = HistFromTree("d_pT", nBins_pt, edges, d_incl, pT, dummy_d_weight);
  TH1D * g_pT = HistFromTree("g_pT", nBins_pt, edges, g_incl, g_pt, g_weight);
  TH1D * p_pT = HistFromTree("p_pT", nBins_pt, edges, p_incl, p_pt, p_weight);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PROJECTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TH2D* etaphi_dat = (TH2D*)ptetaphi_dat->Project3D("zy"); TH2D* etaphi_py =(TH2D*)ptetaphi_py->Project3D("zy"); TH2D* etaphi_ge = (TH2D*)ptetaphi_ge->Project3D("zy");
  const int nHists_etaphi = 1;
  double proj_ranges[nHists_etaphi + 1] = {0,1000};
  
  vector<TH1D*> eta_dat_arr = Projection2D(etaphi_dat, nHists_etaphi, proj_ranges, "y");
  vector<TH1D*> phi_dat_arr = Projection2D(etaphi_dat, nHists_etaphi, proj_ranges, "x");
  vector<TH1D*> eta_ge_arr = Projection2D(etaphi_ge, nHists_etaphi, proj_ranges, "y");
  vector<TH1D*> phi_ge_arr = Projection2D(etaphi_ge, nHists_etaphi, proj_ranges, "x");
  vector<TH1D*> eta_py_arr = Projection2D(etaphi_py, nHists_etaphi, proj_ranges, "y");
  vector<TH1D*> phi_py_arr = Projection2D(etaphi_py, nHists_etaphi, proj_ranges, "x");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                 
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TLEGEND~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TLegend * tpt = TitleLegend(0.244, 0.15, 0.7, 0.35);
  tpt->AddEntry((TObject*)0,(flag1 + " " + flag2 + " " + flag3).c_str(), "");
  tpt->AddEntry(p_pT, "PYTHIA6", "p");
  tpt->AddEntry(g_pT, "PYTHIA6+GEANT", "p");
  tpt->AddEntry(d_pT, "STAR data", "p");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRETTIFYING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                  
  Prettify1D(p_pT, kBlue, kOpenCircle, 1, kBlue, "p^{jet}_{T} [GeV/c]", "prob.", 0, 100, 0, -1);
  Prettify1D(g_pT, kRed, kOpenSquare, 1, kRed, "p^{jet}_{T} [GeV/c]", "prob.", 0, 100, 0, -1);
  Prettify1D(d_pT, kBlack, kFullStar, 1, kBlack, "p^{jet}_{T} [GeV/c]", "prob.", 0, 100, 0, -1);
  
  for (int i = 0; i < nHists_etaphi; ++ i) {
    Prettify1D(eta_dat_arr[i], kBlack, kFullStar, 1, kBlack, "#eta^{jet}", "prob.", -1, 1, 0, 0.05);
    Prettify1D(phi_dat_arr[i], kBlack, kFullStar, 1, kBlack, "#phi^{jet}", "prob.", -1, -1, 0, 0.035);
    Prettify1D(eta_ge_arr[i], kRed, kOpenSquare, 1, kRed, "#eta^{jet}", "prob.", -1, 1, 0, 0.05);
    Prettify1D(phi_ge_arr[i], kRed, kOpenSquare, 1, kRed, "#phi^{jet}", "prob.", -1, -1, 0, 0.035);
    Prettify1D(eta_py_arr[i], kBlue, kOpenCircle, 1, kBlue, "#eta^{jet}", "prob.", -1, 1, 0, 0.05);
    Prettify1D(phi_py_arr[i], kBlue, kOpenCircle, 1, kBlue, "#phi^{jet}", "prob.", -1, -1, 0, 0.035);
  }
  Prettify1D(m_py, kBlue, kOpenCircle, 1, kBlue, "M^{jet} [GeV/c^{2}]", "prob.", 0, 16, -1, -1);
  Prettify1D(m_ge, kRed, kOpenSquare, 1, kRed, "M^{jet} [GeV/c^{2}]", "prob.", 0, 16, -1, -1);
  Prettify1D(m_dat, kBlack, kFullStar, 1, kBlack, "M^{jet} [GeV/c^{2}]", "prob.", 0, 16, -1, -1);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                   
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DRAWING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                   
  cpt->cd(); p_pT->Draw("same"); g_pT->Draw("same"); d_pT->Draw("same"); tpt->Draw("same");
  ceta->cd(); eta_py_arr[0]->Draw("same"); eta_ge_arr[0]->Draw("same"); eta_dat_arr[0]->Draw("same"); tpt->Draw("same");
  cphi->cd(); phi_py_arr[0]->Draw("same"); phi_ge_arr[0]->Draw("same"); phi_dat_arr[0]->Draw("same"); tpt->Draw("same");
  cm->cd(); m_py->Draw("same"); m_ge->Draw("same"); m_dat->Draw("same"); tpt->Draw("same");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                   
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SAVING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                   
  
    cpt->SaveAs((out + "pt_" + flag1 + "_" + flag2 + "_" + flag3 + filetype).c_str());
    ceta->SaveAs((out + "eta_" + flag1 + "_" + flag2 + "_" + flag3 + filetype).c_str());
    cphi->SaveAs((out + "phi_" + flag1 + "_" + flag2 + "_" + flag3 + filetype).c_str());
    cm->SaveAs((out + "m_" + flag1 + "_" + flag2 + "_" + flag3 + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//    
  
  return;
}

//Plots jet observables! M, pT, zg, Rg
void Observables(TFile *matchFile, TFile *simFile, TFile *dataFile, const string out, const string filetype) {
  //plots the mass for various pT ranges
  ObservablePtSlices((TH2D*) dataFile->Get("m_v_pt_full_incl_jet"), (TH2D*) simFile->Get("m_v_pt_full_incl_jet_py"), (TH2D*) simFile->Get("m_v_pt_full_incl_jet_ge"), out, filetype, "M_{jet} [GeV/c^{2}]", 0, 0);
  //plots the ch fraction for various pT ranges
  ObservablePtSlices((TH2D*) dataFile->Get("ch_frac_v_pt"), (TH2D*) simFile->Get("ch_frac_v_pt_py"), (TH2D*) simFile->Get("ch_frac_v_pt_ge"), out, filetype, "ch. frac.", 0, 0);
  THnSparseD* sd_d = (THnSparseD*) dataFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd");
  THnSparseD* sd_py = (THnSparseD*) simFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd_py");
  THnSparseD* sd_ge = (THnSparseD*) simFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd_ge");
  //plots the zg for various pT ranges
  ObservablePtSlices((TH2D*) (sd_d->Projection(3,0)), (TH2D*) (sd_py->Projection(3,0)), (TH2D*) (sd_ge->Projection(3,0)), out, filetype, "z_{g}", 0, 0);
  //plots the rg for various pT ranges
  ObservablePtSlices((TH2D*) (sd_d->Projection(3,2)), (TH2D*) (sd_py->Projection(3,2)), (TH2D*) (sd_ge->Projection(3,2)), out, filetype, "R_{g}", 0, 0);
  return;
}

//Plots jet observable resolutions in slices of pT! M, pT, zg, Rg
void Slices(TFile *matchFile, const string out, const string filetype) {
  TH2D * h_dummy = new TH2D("h_dummy","",1,0,1,1,0,1);
  //plots slices in pT of jet observable resolutions
  ObservablePtSlices (h_dummy, h_dummy, (TH2D*) matchFile->Get("ratioMvPyPt"), out, filetype, "M_{jet}^{det} / M_{jet}^{gen}", 1, 0);
  ObservablePtSlices (h_dummy, h_dummy, (TH2D*) matchFile->Get("ratioPtvPyPt"), out, filetype, "p_{T}^{det-jet} / p_{T}^{gen-jet}", 1, 0);
  ObservablePtSlices (h_dummy, h_dummy, (TH2D*) matchFile->Get("ratioZgvPyPt"), out, filetype, "z_{g}^{det} / z_{g}^{gen}", 1, 1); //NOTE: sliced in bins of GROOMED pT
  ObservablePtSlices (h_dummy, h_dummy, (TH2D*) matchFile->Get("ratioRgvPyPt"), out, filetype, "R_{g}^{det} / R_{g}^{gen}", 1, 1); //NOTE: sliced in bins of GROOMED pT
  return;
}

void Resolutions(TFile *matchFile, const string out, const string filetype) {
  //plots jet observable resolutions v. gen-level pT and associated mean & RMS for each pT bin
  Resolution((TH2D*) matchFile->Get("deltaMvPyPt"), out, filetype, "#Delta M_{jet} (Det - Gen) / M^{gen}_{jet}", 0);
  Resolution((TH2D*) matchFile->Get("deltaPtvPyPt"), out, filetype, "#Delta p_{T}^{jet} (Det - Gen) / p_{T}^{gen-jet}", 0);
  Resolution((TH2D*) matchFile->Get("deltaZgvPyPt"), out, filetype, "#Delta z_{g} (Det - Gen)", 1);
  Resolution((TH2D*) matchFile->Get("deltaRgvPyPt"), out, filetype, "#Delta R_{g} (Det - Gen)", 1);
  return;
}

//bins pT finely or coarsely for the subsequent call to "Response()"
vector<TH1D*> PtBinCorrectly(TFile* simFile) {
  double pT, g_pt, g_weight, p_pt, p_weight;
  
  TTree* py = (TTree*) simFile->Get("py_inclTree");
  TTree* ge = (TTree*) simFile->Get("ge_inclTree");
  
  py->SetBranchAddress("Pt",&p_pt); ge->SetBranchAddress("Pt",&g_pt);
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);
  
  TH1D* py_pt = HistFromTree(py->GetName(), 80, 0, 80, py, p_pt, p_weight);
  TH1D* ge_pt = HistFromTree(ge->GetName(), 80, 0, 80, ge, g_pt, g_weight);
  TH1D* py_ptcoarse = HistFromTree("py_ptcoarse", 15, 5, 80, py, p_pt, p_weight);                                                              
  TH1D* ge_ptcoarse = HistFromTree("ge_ptcoarse", 13, 15, 80, ge, g_pt, g_weight);
  
  vector<TH1D*> hists = {py_pt, ge_pt, py_ptcoarse, ge_ptcoarse};
  
  return hists;
}

void Responses(TFile* matchFile, TFile* simFile, const THnSparseD* sd_py, const THnSparseD* sd_ge, vector<TH1D*> pts, const string out, const string filetype) {
  Response(matchFile, "pt_response", pts[0], pts[1], "p_{T}^{jet} [GeV/c]", out, filetype);
  Response(matchFile, "pt_res_coarse", pts[2], pts[3], "p_{T}^{jet} [GeV/c]", out, filetype);
  Response(matchFile, "m_response", (TH1D*) simFile->Get("m_full_incl_jet_py"), (TH1D*) simFile->Get("m_full_incl_jet_ge"), "M^{jet} [GeV/c^{2}]", out, filetype);
  Response(matchFile, "zg_response", (TH1D*) (sd_py->Projection(0)), (TH1D*) (sd_ge->Projection(0)), "z_{g}", out, filetype);
  Response(matchFile, "rg_response", (TH1D*) (sd_py->Projection(2)), (TH1D*) (sd_ge->Projection(2)), "R_{g}", out, filetype);
  
  return;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void Y12InitialAnalysis () {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILES&STYLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~//    
  
  string dir = "~/jetmass_temp/";
  string simin = "out/sim/";
  string datain = "out/data/";
  string matchin = "out/matching/";
  string file = "full.root";
  string out = "~/jetmass_temp/plots/Y12ResolutionsandResponses/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";
  string flag3 = "jet";
  
  TFile* simFile = new TFile( (dir + simin + file).c_str(), "READ");
  TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  TFile* matchFile = new TFile( (dir + matchin + file).c_str(), "READ");
  
  THnSparseD* sd_py = (THnSparseD*) simFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd_py");
  THnSparseD* sd_ge = (THnSparseD*) simFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd_ge");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  //plots pT, eta, phi, and mass
  JetQA(simFile, dataFile, flag1, flag2, flag3, out, filetype);
  //plots the pT response matrix and the simulation pT spectra for fine-grained and coarse-grained binning
  vector<TH1D*> pts = PtBinCorrectly(simFile);
  Responses(matchFile, simFile, sd_py, sd_ge, pts, out, filetype);
  
  Observables(matchFile, simFile, dataFile, out, filetype); //plots jet observables: M, pT, zg, Rg
  Slices(matchFile, out, filetype); //plots the resolution of the observables in bins of pT
  Resolutions(matchFile, out, filetype); //plots 2D observable resolution with pT and associated means & RMSs
  
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
