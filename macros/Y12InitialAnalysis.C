//This macro is a compendium of various plotting routines on the ppY12 dataset.
//For the most part, it reproduces plots from Raghav's internal group update on 11/30/2017 and his JetCorr update on 2/27/2018.

#include <iostream>
#include <string>
#include <math.h>
//#include "RooUnfold.h"
#include "Plots.h"

using namespace std;

//bins simulation observable in desired manner from tree
vector<TH1D*> PtBinCorrectly(TFile* dataFile, TFile* pyFile, TFile* geFile, const int nBins_py, const double lo_py, const double hi_py, const int nBins_ge, const double lo_ge, const double hi_ge, const string obs, const string append) {
  vector<double> *d_obs = 0; vector<double> *g_obs = 0; vector<double> *p_obs = 0;
  double dummy_d_weight = 1; double g_weight, p_weight;

  TTree* py = (TTree*) pyFile->Get("event");
  TTree* ge = (TTree*) geFile->Get("event");
  TTree* d = (TTree*) dataFile->Get("event");
  
  py->SetBranchAddress(obs.c_str(),&p_obs); ge->SetBranchAddress(obs.c_str(),&g_obs); d->SetBranchAddress(obs.c_str(),&d_obs);
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);
  
  TH1D * py_obs = HistFromTree(("py_" + obs + "_" + append).c_str(), nBins_py, lo_py, hi_py, py, p_obs, p_weight);
  TH1D * ge_obs = HistFromTree(("ge_" + obs + "_" + append).c_str(), nBins_ge, lo_ge, hi_ge, ge, g_obs, g_weight);
  TH1D * data_obs = HistFromTree(("data_" + obs + "_" + append).c_str(), nBins_ge, lo_ge, hi_ge, d, d_obs, dummy_d_weight);
  
  vector<TH1D*> hists = {py_obs, ge_obs, data_obs};
  return hists;
}


//bins simulation observable in desired manner from tree
vector<TH1D*> PtBinCorrectly(TFile* dataFile, TFile* pyFile, TFile* geFile, const int nBins_py, double * edges_py, const int nBins_ge, double * edges_ge, const std::string obs, const string append) {
  vector<double> *d_obs = 0; vector<double> *g_obs = 0; vector<double> *p_obs = 0;
  double dummy_d_weight = 1; double g_weight, p_weight;

  TTree* py = (TTree*) pyFile->Get("event");
  TTree* ge = (TTree*) geFile->Get("event");
  TTree* d = (TTree*) dataFile->Get("event");
  
  py->SetBranchAddress(obs.c_str(),&p_obs); ge->SetBranchAddress(obs.c_str(),&g_obs); d->SetBranchAddress(obs.c_str(),&d_obs);
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);
  
  TH1D * py_obs = HistFromTree(("py_" + obs).c_str(), nBins_py, edges_py, py, p_obs, p_weight);
  TH1D * ge_obs = HistFromTree(("ge_" + obs).c_str(), nBins_ge, edges_ge, ge, g_obs, g_weight);
  TH1D * data_obs = HistFromTree(("data_" + obs).c_str(), nBins_ge, edges_ge, d, d_obs, dummy_d_weight);
  
  vector<TH1D*> hists = {py_obs, ge_obs, data_obs};
  return hists;
}



//bins simulation observable in desired manner from tree
vector<TH2D*> PtBinCorrectly(TFile* dataFile, TFile* pyFile, TFile* geFile, const int xnBins_py, const double xlo_py, const double xhi_py, const int ynBins_py, const double ylo_py, const double yhi_py, const int xnBins_ge, const double xlo_ge, const double xhi_ge, const int ynBins_ge, const double ylo_ge, const double yhi_ge, const std::string xobs, const std::string yobs) {
  cout << "AAAA" << endl;
  vector<double> *d_obsx = 0; vector<double> *g_obsx = 0; vector<double> *p_obsx = 0;
  vector<double> *d_obsy = 0; vector<double> *g_obsy = 0; vector<double> *p_obsy = 0;
  cout << "BBBB" << endl;
  double g_weight, p_weight;
  double dummy_d_weight = 1;
  cout << "CCCC" << endl;
  TTree* py = (TTree*) pyFile->Get("event");
  TTree* ge = (TTree*) geFile->Get("event");
  TTree* d = (TTree*) dataFile->Get("event");
  cout << "DDDD" << endl;
  py->SetBranchAddress(xobs.c_str(),&p_obsx); ge->SetBranchAddress(xobs.c_str(),&g_obsx); d->SetBranchAddress(xobs.c_str(),&d_obsx);
  py->SetBranchAddress(yobs.c_str(),&p_obsy); ge->SetBranchAddress(yobs.c_str(),&g_obsy); d->SetBranchAddress(yobs.c_str(),&d_obsy);  
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);
  cout << "EEEE" << endl;
  TH2D * py_obs = HistFromTree(("py_" + xobs + "_" + yobs).c_str(), xnBins_py, xlo_py, xhi_py, ynBins_py, ylo_py, yhi_py, py, p_obsx, p_obsy, p_weight);
  TH2D * ge_obs = HistFromTree(("ge_" + xobs + "_" + yobs).c_str(), xnBins_ge, xlo_ge, xhi_ge, ynBins_ge, ylo_ge, yhi_ge, ge, g_obsx, g_obsy, g_weight);
  TH2D * d_obs = HistFromTree(("d_" + xobs + "_" + yobs).c_str(), xnBins_ge, xlo_ge, xhi_ge, ynBins_ge, ylo_ge, yhi_ge, d, d_obsx, d_obsy, dummy_d_weight);
  cout << "FFFF" << endl;
  vector<TH2D*> hists = {py_obs, ge_obs, d_obs};
  cout << "GGGG" << endl;
  return hists;
}


//bins simulation observable in desired manner from tree
vector<TH2D*> PtBinCorrectly(TFile* dataFile, TFile* pyFile, TFile* geFile, const int xnBins_py, const double xlo_py, const double xhi_py, const int ynBins_py, double* yedges_py, const int xnBins_ge, const double xlo_ge, const double xhi_ge, const int ynBins_ge, double *yedges_ge, const std::string xobs, const std::string yobs) {
  vector<double> *d_obsx = 0; vector<double> *g_obsx = 0; vector<double> *p_obsx = 0;
  vector<double> *d_obsy = 0; vector<double> *g_obsy = 0; vector<double> *p_obsy = 0;
   
  double g_weight, p_weight;
  double dummy_d_weight = 1;
  
  TTree* py = (TTree*) pyFile->Get("event");
  TTree* ge = (TTree*) geFile->Get("event");
  TTree* d = (TTree*) dataFile->Get("event");
  
  py->SetBranchAddress(xobs.c_str(),&p_obsx); ge->SetBranchAddress(xobs.c_str(),&g_obsx); d->SetBranchAddress(xobs.c_str(),&d_obsx);
  py->SetBranchAddress(yobs.c_str(),&p_obsy); ge->SetBranchAddress(yobs.c_str(),&g_obsy); d->SetBranchAddress(yobs.c_str(),&d_obsy);  
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);
  
  TH2D * py_obs = HistFromTree(("py_" + xobs + "_" + yobs).c_str(), xnBins_py, xlo_py, xhi_py, ynBins_py, yedges_py, py, p_obsx, p_obsy, p_weight);
  TH2D * ge_obs = HistFromTree(("ge_" + xobs + "_" + yobs).c_str(), xnBins_ge, xlo_ge, xhi_ge, ynBins_ge, yedges_ge, ge, g_obsx, g_obsy, g_weight);
  TH2D * d_obs = HistFromTree(("d_" + xobs + "_" + yobs).c_str(), xnBins_ge, xlo_ge, xhi_ge, ynBins_ge, yedges_ge, d, d_obsx, d_obsy, dummy_d_weight);
  
  vector<TH2D*> hists = {py_obs, ge_obs, d_obs};
  
  return hists;
}


//bins simulation observable in desired manner from tree
TH2D* PtBinCorrectly(TFile* matchFile, const int xnBins, const double xlo, const double xhi, const int ynBins, double* yedges, const std::string xobs, const std::string yobs) {
  vector<double> *obsx = 0; vector<double> *obsy = 0;
   
  double wt;
  
  TTree* m = (TTree*) matchFile->Get("event");
  
  m->SetBranchAddress(xobs.c_str(),&obsx); m->SetBranchAddress(yobs.c_str(),&obsy);  
  m->SetBranchAddress("weight", &wt);
  
  TH2D * m_obs = HistFromTree(("m_" + xobs + "_" + yobs).c_str(), xnBins, xlo, xhi, ynBins, yedges, m, obsx, obsy, wt);
  
  return m_obs;
}

//bins simulation observable in desired manner from tree
TH2D* PtBinCorrectly(TFile* matchFile, const int xnBins, const double xlo, const double xhi, const int ynBins, const double ylo, const double yhi, const std::string xobs, const std::string yobs) {
  vector<double> *obsx = 0; vector<double> *obsy = 0;
   
  double wt;
  
  TTree* m = (TTree*) matchFile->Get("event");
  
  m->SetBranchAddress(xobs.c_str(),&obsx); m->SetBranchAddress(yobs.c_str(),&obsy);  
  m->SetBranchAddress("weight", &wt);
  
  TH2D * m_obs = HistFromTree(("m_" + xobs + "_" + yobs).c_str(), xnBins, xlo, xhi, ynBins, ylo, yhi, m, obsx, obsy, wt);
  
  return m_obs;
}


void JetQA (TFile* pyFile, TFile* geFile, TFile* dataFile, const string flag1, const string flag2, const string out, const string filetype) { 
  //~~~~~~~~~~~~~~CANVASES,TREES,HISTOGRAMS~~~~~~~~~~~~~~//
  TCanvas *cpt = MakeCanvas("cpt", "y",800,800);
  TCanvas *ceta = MakeCanvas("ceta","0",800,800);
  TCanvas *cphi = MakeCanvas("cphi","0",800,800);
  TCanvas *cm = MakeCanvas("cm","0",800,800);
 
  //making variable bin size histogram
  const int nBins_pt = 8;
  double edges[nBins_pt + 1] = {5,10,15,20,25,30,40,60,100};

  vector<TH1D*> pts = PtBinCorrectly(dataFile, pyFile, geFile, nBins_pt, edges, nBins_pt, edges, "Pt", "initial");
  vector<TH1D*> etas = PtBinCorrectly(dataFile, pyFile, geFile, 50, -1, 1, 50, -1, 1, "Eta", "initial");
  vector<TH1D*> phis = PtBinCorrectly(dataFile, pyFile, geFile, 50, 0, 2*M_PI, 50, 0, 2*M_PI, "Phi", "initial");
  vector<TH1D*> ms = PtBinCorrectly(dataFile, pyFile, geFile, 20, 0, 10, 20, 0, 10, "M", "initial");

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
void Slices(TFile *matchFile, const string out, const string filetype) {
  TH2D * h_dummy = new TH2D("h_dummy","",1,0,1,1,0,1);
  //plots slices in pT of jet observable resolutions
  const int nBins = 7;
  double edges[nBins + 1] = {5,10,15,20,25,30,40,60};
  TH2D* ratM = PtBinCorrectly(matchFile, 51, 0, 2, nBins, edges, "ratioM", "pyPt");
  ObservablePtSlices (h_dummy, h_dummy, ratM, out, filetype, "M_{jet}^{det} / M_{jet}^{gen}", 1, 0);
  TH2D* ratPt = PtBinCorrectly(matchFile, 51, 0, 2, nBins, edges, "ratioPt", "pyPt");
  ObservablePtSlices (h_dummy, h_dummy, ratPt, out, filetype, "p_{T}^{det-jet} / p_{T}^{gen-jet}", 1, 0);
  TH2D* ratZg = PtBinCorrectly(matchFile, 51, 0, 2, nBins, edges, "ratioZg", "pyPtg");
  ObservablePtSlices (h_dummy, h_dummy, ratZg, out, filetype, "z_{g}^{det} / z_{g}^{gen}", 1, 1); //NOTE: sliced in bins of GROOMED pT
  TH2D* ratRg = PtBinCorrectly(matchFile, 51, 0, 2, nBins, edges, "ratioRg", "pyPtg");
  ObservablePtSlices (h_dummy, h_dummy, ratRg, out, filetype, "R_{g}^{det} / R_{g}^{gen}", 1, 1); //NOTE: sliced in bins of GROOMED pT
  return;
}

void Resolutions(TFile *matchFile, const string out, const string filetype) {
  //plots jet observable resolutions v. gen-level pT and associated mean & RMS for each pT bin
  TH2D* delM = PtBinCorrectly(matchFile, 11,5,60, 220,-1,1, "pyPt", "deltaM");
  Resolution(delM, out, filetype, "#Delta M_{jet} (Det - Gen) / M^{gen}_{jet}", 0);
  TH2D* delPt = PtBinCorrectly(matchFile, 11, 5, 60, 220, -1, 1, "pyPt", "deltaPt");
  Resolution(delPt, out, filetype, "#Delta p_{T}^{jet} (Det - Gen) / p_{T}^{gen-jet}", 0);
  TH2D* delZg = PtBinCorrectly(matchFile, 11, 5, 60, 220, -1, 1, "pyPt", "deltaZg");
  Resolution(delZg, out, filetype, "#Delta z_{g} (Det - Gen)", 1);
  TH2D* delRg = PtBinCorrectly(matchFile, 11, 5, 60, 220, -1, 1, "pyPt", "deltaRg");
  Resolution(delRg, out, filetype, "#Delta R_{g} (Det - Gen)", 1);
  return;
}



//Plots jet observables! M, pT, zg, Rg
void Observables(TFile *matchFile, TFile *pyFile, TFile *geFile, TFile *dataFile, const string out, const string filetype) {
  //plots the mass for various pT ranges
  cout << "AAA" << endl;
  vector<TH2D*> mpt = PtBinCorrectly(dataFile, pyFile, geFile, 20, 0, 10, 11, 5, 60, 20, 0, 10, 11, 5, 60, "M", "Pt");
  cout << "BBB" << endl;
  ObservablePtSlices(mpt[2], mpt[0], mpt[1], out, filetype, "M_{jet} [GeV/c^{2}]", 0, 0);
  //plots the ch fraction for various pT ranges
  cout << "CCC" << endl;
  vector<TH2D*> chfracpt = PtBinCorrectly(dataFile, pyFile, geFile, 10, 0, 1, 11, 5, 60, 10, 0, 1, 11, 5, 60, "ch_frac", "Pt");
  cout << "DDD" << endl;
  ObservablePtSlices(chfracpt[2], chfracpt[0], chfracpt[1], out, filetype, "ch. frac.", 0, 0);
  //plots the zg for various pTg ranges
  cout << "EEE" << endl;
  vector<TH2D*> zgptg = PtBinCorrectly(dataFile, pyFile, geFile, 10, 0, 0.5, 11, 5, 60, 10, 0, 0.5, 11, 5, 60, "zg", "ptg");
  cout << "FFF" << endl;
  ObservablePtSlices(zgptg[2], zgptg[0], zgptg[1], out, filetype, "z_{g}", 0, 0);
  //plots the rg for various pTg ranges
  cout << "GGG" << endl;
  vector<TH2D*> rgptg = PtBinCorrectly(dataFile, pyFile, geFile, 10, 0, 0.5, 11, 5, 60, 10, 0, 0.5, 11, 5, 60, "rg", "ptg");
  cout << "HHH" << endl;
  ObservablePtSlices(rgptg[2], rgptg[0], rgptg[1], out, filetype, "R_{g}", 0, 0);
  cout << "III" << endl;
  return;
}


void Responses(TFile* matchFile, TFile* pyFile, TFile* geFile, TFile* dataFile, const string out, const string filetype) {
  //  std::vector<TH1D*> ptfine = PtBinCorrectly(pyFile, geFile, 80, 0, 80, 80, 0, 80, "Pt");
  // Response(matchFile, "pt_response", ptfine[0], ptfine[1], "p_{T}^{jet} [GeV/c]", out, filetype);
  std::cout << "1" << std::endl;
  std::vector<TH1D*> ptcoarse = PtBinCorrectly(dataFile, pyFile, geFile, 15, 5, 80, 13, 15, 80, "Pt", "for_response");
  std::cout << "2" << std::endl;
  Response(matchFile, pyFile, geFile, "pt_res_coarse", ptcoarse[0], ptcoarse[1], "p_{T}^{j} [GeV/c]", out, filetype);
  std::cout << "3" << std::endl;
  std::vector<TH1D*> ms = PtBinCorrectly(dataFile, pyFile, geFile, 20, 0, 10, 20, 0, 10, "M", "for_response");
  std::cout << "4" << std::endl;
  Response(matchFile, pyFile, geFile, "m_response", ms[0], ms[1], "M^{j} [GeV/c^{2}]", out, filetype);
  std::cout << "5" << std::endl;
  std::vector<TH1D*> zgs = PtBinCorrectly(dataFile, pyFile, geFile, 20, 0, 1, 20, 0, 1, "zg", "for_response");
  Response(matchFile, pyFile, geFile, "zg_response", zgs[0], zgs[1], "z_{g}", out, filetype);
  std::vector<TH1D*> rgs = PtBinCorrectly(dataFile, pyFile, geFile, 20, 0, 1, 20, 0, 1, "rg", "for_response");
  Response(matchFile, pyFile, geFile, "rg_response", rgs[0], rgs[1], "R_{g}", out, filetype);
  
  return;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

void Y12InitialAnalysis () {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILES&STYLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~//    
  
  string dir = "~/jetmass_local/";
  string gein = "out/sim/ge/";
  string pyin = "out/sim/py/";
  string datain = "out/data/";
  string matchin = "out/matching/";
  string file = "full.root";
  string out = "~/jetmass_local/plots/Y12ResolutionsandResponses/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";
  
  TFile* pyFile = new TFile( (dir + pyin + file).c_str(), "READ");
  TFile* geFile = new TFile( (dir + gein + file).c_str(), "READ");
  TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  TFile* matchFile = new TFile( (dir + matchin + file).c_str(), "READ");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  cout << "AA" << endl;
  //plots pT, eta, phi, and mass
    JetQA(pyFile, geFile, dataFile, flag1, flag2, out, filetype);
  //plots the pT response matrix and the simulation pT spectra for fine-grained and coarse-grained binning
  cout << "BB" << endl;
  Responses(matchFile, pyFile, geFile, dataFile, out, filetype);
  cout << "CC" << endl;
  Observables(matchFile, pyFile, geFile, dataFile, out, filetype); //plots jet observables: M, pT, zg, Rg
  cout << "DD" << endl;
  Slices(matchFile, out, filetype); //plots the resolution of the observables in bins of pT
  cout << "EE" << endl;
  Resolutions(matchFile, out, filetype); //plots 2D observable resolution with pT and associated means & RMSs
  cout << "FF" << endl;
  return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
