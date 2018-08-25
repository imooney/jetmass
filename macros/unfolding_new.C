#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "Plots_new.h"

using namespace std;

void UnfoldedObs(TFile* matchFile, TH1D* raw, TH1D* gen, TH1D* det, const string resName, const string xTitle, const string log, const string out, const string filetype, const string flag) {
  string cName = (string) "cu" + (string) raw->GetName();
  TCanvas * cu = MakeCanvas(cName.c_str(),log,800,800);

  RooUnfoldResponse *res = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  
  /*0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests*/
  
  RooUnfoldBayes *unfolded_4iter = new RooUnfoldBayes(res, raw, 4, false, ((string) "unfolded_4iter" + (string) raw->GetName()).c_str(),"");
  TH1D * reco2 = (TH1D*) unfolded_4iter->Hreco((RooUnfold::ErrorTreatment) 1);
  double lox = -2; double hix = -2; double loy =-2; double hiy = -2;
  if (log == "Y") {loy = -1; hiy = -1;} else {loy = 0; hiy = 0.5;}
  if (flag == "mass") {lox = -1; hix = -1;} else {lox = 0; hix = 0.5;}
  
  //  MakeCrossSection(pyFile, geFile, dataFile, gen, det, raw);
  
  Prettify1D(raw, kBlack, kFullStar, 2, kBlack, xTitle, /*("1/N_{j} dN_{j}/d" + xTitle.substr(0,5)).c_str()*/"prob.",lox,hix,loy, hiy);
  Prettify1DwLineStyle(gen, kGreen, kDashed, 5, xTitle, /*("1/N_{j} dN_{j}/d" + xTitle.substr(0,5)).c_str()*/"prob.",lox,hix,loy, hiy);
  Prettify1D(det, kBlue, kOpenCircle, 2, kBlue, xTitle, /*("1/N_{j} dN_{j}/d" + xTitle.substr(0,5)).c_str()*/"prob.",lox,hix,loy, hiy);
  Prettify1D(reco2, kRed, kFullStar, 2, kRed, xTitle, /*("1/N_{j} dN_{j}/d" + xTitle.substr(0,5)).c_str()*/"prob.",lox,hix,loy, hiy);
  reco2->SetTitle("");
  
  TLegend * tu = TitleLegend(0.44,0.57,0.84,0.87);
  tu->AddEntry(gen,"PYTHIA6","l");
  tu->AddEntry(det,"PYTHIA6+GEANT","p");
  tu->AddEntry(raw,"Raw data", "p");
  tu->AddEntry(reco2,"Unfolded data (4 iter)","p");
  
  reco2->Draw(); raw->Draw("same"); gen->Draw("C,same"); det->Draw("same"); tu->Draw("same");
  
  cu->SaveAs((out + "unfolded_" + raw->GetName() + filetype).c_str());
 
  return;
}

void Unfold4D(TFile *matchFile, TH2D* raw, TH2D* gen, TH2D* det, const string log, const string out, const string filetype, const string xTitle, const string resName, const string flag) {
  string cName = (string) "c4" + (string) raw->GetName();
  TCanvas *c4 = MakeCanvas(cName, log, 800,800);
  
  RooUnfoldResponse *res4D = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  cout << resName.c_str() << endl; cout << raw->GetName() << endl; cout << res4D->GetName() << endl;
  RooUnfoldBayes *unfold4D_4iter = new RooUnfoldBayes(res4D, raw, 4, false, ((string) "unfold4D_4iter" + (string) raw->GetName()).c_str(),"");
  TH2D *reco = (TH2D*) unfold4D_4iter->Hreco((RooUnfold::ErrorTreatment) 1);
  
  const unsigned nBins = 1;
  double ranges[nBins+1] = {0,999};
  vector<TH1D*> recoXs = Projection2D(reco, nBins, ranges, "x"); 
  vector<TH1D*> rawXs = Projection2D(raw, nBins, ranges, "x"); 
  vector<TH1D*> genXs = Projection2D(gen, nBins, ranges, "x");
  vector<TH1D*> detXs = Projection2D(det, nBins, ranges, "x");
  
  double lox, hix, loy, hiy;
  
  if (log == "Y") {loy = -1; hiy = -1;} else {loy = 0; hiy = 0.5;}
  if (flag == "mass") {lox = -1; hix = -1;} else {lox = 0; hix = 0.5;}
  
  Prettify1D(rawXs[0], kBlack, kFullStar, 2, kBlack, xTitle, "prob.",lox, hix,loy,hiy);
  Prettify1DwLineStyle(genXs[0], kGreen, kDashed, 5, xTitle, "prob.",lox, hix,loy,hiy);
  Prettify1D(detXs[0], kBlue, kOpenCircle, 2, kBlue, xTitle, "prob.",lox, hix,loy,hiy);
  Prettify1D(recoXs[0], kRed, kFullStar, 2, kRed, xTitle, "prob.",lox, hix,loy,hiy);
  recoXs[0]->SetTitle("");

  TLegend * tu = TitleLegend(0.44,0.57,0.84,0.87);
  tu->AddEntry(genXs[0],"PYTHIA6","l");
  tu->AddEntry(detXs[0],"PYTHIA6+GEANT","p");
  tu->AddEntry(rawXs[0],"Raw data", "p");
  tu->AddEntry(recoXs[0],"Unfolded data (4 iter)","p");
  
  recoXs[0]->Draw(); rawXs[0]->Draw("same"); genXs[0]->Draw("C,same"); detXs[0]->Draw("same"); tu->Draw("same");

  c4->SaveAs((out + "4Dunfolded_" + raw->GetName() + filetype).c_str());
  
  return;
}

void SliceUnfolded4D(TFile *matchFile, TH2D* raw, TH2D* gen, TH2D* det, const string log, const string out, const string filetype, const string xTitle, const string resName, const string flag) {
  string cName = (string) "cslice" + (string) raw->GetName();
  TCanvas *cslice = MakeCanvas(cName.c_str(), "0", 800,600);
  DivideCanvas(cslice,log,3,2);
  RooUnfoldResponse *res4D = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  
  RooUnfoldBayes *unfold4D_4iter = new RooUnfoldBayes(res4D, raw, 4, false, ((string) "unfold4D_4iter" + (string) raw->GetName()).c_str(),"");
  TH2D *reco = (TH2D*) unfold4D_4iter->Hreco((RooUnfold::ErrorTreatment) 1);
  
  const unsigned nBins = 5;
  double ranges[nBins+1] = {2,3,4,5,7,11};
  double corresp_pts[nBins+1] = {15,20,25,30,40,60};
  double ranges_d[nBins+1] = {0,1,2,3,5,9};
  vector<TH1D*> recoXs = Projection2D(reco, nBins, ranges, "x"); 
  vector<TH1D*> rawXs = Projection2D(raw, nBins, ranges_d, "x");
  vector<TH1D*> genXs = Projection2D(gen, nBins, ranges, "x");
  vector<TH1D*> detXs = Projection2D(det, nBins, ranges_d, "x");
  
  double lox, hix, loy, hiy;
  if (log == "Y") {loy = -1; hiy = -1;} else {loy = 0; hiy = 0.5;}
  if (flag == "mass") {lox = -1; hix = -1;} else {lox = 0; hix = 0.5;}

  for(int i = 0; i < nBins; ++ i) { 
    Prettify1D(rawXs[i], kBlack, kFullStar, 2, kBlack, xTitle, "prob.",lox, hix,loy,hiy);
    Prettify1DwLineStyle(genXs[i], kGreen, kDashed, 5, xTitle, "prob.",lox, hix,loy,hiy);
    Prettify1D(detXs[i], kBlue, kOpenCircle, 2, kBlue, xTitle, "prob.",lox, hix,loy,hiy);
    Prettify1D(recoXs[i], kRed, kFullStar, 2, kRed, xTitle, "prob.",lox, hix,loy,hiy);
    recoXs[i]->SetTitle("");
  }
  TLegend * tu = new TLegend(0.1,0.15,0.8,0.45); tu->SetBorderSize(0);
  tu->AddEntry(genXs[0],"PYTHIA6","l");
  tu->AddEntry(detXs[0],"PYTHIA6+GEANT","p");
  tu->AddEntry(rawXs[0],"Raw data", "p");
  tu->AddEntry(recoXs[0],"Unfolded data (4 iter)","p");
  
  TLegend *tslices[nBins];
  for (int i = 0; i < nBins; ++ i) {
    tslices[i] = SliceLegend(((to_string(corresp_pts[i])).substr(0,2) + " < p_{T}^{unfolded jet} < " + (to_string(corresp_pts[i + 1])).substr(0,2) + " GeV/c").c_str(), 0.13,0.8,0.9,0.95);
  }  
  cslice->cd(1); TLatex* title = PanelTitle(); tu->Draw();
  for (int i = 0; i < nBins; ++ i) {
    cslice->cd(i+2); recoXs[i]->Draw(); rawXs[i]->Draw("same"); genXs[i]->Draw("C,same"); detXs[i]->Draw("same"); tslices[i]->Draw("same");
  }
  
  cslice->SaveAs((out + "4Dunfolded_slices_" + raw->GetName() + filetype).c_str());
  
  return;
}

void AllUnfolds(TFile *matchFile, TFile *inFile, const string out, const string filetype) {
    
  vector<TH1D*> pts = {(TH1D*) inFile->Get("pt_d"), (TH1D*) inFile->Get("pt_p"), (TH1D*) inFile->Get("pt_g")};
  vector<TH1D*> ms = {(TH1D*) inFile->Get("m_d"), (TH1D*) inFile->Get("m_p"), (TH1D*) inFile->Get("m_g")};
  vector<TH1D*> zgs = {(TH1D*) inFile->Get("zg_d"), (TH1D*) inFile->Get("zg_p"), (TH1D*) inFile->Get("zg_g")};
  vector<TH1D*> rgs = {(TH1D*) inFile->Get("rg_d"), (TH1D*) inFile->Get("rg_p"), (TH1D*) inFile->Get("rg_g")};
  vector<TH1D*> ptgs = {(TH1D*) inFile->Get("ptg_d"), (TH1D*) inFile->Get("ptg_p"), (TH1D*) inFile->Get("ptg_g")};
  vector<TH1D*> mgs = {(TH1D*) inFile->Get("mg_d"), (TH1D*) inFile->Get("mg_p"), (TH1D*) inFile->Get("mg_g")};
  
  UnfoldedObs(matchFile, pts[0], pts[1], pts[2], "pt_res_coarse", "p_{T}^{j} [GeV/c]", "Y", out, filetype, "mass");
  UnfoldedObs(matchFile, ms[0], ms[1], ms[2], "m_response", "M^{j} [GeV/c^{2}]", "0", out, filetype, "mass");
  UnfoldedObs(matchFile, zgs[0] , zgs[1], zgs[2], "zg_response", "z_{g}", "0", out, filetype, "");
  UnfoldedObs(matchFile, rgs[0], rgs[1], rgs[2], "rg_response", "R_{g}", "0", out, filetype,"");
  UnfoldedObs(matchFile, ptgs[0], ptgs[1], ptgs[2], "ptg_response", "p_{T,g}", "Y", out, filetype,"");
  UnfoldedObs(matchFile, mgs[0], mgs[1], mgs[2], "mg_response", "M_{g}", "0", out, filetype,"mass");
  
  return;
}

void Draw4DResponse(TFile* matchFile, const std::string resName, const std::string out, const std::string filetype, const int nBinsx, const int nBinsy) {
  TCanvas * c4res = MakeCanvas(("c4res_" + resName).c_str(), "z", 800,800);
  RooUnfoldResponse *res4D = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  TH2D* res_matrix = (TH2D*) res4D->Hresponse();
  Prettify2D(res_matrix, "p^{meas.}_{T} [GeV/c]", "p^{truth}_{T} [GeV/c]", -1,-1, -1,-1, 10e-11, 10e-3);
  
  //  res_matrix->GetXaxis()->SetBinLabel(180, "60");
  for (int i = 1; i < (9*nBinsx + 1); i += nBinsx) {
    if (nBinsx == 20) {res_matrix->GetXaxis()->SetBinLabel(i,(std::to_string( (int) ((i-1)/(double) 4) + 15)).c_str()); }
    if (nBinsx == 9) {res_matrix->GetXaxis()->SetBinLabel(i,(std::to_string( (int) (5*(i-1)/(double) 9) + 15)).c_str()); }
  } //0 -> 15, 20 -> 20, 40 -> 25, 60 -> 30, 80 -> 35, 100 -> 40, 120 -> 45, 140 -> 50, 160 ->55, 180 -> 60
  // 1 => (1 - 1) / x + 15 ... 21=> (21 - 1) / x + 15 = 20
  //res_matrix->GetYaxis()->SetNoAlphanumeric();
  for (int i = 1; i <= (14*nBinsy + 1); i += 2*nBinsy) {
    //res_matrix->GetYaxis()->ChangeLabel(i+1,-1,-1,-1,-1,-1,(std::to_string((i)+5)).c_str());
    if (nBinsy == 20) {res_matrix->GetYaxis()->SetBinLabel(i, (std::to_string( (int) ((i-1)/(double) 4) + 5) ).c_str());} 
    if (nBinsy == 15) {res_matrix->GetYaxis()->SetBinLabel(i,(std::to_string( (int) ((i-1)/(double) 3) + 5) ).c_str());}
  } // 1 => (1 - 1) / 4 + 5... 21 => (21 - 1) / 4 + 5
  res_matrix->SetTitle("");
  
  res_matrix->Draw("colz");

  //  c4res->SaveAs((out+"4D_"+resName+filetype).c_str());

  return;
}

void unfolding_new () {
  string dir = "~/jetmass/";
  string matchin = "out/matching/";
  string pyin = "out/sim/py/";
  string gein = "out/sim/ge/";
  string datain = "out/data/";
  string in = "macros/hists/";
  string infile = "hists.root";
  string file = "full.root";
  string out = "~/jetmass/plots/unfolding/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";

  TFile* matchFile = new TFile( (dir + matchin + file).c_str(), "READ");
  //TFile* pyFile = new TFile( (dir + pyin +file).c_str(), "READ");
  //TFile* geFile = new TFile( (dir + gein +file).c_str(), "READ");
  //TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  TFile* inFile = new TFile( (dir + in + infile).c_str(), "READ");
    
  //plots the raw and unfolded data for all observables (pT, m, zg, Rg) along with truth & reco
  AllUnfolds(matchFile, inFile, out, filetype);
  
  vector<TH2D*> m2d = {(TH2D*) inFile->Get("m_v_pt_d"), (TH2D*) inFile->Get("m_v_pt_p"), (TH2D*) inFile->Get("m_v_pt_g")};
  vector<TH2D*> zg2d = {(TH2D*) inFile->Get("zg_v_pt_d"), (TH2D*) inFile->Get("zg_v_pt_p"), (TH2D*) inFile->Get("zg_v_pt_g")};
  vector<TH2D*> rg2d = {(TH2D*) inFile->Get("rg_v_pt_d"), (TH2D*) inFile->Get("rg_v_pt_p"), (TH2D*) inFile->Get("rg_v_pt_g")};
  vector<TH2D*> ptg2d = {(TH2D*) inFile->Get("ptg_v_pt_d"), (TH2D*) inFile->Get("ptg_v_pt_p"), (TH2D*) inFile->Get("ptg_v_pt_g")};
  vector<TH2D*> mg2d = {(TH2D*) inFile->Get("mg_v_pt_d"), (TH2D*) inFile->Get("mg_v_pt_p"), (TH2D*) inFile->Get("mg_v_pt_g")};
  
  
  Unfold4D(matchFile, m2d[0], m2d[1], m2d[2], "0", out, filetype, "M_{jet} [GeV/c^{2}]", "pt_m_response", "mass");
  Unfold4D(matchFile, zg2d[0], zg2d[1], zg2d[2], "0", out, filetype, "z_{g}", "pt_zg_response", "");
  Unfold4D(matchFile, rg2d[0], rg2d[1], rg2d[2], "0", out, filetype, "R_{g}", "pt_rg_response", "");
  Unfold4D(matchFile, ptg2d[0], ptg2d[1], ptg2d[2], "Y", out, filetype, "p_{T,g}", "pt_ptg_response", "");
  Unfold4D(matchFile, mg2d[0], mg2d[1], mg2d[2], "0", out, filetype, "M_{g}", "pt_mg_response", "mass");
  
  SliceUnfolded4D(matchFile, m2d[0], m2d[1], m2d[2], "0", out, filetype, "M_{jet} [GeV/c^{2}]", "pt_m_response", "mass");
  SliceUnfolded4D(matchFile, zg2d[0], zg2d[1], zg2d[2], "0", out, filetype, "z_{g}", "pt_zg_response", "");
  SliceUnfolded4D(matchFile, rg2d[0], rg2d[1], rg2d[2], "0", out, filetype, "R_{g}", "pt_rg_response", "");
  SliceUnfolded4D(matchFile, ptg2d[0], ptg2d[1], ptg2d[2], "Y", out, filetype, "p_{T,g}", "pt_ptg_response", "");
  SliceUnfolded4D(matchFile, mg2d[0], mg2d[1], mg2d[2], "0", out, filetype, "R_{g}", "pt_mg_response", "mass");
  /*
  Draw4DResponse(matchFile, "pt_m_response", out, filetype, 20, 20);
  Draw4DResponse(matchFile, "pt_zg_response", out, filetype, 20, 20);
  Draw4DResponse(matchFile, "pt_rg_response", out, filetype, 20, 20);
  Draw4DResponse(matchFile, "pt_ptg_response", out, filetype, 9, 15);
  Draw4DResponse(matchFile, "pt_mg_response", out, filetype, 20, 20);
  */
  
  return;
}
