#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "Plots.h"

using namespace std;

//bins pT finely or coarsely for the subsequent call to "Response()"                                                                                   
vector<TH1D*> PtBinCorrectly(TFile *dataFile, TFile* simFile) {
  double pT, g_pt, g_weight, p_pt, p_weight;

  TTree* d = (TTree*) dataFile->Get("incl");
  TTree* py = (TTree*) simFile->Get("py_inclTree");
  TTree* ge = (TTree*) simFile->Get("ge_inclTree");

  d->SetBranchAddress("Pt", &pT); double dummy_weight = 1.0;
  py->SetBranchAddress("Pt",&p_pt); ge->SetBranchAddress("Pt",&g_pt);
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);
  
  TH1D* mes_incl_pT = HistFromTree("mes_incl_pT", 9, 15, 60, d, pT, dummy_weight); 
  TH1D* py_ptcoarse = HistFromTree("py_ptcoarse", 15, 5, 80, py, p_pt, p_weight);
  TH1D* ge_ptcoarse = HistFromTree("ge_ptcoarse", 13, 15, 80, ge, g_pt, g_weight);
  
  vector<TH1D*> hists = {mes_incl_pT, py_ptcoarse, ge_ptcoarse};

  return hists;
}

void UnfoldedObs(TFile* matchFile, TH1D* raw, TH1D* gen, TH1D* det, const string resName, const string xTitle, const string log, const string out, const string filetype) {
  string cName = (string) "cu" + (string) raw->GetName();
  TCanvas * cu = MakeCanvas(cName.c_str(),log,800,800);

  RooUnfoldResponse *res = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  
  /*0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests*/
  
  RooUnfoldBayes *unfolded_2iter = new RooUnfoldBayes(res, raw, 2, false, ((string) "unfolded_2iter" + (string) raw->GetName()).c_str(),"");
  TH1D * reco2 = (TH1D*) unfolded_2iter->Hreco((RooUnfold::ErrorTreatment) 1);
  double lowy =-2; double highy = -2;
  if (log == "Y") {lowy = -1; highy = -1;} else {lowy = 0; highy = 0.5;}
  Prettify1D(raw, kBlack, kFullStar, 2, kBlack, xTitle, "arb.",-1,-1,lowy, highy);
  Prettify1DwLineStyle(gen, kGreen, kDashed, 5, xTitle, "arb.",-1,-1,lowy, highy);
  Prettify1D(det, kBlue, kOpenCircle, 2, kBlue, xTitle, "arb.",-1,-1,lowy, highy);
  Prettify1D(reco2, kRed, kFullStar, 2, kRed, xTitle, "arb.",-1,-1,lowy, highy);
  reco2->SetTitle("");
  
  TLegend * tu = TitleLegend(0.44,0.57,0.84,0.87);
  tu->AddEntry(gen,"PYTHIA6","l");
  tu->AddEntry(det,"PYTHIA6+GEANT","p");
  tu->AddEntry(raw,"Raw data", "p");
  tu->AddEntry(reco2,"Unfolded data (2 iter)","p");
  
  reco2->Draw(); raw->Draw("same"); gen->Draw("C,same"); det->Draw("same"); tu->Draw("same");
  
  cu->SaveAs((out + "unfolded_" + raw->GetName() + filetype).c_str());
 
  return;
}

void Unfold4D(TFile *matchFile, TH2D* raw, TH2D* gen, TH2D* det, const string log, const string out, const string filetype, const string xTitle, const string resName) {
  string cName = (string) "c4" + (string) raw->GetName();
  TCanvas *c4 = MakeCanvas(cName, log, 800,800);
  
  RooUnfoldResponse *res4D = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  cout << resName.c_str() << endl; cout << raw->GetName() << endl; cout << res4D->GetName() << endl;
  RooUnfoldBayes *unfold4D_2iter = new RooUnfoldBayes(res4D, raw, 2, false, ((string) "unfold4D_2iter" + (string) raw->GetName()).c_str(),"");
  TH2D *reco = (TH2D*) unfold4D_2iter->Hreco((RooUnfold::ErrorTreatment) 1);
  
  const unsigned nBins = 1;
  double ranges[nBins+1] = {0,999};
  vector<TH1D*> recoXs = Projection2D(reco, nBins, ranges, "x"); 
  vector<TH1D*> rawXs = Projection2D(raw, nBins, ranges, "x"); 
  vector<TH1D*> genXs = Projection2D(gen, nBins, ranges, "x");
  vector<TH1D*> detXs = Projection2D(det, nBins, ranges, "x");

  Prettify1D(rawXs[0], kBlack, kFullStar, 2, kBlack, xTitle, "arb.",-1,-1,0,0.5);
  Prettify1DwLineStyle(genXs[0], kGreen, kDashed, 5, xTitle, "arb.",-1,-1,0,0.5);
  Prettify1D(detXs[0], kBlue, kOpenCircle, 2, kBlue, xTitle, "arb.",-1,-1,0,0.5);
  Prettify1D(recoXs[0], kRed, kFullStar, 2, kRed, xTitle, "arb.",-1,-1,0,0.5);
  recoXs[0]->SetTitle("");

  TLegend * tu = TitleLegend(0.44,0.57,0.84,0.87);
  tu->AddEntry(genXs[0],"PYTHIA6","l");
  tu->AddEntry(detXs[0],"PYTHIA6+GEANT","p");
  tu->AddEntry(rawXs[0],"Raw data", "p");
  tu->AddEntry(recoXs[0],"Unfolded data (2 iter)","p");
  
  recoXs[0]->Draw(); rawXs[0]->Draw("same"); genXs[0]->Draw("C,same"); detXs[0]->Draw("same"); tu->Draw("same");

  c4->SaveAs((out + "4Dunfolded_" + raw->GetName() + filetype).c_str());
  
  return;
}

void SliceUnfolded4D(TFile *matchFile, TH2D* raw, TH2D* gen, TH2D* det, const string log, const string out, const string filetype, const string xTitle, const string resName) {
  string cName = (string) "cslice" + (string) raw->GetName();
  TCanvas *cslice = MakeCanvas(cName.c_str(), log, 800,600);
  DivideCanvas(cslice,"0",3,2);
  RooUnfoldResponse *res4D = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  
  RooUnfoldBayes *unfold4D_2iter = new RooUnfoldBayes(res4D, raw, 2, false, ((string) "unfold4D_2iter" + (string) raw->GetName()).c_str(),"");
  TH2D *reco = (TH2D*) unfold4D_2iter->Hreco((RooUnfold::ErrorTreatment) 1);
  
  const unsigned nBins = 5;
  double ranges[nBins+1] = {2,3,4,5,7,11};
  double corresp_pts[nBins+1] = {15,20,25,30,40,60};
  double ranges_d[nBins+1] = {0,1,2,3,5,9};
  vector<TH1D*> recoXs = Projection2D(reco, nBins, ranges, "x"); 
  vector<TH1D*> rawXs = Projection2D(raw, nBins, ranges_d, "x");
  vector<TH1D*> genXs = Projection2D(gen, nBins, ranges, "x");
  vector<TH1D*> detXs = Projection2D(det, nBins, ranges, "x");
  
  for(int i = 0; i < nBins; ++ i) { 
    Prettify1D(rawXs[i], kBlack, kFullStar, 2, kBlack, xTitle, "arb.",-1,-1,0,0.5);
    Prettify1DwLineStyle(genXs[i], kGreen, kDashed, 5, xTitle, "arb.",-1,-1,0,0.5);
    Prettify1D(detXs[i], kBlue, kOpenCircle, 2, kBlue, xTitle, "arb.",-1,-1,0,0.5);
    Prettify1D(recoXs[i], kRed, kFullStar, 2, kRed, xTitle, "arb.",-1,-1,0,0.5);
    recoXs[i]->SetTitle("");
  }
  TLegend * tu = new TLegend(0.1,0.15,0.8,0.45); tu->SetBorderSize(0);
  tu->AddEntry(genXs[0],"PYTHIA6","l");
  tu->AddEntry(detXs[0],"PYTHIA6+GEANT","p");
  tu->AddEntry(rawXs[0],"Raw data", "p");
  tu->AddEntry(recoXs[0],"Unfolded data (2 iter)","p");
  
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

void AllUnfolds(TFile *matchFile, TFile *simFile, TFile *dataFile, const string out, const string filetype) {
  vector<TH1D*> pts = PtBinCorrectly(dataFile, simFile);
  THnSparseD* dataSD = (THnSparseD*) dataFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd");
  THnSparseD* genSD = (THnSparseD*) simFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd_py");
  THnSparseD* recoSD = (THnSparseD*) simFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd_ge");
 
  UnfoldedObs(matchFile, pts[0], pts[1], pts[2], "pt_res_coarse", "p^{jet}_{T} [GeV/c]", "Y", out, filetype);
  UnfoldedObs(matchFile, (TH1D*) dataFile->Get("m_full_incl_jet"), (TH1D*) simFile->Get("m_full_incl_jet_py"), (TH1D*) simFile->Get("m_full_incl_jet_ge"), "m_response", "M^{jet} [GeV/c^{2}]", "0", out, filetype);
  UnfoldedObs(matchFile, (TH1D*) dataSD->Projection(0), (TH1D*) genSD->Projection(0), (TH1D*) recoSD->Projection(0), "zg_response", "z_{g}", "0", out, filetype);
  UnfoldedObs(matchFile, (TH1D*) dataSD->Projection(2), (TH1D*) genSD->Projection(2), (TH1D*) recoSD->Projection(2), "rg_response", "R_{g}", "0", out, filetype);
 
  return;
}

//cuts off the low pT range of a 2D with another observable
void rebin2(TH1 *h, Int_t ngx, Int_t ngy)
{  
  //make a clone of h
  TH1 *hold = (TH1*)h->Clone();
  hold->SetDirectory(0);

  Int_t  nbinsx = hold->GetXaxis()->GetNbins();
  Int_t  nbinsy = hold->GetYaxis()->GetNbins();
  Float_t xmin  = hold->GetXaxis()->GetXmin();
  Float_t xmax  = hold->GetXaxis()->GetXmax();
  Float_t ymin  = hold->GetYaxis()->GetXmin() + 10;
  Float_t ymax  = hold->GetYaxis()->GetXmax();
  Int_t nx = 20;
  Int_t ny = 9;
  Double_t cu;
  Float_t bx,by;
  Int_t ix,iy,ibin,bin,binx,biny;
  h->SetBins (nx,xmin,xmax,ny,ymin,ymax);
  //loop on all bins to reset contents and errors
  for (biny=1;biny<=nbinsy;biny++) {
    for (binx=1;binx<=nbinsx;binx++) {
      ibin = h->GetBin(binx,biny);
      h->SetBinContent(ibin,0);
    }
  }
  //loop on all bins and refill
  for (biny=1;biny<=nbinsy;biny++) {
    by  = hold->GetYaxis()->GetBinCenter(biny);
    iy  = h->GetYaxis()->FindBin(by);
    for (binx=1;binx<=nbinsx;binx++) {
      bx = hold->GetXaxis()->GetBinCenter(binx);
      ix  = h->GetXaxis()->FindBin(bx);
      bin = hold->GetBin(binx,biny);
      ibin= h->GetBin(ix,iy);
      cu  = hold->GetBinContent(bin);
      h->AddBinContent(ibin,cu);
      h->SetBinError(ibin, (double) sqrt(h->GetBinContent(ibin)));
    }
  }
  h->SetEntries(hold->GetEntries());
  h->Sumw2();
  delete hold;          
}

void unfolding () {
  string dir = "~/jetmass_temp/";
  string matchin = "out/matching/";
  string simin = "out/sim/";
  string datain = "out/data/";
  string file = "full.root";
  string out = "~/jetmass_temp/plots/unfolding/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";
  string flag3 = "jet";

  TFile* matchFile = new TFile( (dir + matchin + file).c_str(), "READ");
  TFile* simFile = new TFile( (dir + simin +file).c_str(), "READ");
  TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  
  TH2D* m_pt_d = (TH2D*) dataFile->Get("m_v_pt_rebin_full_incl_jet");
  TH2D* m_pt_py = (TH2D*) simFile->Get("m_v_pt_full_incl_jet_py");
  TH2D* m_pt_ge = (TH2D*) simFile->Get("m_v_pt_full_incl_jet_ge");
  
  THnSparseD* dat = (THnSparseD*) dataFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd");
  THnSparseD* py = (THnSparseD*) simFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd_py");
  THnSparseD* ge = (THnSparseD*) simFile->Get("zg_mg_thetag_ptg_pt_full_incl_sd_ge");
  
  //plots the raw and unfolded data for all observables (pT, m, zg, Rg) along with truth & reco
  AllUnfolds(matchFile, simFile, dataFile, out, filetype);

  Unfold4D(matchFile, m_pt_d, m_pt_py, m_pt_ge, "0", out, filetype, "M_{jet} [GeV/c^{2}]", "pt_m_response");
  TH2D* rebinned_data_zg = (TH2D*) dat->Projection(3,0);
  TH2D* rebinned_data_rg = (TH2D*) dat->Projection(3,2);
  TH2D* py_zg = (TH2D*) py->Projection(3,0); TH2D* ge_zg = (TH2D*) ge->Projection(3,0);
  TH2D* py_rg = (TH2D*) py->Projection(3,2); TH2D* ge_rg = (TH2D*) ge->Projection(3,2);
  rebin2(rebinned_data_zg, 1, 1); rebin2(rebinned_data_rg, 1, 1);

  Unfold4D(matchFile, rebinned_data_zg, py_zg, ge_zg, "0", out, filetype, "z_{g}", "ptg_zg_response");
  Unfold4D(matchFile, rebinned_data_rg, py_rg, ge_rg, "0", out, filetype, "R_{g}", "ptg_rg_response");
  SliceUnfolded4D(matchFile, m_pt_d, m_pt_py, m_pt_ge, "0", out, filetype, "M_{jet} [GeV/c^{2}]", "pt_m_response");
  SliceUnfolded4D(matchFile, rebinned_data_zg, py_zg, ge_zg, "0", out, filetype, "z_{g}", "ptg_zg_response");
  SliceUnfolded4D(matchFile, rebinned_data_rg, py_rg, ge_rg, "0", out, filetype, "R_{g}", "ptg_rg_response");

  return;
}
