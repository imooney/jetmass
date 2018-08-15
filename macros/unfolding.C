#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "Plots.h"

using namespace std;

//bins pT finely or coarsely for the subsequent call to "Response()"                                                                
vector<TH1D*> PtBinCorrectly(TFile *dataFile, TFile* pyFile, TFile* geFile) {
  vector<double> * pT = 0; vector<double> *g_pt = 0; vector<double> * p_pt = 0;
  double p_weight, g_weight;

  TTree* d = (TTree*) dataFile->Get("event");
  TTree* py = (TTree*) pyFile->Get("event");
  TTree* ge = (TTree*) geFile->Get("event");

  d->SetBranchAddress("Pt", &pT); double dummy_weight = 1.0;
  py->SetBranchAddress("Pt",&p_pt); ge->SetBranchAddress("Pt",&g_pt);
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);
  
  TH1D* mes_incl_pT = HistFromTree("mes_incl_pT", 9, 15, 60, d, pT, dummy_weight); 
  TH1D* py_ptcoarse = HistFromTree("py_ptcoarse", 15, 5, 80, py, p_pt, p_weight);
  TH1D* ge_ptcoarse = HistFromTree("ge_ptcoarse", 13, 15, 80, ge, g_pt, g_weight);
  
  vector<TH1D*> hists = {mes_incl_pT, py_ptcoarse, ge_ptcoarse};

  return hists;
}

//bins simulation observable in desired manner from tree                                                                                                                          
vector<TH1D*> PtBinCorrectly(TFile* dataFile, TFile* pyFile, TFile* geFile, const int nBins_py, const double lo_py, const double hi_py, const int nBins_ge, const double lo_ge, const double hi_ge, const std::string obs) {
  cout << "i" << endl;
  vector<double> *d_obs = 0; vector<double> *g_obs = 0; vector<double> *p_obs = 0;
  double dummy_d_weight = 1.0; double g_weight, p_weight;
  cout << "ii" << endl;
  TTree* py = (TTree*) pyFile->Get("event");
  TTree* ge = (TTree*) geFile->Get("event");
  TTree* d = (TTree*) dataFile->Get("event");

  py->SetBranchAddress(obs.c_str(),&p_obs); ge->SetBranchAddress(obs.c_str(),&g_obs); d->SetBranchAddress(obs.c_str(),&d_obs);
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);
  cout << "iii" << endl;
  TH1D * py_obs = HistFromTree(("py_" + obs).c_str(), nBins_py, lo_py, hi_py, py, p_obs, p_weight);
  cout << "iiia" << endl;
  TH1D * ge_obs = HistFromTree(("ge_" + obs).c_str(), nBins_ge, lo_ge, hi_ge, ge, g_obs, g_weight);
  cout << "iiib" << endl;
  TH1D * data_obs = HistFromTree(("data_" + obs).c_str(), nBins_ge, lo_ge, hi_ge, d, d_obs, dummy_d_weight);
  cout << "iv" << endl;
  
  vector<TH1D*> hists = {data_obs, py_obs, ge_obs};
  cout << "v" << endl;
  return hists;
}


//bins simulation observable in desired manner from tree                                                                                                                          
vector<TH2D*> PtBinCorrectly(TFile* dataFile, TFile* pyFile, TFile* geFile, const int nBins_py, const double lo_py, const double hi_py, const int nBins_d, const double lo_d, const double hi_d, const int nBins_x, const double lo_x, const double hi_x, const std::string obsx, const std::string obsy) {
  vector<double> *d_obsx = 0; vector<double> *g_obsx = 0; vector<double> *p_obsx = 0;
  vector<double> *d_obsy = 0; vector<double> *g_obsy = 0; vector<double> *p_obsy = 0; 
  double dummy_d_weight = 1.0; double g_weight, p_weight;

  TTree* py = (TTree*) pyFile->Get("event");
  TTree* ge = (TTree*) geFile->Get("event");
  TTree* d = (TTree*) dataFile->Get("event");

  py->SetBranchAddress(obsx.c_str(),&p_obsx); ge->SetBranchAddress(obsx.c_str(),&g_obsx); d->SetBranchAddress(obsx.c_str(),&d_obsx);
  py->SetBranchAddress(obsy.c_str(),&p_obsy); ge->SetBranchAddress(obsy.c_str(),&g_obsy); d->SetBranchAddress(obsy.c_str(),&d_obsy); 
  py->SetBranchAddress("weight", &p_weight); ge->SetBranchAddress("weight",&g_weight);

  TH2D * py_obs = HistFromTree(("py4_" + obsx).c_str(), nBins_x, lo_x, hi_x, nBins_py, lo_py, hi_py, py, p_obsx, p_obsy, p_weight);
  TH2D * ge_obs = HistFromTree(("ge4_" + obsx).c_str(), nBins_x, lo_x, hi_x, nBins_py, lo_py, hi_py, ge, g_obsx, g_obsy, g_weight);
  TH2D * data_obs = HistFromTree(("data4_" + obsx).c_str(), nBins_x, lo_x, hi_x, nBins_d, lo_d, hi_d, d, d_obsx, d_obsy, dummy_d_weight);

  vector<TH2D*> hists = {data_obs, py_obs, ge_obs};
  return hists;
}


void UnfoldedObs(TFile* matchFile, TFile* pyFile, TFile* geFile, TFile* dataFile, TH1D* raw, TH1D* gen, TH1D* det, const string resName, const string xTitle, const string log, const string out, const string filetype, const string flag) {
  string cName = (string) "cu" + (string) raw->GetName();
  TCanvas * cu = MakeCanvas(cName.c_str(),log,800,800);

  RooUnfoldResponse *res = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  
  /*0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests*/
  
  RooUnfoldBayes *unfolded_4iter = new RooUnfoldBayes(res, raw, 4, false, ((string) "unfolded_4iter" + (string) raw->GetName()).c_str(),"");
  TH1D * reco2 = (TH1D*) unfolded_4iter->Hreco((RooUnfold::ErrorTreatment) 1);
  double lowx = -2; double hix = -2; double lowy =-2; double highy = -2;
  if (log == "Y") {lowy = -1; highy = -1;} else {lowy = 0; highy = 0.5;}
  if (flag == "mass") {lowx = -1; hix = -1;} else {lowx = 0; hix = 0.5;}
  
  //  MakeCrossSection(pyFile, geFile, dataFile, gen, det, raw);
  
  Prettify1D(raw, kBlack, kFullStar, 2, kBlack, xTitle, /*("1/N_{j} dN_{j}/d" + xTitle.substr(0,5)).c_str()*/"prob.",lowx,hix,lowy, highy);
  Prettify1DwLineStyle(gen, kGreen, kDashed, 5, xTitle, /*("1/N_{j} dN_{j}/d" + xTitle.substr(0,5)).c_str()*/"prob.",lowx,hix,lowy, highy);
  Prettify1D(det, kBlue, kOpenCircle, 2, kBlue, xTitle, /*("1/N_{j} dN_{j}/d" + xTitle.substr(0,5)).c_str()*/"prob.",lowx,hix,lowy, highy);
  Prettify1D(reco2, kRed, kFullStar, 2, kRed, xTitle, /*("1/N_{j} dN_{j}/d" + xTitle.substr(0,5)).c_str()*/"prob.",lowx,hix,lowy, highy);
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
  
  double lox, hix;
  
  if (flag == "mass") {lox = -1; hix = -1;} else {lox = 0; hix = 0.5;}
  Prettify1D(rawXs[0], kBlack, kFullStar, 2, kBlack, xTitle, "prob.",lox, hix,0,0.5);
  Prettify1DwLineStyle(genXs[0], kGreen, kDashed, 5, xTitle, "prob.",lox, hix,0,0.5);
  Prettify1D(detXs[0], kBlue, kOpenCircle, 2, kBlue, xTitle, "prob.",lox, hix,0,0.5);
  Prettify1D(recoXs[0], kRed, kFullStar, 2, kRed, xTitle, "prob.",lox, hix,0,0.5);
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
  TCanvas *cslice = MakeCanvas(cName.c_str(), log, 800,600);
  DivideCanvas(cslice,"0",3,2);
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
  vector<TH1D*> detXs = Projection2D(det, nBins, ranges, "x");
  
  double lox, hix;
  if (flag == "mass") {lox = -1; hix = -1;} else {lox = 0; hix = 0.5;}
  for(int i = 0; i < nBins; ++ i) { 
    Prettify1D(rawXs[i], kBlack, kFullStar, 2, kBlack, xTitle, "prob.",lox, hix,0,0.5);
    Prettify1DwLineStyle(genXs[i], kGreen, kDashed, 5, xTitle, "prob.",lox, hix,0,0.5);
    Prettify1D(detXs[i], kBlue, kOpenCircle, 2, kBlue, xTitle, "prob.",lox, hix,0,0.5);
    Prettify1D(recoXs[i], kRed, kFullStar, 2, kRed, xTitle, "prob.",lox, hix,0,0.5);
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

void AllUnfolds(TFile *matchFile, TFile *pyFile, TFile *geFile, TFile *dataFile, const string out, const string filetype) {

  vector<TH1D*> pts = PtBinCorrectly(dataFile, pyFile, geFile);
  cout << "w" <<endl;
  UnfoldedObs(matchFile, pyFile, geFile, dataFile, pts[0], pts[1], pts[2], "pt_res_coarse", "p_{T}^{j} [GeV/c]", "Y", out, filetype, "mass");
  cout << "x" << endl;
  vector<TH1D*> ms = PtBinCorrectly(dataFile, pyFile, geFile, 20, 0, 10, 20, 0, 10, "M");
  cout << "y" << endl;
  UnfoldedObs(matchFile, pyFile, geFile, dataFile, ms[0], ms[1], ms[2], "m_response", "M^{j} [GeV/c^{2}]", "0", out, filetype, "mass");
  vector<TH1D*> zgs = PtBinCorrectly(dataFile, pyFile, geFile, 20, 0.001, 1.001, 20, 0.0/01, 1.001, "zg");
  UnfoldedObs(matchFile, pyFile, geFile, dataFile, zgs[0] , zgs[1], zgs[2], "zg_response", "z_{g}", "0", out, filetype, "");
  vector<TH1D*> rgs = PtBinCorrectly(dataFile, pyFile, geFile, 20, 0.001, 1.001, 20, 0.001, 1.001, "rg");
  UnfoldedObs(matchFile, pyFile, geFile, dataFile, rgs[0], rgs[1], rgs[2], "rg_response", "R_{g}", "0", out, filetype,"");
  
  return;
}
/*
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
*/
void Draw4DResponse(TFile* matchFile, const std::string resName, const std::string out, const std::string filetype) {
  TCanvas * c4res = MakeCanvas(("c4res_" + resName).c_str(), "z", 800,800);
  RooUnfoldResponse *res4D = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  TH2D* res_matrix = (TH2D*) res4D->Hresponse();
  Prettify2D(res_matrix, "p^{meas.}_{T} [GeV/c]", "p^{truth}_{T} [GeV/c]", -1,-1, -1,-1, 10e-14, 10e-1);
  
  //  res_matrix->GetXaxis()->SetBinLabel(180, "60");
  for (int i = 1; i < 181; i += 20) {
    res_matrix->GetXaxis()->SetBinLabel(i,(std::to_string( (int) ((i-1)/(double) 4) + 15)).c_str());
  } //0 -> 15, 20 -> 20, 40 -> 25, 60 -> 30, 80 -> 35, 100 -> 40, 120 -> 45, 140 -> 50, 160 ->55, 180 -> 60
  // 1 => (1 - 1) / x + 15 ... 21=> (21 - 1) / x + 15 = 20
  //res_matrix->GetYaxis()->SetNoAlphanumeric();
  for (int i = 1; i <= 281; i += 40) {
    //res_matrix->GetYaxis()->ChangeLabel(i+1,-1,-1,-1,-1,-1,(std::to_string((i)+5)).c_str());
    res_matrix->GetYaxis()->SetBinLabel(i, (std::to_string( (int) ((i-1)/(double) 4) + 5) ).c_str()); 
  } // 1 => (1 - 1) / 4 + 5... 21 => (21 - 1) / 4 + 5
  res_matrix->SetTitle("");
  
  res_matrix->Draw("colz");

  //  c4res->SaveAs((out+"4D_"+resName+filetype).c_str());

  return;
}

void unfolding () {
  string dir = "~/jetmass_local/";
  string matchin = "out/matching/";
  string pyin = "out/sim/py/";
  string gein = "out/sim/ge/";
  string datain = "out/data/";
  string file = "full.root";
  string out = "~/jetmass_local/plots/unfolding/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";

  TFile* matchFile = new TFile( (dir + matchin + file).c_str(), "READ");
  TFile* pyFile = new TFile( (dir + pyin +file).c_str(), "READ");
  TFile* geFile = new TFile( (dir + gein +file).c_str(), "READ");
  TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  
  //plots the raw and unfolded data for all observables (pT, m, zg, Rg) along with truth & reco
  /* AllUnfolds(matchFile, pyFile, geFile, dataFile, out, filetype);
  
  vector<TH2D*> m2d = PtBinCorrectly(dataFile, pyFile, geFile, 11, 5, 60, 9, 15, 60, 20, 0, 10, "M", "Pt");
  Unfold4D(matchFile, m2d[0], m2d[1], m2d[2], "0", out, filetype, "M_{jet} [GeV/c^{2}]", "pt_m_response", "mass");
  vector<TH2D*> zg2d = PtBinCorrectly(dataFile, pyFile, geFile, 11, 5, 60, 9, 15, 60, 20, 0.001, 1.001, "zg", "ptg"); 
  Unfold4D(matchFile, zg2d[0], zg2d[1], zg2d[2], "0", out, filetype, "z_{g}", "ptg_zg_response", "");
  vector<TH2D*> rg2d = PtBinCorrectly(dataFile, pyFile, geFile, 11, 5, 60, 9, 15, 60, 20, 0.001, 1.001, "rg", "ptg");
  Unfold4D(matchFile, rg2d[0], rg2d[1], rg2d[2], "0", out, filetype, "R_{g}", "ptg_rg_response", "");
  
  SliceUnfolded4D(matchFile, m2d[0], m2d[1], m2d[2], "0", out, filetype, "M_{jet} [GeV/c^{2}]", "pt_m_response", "mass");
  SliceUnfolded4D(matchFile, zg2d[0], zg2d[1], zg2d[2], "0", out, filetype, "z_{g}", "ptg_zg_response", "");
  SliceUnfolded4D(matchFile, rg2d[0], rg2d[1], rg2d[2], "0", out, filetype, "R_{g}", "ptg_rg_response", "");
*/  
  Draw4DResponse(matchFile, "pt_m_response", out, filetype);
  Draw4DResponse(matchFile, "ptg_zg_response", out, filetype);
  Draw4DResponse(matchFile, "ptg_rg_response", out, filetype);
  
  return;
}
