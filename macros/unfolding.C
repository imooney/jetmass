#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "Plots.h"

using namespace std;

void UnfoldedObs(TFile* matchFile, TH1D* raw, TH1D* gen, TH1D* det/*, TH1D* raw15*/, const string resName, const string xTitle, const string log, const string out, const string filetype, const string flag) {
  string cName = (string) "cu" + (string) raw->GetName();
  TCanvas * cu = MakeCanvas(cName.c_str(),log,800,800);

  RooUnfoldResponse *res = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  
  /*0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests*/
  
  std::string yTitle = "";
  /*
  if (flag == "pt") {
    yTitle = "arb.";
    raw->Scale(1/(double)raw->Integral()); det->Scale(1/(double)det->Integral());
    double rat_scale = raw->GetBinContent(3) / (double) gen->GetBinContent(3);
    gen->Scale(rat_scale); 
  } 
  else {*/yTitle = ("1/N_{j} dN_{j}/d" + xTitle.substr(0,xTitle.find(' '))).c_str();//}

  RooUnfoldBayes *unfolded_4iter = new RooUnfoldBayes(res, raw, 4, false, ((string) "unfolded_4iter" + (string) raw->GetName()).c_str(),"");
  TH1D * reco2 = (TH1D*) unfolded_4iter->Hreco((RooUnfold::ErrorTreatment) 2);
  double lox = -2; double hix = -2; double loy = -2; double hiy = -2;
  if (log != "Y") {loy = 0; hiy = 7;}
  if (flag == "mass") {lox = -1; hix = -1; loy = 0; hiy = 2;} else if (log == "Y") {lox = 15; hix = 60; loy = 4e-7; hiy = 400;} else {lox = 0; hix = 0.5;}
  
  Prettify1D(raw, kBlack, kFullStar, 2, kBlack, xTitle, yTitle,lox,hix,loy, hiy);
  Prettify1DwLineStyle(gen, kGreen, kDashed, 5, xTitle, yTitle,lox,hix,loy, hiy);
  Prettify1D(det, kBlue, kOpenCircle, 2, kBlue, xTitle, yTitle,lox,hix,loy, hiy);
  Prettify1D(reco2, kRed, kFullStar, 2, kRed, xTitle, yTitle,lox,hix,loy, hiy);
  reco2->SetTitle("");
  if (flag == "pt") {
    TAxis *reco_axis = reco2->GetXaxis();
    int bmin = reco_axis->FindBin(15);
    int bmax = reco_axis->FindBin(60);
    double binwidth = (reco_axis->GetXmax() - reco_axis->GetXmin()) / (double)reco_axis->GetNbins();    
    reco2->Scale(1/(double)reco2->Integral(bmin,bmax));
    gen->Scale(1/(double)gen->Integral(bmin,bmax));
    gen->Scale(1/(double)binwidth);
    reco2->Scale(1/(double)binwidth);
  }
  TLegend *tu;
  if (xTitle.find("R_{g}") != std::string::npos) {tu = TitleLegend(0.15,0.57,0.51,0.87);}
  else {tu = TitleLegend(0.48,0.57,0.84,0.87);}
  tu->AddEntry(gen,"PYTHIA6","l");
  tu->AddEntry(det,"PYTHIA6+GEANT","p");
  tu->AddEntry(raw,"Raw data", "p");
  tu->AddEntry(reco2,"Unfolded data (4 iter)","p");
  
  reco2->Draw(); raw->Draw("same"); gen->Draw("C,same"); det->Draw("same"); tu->Draw("same");

  cout << (gen->GetBinContent(3) / (double) reco2->GetBinContent(3)) << endl;
  cout << (gen->GetBinContent(4) / (double) reco2->GetBinContent(4)) << endl;
  cout << (gen->GetBinContent(5) / (double) reco2->GetBinContent(5)) << endl;
  cout << (gen->GetBinContent(6) / (double) reco2->GetBinContent(6)) << endl;
  
 
  cu->SaveAs((out + "unfolded_" + raw->GetName() + filetype).c_str());
 
  return;
}

void Unfold4D(TFile *matchFile, TH2D* raw, TH2D* gen, TH2D* det, const string log, const string out, const string filetype, const string xTitle, const string resName, const string flag) {
  string cName = (string) "c4" + (string) raw->GetName();
  TCanvas *c4 = MakeCanvas(cName, log, 800,800);

  std::string yTitle = "";
  if (flag == "pt") {
    yTitle = "arb.";
    raw->Scale(1/(double)raw->Integral()); det->Scale(1/(double)det->Integral());
  }
  else {yTitle = ("1/N_{j} dN_{j}/d" + xTitle.substr(0,xTitle.find(' '))).c_str();}
  
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

  if (flag == "pt") {
    rawXs[0]->Scale(1/(double)rawXs[0]->Integral()); detXs[0]->Scale(1/(double)detXs[0]->Integral());
    double rat1D_scale = rawXs[0]->GetBinContent(1) / (double) genXs[0]->GetBinContent(3);
    genXs[0]->Scale(rat1D_scale);
  }

  double lox = -2; double hix = -2; double loy = -2; double hiy = -2;

  if (log != "Y") {loy = 0; hiy = 7;}
  if (flag == "mass") {lox = -1; hix = -1; loy = 0; hiy = 2;} else if (log == "Y") {loy = 4e-7; hiy = 400;} else {lox = 0; hix = 0.5;}

  
  Prettify1D(rawXs[0], kBlack, kFullStar, 2, kBlack, xTitle, yTitle,lox, hix,loy,hiy);
  Prettify1DwLineStyle(genXs[0], kGreen, kDashed, 5, xTitle, yTitle,lox, hix,loy,hiy);
  Prettify1D(detXs[0], kBlue, kOpenCircle, 2, kBlue, xTitle, yTitle,lox, hix,loy,hiy);
  Prettify1D(recoXs[0], kRed, kFullStar, 2, kRed, xTitle, yTitle,lox, hix,loy,hiy);
  recoXs[0]->SetTitle("");

  TLegend *tu;
  if (xTitle.find("R_{g}") != std::string::npos) {tu = TitleLegend(0.15,0.57,0.51,0.87);}
  else {tu = TitleLegend(0.48,0.57,0.84,0.87);}
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
  
  std::string yTitle = "";
  if (flag == "pt") {
    yTitle = "arb.";
    raw->Scale(1/(double)raw->Integral()); det->Scale(1/(double)det->Integral());
  }
  else {yTitle = ("1/N_{j} dN_{j}/d" + xTitle.substr(0,xTitle.find(' '))).c_str();}

  RooUnfoldBayes *unfold4D_4iter = new RooUnfoldBayes(res4D, raw, 4, false, ((string) "unfold4D_4iter" + (string) raw->GetName()).c_str(),"");
  TH2D *reco = (TH2D*) unfold4D_4iter->Hreco((RooUnfold::ErrorTreatment) 2);
  
  const unsigned nBins = 5;
  double ranges[nBins+1] = {2,3,4,5,7,11};
  double corresp_pts[nBins+1] = {15,20,25,30,40,60};
  double ranges_d[nBins+1] = {0,1,2,3,5,9};
  vector<TH1D*> recoXs = Projection2D(reco, nBins, ranges, "x"); 
  vector<TH1D*> rawXs = Projection2D(raw, nBins, ranges_d, "x");
  vector<TH1D*> genXs = Projection2D(gen, nBins, ranges, "x");
  vector<TH1D*> detXs = Projection2D(det, nBins, ranges_d, "x");
  
  if (flag == "pt") {
    for (int i = 0; i < 2; ++ i) {
      rawXs[i]->Scale(1/(double)rawXs[i]->Integral()); detXs[i]->Scale(1/(double)detXs[i]->Integral());
      double rat1D_scale = rawXs[i]->GetBinContent(1) / (double) genXs[i]->GetBinContent(3);
      genXs[i]->Scale(rat1D_scale);
    }
    for (int i = 2; i < genXs.size(); ++ i) {
      rawXs[i]->Scale(1/(double)rawXs[i]->Integral()); detXs[i]->Scale(1/(double)detXs[i]->Integral());
      double rat1D_scale = rawXs[i]->GetBinContent(3) / (double) genXs[i]->GetBinContent(5);
      genXs[i]->Scale(rat1D_scale);
    }
  }
  /*
  double lox, hix, loy, hiy;
  if (log == "Y") {loy = 0.0001; hiy = 50;} else {loy = 0; hiy = 0.5;}
  if (flag == "mass") {lox = -1; hix = -1;} else {lox = 0; hix = 0.5;}
*/
  double lox = -2; double hix = -2; double loy = -2; double hiy = -2;
  
  if (log != "Y") {loy = 0; hiy = 7;}
  if (flag == "mass") {lox = -1; hix = -1; loy = 0; hiy = 2;} else if (log == "Y") {loy = 0.001; hiy = 50;} else {lox = 0; hix = 0.5;}
  

  for(int i = 0; i < nBins; ++ i) { 
    Prettify1D(rawXs[i], kBlack, kFullStar, 2, kBlack, xTitle, yTitle,lox, hix,loy,hiy);
    Prettify1DwLineStyle(genXs[i], kGreen, kDashed, 5, xTitle, yTitle,lox, hix,loy,hiy);
    Prettify1D(detXs[i], kBlue, kOpenCircle, 2, kBlue, xTitle, yTitle,lox, hix,loy,hiy);
    Prettify1D(recoXs[i], kRed, kFullStar, 2, kRed, xTitle, yTitle,lox, hix,loy,hiy);
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
  //vector<TH1D*> pts_15 = {(TH1D*) inFile->Get("pt_d_pt15"), (TH1D*) inFile->Get("pt_p_pt15"), (TH1D*) inFile->Get("pt_g_pt15")};
  
  vector<TH1D*> ms = {(TH1D*) inFile->Get("m_d"), (TH1D*) inFile->Get("m_p"), (TH1D*) inFile->Get("m_g")};
  vector<TH1D*> zgs = {(TH1D*) inFile->Get("zg_d"), (TH1D*) inFile->Get("zg_p"), (TH1D*) inFile->Get("zg_g")};
  vector<TH1D*> rgs = {(TH1D*) inFile->Get("rg_d"), (TH1D*) inFile->Get("rg_p"), (TH1D*) inFile->Get("rg_g")};
  vector<TH1D*> ptgs = {(TH1D*) inFile->Get("ptg_d"), (TH1D*) inFile->Get("ptg_p"), (TH1D*) inFile->Get("ptg_g")};
  vector<TH1D*> mgs = {(TH1D*) inFile->Get("mg_d"), (TH1D*) inFile->Get("mg_p"), (TH1D*) inFile->Get("mg_g")};
  
  UnfoldedObs(matchFile, pts[0], pts[1], pts[2], "pt_res_coarse", "p_{T} [GeV/c]", "Y", out, filetype, "pt");
  UnfoldedObs(matchFile, ms[0], ms[1], ms[2], "m_response", "M [GeV/c^{2}]", "0", out, filetype, "mass");
  UnfoldedObs(matchFile, zgs[0] , zgs[1], zgs[2], "zg_response", "z_{g}", "0", out, filetype, "");
  UnfoldedObs(matchFile, rgs[0], rgs[1], rgs[2], "rg_response", "R_{g}", "0", out, filetype,"");
  UnfoldedObs(matchFile, ptgs[0], ptgs[1], ptgs[2], "ptg_response", "p_{T,g} [GeV/c]", "Y", out, filetype,"pt");
  UnfoldedObs(matchFile, mgs[0], mgs[1], mgs[2], "mg_response", "M_{g} [GeV/c^{2}]", "0", out, filetype,"mass");
 
  return;
}

void Draw4DResponse(TFile* matchFile, const std::string resName, const std::string out, const std::string filetype, const int nBinsx, const int nBinsy) {
  TCanvas * c4res = MakeCanvas(("c4res_" + resName).c_str(), "z", 800,800);
  RooUnfoldResponse *res4D = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  TH2D* res_matrix = (TH2D*) res4D->Hresponse();
  Prettify2D(res_matrix, "p^{meas.}_{T} [GeV/c]", "p^{truth}_{T} [GeV/c]", -1,-1, -1,-1, 10e-11, 10e-6);
  
  //  res_matrix->GetXaxis()->SetBinLabel(180, "60");
  for (int i = 1; i < (9*nBinsx + 1); i += nBinsx) {
    if (nBinsx == 20) {res_matrix->GetXaxis()->SetBinLabel(i,(std::to_string( (int) ((i-1)/(double) 4) + 15)).c_str()); }
    if (nBinsx == 9) {res_matrix->GetXaxis()->SetBinLabel(i,(std::to_string( (int) (5*(i-1)/(double) 9) + 15)).c_str()); }
  }

  for (int i = 1; i <= (14*nBinsy + 1); i += 2*nBinsy) {
    if (nBinsy == 20) {res_matrix->GetYaxis()->SetBinLabel(i, (std::to_string( (int) ((i-1)/(double) 4) + 5) ).c_str());} 
    if (nBinsy == 15) {res_matrix->GetYaxis()->SetBinLabel(i,(std::to_string( (int) ((i-1)/(double) 3) + 5) ).c_str());}
  }
  res_matrix->SetTitle("");
  
  res_matrix->Draw("colz");

  c4res->SaveAs((out+"4D_"+resName+filetype).c_str());

  return;
}

void MCClosure1D(TFile *File, const std::string xTitle, const std::string resName,const std::string out,const std::string filetype, TH1D* detOdd, TH1D* genOdd, TH1D* genEven) {
  TCanvas* cclos = MakeCanvas("cclos","0",800,800); cclos->cd();
  
  std::string resOdd = (resName + "_odd").c_str(); std::string resEven = (resName + "_even").c_str();
  //  RooUnfoldResponse *res_odd = (RooUnfoldResponse*) File->Get(resOdd.c_str());
  RooUnfoldResponse *res_even = (RooUnfoldResponse*) File->Get(resEven.c_str());

  RooUnfoldBayes *unfolded_opp = new RooUnfoldBayes(res_even, detOdd, 4, false, ("unfolded_4iter" + (string) detOdd->GetName() + "_opp").c_str(),"");
  //RooUnfoldBayes *unfolded_same = new RooUnfoldBayes(res_odd, detOdd, 4, false, ("unfolded_4iter" + (string) detOdd->GetName() + "_same").c_str(),"");
  
  TH1D *reco_opp = (TH1D*) unfolded_opp->Hreco((RooUnfold::ErrorTreatment) 3);
  //TH1D *reco_same = (TH1D*) unfolded_same->Hreco((RooUnfold::ErrorTreatment) 3);

  //  reco_opp->Scale(1/(double)reco_opp->Integral()); reco_same->Scale(1/(double)reco_same->Integral());
  //genEven->Scale(1/(double)genEven->Integral()); genOdd->Scale(1/(double)genOdd->Integral());
  //detOdd->Scale(1/(double)detOdd->Integral());
  /*  
  for (int i =0; i < 11; ++ i) {
    cout << genOdd->GetBinContent(i) << endl;
  }
  */
  reco_opp->Divide(genEven); //reco_same->Divide(genOdd);
  /*  
  detOdd->SetMarkerStyle(kOpenTriangleUp);
  genOdd->SetMarkerStyle(kOpenCircle); genEven->SetMarkerStyle(kOpenCircle); reco_opp->SetMarkerStyle(kOpenSquare); reco_same->SetMarkerStyle(kOpenSquare);
  reco_opp->SetMarkerColor(kRed); reco_same->SetMarkerColor(kBlue); genEven->SetMarkerColor(kRed); genOdd->SetMarkerColor(kBlue);
  genOdd->Draw(); genEven->Draw("same"); reco_opp->Draw("same"); reco_same->Draw("same"); detOdd->Draw("same");
  */
  //reco_opp->SetMarkerStyle(kOpenSquare); reco_same->SetMarkerStyle(kFullCircle);
  
  double loy = 0.5; double hiy = 1.5; double lox = 0; double hix = 60;
  
  Prettify1D(reco_opp, kBlue, kOpenSquare, 1, kBlue, xTitle, "arb.",lox,hix,loy,hiy);
  //Prettify1D(reco_same, kViolet, kFullCircle, 1, kViolet, xTitle, "arb.", lox,hix,loy,hiy);
  reco_opp->SetTitle("STAR Simulation"); //reco_same->SetTitle("");
  
  TLine *one = new TLine(0,1,60,1); TLine *low = new TLine(0,0.95,60,0.95); TLine *high = new TLine(0,1.05,60,1.05); 
  one->SetLineStyle(kDashed); low->SetLineStyle(kDashed); high->SetLineStyle(kDashed);
  
  TLegend *t = TitleLegend(0.15,0.15, 0.45, 0.35);
  t->AddEntry((TObject*)0,"Bayesian unfolding closure test","");
  TLegend *l = new TLegend(0.15,0.7,0.45, 0.85); l->SetBorderSize(0);
  l->AddEntry(reco_opp,"Opp. side (4 iter.)","p");
  //l->AddEntry(reco_same,"Same side (4 iter.)","p");
  
  int nBins = reco_opp->GetNbinsX();
  double ygen[nBins]; double ydet[nBins]; double yunf[nBins]; double xpart[nBins]; double xdet[nBins];
  for (int i = 1; i <= nBins; ++ i) {
    yunf[i-1] = reco_opp->GetBinError(i) / (double) reco_opp->GetBinContent(i);
    ygen[i-1] = genOdd->GetBinError(i) / (double) genOdd->GetBinContent(i);
    ydet[i-1] = detOdd->GetBinError(i) / (double) detOdd->GetBinContent(i);
    xpart[i-1] = 2.5 + 5*i;
    xdet[i-1] = 12.5 + 5*i;
  }
  
  TGraph *error_unf = new TGraph(nBins, xpart, yunf);
  TGraph *error_gen = new TGraph(nBins, xpart, ygen);
  TGraph *error_det = new TGraph(nBins, xdet, ydet);

  error_unf->SetMarkerStyle(kOpenCircle); error_unf->SetMarkerColor(kBlue);
  error_gen->SetMarkerStyle(kOpenCircle); error_gen->SetMarkerColor(kRed);
  error_det->SetMarkerStyle(kOpenCircle); error_det->SetMarkerColor(kGreen);

  TLegend * e = new TLegend(0.6,0.2,0.9,0.4); e->SetBorderSize(0);
  e->AddEntry(error_unf,"Unfolded","p");
  e->AddEntry(error_gen,"Gen-level","p");
  e->AddEntry(error_det,"Det-level","p");
  
  /*  
  error_unf->GetXaxis()->SetTitle("p_{T} [GeV/c]"); error_unf->GetYaxis()->SetTitle("Relative stat. error (bin err. / content)");
  error_unf->GetYaxis()->SetRangeUser(1e-3,5); error_unf->SetTitle(""); cclos->SetLogy();
  
  error_unf->Draw("AP");
  error_gen->Draw("P,same");
  error_det->Draw("P,same");
  e->Draw("same");
  */
  reco_opp->Draw(); /*reco_same->Draw("same");*/ one->Draw("same"); low->Draw("same"); high->Draw("same");
  t->Draw("same"); l->Draw("same");
  
  cclos->SaveAs((out+ "1Dclosure_full_"+ resName +filetype).c_str());
  
  return;
}

void MCClosure2DFull(TFile *File, const std::string xTitle, const std::string resName,const std::string out,const std::string filetype, TH2D* det, TH2D* gen_opp) {
  TCanvas* cclos2D = MakeCanvas(("cclos2D" + resName).c_str(),"0",800,800); cclos2D->cd();
  
  RooUnfoldResponse *res = (RooUnfoldResponse*) File->Get(resName.c_str());
  
  RooUnfoldBayes *unfolded_opp = new RooUnfoldBayes(res, det, 4, false, ("unfolded_4iter" + (string) det->GetName() + "_opp").c_str(),"");
  
  TH2D *reco_opp = (TH2D*) unfolded_opp->Hreco((RooUnfold::ErrorTreatment) 3);
  
  const unsigned nBins = 1;
  double ranges[nBins+1] = {0,999};
  vector<TH1D*> recoYs = Projection2D(reco_opp, nBins, ranges, "y");
  vector<TH1D*> genYs = Projection2D(gen_opp, nBins, ranges, "y");
  vector<TH1D*> detYs = Projection2D(det, nBins, ranges, "y");
  
  recoYs[0]->Divide(genYs[0]);
  
  double loy = 0.5; double hiy = 1.5; double lox = 0; double hix = 60;
  
  Prettify1D(recoYs[0], kBlue, kOpenSquare, 1, kBlue, xTitle, "arb.",lox,hix,loy,hiy);
  recoYs[0]->SetTitle("STAR Simulation");
  
  TLine *one = new TLine(0,1,60,1); TLine *low = new TLine(0,0.95,60,0.95); TLine *high = new TLine(0,1.05,60,1.05); 
  one->SetLineStyle(kDashed); low->SetLineStyle(kDashed); high->SetLineStyle(kDashed);
   
  TLegend *t = TitleLegend(0.15,0.15, 0.45, 0.35);
  t->AddEntry((TObject*)0,"2D Bayesian unfolding closure test","");
  TLegend *l = new TLegend(0.15,0.7,0.45, 0.85); l->SetBorderSize(0);
  l->AddEntry(recoYs[0],"Opp. side (4 iter.)","p");

  recoYs[0]->Draw(); one->Draw("same"); low->Draw("same"); high->Draw("same");
  t->Draw("same"); l->Draw("same");
  
  cclos2D->SaveAs((out+ "2Dclosure_full_"+ resName +filetype).c_str());
  
  return;
}

void MCClosure2DSlices(TFile *File, const std::string xTitle, const std::string resName,const std::string out,const std::string filetype, TH2D* det, TH2D* gen_opp) {
  TCanvas* cclosSlice = MakeCanvas(("cclosSlice" + resName).c_str(),"0",800,800); cclosSlice->cd();
  DivideCanvas(cclosSlice,"0",3,2);
  
  RooUnfoldResponse *res = (RooUnfoldResponse*) File->Get(resName.c_str());
  
  RooUnfoldBayes *unfolded_opp = new RooUnfoldBayes(res, det, 4, false, ("unfolded_4iter" + (string) det->GetName() + "_opp").c_str(),"");
  
  TH2D *reco_opp = (TH2D*) unfolded_opp->Hreco((RooUnfold::ErrorTreatment) 3);
  
  const unsigned nBins = 5;
  double ranges[nBins+1] = {2,3,4,5,7,11};
  double corresp_pts[nBins+1] = {15,20,25,30,40,60};
  //double ranges_d[nBins+1] = {0,1,2,3,5,9};
  vector<TH1D*> recoXs = Projection2D(reco_opp, nBins, ranges, "x");
  vector<TH1D*> genXs = Projection2D(gen_opp, nBins, ranges, "x");
  //vector<TH1D*> detYs = Projection2D(det, nBins, ranges_d, "y");

  for (int i = 0; i < recoXs.size(); ++ i) {
    recoXs[i]->Divide(genXs[i]);
  }
  
  double loy = 0.5; double hiy = 1.5; double lox = 0; double hix = 0.5; 
  if (xTitle.find("p_{T}") != std::string::npos) {hix *= 120;} else if (xTitle.find("M") != std::string::npos) {hix *= 20;}
  
  for (int i = 0; i < recoXs.size(); ++ i) {
    Prettify1D(recoXs[i], kBlue, kOpenSquare, 1, kBlue, xTitle, "arb.",lox,hix,loy,hiy);  
  }
  recoXs[0]->SetTitle("STAR Simulation");

  TLine *one = new TLine(0,1,60,1); TLine *low = new TLine(0,0.95,60,0.95); TLine *high = new TLine(0,1.05,60,1.05); 
  one->SetLineStyle(kDashed); low->SetLineStyle(kDashed); high->SetLineStyle(kDashed);
   
  TLegend *t = TitleLegend(0.15,0.15, 0.45, 0.35);
  t->AddEntry((TObject*)0,"2D Bayesian unfolding closure test","");
  TLegend *l = new TLegend(0.15,0.7,0.45, 0.85); l->SetBorderSize(0);
  l->AddEntry(recoXs[0],"Opp. side (4 iter.)","p");
  
  TLegend *tslices[nBins];
  for (int i = 0; i < nBins; ++ i) {
    tslices[i] = SliceLegend(((to_string(corresp_pts[i])).substr(0,2) + " < p_{T}^{unfolded jet} < " + (to_string(corresp_pts[i + 1])).substr(0,2) + " GeV/c").c_str(), 0.13,0.8,0.9,0.95);
  }

  cclosSlice->cd(1); t->Draw(); l->Draw("same");
  for (int i = 0; i < recoXs.size(); ++ i) {
    cclosSlice->cd(i+2);
    recoXs[i]->Draw(); one->Draw("same"); low->Draw("same"); high->Draw("same");
    tslices[i]->Draw("same");
  }
  
  cclosSlice->SaveAs((out+ "2Dclosure_slices_"+ resName +filetype).c_str());
  
  return;
}


void unfolding () {
  string dir = "~/jetmass/";
  string matchin = "out/matching/";
  string closein = "out/closure/";
  string pyin = "out/sim/py/";
  string gein = "out/sim/ge/";
  string datain = "out/data/";
  string in = "macros/hists/";
  string infile = "hists.root";
  string file = "full.root";
  string out = "~/jetmass/plots/unfolding/";
  string closeout = "~/jetmass/plots/closure/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";

  TFile* matchFile = new TFile( (dir + matchin + file).c_str(), "READ");
  //TFile* pyFile = new TFile( (dir + pyin +file).c_str(), "READ");
  //TFile* geFile = new TFile( (dir + gein +file).c_str(), "READ");
  //TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  TFile* inFile = new TFile( (dir + in + infile).c_str(), "READ");
  TFile* closureFile = new TFile( (dir + closein + file).c_str(), "READ");
  //TFile* raghavFile = new TFile("~/jetmass/GeantToPythia_Match_ppRun12_200GeV_NoEff_NoBg_JP2_R04.root","READ");

  //plots the raw and unfolded data for all observables (pT, m, zg, Rg) along with truth & reco
  //  AllUnfolds(matchFile, inFile, out, filetype);
  
  vector<TH2D*> m2d = {(TH2D*) inFile->Get("m_v_pt_d"), (TH2D*) inFile->Get("m_v_pt_p"), (TH2D*) inFile->Get("m_v_pt_g")};
  vector<TH2D*> zg2d = {(TH2D*) inFile->Get("zg_v_pt_d"), (TH2D*) inFile->Get("zg_v_pt_p"), (TH2D*) inFile->Get("zg_v_pt_g")};
  vector<TH2D*> rg2d = {(TH2D*) inFile->Get("rg_v_pt_d"), (TH2D*) inFile->Get("rg_v_pt_p"), (TH2D*) inFile->Get("rg_v_pt_g")};
  vector<TH2D*> ptg2d = {(TH2D*) inFile->Get("ptg_v_pt_d"), (TH2D*) inFile->Get("ptg_v_pt_p"), (TH2D*) inFile->Get("ptg_v_pt_g")};
  vector<TH2D*> mg2d = {(TH2D*) inFile->Get("mg_v_pt_d"), (TH2D*) inFile->Get("mg_v_pt_p"), (TH2D*) inFile->Get("mg_v_pt_g")};

  /*  
  Unfold4D(matchFile, m2d[0], m2d[1], m2d[2], "0", out, filetype, "M [GeV/c^{2}]", "pt_m_response", "mass");
  Unfold4D(matchFile, zg2d[0], zg2d[1], zg2d[2], "0", out, filetype, "z_{g}", "pt_zg_response", "");
  Unfold4D(matchFile, rg2d[0], rg2d[1], rg2d[2], "0", out, filetype, "R_{g}", "pt_rg_response", "");
  Unfold4D(matchFile, ptg2d[0], ptg2d[1], ptg2d[2], "Y", out, filetype, "p_{T,g} [GeV/c]", "pt_ptg_response", "pt");
  Unfold4D(matchFile, mg2d[0], mg2d[1], mg2d[2], "0", out, filetype, "M_{g} [GeV/c^{2}]", "pt_mg_response", "mass");
  */
  
  SliceUnfolded4D(matchFile, m2d[0], m2d[1], m2d[2], "0", out, filetype, "M_{jet} [GeV/c^{2}]", "pt_m_response", "mass");
  SliceUnfolded4D(matchFile, zg2d[0], zg2d[1], zg2d[2], "0", out, filetype, "z_{g}", "pt_zg_response", "");
  SliceUnfolded4D(matchFile, rg2d[0], rg2d[1], rg2d[2], "0", out, filetype, "R_{g}", "pt_rg_response", "");
  SliceUnfolded4D(matchFile, ptg2d[0], ptg2d[1], ptg2d[2], "Y", out, filetype, "p_{T,g}", "pt_ptg_response", "pt");
  SliceUnfolded4D(matchFile, mg2d[0], mg2d[1], mg2d[2], "0", out, filetype, "R_{g}", "pt_mg_response", "mass");
  
  /*
  Draw4DResponse(matchFile, "pt_m_response", out, filetype, 20, 20);
  Draw4DResponse(matchFile, "pt_zg_response", out, filetype, 20, 20);
  Draw4DResponse(matchFile, "pt_rg_response", out, filetype, 20, 20);
  Draw4DResponse(matchFile, "pt_ptg_response", out, filetype, 9, 15);
  Draw4DResponse(matchFile, "pt_mg_response", out, filetype, 20, 20);
  */
  
  TH1D* genPtSame = (TH1D*) closureFile->Get("pt_gen_odd"); TH1D* detPt = (TH1D*) closureFile->Get("pt_det_odd");
  TH1D* genPtOpp = (TH1D*) closureFile->Get("pt_gen_even");
  
  TH2D* genPtM = (TH2D*) closureFile->Get("pt_m_gen"); TH2D* detPtM = (TH2D*) closureFile->Get("pt_m_det");
  TH2D* genPtZg = (TH2D*) closureFile->Get("pt_zg_gen"); TH2D* detPtZg = (TH2D*) closureFile->Get("pt_zg_det");
  TH2D* genPtRg = (TH2D*) closureFile->Get("pt_rg_gen"); TH2D* detPtRg = (TH2D*) closureFile->Get("pt_rg_det");
  
  
  //MCClosure1D(closureFile,"p_{T} [GeV/c]", "pt",closeout,filetype, detPt, genPtSame, genPtOpp);
  /*  
  MCClosure2DFull(closureFile,"p_{T} [GeV/c]", "pt_m_response",closeout,filetype, detPtM, genPtM);
  MCClosure2DFull(closureFile,"p_{T} [GeV/c]", "pt_zg_response",closeout,filetype, detPtZg, genPtZg);
  MCClosure2DFull(closureFile,"p_{T} [GeV/c]", "pt_rg_response",closeout,filetype, detPtRg, genPtRg);
  */
  /*
  MCClosure2DSlices(closureFile,"M [GeV/c^{2}]", "pt_m_response",closeout,filetype, detPtM, genPtM);
  MCClosure2DSlices(closureFile,"z_{g}", "pt_zg_response",closeout,filetype, detPtZg, genPtZg);
  MCClosure2DSlices(closureFile,"R_{g}", "pt_rg_response",closeout,filetype, detPtRg, genPtRg);
  */
  return;
}
