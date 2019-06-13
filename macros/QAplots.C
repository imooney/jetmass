
#include "Plots.h"

using namespace std;

void ratioplot(TH1D* pp, TH1D* pAu, TH1D* pAuBBCMB, const int nevts_pp, const int nevts_pAu, const int nevts_pAuBBCMB, const string uniqid) {
  pp->SetTitle(""); pAu->SetTitle(""); pAuBBCMB->SetTitle("");
  
  pp->Scale(1/(double)nevts_pp);
  pAu->Scale(1/(double)nevts_pAu);
  pAuBBCMB->Scale(1/(double)nevts_pAuBBCMB);

  pp->SetMarkerStyle(kOpenCircle);
  pAu->SetMarkerStyle(kOpenCircle);
  pAuBBCMB->SetMarkerStyle(kOpenCircle);
  //pp->SetMarkerSize(2);
  //pAu->SetMarkerSize(2);
  pp->SetMarkerColor(kRed);
  pAu->SetMarkerColor(kMagenta);
  pAuBBCMB->SetMarkerColor(kBlue);
  pp->SetLineColor(kRed);
  pAu->SetLineColor(kMagenta);
  pAuBBCMB->SetLineColor(kBlue);

  pp->GetXaxis()->SetRangeUser(0,30);
  pAu->GetXaxis()->SetRangeUser(0,30);
  pAuBBCMB->GetXaxis()->SetRangeUser(0,30);

  const double meanpp = (double) pp->GetMean();
  const double meanpAu = (double) pAu->GetMean();
  const double meanpAuBBCMB = (double) pAuBBCMB->GetMean();

  TCanvas *c = new TCanvas(("c"+uniqid).c_str(),("c"+uniqid).c_str(),800,800);

  TLegend *t = new TLegend(0.5,0.6,0.7,0.8); t->SetBorderSize(0);
  t->AddEntry(pp,"ppJP2", "p");
  t->AddEntry(pAu,"pAuJP2", "p");
  t->AddEntry(pAuBBCMB,"pAuBBCMB","p");

  TLine *lpp = new TLine(meanpp,0,meanpp,1e10);
  TLine *lpAu = new TLine(meanpAu,0,meanpAu,1e10);
  TLine *lpAuBBCMB = new TLine(meanpAuBBCMB,0,meanpAuBBCMB,1e10);


  lpp->SetLineStyle(kDashed); lpAu->SetLineStyle(kDashed); lpAuBBCMB->SetLineStyle(kDashed);
  lpp->SetLineColor(kRed); lpAu->SetLineColor(kMagenta); lpAuBBCMB->SetLineColor(kBlue);

  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  pp->SetStats(0);          // No statistics on upper plot
  pp->Draw();               // Draw pp
  pAu->Draw("same");         // Draw h2 on top of pp
  pAuBBCMB->Draw("same");
  t->Draw("same");
  lpp->Draw("same");
  lpAu->Draw("same");
  lpAuBBCMB->Draw("same");
  // Do not draw the Y axis label on the upper plot and redraw a small
  // axis instead, in order to avoid the first label (0) to be clipped.
  /*pp->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis->SetLabelSize(15);
  axis->Draw();
  */
  // lower plot will be in pad
  c->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  
  // Define the ratio plot
  TH1D *rat = (TH1D*)pp->Clone("rat");
  rat->SetLineColor(kRed);
  // rat->SetMinimum(0);  // Define Y ..
  //rat->SetMaximum(2); // .. range
  rat->Sumw2();
  rat->SetStats(0);      // No statistics on lower plot
  rat->Divide(pAu);
  rat->SetMarkerStyle(kOpenCircle);
  rat->Draw("ep");       // Draw the ratio plot
  
  TH1D *ratpAu = (TH1D*)pAuBBCMB->Clone("ratpAu");
  ratpAu->SetLineColor(kBlue);
  ratpAu->Sumw2();
  ratpAu->SetStats(0);
  ratpAu->Divide(pAu);
  ratpAu->SetMarkerStyle(kOpenCircle);
  ratpAu->Draw("epsame");

  // Y axis pp plot settings
  /*
  pp->GetYaxis()->SetTitleSize(20);
  pp->GetYaxis()->SetTitleFont(43);
  pp->GetYaxis()->SetTitleOffset(1.55);
  */

  // Ratio plot (rat) settings
  rat->SetTitle(""); // Remove the ratio title
  
  // Y axis ratio plot settings
  rat->GetYaxis()->SetTitle("pp/pAu & pAu BBCMB/JP2 ");
  /* rat->GetYaxis()->SetNdivisions(505);
  rat->GetYaxis()->SetTitleSize(20);
  rat->GetYaxis()->SetTitleFont(43);
  rat->GetYaxis()->SetTitleOffset(1.55);
  rat->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  rat->GetYaxis()->SetLabelSize(15);
  */
  /*
  // X axis ratio plot settings
  rat->GetXaxis()->SetTitleSize(20);
  rat->GetXaxis()->SetTitleFont(43);
  rat->GetXaxis()->SetTitleOffset(4.);
  rat->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  rat->GetXaxis()->SetLabelSize(15);
  */
}

void plot(TH1D* pp, TH1D* pAu, TH1D* pAuBBCMB, const int nevts_pp, const int nevts_pAu, const int nevts_pAuBBCMB, const string uniqid) {
  pp->SetTitle(""); pAu->SetTitle(""); pAuBBCMB->SetTitle("");
  
  pp->Scale(1/(double)nevts_pp);
  pAu->Scale(1/(double)nevts_pAu);
  pAuBBCMB->Scale(1/(double)nevts_pAuBBCMB);

  pp->SetMarkerStyle(kOpenCircle);
  pAu->SetMarkerStyle(kOpenCircle);
  pAuBBCMB->SetMarkerStyle(kOpenCircle);
  //pp->SetMarkerSize(2);
  //pAu->SetMarkerSize(2);
  pp->SetMarkerColor(kRed);
  pAu->SetMarkerColor(kMagenta);
  pAuBBCMB->SetMarkerColor(kBlue);
  pp->SetLineColor(kRed);
  pAu->SetLineColor(kMagenta);
  pAuBBCMB->SetLineColor(kBlue);
  
  const double meanpp = (double) pp->GetMean();
  const double meanpAu = (double) pAu->GetMean();
  const double meanpAuBBCMB = (double) pAuBBCMB->GetMean();

  TCanvas *c = new TCanvas(("c"+uniqid).c_str(),("c"+uniqid).c_str(),800,800);
  
  TLegend *t = new TLegend(0.5,0.6,0.7,0.8); t->SetBorderSize(0);
  t->AddEntry(pp,"ppJP2", "p");
  t->AddEntry(pAu,"pAuJP2", "p");
  t->AddEntry(pAuBBCMB,"pAuBBCMB","p");

  TLine *lpp = new TLine(meanpp,0,meanpp,1e10);
  TLine *lpAu = new TLine(meanpAu,0,meanpAu,1e10);
  TLine *lpAuBBCMB = new TLine(meanpAuBBCMB,0,meanpAuBBCMB,1e10);
  
  lpp->SetLineStyle(kDashed); lpAu->SetLineStyle(kDashed); lpAuBBCMB->SetLineStyle(kDashed);
  lpp->SetLineColor(kRed); lpAu->SetLineColor(kMagenta); lpAuBBCMB->SetLineColor(kBlue);

  pp->Draw();
  pAu->Draw("same");
  pAuBBCMB->Draw("same");
  t->Draw("same");
  lpp->Draw("same");
  lpAu->Draw("same");
  lpAuBBCMB->Draw("same");

  return;
}

void projections (TH2D* h, const string uniqid, const vector<int> bins, const string axis, const string projecting_over, vector<double> proj_integrals) {
  TCanvas *c = new TCanvas(("c_proj_"+uniqid).c_str(),("c_proj_"+uniqid).c_str(),800,800);
  
  const int nBins = bins.size() - 1;
  
  vector<TH1D*> hists;
  
  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4); leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  
  //for BBC coinc. legend:
  vector<string> mhz_strings = {"0","0.25","0.5","0.75","1","1.25","1.5"};
  
  for (int i = 0; i < nBins; ++ i) {
    if (axis == "x") {
      hists.push_back(h->ProjectionX(("proj_"+to_string(i)+uniqid).c_str(),h->GetYaxis()->FindBin(bins[i]),h->GetYaxis()->FindBin(bins[i+1])));
    }
    else if (axis == "y") {
      hists.push_back(h->ProjectionY(("proj_"+to_string(i)+uniqid).c_str(),h->GetXaxis()->FindBin(bins[i]),h->GetXaxis()->FindBin(bins[i+1])));
    }
    if (projecting_over == "BBC coinc. [MHz]") {
      leg->AddEntry(hists[i],(mhz_strings[i]+" < "+projecting_over+" < "+mhz_strings[i+1]).c_str(),"p");
    }
    else {
      leg->AddEntry(hists[i],(to_string(bins[i])+" < "+projecting_over+" < "+to_string(bins[i+1])).c_str(),"p");
    }
  }
  
  for (int i = 0; i < nBins; ++ i) {
    if (hists[i]->Integral()) {
      cout << "evt scale " << proj_integrals[i] << " and trk/tow scale " << hists[i]->Integral() << endl;
      if (axis == "y") {
	cout << "for " << hists[i]->GetName() << " evt scale " << proj_integrals[i] << endl;
	hists[i]->Scale(1/(double)proj_integrals[i]);
      }
      else {
	cout << "for " << hists[i]->GetName() << " trk/tow scale " << hists[i]->Integral() << endl;
	hists[i]->Scale(1/(double)hists[i]->Integral());
      }
    }
    hists[i]->SetMarkerColor(i+1); hists[i]->SetLineColor(i+1);
    if (i >= 4) {hists[i]->SetMarkerColor(i+2); hists[i]->SetLineColor(i+2);}
    hists[i]->SetMarkerStyle(20+i);
    //hists[i]->GetYaxis()->SetRangeUser(0,0.012);
    //hists[i]->GetYaxis()->SetTitle("prob.");
    hists[i]->SetTitle("");
    if (hists[i]->GetEntries() > 20) {hists[i]->Draw("same");}
  }
  leg->Draw("same");
}

void plot2D(TH2D* h, const string uniqid, vector<double> proj_integrals) {
  /*
  TCanvas *c = new TCanvas(("c"+uniqid).c_str(),("c"+uniqid).c_str(),800,800);
  
  h->Draw("colz");
  */
  vector<int> bins;
  string axis; string projecting_over;
  //change the numbers later!
  if (uniqid.find("3") != string::npos) {
    axis = "x"; projecting_over = "N_{trk}"; bins = {0,25,50,75,100,125}; 
  }
  else if (uniqid.find("0") != string::npos || uniqid.find("1") != string::npos || uniqid.find("2") != string::npos || uniqid.find("4") != string::npos || uniqid.find("5") != string::npos || uniqid.find("6") != string::npos || uniqid.find("7") != string::npos || uniqid.find("9") != string::npos || uniqid.find("10") != string::npos) {
    axis = "y"; projecting_over = "BBC coinc. [MHz]"; bins = {0,(int) 2.5e5,(int) 5e5,(int) 7.5e5,(int) 1e6,(int) 1.25e6,(int) 1.5e6};
  }
  else {
    axis = "x"; projecting_over = "N/A"; bins = {0};
  }
  if (projecting_over != "N/A") {
    projections(h, uniqid, bins,axis,projecting_over, proj_integrals);
  }
  return;
}

void plotGr(TGraph* g, const string uniqid) {  
  g->GetXaxis()->SetTitle("Tower ID");
  g->SetTitle("");

  
  if (uniqid.find("0") != string::npos) {
    g->GetYaxis()->SetRangeUser(0,2.7);
    g->GetYaxis()->SetTitle("<E_{T}>");
  }
  if (uniqid.find("1") != string::npos) {
    g->GetYaxis()->SetRangeUser(0,24);
    g->GetYaxis()->SetTitle("<(E_{T} > 2 GeV)>");
  }
  g->SetMarkerStyle(kOpenCircle);
  g->SetMarkerSize(0.2);
  g->SetMarkerColor(kBlue+2);
  g->SetLineColor(kBlue+2);

  //RMS:
  double sumx2 = 0;//(double) g->GetRMS(2);
  double sumx = 0;
  double x = 0.0; double y = 0.0;
  for (int i = 0; i < g->GetN(); ++ i) {
    g->GetPoint(i, x, y);
    //cout << "y: " << y << endl;
    if (y != y) continue; //skips empty towers since they mess up the calculation
    sumx += y;
    sumx2 += y*y;
    //cout << "sumx2: " << sumx2 << endl;
  }
  double mean = sumx / (double) g->GetN();
  //cout << "mean: " << mean << endl;
  double sigma = TMath::Abs((sumx2/(double) g->GetN()) - mean*mean);
  //cout << "sigma: " << sigma << endl;
  const double three_sigma = 3*sigma;
  
  TCanvas *c = new TCanvas(("c"+uniqid).c_str(),("c"+uniqid).c_str(),800,800);

  double lowerbound = 0;
  //  if (uniqid.find("1") != string::npos) { lowerbound = 2; } //for tower Et average with towers above 2 GeV
  if ((mean - three_sigma) > lowerbound) {lowerbound = mean - three_sigma;}
  
  TLine *lplus = new TLine(0,mean+three_sigma,4800,mean+three_sigma);
  TLine *lmean = new TLine(0,mean,4800,mean);
  TLine *lminus = new TLine(0,lowerbound,4800,lowerbound);
  
  cout << "mean: " << mean << " 3sigma: " << three_sigma << endl;
  cout << "mean+3sigma: " << mean+three_sigma << endl;
  cout << "mean-3sigma: " << lowerbound << endl;
  
  lplus->SetLineStyle(kSolid); lminus->SetLineStyle(kSolid); lmean->SetLineStyle(kSolid);
  lplus->SetLineColor(kRed); lminus->SetLineColor(kRed); lmean->SetLineColor(kRed);

  g->Draw("AP");
  lplus->Draw("same");
  lminus->Draw("same");
  lmean->Draw("same");

  return;
}

void QAplots () {
  TFile *f = new TFile("~/jetmass/macros/QA_pp_pAu_w_MB_and_vpdvzcut_3.root","READ");

  vector<TH1D*> pp_vec = {(TH1D*) f->Get("hn_trks_ppJP2"),(TH1D*) f->Get("hn_tows_ppJP2"),(TH1D*) f->Get("hbbc_coinc_ppJP2"),(TH1D*) f->Get("hevt_vtx_ppJP2"),(TH1D*) f->Get("htrackPt_ppJP2"),(TH1D*) f->Get("htrackEta_ppJP2"),(TH1D*) f->Get("htrackPhi_ppJP2"),(TH1D*) f->Get("htrackDCA_ppJP2"),(TH1D*) f->Get("htowerEt_ppJP2"),(TH1D*) f->Get("htowerEta_ppJP2"),(TH1D*) f->Get("htowerPhi_ppJP2"),(TH1D*) f->Get("htowerId_ppJP2"),(TH1D*) f->Get("htrackNhits_ppJP2"),(TH1D*) f->Get("htrackNhitsposs_ppJP2"),(TH1D*) f->Get("htrackNhitsratio_ppJP2"),(TH1D*) f->Get("hn_globals_ppJP2")};

  vector<TH1D*> pAu_vec = {(TH1D*) f->Get("hn_trks_pAuJP2"),(TH1D*) f->Get("hn_tows_pAuJP2"),(TH1D*) f->Get("hbbc_coinc_pAuJP2"),(TH1D*) f->Get("hevt_vtx_pAuJP2"),(TH1D*) f->Get("htrackPt_pAuJP2"),(TH1D*) f->Get("htrackEta_pAuJP2"),(TH1D*) f->Get("htrackPhi_pAuJP2"),(TH1D*) f->Get("htrackDCA_pAuJP2"),(TH1D*) f->Get("htowerEt_pAuJP2"),(TH1D*) f->Get("htowerEta_pAuJP2"),(TH1D*) f->Get("htowerPhi_pAuJP2"),(TH1D*) f->Get("htowerId_pAuJP2"),(TH1D*) f->Get("htrackNhits_pAuJP2"),(TH1D*) f->Get("htrackNhitsposs_pAuJP2"),(TH1D*) f->Get("htrackNhitsratio_pAuJP2"),(TH1D*) f->Get("hn_globals_pAuJP2")};

    vector<TH1D*> pAuBBCMB_vec = {(TH1D*) f->Get("hn_trks_pAuBBCMB"),(TH1D*) f->Get("hn_tows_pAuBBCMB"),(TH1D*) f->Get("hbbc_coinc_pAuBBCMB"),(TH1D*) f->Get("hevt_vtx_pAuBBCMB"),(TH1D*) f->Get("htrackPt_pAuBBCMB"),(TH1D*) f->Get("htrackEta_pAuBBCMB"),(TH1D*) f->Get("htrackPhi_pAuBBCMB"),(TH1D*) f->Get("htrackDCA_pAuBBCMB"),(TH1D*) f->Get("htowerEt_pAuBBCMB"),(TH1D*) f->Get("htowerEta_pAuBBCMB"),(TH1D*) f->Get("htowerPhi_pAuBBCMB"),(TH1D*) f->Get("htowerId_pAuBBCMB"),(TH1D*) f->Get("htrackNhits_pAuBBCMB"),(TH1D*) f->Get("htrackNhitsposs_pAuBBCMB"),(TH1D*) f->Get("htrackNhitsratio_pAuBBCMB"),(TH1D*) f->Get("hn_globals_pAuBBCMB")};


    vector<TH2D*> pp2D_vec = {(TH2D*) f->Get("hbbc_coinc_evt_vtx_ppJP2"),(TH2D*) f->Get("hbbc_coinc_n_trks_ppJP2"),(TH2D*) f->Get("hbbc_coinc_n_tows_ppJP2"),(TH2D*) f->Get("htrackDCA_n_trks_ppJP2"),(TH2D*) f->Get("hbbc_coinc_trackDCA_ppJP2"),(TH2D*) f->Get("hbbc_coinc_trackPt_ppJP2"),(TH2D*) f->Get("hbbc_coinc_towerEt_ppJP2"),(TH2D*) f->Get("htrackEta_Phi_ppJP2"),(TH2D*) f->Get("htowerEta_Phi_ppJP2"),(TH2D*) f->Get("hbbc_coinc_n_vertices_ppJP2"),(TH2D*) f->Get("hbbc_coinc_n_globals_ppJP2")};

  vector<TH2D*> pAu2D_vec = {(TH2D*) f->Get("hbbc_coinc_evt_vtx_pAuJP2"),(TH2D*) f->Get("hbbc_coinc_n_trks_pAuJP2"),(TH2D*) f->Get("hbbc_coinc_n_tows_pAuJP2"),(TH2D*) f->Get("htrackDCA_n_trks_pAuJP2"),(TH2D*) f->Get("hbbc_coinc_trackDCA_pAuJP2"),(TH2D*) f->Get("hbbc_coinc_trackPt_pAuJP2"),(TH2D*) f->Get("hbbc_coinc_towerEt_pAuJP2"),(TH2D*) f->Get("htrackEta_Phi_pAuJP2"),(TH2D*) f->Get("htowerEta_Phi_pAuJP2"),(TH2D*) f->Get("hbbc_coinc_n_vertices_pAuJP2"),(TH2D*) f->Get("hbbc_coinc_n_globals_pAuJP2")};
  
   vector<TH2D*> pAuBBCMB2D_vec = {(TH2D*) f->Get("hbbc_coinc_evt_vtx_pAuBBCMB"),(TH2D*) f->Get("hbbc_coinc_n_trks_pAuBBCMB"),(TH2D*) f->Get("hbbc_coinc_n_tows_pAuBBCMB"),(TH2D*) f->Get("htrackDCA_n_trks_pAuBBCMB"),(TH2D*) f->Get("hbbc_coinc_trackDCA_pAuBBCMB"),(TH2D*) f->Get("hbbc_coinc_trackPt_pAuBBCMB"),(TH2D*) f->Get("hbbc_coinc_towerEt_pAuBBCMB"),(TH2D*) f->Get("htrackEta_Phi_pAuBBCMB"),(TH2D*) f->Get("htowerEta_Phi_pAuBBCMB"),(TH2D*) f->Get("hbbc_coinc_n_vertices_pAuBBCMB"),(TH2D*) f->Get("hbbc_coinc_n_globals_pAuBBCMB")};
  
  vector<TGraph*> ppgr_vec = {(TGraph*) f->Get("htowerId_meanEt_ppJP2"),(TGraph*) f->Get("htowerId_meanEtg2GeV_ppJP2")};

  vector<TGraph*> pAugr_vec = {(TGraph*) f->Get("htowerId_meanEt_pAuJP2"),(TGraph*) f->Get("htowerId_meanEtg2GeV_pAuJP2")};

  vector<TGraph*> pAuBBCMBgr_vec = {(TGraph*) f->Get("htowerId_meanEt_pAuBBCMB"),(TGraph*) f->Get("htowerId_meanEtg2GeV_pAuBBCMB")};

  const int nevts_pp = pp_vec[0]->Integral();
  const int nevts_pAu = pAu_vec[0]->Integral();
  const int nevts_pAuBBCMB = pAuBBCMB_vec[0]->Integral();


  //for scaling track/tower info by n_events in a given slice of BBC coincidence:
  TH2D* hppJP2 = (TH2D*) f->Get("hbbc_coinc_evt_vtx_ppJP2");
  TH2D* hpAuJP2 = (TH2D*) f->Get("hbbc_coinc_evt_vtx_pAuJP2");
  TH2D* hpAuBBCMB = (TH2D*) f->Get("hbbc_coinc_evt_vtx_pAuBBCMB");

  const vector<int> bins_bbccoinc = {0,(int) 2.5e5,(int) 5e5,(int) 7.5e5,(int) 1e6,(int) 1.25e6,(int) 1.5e6};
  const int nBins_bbccoinc = 6;

  vector<TH1D*> projspp; vector<TH1D*> projspAu; vector<TH1D*> projspAuBBCMB;
  vector<double> intpp; vector<double> intpAu; vector<double> intpAuBBCMB;
  for (int i = 0; i < nBins_bbccoinc; ++ i) {
    projspp.push_back(hppJP2->ProjectionY(("proj_"+to_string(i)+"pp").c_str(),hppJP2->GetXaxis()->FindBin(bins_bbccoinc[i]),hppJP2->GetXaxis()->FindBin(bins_bbccoinc[i+1])));
    projspAu.push_back(hpAuJP2->ProjectionY(("proj_"+to_string(i)+"pAu").c_str(),hpAuJP2->GetXaxis()->FindBin(bins_bbccoinc[i]),hpAuJP2->GetXaxis()->FindBin(bins_bbccoinc[i+1])));
    projspAuBBCMB.push_back(hpAuBBCMB->ProjectionY(("proj_"+to_string(i)+"pAuBBCMB").c_str(),hpAuBBCMB->GetXaxis()->FindBin(bins_bbccoinc[i]),hpAuBBCMB->GetXaxis()->FindBin(bins_bbccoinc[i+1])));
    
    intpp.push_back(projspp[i]->Integral());
    intpAu.push_back(projspAu[i]->Integral());
    intpAuBBCMB.push_back(projspAuBBCMB[i]->Integral());

  }
  
  /*
  for (int i = pp_vec.size() - 4; i < pp_vec.size(); ++ i) {
    if (i == 4 || i == 8) {
      // ratioplot(pp_vec[i],pAu_vec[i],nevts_pp,nevts_pAu,to_string(i));
      ratioplot(pp_vec[i],pAu_vec[i],pAuBBCMB_vec[i],nevts_pp,nevts_pAu,nevts_pAuBBCMB,to_string(i));
    }
    if (i == 12 || i == 13 || i == 14) {
      plot(pp_vec[i],pAu_vec[i],pAuBBCMB_vec[i],pp_vec[i]->Integral(),pAu_vec[i]->Integral(),pAuBBCMB_vec[i]->Integral(),to_string(i));
    }
    else {
      plot(pp_vec[i],pAu_vec[i],pAuBBCMB_vec[i],nevts_pp,nevts_pAu,nevts_pAuBBCMB,to_string(i));
    }
  }
  
    
  for (int i = pp2D_vec.size() - 1; i < pp2D_vec.size(); ++ i) {
    plot2D(pp2D_vec[i],(to_string(i)+"pp2D").c_str(), intpp);
    plot2D(pAu2D_vec[i],(to_string(i)+"pAu2D").c_str(), intpAu);
    plot2D(pAuBBCMB2D_vec[i],(to_string(i)+"pAuBBCMB2D").c_str(), intpAuBBCMB);
  }
  */
  for (int i = 0; i < ppgr_vec.size(); ++ i) {
    //    plotGr(ppgr_vec[i],(to_string(i)+"ppGr").c_str());
    //plotGr(pAugr_vec[i],(to_string(i)+"pAuGr").c_str());
    plotGr(pAuBBCMBgr_vec[i],(to_string(i)+"pAuBBCMBGr").c_str());
  }
  
  return;
}

