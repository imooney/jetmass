
#include "Plots.h"

using namespace std;

void plot_from3D(int sim, string obs, string xtitle, string ytitle, string xunit, string yunit, TH3D* hist, vector<double> xaxes, vector<double> yaxes, vector<double> zaxes) {
  const int nBins = 5;
  TAxis* zax = (TAxis*) hist->GetZaxis();
  int bins[nBins+1] = {(int)zax->FindBin(15),(int)zax->FindBin(20),(int)zax->FindBin(25),(int)zax->FindBin(30),(int)zax->FindBin(40),(int)zax->FindBin(60)};
  string pts[nBins+1] = {"15","20","25","30","40","60"};
  
  vector<TH2D*> hist2Ds = Projection3D (hist, nBins, bins, "yx");

  for (int i = 0; i < nBins; ++ i) {
    Prettify2D (hist2Ds[i], (xtitle+ " " + xunit).c_str(), (ytitle+" "+yunit).c_str(), xaxes[0],xaxes[1], yaxes[0], yaxes[1], zaxes[0], zaxes[1]);
  }
    
  TH2D* hdummy = new TH2D(("hdummy"+obs+to_string(sim)).c_str(),"",1,xaxes[0],xaxes[1],1,yaxes[0],yaxes[1]); hdummy->GetZaxis()->SetRangeUser(zaxes[0],zaxes[1]);
  
  TLatex *p = new TLatex(); TLatex *slice = new TLatex();
  TCanvas *c = new TCanvas (("c"+obs+to_string(sim)).c_str(),("c"+obs+to_string(sim)).c_str(),1400,1000); c->cd();
  DivideCanvas(c,"0",3,2);
  
  c->cd(1); hdummy->Draw(); p = PanelTitle();
  for (int i = 0; i < nBins; ++ i) {
    c->cd(i+2);
    hist2Ds[i]->Draw("colz");
    slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
  
  string fileappend;
  if (sim == 0) { fileappend = "py8dec"; }
  if (sim == 1) { fileappend = "py8undec"; }
  if (sim == 2) { fileappend = "py6dec"; }
  if (sim == 3) { fileappend = "py6undec"; }
  if (sim == 4) { fileappend = "hw7dec"; }
  if (sim == 5) { fileappend = "hw7undec"; }

  c->SaveAs(("~/jetmass/plots/obs_correlations/"+obs+"_"+fileappend+".pdf").c_str());   
  
  return;
}

void ObservableCorrelations() {
  
  TFile *f = new TFile("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved.root","READ");
  
  vector<TH3D*> p8on = {(TH3D*)f->Get("mvmgvpt_p8on"),(TH3D*)f->Get("mvzgvpt_p8on"),(TH3D*)f->Get("mvrgvpt_p8on"),
		(TH3D*)f->Get("mgvzgvpt_p8on"),(TH3D*)f->Get("mgvrgvpt_p8on"),(TH3D*)f->Get("zgvrgvpt_p8on")};
  vector<TH3D*> p8off = {(TH3D*)f->Get("mvmgvpt_p8off"),(TH3D*)f->Get("mvzgvpt_p8off"),(TH3D*)f->Get("mvrgvpt_p8off"),
		 (TH3D*)f->Get("mgvzgvpt_p8off"),(TH3D*)f->Get("mgvrgvpt_p8off"),(TH3D*)f->Get("zgvrgvpt_p8off")};
  vector<TH3D*> p6on = {(TH3D*)f->Get("mvmgvpt_p6on"),(TH3D*)f->Get("mvzgvpt_p6on"),(TH3D*)f->Get("mvrgvpt_p6on"),
		(TH3D*)f->Get("mgvzgvpt_p6on"),(TH3D*)f->Get("mgvrgvpt_p6on"),(TH3D*)f->Get("zgvrgvpt_p6on")};
  vector<TH3D*> p6off = {(TH3D*)f->Get("mvmgvpt_p6off"),(TH3D*)f->Get("mvzgvpt_p6off"),(TH3D*)f->Get("mvrgvpt_p6off"),
		 (TH3D*)f->Get("mgvzgvpt_p6off"),(TH3D*)f->Get("mgvrgvpt_p6off"),(TH3D*)f->Get("zgvrgvpt_p6off")};
  vector<TH3D*> h7on = {(TH3D*)f->Get("mvmgvpt_h7on"),(TH3D*)f->Get("mvzgvpt_h7on"),(TH3D*)f->Get("mvrgvpt_h7on"),
		(TH3D*)f->Get("mgvzgvpt_h7on"),(TH3D*)f->Get("mgvrgvpt_h7on"),(TH3D*)f->Get("zgvrgvpt_h7on")};
  vector<TH3D*> h7off = {(TH3D*)f->Get("mvmgvpt_h7off"),(TH3D*)f->Get("mvzgvpt_h7off"),(TH3D*)f->Get("mvrgvpt_h7off"),
		 (TH3D*)f->Get("mgvzgvpt_h7off"),(TH3D*)f->Get("mgvrgvpt_h7off"),(TH3D*)f->Get("zgvrgvpt_h7off")};

  vector<vector<TH3D*> > all = {p8on, p8off, p6on, p6off, h7on, h7off};
  
  vector<double> m_axes = {0,10};
  vector<double> zr_axes = {0,0.5};
  vector<double> z_range = {-1,-1};
  
  for (int i = 0; i < all.size(); ++ i) {
    plot_from3D(i,"mvmg", "M", "M_{g}", "[GeV/c^{2}]", "[GeV/c^{2}]", all[i][0], m_axes, m_axes, z_range);
  }
  for (int i = 0; i < all.size(); ++ i) {
    plot_from3D(i,"mvzg", "M", "z_{g}", "[GeV/c^{2}]", "", all[i][1], m_axes, zr_axes, z_range);
  }
  for (int i = 0; i < all.size(); ++ i) {
    plot_from3D(i,"mvrg", "M", "R_{g}", "[GeV/c^{2}]", "", all[i][2], m_axes, zr_axes, z_range);
  }
  for (int i = 0; i < all.size(); ++ i) {
    plot_from3D(i,"mgvzg", "M_{g}", "z_{g}", "[GeV/c^{2}]", "", all[i][3], m_axes, zr_axes, z_range);
  }
  for (int i = 0; i < all.size(); ++ i) {
    plot_from3D(i,"mgvrg", "M_{g}", "R_{g}", "[GeV/c^{2}]", "", all[i][4], m_axes, zr_axes, z_range);
  }
  for (int i = 0; i < all.size(); ++ i) {
    plot_from3D(i,"zgvrg", "z_{g}", "R_{g}", "", "", all[i][5], zr_axes, zr_axes, z_range);
  }
  
  return;
}
