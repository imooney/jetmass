#include "Plots.h"

using namespace std;

void compareTauM () {
  
  TFile *t = new TFile ("~/jetmass/production/macros/hists/hists.root","READ");
 
  //note: typically defined with 2 - a instead of a, so tau2 = mass, tau 1 = girth, etc.
  TH3D *tau1 = (TH3D*) t->Get("jetM_v_a1_v_pt_dec");
  TH3D *tau05 = (TH3D*) t->Get("jetM_v_a05_v_pt_dec");
  TH3D *tau0 = (TH3D*) t->Get("jetM_v_a0_v_pt_dec");
  TH3D *tau_05 = (TH3D*) t->Get("jetM_v_a_05_v_pt_dec");
  TH3D *tau_1 = (TH3D*) t->Get("jetM_v_a_1_v_pt_dec");

  vector<TH2D*> tau1_v_m;
  vector<TH2D*> tau05_v_m;
  vector<TH2D*> tau0_v_m;
  vector<TH2D*> tau_05_v_m;
  vector<TH2D*> tau_1_v_m;
  
  const int nBins_pt = 5;
  int ranges_pt[nBins_pt+1] = {3,4,5,6,8,12};
  string pts[nBins_pt+1] = {"15","20","25","30","40","60"};
  
  tau1_v_m = Projection3D (tau1, nBins_pt,ranges_pt, "yx");
  tau05_v_m = Projection3D (tau05, nBins_pt,ranges_pt, "yx");
  tau0_v_m = Projection3D (tau0, nBins_pt,ranges_pt, "yx");
  tau_05_v_m = Projection3D (tau_05, nBins_pt,ranges_pt, "yx");
  tau_1_v_m = Projection3D (tau_1, nBins_pt,ranges_pt, "yx");
  
  TLatex *p = new TLatex(); TLatex *slice = new TLatex();
  
  TCanvas* c1 = MakeCanvas("c1","0",1400,1000);
  DivideCanvas(c1,"z",3,2);
  TCanvas* c05 = MakeCanvas("c05","0",1400,1000);
  DivideCanvas(c05,"z",3,2);
  TCanvas* c0 = MakeCanvas("c0","0",1400,1000);
  DivideCanvas(c0,"z",3,2);
  TCanvas* c_05 = MakeCanvas("c_05","0",1400,1000);
  DivideCanvas(c_05,"z",3,2);
  TCanvas* c_1 = MakeCanvas("c_1","0",1400,1000);
  DivideCanvas(c_1,"z",3,2);
  
  TH2D* hdummy1 = new TH2D("hdummy1","",30,-8,0,30,-8,0);
  TH2D* hdummy05 = new TH2D("hdummy05","",30,-8,0,30,-8,0);
  TH2D* hdummy0 = new TH2D("hdummy0","",30,-8,0,30,-8,0);
  TH2D* hdummy_05 = new TH2D("hdummy_05","",30,-8,0,30,-8,0);
  TH2D* hdummy_1 = new TH2D("hdummy_1","",30,-8,0,30,-8,0);

  TLine *yx = new TLine(-8,-8,0,0); yx->SetLineStyle(kDashed);
  
  c1->cd(1); hdummy1->Draw("same"); p=PanelTitle();
  for(int i = 0; i < nBins_pt; ++ i) {
    tau1_v_m[i]->GetXaxis()->SetTitle("M^{2}/p_{T}^{2} [c^{-2}]"); tau1_v_m[i]->GetYaxis()->SetTitle("#tau_{1}");
    tau1_v_m[i]->GetZaxis()->SetRangeUser(1e4,pow(10,10-i));
    c1->cd(i+2); tau1_v_m[i]->Draw("colz"); yx->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
  
  c05->cd(1); hdummy05->Draw("same"); p=PanelTitle();
  for(int i = 0; i < nBins_pt; ++ i) {
    tau05_v_m[i]->GetXaxis()->SetTitle("M^{2}/p_{T}^{2} [c^{-2}]"); tau05_v_m[i]->GetYaxis()->SetTitle("#tau_{0.5}");
    tau05_v_m[i]->GetZaxis()->SetRangeUser(1e4,pow(10,10-i));
    c05->cd(i+2); tau05_v_m[i]->Draw("colz"); yx->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }

  c0->cd(1); hdummy0->Draw("same"); p=PanelTitle();
  for(int i = 0; i < nBins_pt; ++ i) {
    tau0_v_m[i]->GetXaxis()->SetTitle("M^{2}/p_{T}^{2} [c^{-2}]"); tau0_v_m[i]->GetYaxis()->SetTitle("#tau_{0}");
    tau0_v_m[i]->GetZaxis()->SetRangeUser(1e4,pow(10,10-i));
    c0->cd(i+2); tau0_v_m[i]->Draw("colz"); yx->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }

  c_05->cd(1); hdummy_05->Draw("same"); p=PanelTitle();
  for(int i = 0; i < nBins_pt; ++ i) {
    tau_05_v_m[i]->GetXaxis()->SetTitle("M^{2}/p_{T}^{2} [c^{-2}]"); tau_05_v_m[i]->GetYaxis()->SetTitle("#tau_{-0.5}");
    tau_05_v_m[i]->GetZaxis()->SetRangeUser(1e4,pow(10,10-i));
    c_05->cd(i+2); tau_05_v_m[i]->Draw("colz"); yx->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
  
  c_1->cd(1); hdummy_1->Draw("same"); p=PanelTitle();
  for(int i = 0; i < nBins_pt; ++ i) {
    tau_1_v_m[i]->GetXaxis()->SetTitle("M^{2}/p_{T}^{2} [c^{-2}]"); tau_1_v_m[i]->GetYaxis()->SetTitle("#tau_{-1}");
    tau_1_v_m[i]->GetZaxis()->SetRangeUser(1e4,pow(10,10-i));
    c_1->cd(i+2); tau_1_v_m[i]->Draw("colz"); yx->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
  
  
  c1->SaveAs("~/jetmass/plots/angularity/tau1_v_m.pdf");
  c05->SaveAs("~/jetmass/plots/angularity/tau05_v_m.pdf");
  c0->SaveAs("~/jetmass/plots/angularity/tau0_v_m.pdf");
  c_05->SaveAs("~/jetmass/plots/angularity/tau_05_v_m.pdf");
  c_1->SaveAs("~/jetmass/plots/angularity/tau_1_v_m.pdf");
  
  return;
}
