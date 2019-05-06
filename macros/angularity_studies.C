#include "Plots.h"

using namespace std;

void MakePlot(TH2D* hist2D0, TH2D* hist2D05, TH2D* hist2D_05, TH2D* hist2D_1, const string type) {
  const unsigned nBins = 5;
  string corresp_pts[nBins+1] = {"15","20","25","30","40","60"};
  double ranges_d[nBins+1] = {1,2,3,4,6,10};
  
  vector<TH1D*> hists0 = Projection2D(hist2D0, nBins, ranges_d, "x");
  vector<TH1D*> hists05 = Projection2D(hist2D05, nBins, ranges_d, "x");
  vector<TH1D*> hists_05 = Projection2D(hist2D_05, nBins, ranges_d, "x");
  vector<TH1D*> hists_1 = Projection2D(hist2D_1, nBins, ranges_d, "x");

  for(int i = 0; i < nBins; ++ i) {
    Prettify1DwLineStyle(hists0[i], kBlack, kSolid, 2, "log_{10}(#tau_{a})", "1/N_{jets} dN/dlog_{10}(#tau_{a})",-8,0,0,1.6);
    hists0[i]->SetFillColor(kBlack); hists0[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(hists05[i], kBlue, kSolid, 2, "log_{10}(#tau_{a})", "1/N_{jets} dN/dlog_{10}(#tau_{a})",-8,0,0,1.6);
    hists05[i]->SetFillColor(kBlue); hists05[i]->SetFillStyle(3395);
    Prettify1DwLineStyle(hists_05[i], kRed, kSolid, 2, "log_{10}(#tau_{a})", "1/N_{jets} dN/dlog_{10}(#tau_{a})",-8,0,0,1.6);
    hists_05[i]->SetFillColor(kRed); hists_05[i]->SetFillStyle(3490);
    Prettify1DwLineStyle(hists_1[i], kOrange, kSolid, 2, "log_{10}(#tau_{a})", "1/N_{jets} dN/dlog_{10}(#tau_{a})",-8,0,0,1.6); 
    hists_1[i]->SetFillColor(kOrange); hists_1[i]->SetFillStyle(3436);
  }
  
  TLegend * t1 = new TLegend(0.4,0.6,0.55,0.8); t1->SetBorderSize(0);
  t1->AddEntry(hists0[0],"a = 0","f");
  t1->AddEntry(hists05[0],"a = 0.5","f");
  TLegend * t2 = new TLegend(0.12,0.65,0.32,0.84); t2->SetBorderSize(0);
  t2->AddEntry(hists_05[0],"a = -0.5", "f");
  t2->AddEntry(hists_1[0],"a = -1","f");
  
  TLatex *slice = new TLatex();

  TH1D* hdummy = new TH1D(("hdummy"+type).c_str(),";log_{10}(#tau_{a});1/N_{jets} dN/dlog_{10}(#tau_{a})",1,-8,0);
  hdummy->GetYaxis()->SetRangeUser(0,1.6);

  TCanvas *ca = new TCanvas(("ca"+type).c_str(),("ca"+type).c_str(),1400,1000);
  DivideCanvas(ca,"0",3,2);
  
  ca->cd(1); hdummy->Draw(); TLatex* title = PanelTitle(); title->DrawLatexNDC(0.2,0.35, "Charged tracks assigned m_{#pi}");
  
  for (int i = 0; i < nBins; ++ i) {
    ca->cd(i+2);
    hists0[i]->Draw("lf2same"); hists05[i]->Draw("lf2same"); hists_05[i]->Draw("lf2same"); hists_1[i]->Draw("lf2same");
    slice->DrawLatexNDC(0.15,0.85,(corresp_pts[i]+" < p_{T} < "+corresp_pts[i+1]+" GeV/c").c_str());
    if (i == 1) {t1->Draw("same");} if (i == 4) {t2->Draw("same");}
  }
  
  ca->SaveAs(("~/jetmass/plots/angularity/angularities_" +type+".pdf").c_str());
  
  return;
}

void angularity_studies () {
  TFile *f = new TFile ("~/jetmass/macros/hists/angularity_hists.root","READ");
  
  TH2D* tau0_v_pt_d = (TH2D*) f->Get("tau0_v_pt_d");
  TH2D* tau05_v_pt_d = (TH2D*) f->Get("tau05_v_pt_d");
  TH2D* tau_05_v_pt_d = (TH2D*) f->Get("tau_05_v_pt_d");
  TH2D* tau_1_v_pt_d = (TH2D*) f->Get("tau_1_v_pt_d");
  TH2D* tau0_g_v_pt_d = (TH2D*) f->Get("tau0_g_v_pt_d");
  TH2D* tau05_g_v_pt_d = (TH2D*) f->Get("tau05_g_v_pt_d");
  TH2D* tau_05_g_v_pt_d = (TH2D*) f->Get("tau_05_g_v_pt_d");
  TH2D* tau_1_g_v_pt_d = (TH2D*) f->Get("tau_1_g_v_pt_d");

  MakePlot(tau0_v_pt_d, tau05_v_pt_d, tau_05_v_pt_d, tau_1_v_pt_d, "ungroomed");
  MakePlot(tau0_g_v_pt_d, tau05_g_v_pt_d, tau_05_g_v_pt_d, tau_1_g_v_pt_d, "groomed");
  
  return;
}
