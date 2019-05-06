#include "Plots.h"

using namespace std;

void poster_plots () {
  
  TFile *f = new TFile("~/jetmass/macros/hists/hists_w_o_bin_drop.root","READ");
  TH2D* det_res = (TH2D*) f->Get("ratioMvPyPt");
  TAxis* yax = (TAxis*) det_res->GetYaxis();
  
  const int nBins = 5;
  double ranges[nBins+1] = {(double) yax->FindBin(15),(double) yax->FindBin(20),(double) yax->FindBin(25),(double) yax->FindBin(30),(double) yax->FindBin(40),(double) yax->FindBin(60)};
  
  vector<TH1D*> res_projs = Projection2D(det_res, nBins, ranges, "x");
  vector<Color_t> cols_tom_blue = {kBlack, kBlue, kGreen+2, kMagenta, kOrange};
  vector<Style_t> marks = {kOpenCircle, kFullCircle, kOpenStar, kFullSquare, kOpenSquare};
  vector<string> pts = {"15", "20", "25", "30", "40", "60"};
  
  TCanvas *cres = new TCanvas("cres","cres",1400,1000); cres->cd();
  TArrow *ar = new TArrow(0.7,1.8,1,1.8,0.03,"<|");
  ar->SetAngle(60);
  ar->SetLineWidth(8);
  ar->SetLineColor(4);
  //ar->SetFillStyle(3008);
  ar->SetFillColor(kBlue);
  //ar->Draw();
  TLatex *trend = new TLatex(); trend->SetTextColor(kRed);
  TLatex *trend2 = new TLatex();
  
  TLegend *tres = new TLegend(0.55,0.3,0.6,0.8); tres->SetBorderSize(0);
  TLine *vert = new TLine(1,0,1,2.5); vert->SetLineStyle(kDashed);
  
  for (int i = 1; i < res_projs.size() - 1; ++ i) {
    tres->AddEntry(res_projs[i],("p_{T} #in ("+pts[i]+","+pts[i+1]+") GeV/c").c_str(),"p");
    Prettify1D(res_projs[i],cols_tom_blue[i],marks[i],2,cols_tom_blue[i],"R = M^{det}_{j} / M^{gen}_{j}", "1/N_{pairs} dN/dR",0,2,0,2.5);
    res_projs[i]->Draw("same");
  }
  gStyle->SetLegendTextSize(0.05);
  tres->Draw("same");
  vert->Draw("same");
  ar->Draw();
  trend2->DrawLatex(0.65,1.9,"mass loss");
  trend->DrawLatex(0.1,2.3,"PREVIEW - Work in Progress");
  cres->SaveAs("~/poster_plots/m_res.pdf");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~systematics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  




  return;
}
