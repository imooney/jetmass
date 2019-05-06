#include "Plots.h"

using namespace std;

void ratio_plots () {
  TFile *f = new TFile("~/jetmass/production/macros/hists/ratio.root","READ");
  
  TH1D* off1520 = (TH1D*) f->Get("hoff1520");
  TH1D* off2025 = (TH1D*) f->Get("hoff2025");
  TH1D* off2530 = (TH1D*) f->Get("hoff2530");
  TH1D* off3040 = (TH1D*) f->Get("hoff3040");
  TH1D* off4060 = (TH1D*) f->Get("hoff4060");
  TH1D* on1520 = (TH1D*) f->Get("hon1520");
  TH1D* on2025 = (TH1D*) f->Get("hon2025");
  TH1D* on2530 = (TH1D*) f->Get("hon2530");
  TH1D* on3040 = (TH1D*) f->Get("hon3040");
  TH1D* on4060 = (TH1D*) f->Get("hon4060");
  TH1D* ratio1520 = (TH1D*) f->Get("hratio1520");
  TH1D* ratio2025 = (TH1D*) f->Get("hratio2025");
  TH1D* ratio2530 = (TH1D*) f->Get("hratio2530");
  TH1D* ratio3040 = (TH1D*) f->Get("hratio3040");
  TH1D* ratio4060 = (TH1D*) f->Get("hratio4060");

  TH1D* pton = (TH1D*) f->Get("pton");
  TH1D* ptoff = (TH1D*) f->Get("ptoff");
  TH1D* ptratio = (TH1D*) f->Get("ptratio");
  
  Prettify1D (off1520, kRed, kFullCircle, 1, kRed, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (off2025, kRed, kFullCircle, 1, kRed, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (off2530, kRed, kFullCircle, 1, kRed, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (off3040, kRed, kFullCircle, 1, kRed, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (off4060, kRed, kFullCircle, 1, kRed, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (on1520, kBlue, kOpenCircle, 1, kBlue, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (on2025, kBlue, kOpenCircle, 1, kBlue, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (on2530, kBlue, kOpenCircle, 1, kBlue, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (on3040, kBlue, kOpenCircle, 1, kBlue, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (on4060, kBlue, kOpenCircle, 1, kBlue, "M^{jet} [GeV/c^{2}]", "1/N^{j} dN/dM^{j}", -1,-1,0,0.5);
  Prettify1D (ratio1520, kBlue, kOpenSquare, 2, kBlue, "M^{jet} [GeV/c^{2}]", "decays off / on", -1,-1,0,2);
  Prettify1D (ratio2025, kBlue, kOpenSquare, 2, kBlue, "M^{jet} [GeV/c^{2}]", "decays off / on", -1,-1,0,2);
  Prettify1D (ratio2530, kBlue, kOpenSquare, 2, kBlue, "M^{jet} [GeV/c^{2}]", "decays off / on", -1,-1,0,2);
  Prettify1D (ratio3040, kBlue, kOpenSquare, 2, kBlue, "M^{jet} [GeV/c^{2}]", "decays off / on", -1,-1,0,2);
  Prettify1D (ratio4060, kBlue, kOpenSquare, 2, kBlue, "M^{jet} [GeV/c^{2}]", "decays off / on", -1,-1,0,2);

  Prettify1D (pton, kBlue, kOpenCircle, 1, kBlue, "p_{T}^{jet} [GeV/c]", "1/N^{j} dN/dp_{T}", -1,-1,-1,-1);
  Prettify1D (ptoff, kRed, kFullCircle, 1, kRed, "p_{T}^{jet} [GeV/c]", "1/N^{j} dN/dp_{T}", -1,-1,-1,-1);
  Prettify1D (ptratio, kBlack, kOpenSquare, 1, kBlack, "p_{T}^{jet} [GeV/c]", "Decays off / on", -1,-1,0,2);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  TCanvas *cpt = MakeCanvas("cpt","0",1000,1000);
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetLogy();
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  ptoff->SetStats(0);          // No statistics on upper plot
  ptoff->Draw();               // Draw h1
  pton->Draw("same");         // Draw h2 on top of h1

  TLegend *t = new TLegend(0.5,0.7,0.87,0.87); t->SetBorderSize(0);
  t->AddEntry(ptoff,"Decays off","p");
  t->AddEntry(pton,"Decays on","p");
  t->Draw("same");
  
  // Do not draw the Y axis label on the upper plot and redraw a small
  // axis instead, in order to avoid the first label (0) to be clipped.
  /*  ptoff->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
  axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis->SetLabelSize(15);
  axis->Draw();
  */
  // lower plot will be in pad
  cpt->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  // Define the ratio plot
  TH1D *h3 = (TH1D*)ptoff->Clone("h3");
  h3->SetLineColor(kBlack);
  h3->SetMinimum(0.8);  // Define Y ..
  h3->SetMaximum(1.2); // .. range
  h3->Sumw2();
  h3->SetStats(0);      // No statistics on lower plot
  h3->Divide(pton);
  h3->SetMarkerStyle(21);
  h3->Draw("ep");       // Draw the ratio plot

  // ptoff settings
  ptoff->SetLineColor(kBlue+1);
  ptoff->SetLineWidth(2);

  // Y axis ptoff plot settings
  ptoff->GetYaxis()->SetTitleSize(20);
  ptoff->GetYaxis()->SetTitleFont(43);
  ptoff->GetYaxis()->SetTitleOffset(1.55);

  // pton settings
  pton->SetLineColor(kRed);
  pton->SetLineWidth(2);

  // Ratio plot (h3) settings
  h3->SetTitle(""); // Remove the ratio title

  // Y axis ratio plot settings
  h3->GetYaxis()->SetTitle("Decays off / on ");
  h3->GetYaxis()->SetNdivisions(505);
  h3->GetYaxis()->SetTitleSize(20);
  h3->GetYaxis()->SetTitleFont(43);
  h3->GetYaxis()->SetTitleOffset(1.55);
  h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h3->GetYaxis()->SetLabelSize(15);

  // X axis ratio plot settings
  h3->GetXaxis()->SetTitleSize(20);
  h3->GetXaxis()->SetTitleFont(43);
  h3->GetXaxis()->SetTitleOffset(4.);
  h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h3->GetXaxis()->SetLabelSize(15);
  
  cpt->SaveAs("../plots/decays_correction/pt_ratio_decays_on_off.pdf");
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SPECTRA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  TCanvas *cspec = MakeCanvas("cspec","0",1200,800);
  DivideCanvas(cspec,"0",3,2);

  TLegend *tslices[5];
  double corresp_pts[5+1] = {15,20,25,30,40,60};
  for (int i = 0; i < 5; ++ i) {
    tslices[i] = SliceLegend(((to_string(corresp_pts[i])).substr(0,2) + " < p_{T}^{jet} < " + (to_string(corresp_pts[i + 1])).substr(0,2) + " GeV/c").c_str(), 0.13,0.8,0.9,0.95);
  }
  
  TLegend *titlepan = TitleLegend(0.15,0.15,0.9,0.9); titlepan->SetBorderSize(0);
  //TLegend *donoff = new TLegend(0.15,0.1,0.9,0.39); donoff->SetBorderSize(0);
  titlepan->AddEntry(off1520,"Decays off","p");
  titlepan->AddEntry(on1520,"Decays on","p");
  
  cspec->cd(1); /*TLatex *pan = PanelTitle();*/ titlepan->Draw("same");
  cspec->cd(2); on1520->Draw(); off1520->Draw("same"); tslices[0]->Draw("same");
  cspec->cd(3); on2025->Draw(); off2025->Draw("same"); tslices[1]->Draw("same");
  cspec->cd(4); on2530->Draw(); off2530->Draw("same"); tslices[2]->Draw("same");
  cspec->cd(5); on3040->Draw(); off3040->Draw("same"); tslices[3]->Draw("same");
  cspec->cd(6); on4060->Draw(); off4060->Draw("same"); tslices[4]->Draw("same");

  cspec->SaveAs("../plots/decays_correction/spectra_decays_on_off.pdf");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RATIO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  string pts[6] = {"15","20","25","30","40","60"};
  
  TCanvas *crat = MakeCanvas("crat","0",1200,800);
  DivideCanvas(crat,"0",3,2);

  TLegend *title = TitleLegend(0.15,0.15,0.9,0.9); title->SetBorderSize(0);

  TLine *one = new TLine(0,1,10,1);
  TLine *plus10 = new TLine(0,1.1,10,1.1);
  TLine *minus10 = new TLine(0,0.9,10,0.9);
  
  one->SetLineStyle(kDashed); plus10->SetLineStyle(kDashed); minus10->SetLineStyle(kDashed);
  
  TH1D* hdummy = new TH1D("hdummy",";M_{jet} [GeV/c^{2}];decays off / on",1,0,10);
  hdummy->GetYaxis()->SetRangeUser(0,2);
  
  TLatex* slice = new TLatex(); TLatex* p = new TLatex();
  
  crat->cd(1); hdummy->Draw(); p = PanelTitle();
  crat->cd(2); ratio1520->Draw(); one->Draw("same"); slice->DrawLatexNDC(0.3,0.2,(pts[0]+" < p_{T} < "+pts[1]+" GeV/c").c_str());
  crat->cd(3); ratio2025->Draw(); one->Draw("same"); slice->DrawLatexNDC(0.3,0.2,(pts[1]+" < p_{T} < "+pts[2]+" GeV/c").c_str());
  crat->cd(4); ratio2530->Draw(); one->Draw("same"); slice->DrawLatexNDC(0.3,0.2,(pts[2]+" < p_{T} < "+pts[3]+" GeV/c").c_str());
  crat->cd(5); ratio3040->Draw(); one->Draw("same"); slice->DrawLatexNDC(0.3,0.2,(pts[3]+" < p_{T} < "+pts[4]+" GeV/c").c_str());
  crat->cd(6); ratio4060->Draw(); one->Draw("same"); slice->DrawLatexNDC(0.3,0.2,(pts[4]+" < p_{T} < "+pts[5]+" GeV/c").c_str());
  
  crat->SaveAs("../plots/decays_correction/ratio_decays_on_off.pdf");
  
  return;
}
