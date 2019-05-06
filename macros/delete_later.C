
#include "Plots.h"

using namespace std;

void delete_later() {
  
  bool mat = 1;

  TH1D* p; TH1D* g_wo;
  
  if (!mat) {
    TFile *fpy = new TFile ("~/jetmass/out/sim/py/full_w_o_bin_drop_R06.root","READ");
    //TFile *fgew = new TFile ("~/jetmass/out/sim/ge/full_w_o_bin_drop_R06.root","READ");
    TFile *fgewo = new TFile ("~/jetmass/out/sim/ge/full_w_o_bin_drop_w_o_pt_cut_R06.root","READ");
    
    
    g_wo = new TH1D("g_wo","",15,5,80);
    //  TH1D* g_w = new TH1D("g_w","",15,5,80);
    p = new TH1D("p","",15,5,80);
    
    TTree *tpy = (TTree*) fpy->Get("event");
    //TTree *tgew = (TTree*) fgew->Get("event");
    TTree *tgewo = (TTree*) fgewo->Get("event");
    
    tpy->Draw("Pt>>p","weight","same");
    //tgew->Draw("Pt>>g_w","weight","same");
    tgewo->Draw("Pt>>g_wo","weight","same");
  }
  
  if (mat) { 
    TFile *fmat_wo = new TFile("~/jetmass/out/matching/full_w_o_bin_drop_w_o_pt_cut_R06.root","READ");
    
    RooUnfoldResponse* res = (RooUnfoldResponse*) fmat_wo->Get("pt_response");
    
    TH2D* hres = (TH2D*) res->Hresponse();
    
    p = (TH1D*) hres->ProjectionY("p",hres->GetXaxis()->FindBin(5),hres->GetXaxis()->FindBin(80));
    g_wo = (TH1D*) hres->ProjectionX("g_wo",hres->GetYaxis()->FindBin(5),hres->GetYaxis()->FindBin(80));
  }
  
  Prettify1D(p,kBlue,kOpenCircle,2,kBlue,"p_{T} [GeV/c]","arb.",5,80,1e-12,1e-1);
  //Prettify1D(g_w,kRed,kOpenCircle,2,kRed,"p_{T} [GeV/c]","arb.",5,80,1e-12,1e-1);
  Prettify1D(g_wo,kRed,kFullCircle,2,kRed,"p_{T} [GeV/c]","arb.",5,80,1e-12,1e-1);

  //p->Scale(/(double) p->GetBinContent(3));
  //g_w->Scale(p->GetBinContent(4) /(double) g_w->GetBinContent(4));
  g_wo->Scale(p->GetBinContent(4)/(double) g_wo->GetBinContent(4));

  TCanvas *c = new TCanvas("c","c",1400,1000); c->cd(); gPad->SetLogy();
  
  TLegend *t = new TLegend(0.5,0.6,0.8,0.8); t->SetBorderSize(0);
  t->AddEntry(p,"Pythia","p");
  t->AddEntry(g_wo,"Geant","p");
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad	
  pad1->SetLogy();
  p->SetStats(0);          // No statistics on upper plot
  p->Draw();               // Draw p
  g_wo->Draw("same");         // Draw g_wo on top of p
  t->Draw("same");
  
  // Do not draw the Y axis label on the upper plot and redraw a small
  // axis instead, in order to avoid the first label (0) to be clipped.
  /*  p->GetYaxis()->SetLabelSize(0.);
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
  TH1F *rat = (TH1F*)p->Clone("rat");
  rat->SetLineColor(kBlack);
  rat->SetMinimum(0);  // Define Y ..
  rat->SetMaximum(2); // .. range
  rat->Sumw2();
  rat->SetStats(0);      // No statistics on lower plot
  rat->Divide(g_wo);
  rat->SetMarkerStyle(21);
  rat->Draw("ep");       // Draw the ratio plot

  TLine *one = new TLine(5,1,80,1); one->SetLineStyle(kDashed);
  one->Draw("same");
  
  // p settings
  p->SetLineColor(kBlue+1);
  p->SetLineWidth(2);

  // Y axis p plot settings
  p->GetYaxis()->SetTitleSize(20);
  p->GetYaxis()->SetTitleFont(43);
  p->GetYaxis()->SetTitleOffset(1.55);

  // g_wo settings
  g_wo->SetLineColor(kRed);
  g_wo->SetLineWidth(2);

  // Ratio plot (rat) settings
  rat->SetTitle(""); // Remove the ratio title

  // Y axis ratio plot settings
  rat->GetYaxis()->SetTitle("gen. / det.");
  rat->GetYaxis()->SetNdivisions(505);
  rat->GetYaxis()->SetTitleSize(20);
  rat->GetYaxis()->SetTitleFont(43);
  rat->GetYaxis()->SetTitleOffset(1.55);
  rat->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  rat->GetYaxis()->SetLabelSize(15);

  // X axis ratio plot settings
  rat->GetXaxis()->SetTitleSize(20);
  rat->GetXaxis()->SetTitleFont(43);
  rat->GetXaxis()->SetTitleOffset(4.);
  rat->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  rat->GetXaxis()->SetLabelSize(15);
  
  //p->Draw(); g_wo->Draw("same"); t->Draw("same");
  if (!mat) {
    c->SaveAs("~/spectra_nodetcut_R06.pdf");
  }
  if (mat) {
    c->SaveAs("~/matched_spectra_nodetcut_R06.pdf");
  }

  return;
}
