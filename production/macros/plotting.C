#include "Plots.h"

using namespace std;

void plotting () {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  const string dir = "~/jetmass/production/Results/";
  const string out = "~/jetmass/production/plots/";
  string flag = "30";
  const string file = ("pythia8_recSD_histograms_R02_max" + flag + ".root").c_str();
  const string filetype = ".pdf";
  
  TFile *f = new TFile( (dir + file).c_str(),"READ");
  
  TH2D *pt_trig_v_rec = (TH2D*) f->Get("hrecpTvtrigpT");
  
  //2D projections!
  TCanvas * crec = MakeCanvas("crec","y",1000,800);
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->SetLogy();
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  //h1->SetStats(0);          // No statistics on upper plot
  //DivideCanvas(crec,"y",3,2);

  const int nBins = 3;
  double ranges[nBins+1] = {1,3,7,22}; double range_full[2] = {1,22}; double range_recs[2] = {1,81};
  string corresponding_pt[nBins+1] = {"9","11","15","30"};
  vector<TH1D*> slices = Projection2D( pt_trig_v_rec, nBins, ranges, "x");
  vector<TH1D*> whole = Projection2D( pt_trig_v_rec, 1, range_full, "x");

  vector<TH1D*> for_trigs = Projection2D(pt_trig_v_rec, 1, range_recs, "y");
  TAxis *axis = for_trigs[0]->GetXaxis();
  int b9 = axis->FindBin(9); int b11 = axis->FindBin(11); int b15 = axis->FindBin(15);
  int b30 = axis->FindBin(30);
  double ntrigs911 = for_trigs[0]->Integral(b9, b11);
  double ntrigs1115 = for_trigs[0]->Integral(b11,b15);
  double ntrigs1530 = for_trigs[0]->Integral(b15,b30); 
  double ntrigstot = ntrigs911 + ntrigs1115 + ntrigs1530;
  
  
  Prettify1D(slices[0], kBlue, kOpenCircle, 1, kBlue, "p^{rec. jet}_{T} [GeV/c]", "1/N_{trig.} dN/dp_{T}", 0,35,1e-4,5);
  Prettify1D(slices[1], kRed, kFullCircle, 1, kRed, "p^{rec. jet}_{T} [GeV/c]", "1/N_{trig.} dN/dp_{T}", 0,35,1e-4,5);
  Prettify1D(slices[2], kOrange, kOpenSquare, 1, kOrange, "p^{rec. jet}_{T} [GeV/c]", "1/N_{trig.} dN/dp_{T}", 0,35,1e-4,5);
  Prettify1D(whole[0], kViolet, kFullSquare, 1, kViolet, "p^{rec. jet}_{T} [GeV/c]", "1/N_{trig.} dN/dp_{T}", 0,35,1e-4,5);
  
  slices[0]->Scale(1/(double)ntrigs911, "width");
  slices[1]->Scale(1/(double)ntrigs1115, "width");
  slices[2]->Scale(1/(double)ntrigs1530, "width");
  whole[0]->Scale(1/(double)ntrigstot, "width");
  //binwidth is 1 everywhere so dont need to scale by it
    
  TLegend *leg = TitleLegend(0.5,0.6,0.88,0.85);
  /*TLegend *leg_slices[nBins + 1];
  for (int i = 0; i < slices.size(); ++ i) {
    leg_slices[i] = SliceLegend((corresponding_pt[i] + " < p^{trig.}_{T} < " + corresponding_pt[i+1] + " GeV/c").c_str(), 0.13,0.8,0.9,0.95);
  }
  leg_slices[nBins] = SliceLegend((corresponding_pt[0] + " < p^{trig.}_{T} < " + corresponding_pt[nBins] + " GeV/c").c_str(), 0.13, 0.8,0.9,0.95);
  */
  TLegend *leg_w_markers = new TLegend(0.5,0.35,0.88,0.59); leg_w_markers->SetBorderSize(0);
  for (int i = 0; i < slices.size(); ++ i) {
    leg_w_markers->AddEntry(slices[i],(corresponding_pt[i] + " < p^{trig.}_{T} < " + corresponding_pt[i+1] + " GeV/c").c_str(),"p");
    slices[i]->GetXaxis()->SetRangeUser(0,35); slices[i]->GetYaxis()->SetRangeUser(1e-4,5);
  }
  leg_w_markers->AddEntry(whole[0],(corresponding_pt[0]+" < p^{trig.}_{T} < " + corresponding_pt[nBins] + " GeV/c").c_str(),"p");
  whole[0]->GetXaxis()->SetRangeUser(0,35); whole[0]->GetYaxis()->SetRangeUser(1e-4,5);
  /*
  crec->cd(1); leg->Draw();
  for (int i = 0; i < slices.size(); ++ i) {
    slices[i]->GetXaxis()->SetRangeUser(0,40); slices[i]->GetYaxis()->SetRangeUser(1e-5,5);
    crec->cd(i+2); slices[i]->Draw(); leg_slices[i]->Draw("same");
  } 
  crec->cd(slices.size() + 2);
  whole[0]->GetXaxis()->SetRangeUser(0,40); whole[0]->GetYaxis()->SetRangeUser(1e-5,5);
  whole[0]->Draw(); leg_slices[nBins]->Draw("same");
  */
  
  slices[0]->Draw(); slices[1]->Draw("same"); slices[2]->Draw("same"); whole[0]->Draw("same"); leg->Draw("same"); leg_w_markers->Draw("same");
  //crec->cd(6); //leg->Draw();
  //crec->SaveAs((out+"recoil_spectra"+filetype).c_str());
  
  // lower plot will be in pad 2
  crec->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.223, 1, 0.35);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  
  TH1D *ratlowhigh = (TH1D*) slices[0]->Clone("ratlowhigh");
  ratlowhigh->SetLineColor(kBlack);
  ratlowhigh->SetMinimum(0);  // Define Y ..
  ratlowhigh->SetMaximum(2); // .. range
  ratlowhigh->Sumw2();
  ratlowhigh->SetStats(0);      // No statistics on lower plot
  ratlowhigh->Divide(slices[2]);
  ratlowhigh->SetMarkerStyle(21);
  ratlowhigh->Draw("ep");       // Draw the ratio plot
  
  TLine* one = new TLine(-5,1,38,1); one->SetLineStyle(kDashed);
  one->Draw("same");
  ratlowhigh->SetTitle("");
  
  // Y axis ratio plot settings
  ratlowhigh->GetYaxis()->SetTitle("low/high ");
  ratlowhigh->GetYaxis()->SetNdivisions(505);
  ratlowhigh->GetYaxis()->SetTitleSize(20);
  ratlowhigh->GetYaxis()->SetTitleFont(43);
  ratlowhigh->GetYaxis()->SetTitleOffset(1.55);
  ratlowhigh->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratlowhigh->GetYaxis()->SetLabelSize(15);

  // X axis ratio plot settings
  ratlowhigh->GetXaxis()->SetTitleSize(20);
  ratlowhigh->GetXaxis()->SetTitleFont(43);
  ratlowhigh->GetXaxis()->SetTitleOffset(4.);
  ratlowhigh->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratlowhigh->GetXaxis()->SetLabelSize(15);
  
  crec->cd();
  TPad *pad3 = new TPad("pad3", "pad3", 0, 0.05, 1, 0.223);
  pad3->SetTopMargin(0);
  pad3->SetBottomMargin(0.35);
  pad3->SetGridx(); // vertical grid
  pad3->Draw();
  pad3->cd();       // pad3 becomes the current pad
  
  TH1D *ratmidhigh = (TH1D*) slices[1]->Clone("ratmidhigh");
  ratmidhigh->SetLineColor(kBlack);
  ratmidhigh->SetMinimum(0);  // Define Y ..
  ratmidhigh->SetMaximum(2); // .. range
  ratmidhigh->Sumw2();
  ratmidhigh->SetStats(0);      // No statistics on lower plot
  ratmidhigh->Divide(slices[2]);
  ratmidhigh->SetMarkerStyle(21);
  ratmidhigh->Draw("ep");       // Draw the ratio plot
  one->Draw("same");

  ratmidhigh->SetTitle("");
  
  // Y axis ratio plot settings
  ratmidhigh->GetYaxis()->SetTitle("mid/high ");
  ratmidhigh->GetYaxis()->SetNdivisions(505);
  ratmidhigh->GetYaxis()->SetTitleSize(20);
  ratmidhigh->GetYaxis()->SetTitleFont(43);
  ratmidhigh->GetYaxis()->SetTitleOffset(1.55);
  ratmidhigh->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratmidhigh->GetYaxis()->SetLabelSize(15);

  // X axis ratio plot settings
  ratmidhigh->GetXaxis()->SetTitleSize(20);
  ratmidhigh->GetXaxis()->SetTitleFont(43);
  ratmidhigh->GetXaxis()->SetTitleOffset(4.);
  ratmidhigh->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratmidhigh->GetXaxis()->SetLabelSize(15);
  
  crec->SaveAs((out + "nihar_check_"+ flag + filetype).c_str());
  
  return;
}
