#include "Plots.h"

using namespace std;

//plot raw mass, geant mass, and pythia6 mass for both groomed and ungroomed cases
//along with ratios between MC and data.
void rawmass(TFile *f, string obs) {

  TH2D* p = (TH2D*) f->Get((obs+"_v_pt_p").c_str());
  TH2D* g = (TH2D*) f->Get((obs+"_v_pt_g").c_str());
  TH2D* d = (TH2D*) f->Get((obs+"_v_pt_d").c_str());
  
  TH1D* px = (TH1D*) p->ProjectionX(("px"+obs).c_str(),p->GetYaxis()->FindBin(25),p->GetYaxis()->FindBin(30));
  TH1D* gx = (TH1D*) g->ProjectionX(("gx"+obs).c_str(),g->GetYaxis()->FindBin(25),g->GetYaxis()->FindBin(30));
  TH1D* dx = (TH1D*) d->ProjectionX(("dx"+obs).c_str(),d->GetYaxis()->FindBin(25),d->GetYaxis()->FindBin(30));
  
  string title = "M";
  if (obs == "mg") { title = "M_{g}"; }
  
  px->Scale(1/(double)px->Integral()); gx->Scale(1/(double)gx->Integral()); dx->Scale(1/(double)dx->Integral());
  
  TH1D* ratpd = (TH1D*) px->Clone(("ratpd"+obs).c_str());
  TH1D* ratgd = (TH1D*) gx->Clone(("ratgd"+obs).c_str());
  TH1D* ratpd_marks = (TH1D*) px->Clone(("ratpd_marks"+obs).c_str());
  
  Prettify1DwLineStyle(px,kBlue,kSolid,3,(title+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+title).c_str(), 0,10,0,0.5);
  Prettify1D(gx,kBlue,kOpenCircle,3,kBlue,(title+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+title).c_str(), 0,10,0,0.5);
  Prettify1D(dx,kBlack,kFullStar,4,kBlack,(title+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+title).c_str(), 0,10,0,0.5);
  
  px->GetYaxis()->SetTitleSize(30);
  px->GetYaxis()->SetTitleFont(43);
  px->GetYaxis()->SetTitleOffset(1.55);
  
  TCanvas *c = new TCanvas(("c"+obs).c_str(),("c"+obs).c_str(),800,1000);
  TPad *pad1 = new TPad((obs+"1").c_str(),"pad1",0,0.3,1,1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  px->SetStats(0); gx->SetStats(0); dx->SetStats(0);
  px->Draw("c"); gx->Draw("same"); dx->Draw("same");
  
  TLatex *l = new TLatex();
  TLegend * leg = new TLegend(0.4,0.4,0.65,0.65); leg->SetBorderSize(0);
  leg->AddEntry(px, "PYTHIA6","l");
  leg->AddEntry(gx, "PYTHIA6+GEANT","p");
  leg->AddEntry(dx, "Raw STAR data","p");
  
  if (obs == "m") {
    l->DrawLatexNDC(0.15,0.8,"pp #sqrt{s_{NN}} = 200 GeV, JP2");
    l->DrawLatexNDC(0.15,0.75,"anti-k_{t} R = 0.4 jets, |#eta| < 1 - R_{jet}");
  }
  else {
    l->DrawLatexNDC(0.2,0.8,"20 < p_{T}^{jet} < 25 GeV/c");
    l->DrawLatexNDC(0.2,0.75,"SoftDrop z_{cut} = 0.1, #beta = 0");
  }
  
  if (obs == "mg") {leg->Draw("same");}

  c->cd();
  
  TPad *pad2 = new TPad((obs+"2").c_str(),"pad2",0,0.05,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();
  
  ratpd->SetMinimum(0); ratpd->SetMaximum(2);
  ratgd->SetMinimum(0); ratgd->SetMaximum(2);
  ratpd_marks->SetMinimum(0); ratpd_marks->SetMaximum(2);
  ratgd->Sumw2(); ratpd_marks->Sumw2();
  ratgd->SetStats(0); ratpd->SetStats(0); ratpd_marks->SetStats(0);
  ratpd->Divide(dx); ratgd->Divide(dx); ratpd_marks->Divide(dx);
  ratpd->Sumw2(0);
  
  ratpd->SetLineColor(kBlue); ratpd->SetLineWidth(3);
  ratpd->GetXaxis()->SetTitle((title+" [GeV/c^{2}]").c_str());
  ratpd->GetYaxis()->SetTitle(("1/N dN/d"+title).c_str());
  ratpd->GetXaxis()->SetRangeUser(0,10);

  ratgd->SetLineColor(kBlue); ratgd->SetMarkerColor(kBlue); ratgd->SetMarkerStyle(kOpenCircle); ratgd->SetMarkerSize(3);
  ratgd->GetXaxis()->SetTitle((title+" [GeV/c^{2}]").c_str());
  ratgd->GetYaxis()->SetTitle(("1/N dN/d"+title).c_str());
  ratgd->GetXaxis()->SetRangeUser(0,10);
  
  ratpd->Draw("][same"); ratgd->Draw("epsame"); ratpd_marks->Draw("esame");

  TLine *one = new TLine(0,1,10,1); one->SetLineStyle(kDashed);
  one->Draw("same");
  
  ratpd->SetTitle(""); ratgd->SetTitle("");
  
  ratpd->GetYaxis()->SetTitle("MC / Data"); ratgd->GetYaxis()->SetTitle("MC / Data");
  
  ratpd->GetYaxis()->SetNdivisions(505);
  ratpd->GetYaxis()->SetTitleSize(30);
  ratpd->GetYaxis()->SetTitleFont(43);
  ratpd->GetYaxis()->SetTitleOffset(1.55);
  ratpd->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratpd->GetYaxis()->SetLabelSize(15);
  
  ratpd->GetXaxis()->SetTitleSize(30);
  ratpd->GetXaxis()->SetTitleFont(43);
  ratpd->GetXaxis()->SetTitleOffset(4.);
  ratpd->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratpd->GetXaxis()->SetLabelSize(15);

  string filename = "mass"; if (obs == "mg") {filename = "mg";}
  c->SaveAs(("~/jetmass/plots/paper/raw_"+filename+".pdf").c_str());
  
}

//plot the pT response and the mean particle-level pT for a given detector-level bin
void response(TFile * fres) {
  //  gStyle->SetPalette(kBlackBody);

  RooUnfoldResponse* res = (RooUnfoldResponse*) fres->Get("pt_res_coarse");
  TH2D* hres = (TH2D*) res->Hresponse();
  
  hres->GetXaxis()->SetTitle("p_{T}^{det jet} [GeV/c]");
  hres->GetYaxis()->SetTitle("p_{T}^{part jet} [GeV/c]");
  hres->GetZaxis()->SetTitle("d#sigma /dp_{T} [mb]");

  hres->GetZaxis()->SetRangeUser(3e-11,2e-5);
  
  hres->GetXaxis()->SetTitleOffset(1.2);
  hres->GetXaxis()->SetTitleSize(0.04);
  hres->GetYaxis()->SetTitleOffset(1.2);
  hres->GetYaxis()->SetTitleSize(0.04);
  hres->GetZaxis()->SetTitleOffset(1.2);
  hres->GetZaxis()->SetTitleSize(0.04);
  
  hres->SetTitle("");

  TCanvas * c = new TCanvas("cres","cres",1000,800); c->SetLogz();
  gPad->SetBottomMargin(0.11);
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  
  hres->Draw("colz");

  TLatex * l = new TLatex();

  l->DrawLatexNDC(0.25,0.8,"pp #sqrt{s_{NN}} = 200 GeV");
  l->DrawLatexNDC(0.25,0.75,"PYTHIA6+GEANT");
  l->DrawLatexNDC(0.25,0.7,"anti-k_{t} R = 0.4 jets, |#eta| < 1 - R_{jet}");

  TLine *diag = new TLine(15,15,60,60); diag->SetLineStyle(kDashed); diag->SetLineWidth(2);
  diag->Draw("same");
  
  vector<double> means = {}; vector<double> rmss = {};
  
  TProfile * p = (TProfile*) hres->ProfileX("p",1,-1/*hres->GetNbinsX()*/, "sdsame");

  p->SetMarkerStyle(kFullStar); p->SetMarkerColor(kBlack); p->SetLineColor(kBlack);
  p->SetMarkerSize(3);

  TLegend *pleg = new TLegend(0.6,0.15,0.75,0.2); pleg->SetBorderSize(0);
  pleg->AddEntry(p,"<p_{T}^{part}> #pm RMS","p");
  pleg->SetTextSize(0.03);

  p->Draw("same");
  pleg->Draw("same");
  
  c->SaveAs("~/jetmass/plots/paper/pt_response.pdf");
  return;
}

void resolutions(TFile* fnodrop, string obs) {

  TH2D* res2D = (TH2D*) fnodrop->Get(("ratio"+obs+"vPyPt").c_str());
  
  TH1D* rat2025 = (TH1D*) res2D->ProjectionX(("ratio_2025"+obs).c_str(),res2D->GetYaxis()->FindBin(20),res2D->GetYaxis()->FindBin(25));
  TH1D* rat2530 = (TH1D*) res2D->ProjectionX(("ratio_2530"+obs).c_str(),res2D->GetYaxis()->FindBin(25),res2D->GetYaxis()->FindBin(30));
  TH1D* rat3040 = (TH1D*) res2D->ProjectionX(("ratio_3040"+obs).c_str(),res2D->GetYaxis()->FindBin(30),res2D->GetYaxis()->FindBin(40));

  rat2025->Scale(1/(double)rat2025->Integral());
  rat2530->Scale(1/(double)rat2530->Integral());
  rat3040->Scale(1/(double)rat3040->Integral());
  
  rat2025->SetMarkerStyle(kOpenCircle); rat2530->SetMarkerStyle(kOpenSquare); rat3040->SetMarkerStyle(kOpenCross); 
  rat2025->SetMarkerColor(kBlue); rat2530->SetMarkerColor(kMagenta); rat3040->SetMarkerColor(kGreen+2);
  rat2025->SetLineColor(kBlue); rat2530->SetLineColor(kMagenta); rat3040->SetLineColor(kGreen+2);
  rat2025->SetMarkerSize(3); rat2530->SetMarkerSize(3); rat3040->SetMarkerSize(3);

  string titlestring = "M^{det} / M^{part}";
  if (obs == "Mg") { titlestring = "M_{g}^{det} / M_{g}^{part}";}
  rat2025->GetXaxis()->SetTitle(titlestring.c_str()); rat2025->GetYaxis()->SetTitle("Normalized integral");
  
  rat2025->GetYaxis()->SetRangeUser(0,0.2); rat2025->GetXaxis()->SetRangeUser(0,2);
  
  rat2025->GetXaxis()->SetTitleOffset(1.2);
  rat2025->GetXaxis()->SetTitleSize(0.04);
  rat2025->GetYaxis()->SetTitleOffset(1.2);
  rat2025->GetYaxis()->SetTitleSize(0.04);

  
  TCanvas *c = new TCanvas(("c"+obs+"res").c_str(),"cmres",1000,800); c->cd();
  gPad->SetBottomMargin(0.11);
  
  rat2025->Draw(); rat2530->Draw("same"); rat3040->Draw("same");
  
  TLegend *tpts = new TLegend(0.51,0.58,0.8,0.85); tpts->SetBorderSize(0);
  tpts->AddEntry(rat2025,"20 < p^{part jet}_{T} < 25 GeV/c","p");
  tpts->AddEntry(rat2530,"25 < p^{part jet}_{T} < 30 GeV/c","p");
  tpts->AddEntry(rat3040,"30 < p^{part jet}_{T} < 40 GeV/c","p");

  tpts->SetTextSize(0.038); 

  TLatex *sd = new TLatex(); sd->SetTextSize(0.038);
  
  if (obs == "M") {  tpts->Draw("same");}
  if (obs == "Mg") { sd->DrawLatexNDC(0.51,0.7,"SoftDrop z_{cut} = 0.1, #beta = 0");}
  TLine *one = new TLine (1,0,1,0.2); one->SetLineStyle(kDashed);
  
  one->Draw("same");
  
  c->SaveAs(("~/jetmass/plots/paper/resolution_"+obs+".pdf").c_str());
  
  return;
}

void paperplots() {
  
  //Files
  TFile * f = new TFile("~/jetmass/macros/hists/hists_R04.root","READ");
  TFile * fnodrop = new TFile("~/jetmass/macros/hists/hists_w_o_bin_drop_R04.root","READ");
  TFile * fres = new TFile("~/jetmass/out/matching/full_w_o_bin_drop_R04.root","READ");
  
  //Function calls
  //  rawmass(f, "m");
  //  rawmass(f, "mg");
  //  response(fres);
  //  resolutions(fnodrop,"M"); 
  //  resolutions(fnodrop,"Mg");

  return;
}
