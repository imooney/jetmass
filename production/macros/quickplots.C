

using namespace std;

void quickplots () {

  TFile *f = new TFile ("~/jetmass/production/Results/pythia8_decays_off_and_noHad_R04_isaac.root", "READ");
  
  TH1D* h2025 = new TH1D("h2025","",14,0,14);
  TH1D* h2530 = new TH1D("h2530","",14,0,14);
  TH1D* h3040 = new TH1D("h3040","",14,0,14);
 
  TH1D* hm02025 = new TH1D("hm02025","",14,0,14);
  TH1D* hm02530 = new TH1D("hm02530","",14,0,14);
  TH1D* hm03040 = new TH1D("hm03040","",14,0,14);

 
  TTree *t = (TTree*) f->Get("PartonTree");
  TTree *tm0 = (TTree*) f->Get("m0PartonTree");
 
  t->Draw("PLm>>h2025","mcweight*(PLpt>20 && PLpt<25)");
  t->Draw("PLm>>h2530","mcweight*(PLpt>25 && PLpt<30)");
  t->Draw("PLm>>h3040","mcweight*(PLpt>30 && PLpt<40)");

  tm0->Draw("m0PLm>>hm02025","mcweight*(m0PLpt>20 && m0PLpt<25)");
  tm0->Draw("m0PLm>>hm02530","mcweight*(m0PLpt>25 && m0PLpt<30)");
  tm0->Draw("m0PLm>>hm03040","mcweight*(m0PLpt>30 && m0PLpt<40)");

  h2025->Scale(1/(double)h2025->Integral());
  h2530->Scale(1/(double)h2530->Integral());
  h3040->Scale(1/(double)h3040->Integral());
  hm02025->Scale(1/(double)hm02025->Integral());
  hm02530->Scale(1/(double)hm02530->Integral());
  hm03040->Scale(1/(double)hm03040->Integral());

  h2025->Divide(hm02025);
  h2530->Divide(hm02530);
  h3040->Divide(hm03040);

  h2025->SetMarkerStyle(kOpenCircle);
  h2530->SetMarkerStyle(kOpenCircle);
  h3040->SetMarkerStyle(kOpenCircle);

  h2025->SetMarkerColor(kBlue);
  h2530->SetMarkerColor(kBlue);
  h3040->SetMarkerColor(kBlue);

  h2025->SetMarkerSize(3);
  h2530->SetMarkerSize(3);
  h3040->SetMarkerSize(3);

  h2025->GetYaxis()->SetRangeUser(0,2);
  h2530->GetYaxis()->SetRangeUser(0,2);
  h3040->GetYaxis()->SetRangeUser(0,2);

  h2025->GetXaxis()->SetTitle("M^{PL}_{jet} [GeV/c^{2}]");
  h2530->GetXaxis()->SetTitle("M^{PL}_{jet} [GeV/c^{2}]");
  h3040->GetXaxis()->SetTitle("M^{PL}_{jet} [GeV/c^{2}]");

  h2025->GetYaxis()->SetTitle("PDG / massless");
  h2530->GetYaxis()->SetTitle("PDG / massless");
  h3040->GetYaxis()->SetTitle("PDG / massless");
  
  TLatex *l = new TLatex();

  TLine *one = new TLine(0,1,14,1); one->SetLineStyle(kDashed);  

  TCanvas *c = new TCanvas("c","c",1200,500);
  c->Divide(3,1,0,0);

  c->cd(1); h2025->Draw(); one->Draw("same"); l->DrawLatexNDC(0.3,0.2,"20 < p^{PL}_{T} < 25 GeV/c");
  c->cd(2); h2530->Draw(); one->Draw("same"); l->DrawLatexNDC(0.3,0.2,"25 < p^{PL}_{T} < 30 GeV/c");
  c->cd(3); h3040->Draw(); one->Draw("same"); l->DrawLatexNDC(0.3,0.2,"30 < p^{PL}_{T} < 40 GeV/c");
  
  c->SaveAs("~/jetmass/production/plots/partonmassassignment.pdf");

  return;
}
