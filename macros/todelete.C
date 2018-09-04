void todelete () {

  string dir = "~/jetmass/";
  string pyin = "out/sim/py/";
  string file = "full.root";
  string out = "~/jetmass/plots/tests/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";
  
  TFile* pyFile = new TFile( (dir + pyin +file).c_str(), "READ");                                                                                
  
  TH1D* pt23 = (TH1D*) pyFile->Get("pt23");
  TH1D* pt34 = (TH1D*) pyFile->Get("pt34");
  TH1D* pt45 = (TH1D*) pyFile->Get("pt45");
  TH1D* pt57 = (TH1D*) pyFile->Get("pt57");
  TH1D* pt79 = (TH1D*) pyFile->Get("pt79");
  TH1D* pt911 = (TH1D*) pyFile->Get("pt911");
  TH1D* pt1115 = (TH1D*) pyFile->Get("pt1115");
  TH1D* pt1520 = (TH1D*) pyFile->Get("pt1520");
  TH1D* pt2025 = (TH1D*) pyFile->Get("pt2025");
  TH1D* pt2535 = (TH1D*) pyFile->Get("pt2535");
  TH1D* pt35plus = (TH1D*) pyFile->Get("pt35plus");
 
  TH1D* pt23w = (TH1D*) pyFile->Get("pt23w");
  TH1D* pt34w = (TH1D*) pyFile->Get("pt34w");
  TH1D* pt45w = (TH1D*) pyFile->Get("pt45w");
  TH1D* pt57w = (TH1D*) pyFile->Get("pt57w");
  TH1D* pt79w = (TH1D*) pyFile->Get("pt79w");
  TH1D* pt911w = (TH1D*) pyFile->Get("pt911w");
  TH1D* pt1115w = (TH1D*) pyFile->Get("pt1115w");
  TH1D* pt1520w = (TH1D*) pyFile->Get("pt1520w");
  TH1D* pt2025w = (TH1D*) pyFile->Get("pt2025w");
  TH1D* pt2535w = (TH1D*) pyFile->Get("pt2535w");
  TH1D* pt35plusw = (TH1D*) pyFile->Get("pt35plusw");

  TH1D* ptall = (TH1D*) pyFile->Get("ptall");
  TH1D* ptallw = (TH1D*) pyFile->Get("ptallw");
  
  pt23->SetMarkerStyle(20); pt34->SetMarkerStyle(21);
  pt45->SetMarkerStyle(22); pt57->SetMarkerStyle(23);
  pt79->SetMarkerStyle(24); pt911->SetMarkerStyle(25);
  pt1115->SetMarkerStyle(26); pt1520->SetMarkerStyle(27);
  pt2025->SetMarkerStyle(28); pt2535->SetMarkerStyle(29);
  pt35plus->SetMarkerStyle(30);

  pt23->SetMarkerColor(1); pt34->SetMarkerColor(2); pt45->SetMarkerColor(3);
  pt57->SetMarkerColor(4); pt79->SetMarkerColor(13); pt911->SetMarkerColor(6);
  pt1115->SetMarkerColor(7); pt1520->SetMarkerColor(8); pt2025->SetMarkerColor(9);
  pt2535->SetMarkerColor(11); pt35plus->SetMarkerColor(12);

  pt23w->SetMarkerStyle(20); pt34w->SetMarkerStyle(21);
  pt45w->SetMarkerStyle(22); pt57w->SetMarkerStyle(23);
  pt79w->SetMarkerStyle(24); pt911w->SetMarkerStyle(25);
  pt1115w->SetMarkerStyle(26); pt1520w->SetMarkerStyle(27);
  pt2025w->SetMarkerStyle(28); pt2535w->SetMarkerStyle(29);
  pt35plusw->SetMarkerStyle(30);

  pt23w->SetMarkerColor(1); pt34w->SetMarkerColor(2); pt45w->SetMarkerColor(3);
  pt57w->SetMarkerColor(4); pt79w->SetMarkerColor(13); pt911w->SetMarkerColor(6);
  pt1115w->SetMarkerColor(7); pt1520w->SetMarkerColor(8); pt2025w->SetMarkerColor(9);
  pt2535w->SetMarkerColor(11); pt35plusw->SetMarkerColor(12);

  ptall->SetMarkerColor(kRed); ptall->SetMarkerStyle(kFullSquare);
  ptallw->SetMarkerColor(kRed); ptallw->SetMarkerStyle(kFullSquare);
  
  pt23->GetYaxis()->SetRangeUser(1e-1,1e6);
  pt23w->GetYaxis()->SetRangeUser(1e-11,1e-1);
  pt34->GetYaxis()->SetRangeUser(1e-1,1e6);
  pt34w->GetYaxis()->SetRangeUser(1e-11,1e-1);
  //gStyle->SetPalette(kSunset);
  
  TCanvas * c = new TCanvas ("c","c",800,800); c->SetLogy();
  TCanvas * cw = new TCanvas ("cw","cw",800,800); cw->SetLogy();
  TLegend * t = new TLegend (0.6,0.6,0.85,0.85); t->SetBorderSize(0);
  t->AddEntry(pt23,"2 - 3", "p");
  t->AddEntry(pt34,"3 - 4", "p");
  t->AddEntry(pt45,"4 - 5", "p");
  t->AddEntry(pt57,"5 - 7", "p");
  t->AddEntry(pt79,"7 - 9", "p");
  t->AddEntry(pt911,"9 - 11", "p");
  t->AddEntry(pt1115,"11 - 15", "p");
  t->AddEntry(pt1520,"15 - 20", "p");
  t->AddEntry(pt2025,"20 - 25", "p");
  t->AddEntry(pt2535,"25 - 35", "p");
  t->AddEntry(pt35plus,"35 +", "p");
  t->AddEntry(ptall,"all pT", "p");
  c->cd();
  pt23->Draw(); pt34->Draw("same"); pt45->Draw("same"); pt57->Draw("same"); pt79->Draw("same"); pt911->Draw("same"); pt1115->Draw("same"); pt1520->Draw("same"); pt2025->Draw("same"); pt2535->Draw("same"); pt35plus->Draw("same"); ptall->Draw("same"); t->Draw("same");
  cw->cd();
  pt23w->Draw(); pt34w->Draw("same"); pt45w->Draw("same"); pt57w->Draw("same"); pt79w->Draw("same"); pt911w->Draw("same"); pt1115w->Draw("same"); pt1520w->Draw("same"); pt2025w->Draw("same"); pt2535w->Draw("same"); pt35plusw->Draw("same"); ptallw->Draw("same"); t->Draw("same");

  c->SaveAs((out + "jetpt_by_pthardbin_unweighted_wo2_3" + filetype).c_str());
  cw->SaveAs((out + "jetpt_by_pthardbin_weighted_wo2_3" + filetype).c_str());
}
