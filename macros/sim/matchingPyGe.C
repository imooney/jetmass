#include <string>
#include <iostream>

using namespace std;

void matchingPyGe () {
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  string input = "~/jetmass/out/sim/full.root";
  string out_path = "~/jetmass/plots/sim/mass_resolution/";
  string filetype = ".pdf";

  TFile *f = new TFile(input.c_str(), "READ");
  TTree *matched = (TTree*) f->Get("py_ge_matchedTree");
  TTree *ch_matched = (TTree*) f->Get("ch_py_ge_matchedTree");
  
  double deltaM, deltaPt, ch_deltaM, ch_deltaPt, weight, ch_weight; // should I be using geant or pythia weight??
  matched->SetBranchAddress("DeltaM", &deltaM);
  matched->SetBranchAddress("DeltaPt", &deltaPt);
  matched->SetBranchAddress("ge_weight", &weight);
  ch_matched->SetBranchAddress("DeltaM", &ch_deltaM);
  ch_matched->SetBranchAddress("DeltaPt", &ch_deltaPt);
  ch_matched->SetBranchAddress("ge_weight", &ch_weight);

  TH1D * fulldeltaM = new TH1D("fulldeltaM","", 30,-15,15);
  TH1D * chdeltaM = new TH1D("chdeltaM","",30,-15,15);
  TH1D * fulldeltaPt = new TH1D("fulldeltaPt","", 140,-70,70);
  TH1D * chdeltaPt = new TH1D("chdeltaPt","",140,-70,70);
  
  TCanvas * cM = new TCanvas("cM","cM",800,800); cM->SetLogy();
  TCanvas * cPt = new TCanvas("cPt","cPt",800,800); cPt->SetLogy();  

  for (int i = 0; i < matched->GetEntries(); ++ i) {
    matched->GetEntry(i);
    fulldeltaM->Fill(deltaM, weight); fulldeltaPt->Fill(deltaPt, weight);
  }

  for (int i = 0; i < ch_matched->GetEntries(); ++ i) {
    ch_matched->GetEntry(i);
    chdeltaM->Fill(ch_deltaM, ch_weight); chdeltaPt->Fill(ch_deltaPt, ch_weight);
  }

  fulldeltaM->SetLineColor(kRed); chdeltaM->SetLineColor(kBlue);
  fulldeltaM->SetMarkerColor(kRed); chdeltaM->SetMarkerColor(kBlue);
  fulldeltaM->SetMarkerStyle(20); chdeltaM->SetMarkerStyle(21);
  fulldeltaM->GetXaxis()->SetTitle("#Delta M [GeV/c^{2}]"); fulldeltaM->GetYaxis()->SetTitle("prob.");
  chdeltaM->GetXaxis()->SetTitle("#Delta M [GeV/c^{2}]"); chdeltaM->GetYaxis()->SetTitle("prob.");
  fulldeltaM->Scale(1/(double) fulldeltaM->Integral()); chdeltaM->Scale(1/(double) chdeltaM->Integral());
  fulldeltaM->GetYaxis()->SetTitleOffset(1.3); chdeltaM->GetYaxis()->SetTitleOffset(1.3);
  
  fulldeltaPt->SetLineColor(kRed); chdeltaPt->SetLineColor(kBlue);
  fulldeltaPt->SetMarkerColor(kRed); chdeltaPt->SetMarkerColor(kBlue);
  fulldeltaPt->SetMarkerStyle(20); chdeltaPt->SetMarkerStyle(21);
  fulldeltaPt->GetXaxis()->SetTitle("#Delta p_{T} [GeV/c]"); fulldeltaPt->GetYaxis()->SetTitle("prob.");
  chdeltaPt->GetXaxis()->SetTitle("#Delta p_{T} [GeV/c]"); chdeltaPt->GetYaxis()->SetTitle("prob.");
  fulldeltaPt->Scale(1/(double) fulldeltaPt->Integral()); chdeltaPt->Scale(1/(double) chdeltaPt->Integral());
  fulldeltaPt->GetYaxis()->SetTitleOffset(1.3); chdeltaPt->GetYaxis()->SetTitleOffset(1.3); 
  
  TLegend *tM = new TLegend(0.55,0.7,0.86,0.8); tM->SetBorderSize(0);
  tM->AddEntry((TObject*)0, "Pythia+Geant - Pythia jet mass resolution","");
  tM->AddEntry((TObject*)0, "Leading Pythia jet matched to Pythia+Geant jet","");
  tM->AddEntry((TObject*)0, "anti-kT, R = 0.4, pp 200 GeV, HT trigger","");
  tM->AddEntry(fulldeltaM,"Ch+Ne","p"); tM->AddEntry(chdeltaM,"Ch","p");

  TLegend *tPt = new TLegend(0.55,0.7,0.86,0.8); tPt->SetBorderSize(0);
  tPt->AddEntry((TObject*)0, "Pythia+Geant - Pythia jet p_{T} resolution","");
  tPt->AddEntry((TObject*)0, "Leading Pythia jet matched to Pythia+Geant jet","");
  tPt->AddEntry((TObject*)0, "anti-kT, R = 0.4, pp 200 GeV, HT trigger","");
  tPt->AddEntry(fulldeltaPt,"Ch+Ne","p"); tPt->AddEntry(chdeltaPt,"Ch","p");

  std::cout << "charged pt mean: " << chdeltaPt->GetMean() << " charged pt RMS: " << chdeltaPt->GetRMS() << std::endl;
  std::cout << "full pt mean: " << fulldeltaPt->GetMean() << " full pt RMS: " << fulldeltaPt->GetRMS() << std::endl;
  std::cout << "charged m mean: " << chdeltaM->GetMean() << " charged m RMS: " << chdeltaM->GetRMS() << std::endl;
  std::cout << "full m mean: " << fulldeltaM->GetMean() << " full m RMS: " << fulldeltaM->GetRMS() << std::endl;

 
  cM->cd(); fulldeltaM->Draw("same"); chdeltaM->Draw("same"); tM->Draw("same");
  cPt->cd(); fulldeltaPt->Draw("same"); chdeltaPt->Draw("same"); tPt->Draw("same");
  /*
  cM->SaveAs((out_path + "deltaM" + filetype).c_str());
  cPt->SaveAs((out_path + "deltaPt" + filetype).c_str());
  */
  return;
}
