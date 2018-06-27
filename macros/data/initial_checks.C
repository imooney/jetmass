//I started this project by plotting pt, eta, phi & mass
//for inclusive and leading charged+neutral jets.
//This is a plotting macro for these quantities

#include <string>
#include <iostream>

using namespace std;

void initial_checks () {
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  string input = "~/jetmass/out/data/full.root";
  string out_path = "~/jetmass/plots/data/initial_checks/";
  string filetype = ".pdf";

  TFile *f = new TFile(input.c_str(), "READ");

  //~~~~~~~~~~~~~~~~~~MASS~~~~~~~~~~~~~~~~~~~~~~~~

  TH1D * m_incl = (TH1D*)f->Get("m_inclusive");
  TH1D * m_lead = (TH1D*)f->Get("m_leading");

  m_incl->SetMarkerStyle(20);
  m_lead->SetMarkerStyle(21);
  m_incl->SetMarkerColor(kBlue);
  m_lead->SetMarkerColor(kRed);
  m_incl->SetLineColor(kBlue);
  m_lead->SetLineColor(kRed);

  m_incl->GetXaxis()->SetTitle("M^{jet}");
  m_incl->GetYaxis()->SetTitle("arb.");

  m_lead->Scale(1/(double)m_lead->Integral());
  m_incl->Scale(1/(double)m_incl->Integral());

  m_incl->GetXaxis()->SetRangeUser(0,11);
  m_lead->GetXaxis()->SetRangeUser(0,11);

  TLegend *tmass = new TLegend(0.7,0.7,0.9,0.8); tmass->SetBorderSize(0);
  tmass->AddEntry(m_lead, "Ch+Ne leading jets", "p");
  tmass->AddEntry(m_incl, "Ch+Ne inclusive jets","p");
  tmass->AddEntry((TObject*)0, "STAR pp Run6 200 GeV", "");
  tmass->AddEntry((TObject*)0, "anti-kt, R = 0.4", "");
  tmass->Draw();

  TCanvas * cmass = new TCanvas("cmass","cmass",800,800); cmass->cd(); cmass->SetLogy();
  m_incl->Draw(); m_lead->Draw("same"); tmass->Draw("same");
  
  cmass->SaveAs((out_path + "mass" + filetype).c_str());
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~3D HISTS~~~~~~~~~~~~~~~~~~~
  
  TH3D * j_PtEtaPhi_incl = (TH3D*)f->Get("PtEtaPhi_incl");
  TH3D * j_PtEtaPhi_lead = (TH3D*)f->Get("PtEtaPhi_lead");
  TH3D * c_PtEtaPhi_incl = (TH3D*)f->Get("cons_PtEtaPhi_incl");
  TH3D * c_PtEtaPhi_lead = (TH3D*)f->Get("cons_PtEtaPhi_lead");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   //marker, axis, scale - incl blue 20, lead red 21

  //~~~~~~~~~~~~~~~~~~~JET PHI~~~~~~~~~~~~~~~~~~~~
  TH1 * j_phi_incl = j_PtEtaPhi_incl->Project3D("z");
  TH1 * j_phi_lead = j_PtEtaPhi_lead->Project3D("z");
  
  j_phi_incl->SetTitle(""); j_phi_lead->SetTitle("");
  
  j_phi_incl->SetMarkerStyle(20);
  j_phi_lead->SetMarkerStyle(21);
  j_phi_incl->SetMarkerColor(kBlue);
  j_phi_lead->SetMarkerColor(kRed);
  j_phi_incl->SetLineColor(kBlue);
  j_phi_lead->SetLineColor(kRed);

  j_phi_incl->GetXaxis()->SetTitle("#phi^{jet}");
  j_phi_incl->GetYaxis()->SetTitle("arb.");
  
  j_phi_incl->Scale(1/(double)j_phi_incl->Integral());
  j_phi_lead->Scale(1/(double)j_phi_lead->Integral());

  j_phi_incl->GetYaxis()->SetRangeUser(0,0.2);
  j_phi_lead->GetYaxis()->SetRangeUser(0,0.2);

  TLegend *tjphi = new TLegend(0.7,0.7,0.9,0.8); tjphi->SetBorderSize(0);
  tjphi->AddEntry(j_phi_lead, "Ch+Ne leading jets", "p");
  tjphi->AddEntry(j_phi_incl, "Ch+Ne inclusive jets","p");

  TCanvas * cjphi = new TCanvas("cjphi","cjphi",800,800); cjphi->cd();
  j_phi_incl->Draw(); j_phi_lead->Draw("same"); tjphi->Draw("same");

  cjphi->SaveAs((out_path + "jet_phi" + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~JET ETA~~~~~~~~~~~~~~~~~~~~
  
  TH1 * j_eta_incl = j_PtEtaPhi_incl->Project3D("y");
  TH1 * j_eta_lead = j_PtEtaPhi_lead->Project3D("y");                                                                                                                                                    

  j_eta_incl->SetTitle(""); j_eta_lead->SetTitle("");

  j_eta_incl->SetMarkerStyle(20);
  j_eta_lead->SetMarkerStyle(21);
  j_eta_incl->SetMarkerColor(kBlue);
  j_eta_lead->SetMarkerColor(kRed);
  j_eta_incl->SetLineColor(kBlue);
  j_eta_lead->SetLineColor(kRed);

  j_eta_incl->GetXaxis()->SetTitle("#eta^{jet}");
  j_eta_incl->GetYaxis()->SetTitle("arb.");

  j_eta_incl->Scale(1/(double)j_eta_incl->Integral());
  j_eta_lead->Scale(1/(double)j_eta_lead->Integral());

  TLegend *tjeta = new TLegend(0.7,0.7,0.9,0.8); tjeta->SetBorderSize(0);
  tjeta->AddEntry(j_eta_lead, "Ch+Ne leading jets", "p");
  tjeta->AddEntry(j_eta_incl, "Ch+Ne inclusive jets","p");

  TCanvas * cjeta = new TCanvas("cjeta","cjeta",800,800); cjeta->cd();
  j_eta_incl->Draw(); j_eta_lead->Draw("same"); tjeta->Draw("same");

  cjeta->SaveAs((out_path + "jet_eta" + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~JET PT~~~~~~~~~~~~~~~~~~~~                                                                                                                                                         

  TH1 * j_pt_incl = j_PtEtaPhi_incl->Project3D("x");
  TH1 * j_pt_lead = j_PtEtaPhi_lead->Project3D("x");

  j_pt_incl->SetTitle(""); j_pt_lead->SetTitle("");

  j_pt_incl->SetMarkerStyle(20);
  j_pt_lead->SetMarkerStyle(21);
  j_pt_incl->SetMarkerColor(kBlue);
  j_pt_lead->SetMarkerColor(kRed);
  j_pt_incl->SetLineColor(kBlue);
  j_pt_lead->SetLineColor(kRed);

  j_pt_incl->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  j_pt_incl->GetYaxis()->SetTitle("arb.");

  j_pt_incl->Scale(1/(double)j_pt_incl->Integral());
  j_pt_lead->Scale(1/(double)j_pt_lead->Integral());

  TLegend *tjpt = new TLegend(0.7,0.7,0.9,0.8); tjpt->SetBorderSize(0);
  tjpt->AddEntry(j_pt_lead, "Ch+Ne leading jets", "p");
  tjpt->AddEntry(j_pt_incl, "Ch+Ne inclusive jets","p");

  TCanvas * cjpt = new TCanvas("cjpt","cjpt",800,800); cjpt->cd(); cjpt->SetLogy();
  j_pt_incl->Draw(); j_pt_lead->Draw("same"); tjpt->Draw("same");

  cjpt->SaveAs((out_path + "jet_pt" + filetype).c_str());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~CONS PHI~~~~~~~~~~~~~~~~~~~~
  TH1 * c_phi_incl = c_PtEtaPhi_incl->Project3D("z");
  TH1 * c_phi_lead = c_PtEtaPhi_lead->Project3D("z");
  
  c_phi_incl->SetTitle(""); c_phi_lead->SetTitle("");
  
  c_phi_incl->SetMarkerStyle(20);
  c_phi_lead->SetMarkerStyle(21);
  c_phi_incl->SetMarkerColor(kBlue);
  c_phi_lead->SetMarkerColor(kRed);
  c_phi_incl->SetLineColor(kBlue);
  c_phi_lead->SetLineColor(kRed);

  c_phi_incl->GetXaxis()->SetTitle("#phi^{cons}");
  c_phi_incl->GetYaxis()->SetTitle("arb.");
  
  c_phi_incl->Scale(1/(double)c_phi_incl->Integral());
  c_phi_lead->Scale(1/(double)c_phi_lead->Integral());

  c_phi_incl->GetYaxis()->SetRangeUser(0,0.2);
  c_phi_lead->GetYaxis()->SetRangeUser(0,0.2);

  TLegend *tcphi = new TLegend(0.7,0.7,0.9,0.8); tcphi->SetBorderSize(0);
  tcphi->AddEntry(c_phi_lead, "Ch+Ne leading jet constituents", "p");
  tcphi->AddEntry(c_phi_incl, "Ch+Ne inclusive jet constituents","p");

  TCanvas * ccphi = new TCanvas("ccphi","ccphi",800,800); ccphi->cd();
  c_phi_incl->Draw(); c_phi_lead->Draw("same"); tcphi->Draw("same");

  ccphi->SaveAs((out_path + "cons_phi" + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~CONS ETA~~~~~~~~~~~~~~~~~~~~
  
  TH1 * c_eta_incl = c_PtEtaPhi_incl->Project3D("y");
  TH1 * c_eta_lead = c_PtEtaPhi_lead->Project3D("y");                                                                                                                                                    

  c_eta_incl->SetTitle(""); c_eta_lead->SetTitle("");

  c_eta_incl->SetMarkerStyle(20);
  c_eta_lead->SetMarkerStyle(21);
  c_eta_incl->SetMarkerColor(kBlue);
  c_eta_lead->SetMarkerColor(kRed);
  c_eta_incl->SetLineColor(kBlue);
  c_eta_lead->SetLineColor(kRed);

  c_eta_incl->GetXaxis()->SetTitle("#eta^{cons}");
  c_eta_incl->GetYaxis()->SetTitle("arb.");

  c_eta_incl->Scale(1/(double)c_eta_incl->Integral());
  c_eta_lead->Scale(1/(double)c_eta_lead->Integral());

  TLegend *tceta = new TLegend(0.7,0.7,0.9,0.8); tceta->SetBorderSize(0);
  tceta->AddEntry(c_eta_lead, "Ch+Ne leading jet constituents", "p");
  tceta->AddEntry(c_eta_incl, "Ch+Ne inclusive jet constituents","p");

  TCanvas * cceta = new TCanvas("cceta","cceta",800,800); cceta->cd();
  c_eta_incl->Draw(); c_eta_lead->Draw("same"); tceta->Draw("same");

  cceta->SaveAs((out_path + "cons_eta" + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~CONS PT~~~~~~~~~~~~~~~~~~~~                                                                                                                                                       

  TH1 * c_pt_incl = c_PtEtaPhi_incl->Project3D("x");
  TH1 * c_pt_lead = c_PtEtaPhi_lead->Project3D("x");

  c_pt_incl->SetTitle(""); c_pt_lead->SetTitle("");

  c_pt_incl->SetMarkerStyle(20);
  c_pt_lead->SetMarkerStyle(21);
  c_pt_incl->SetMarkerColor(kBlue);
  c_pt_lead->SetMarkerColor(kRed);
  c_pt_incl->SetLineColor(kBlue);
  c_pt_lead->SetLineColor(kRed);

  c_pt_incl->GetXaxis()->SetTitle("p_{T}^{cons} [GeV]");
  c_pt_incl->GetYaxis()->SetTitle("arb.");

  c_pt_incl->Scale(1/(double)c_pt_incl->Integral());
  c_pt_lead->Scale(1/(double)c_pt_lead->Integral());

  TLegend *tcpt = new TLegend(0.7,0.7,0.9,0.8); tcpt->SetBorderSize(0);
  tcpt->AddEntry(c_pt_lead, "Ch+Ne leading jet constituents", "p");
  tcpt->AddEntry(c_pt_incl, "Ch+Ne inclusive jet constituents","p");

  TCanvas * ccpt = new TCanvas("ccpt","ccpt",800,800); ccpt->cd(); ccpt->SetLogy();
  c_pt_incl->Draw(); c_pt_lead->Draw("same"); tcpt->Draw("same");

  ccpt->SaveAs((out_path + "cons_pt" + filetype).c_str());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  return;
}
