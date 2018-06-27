//I started this project by plotting pt, eta, phi & mass
//for inclusive and leading charged+neutral jets.
//This is a plotting macro for these quantities

#include <string>
#include <iostream>

using namespace std;

void initial_checks () {
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  string input = "~/jetmass/out/sim/full.root";
  string out_path = "~/jetmass/plots/sim/initial_checks/";
  string filetype = ".pdf";

  TFile *f = new TFile(input.c_str(), "READ");

  //~~~~~~~~~~~~~~~~~~MASS~~~~~~~~~~~~~~~~~~~~~~~~

  TH1D * p_m_incl = (TH1D*)f->Get("p_m_incl");
  TH1D * p_m_lead = (TH1D*)f->Get("p_m_lead");
  TH1D * g_m_incl = (TH1D*)f->Get("g_m_incl");
  TH1D * g_m_lead = (TH1D*)f->Get("g_m_lead");

  p_m_incl->SetMarkerStyle(20);
  p_m_lead->SetMarkerStyle(20);
  p_m_incl->SetMarkerColor(kBlue);
  p_m_lead->SetMarkerColor(kRed);
  p_m_incl->SetLineColor(kBlue);
  p_m_lead->SetLineColor(kRed);
  g_m_incl->SetMarkerStyle(21);
  g_m_lead->SetMarkerStyle(21);
  g_m_incl->SetMarkerColor(kBlue);
  g_m_lead->SetMarkerColor(kRed);
  g_m_incl->SetLineColor(kBlue);
  g_m_lead->SetLineColor(kRed);

  p_m_incl->GetXaxis()->SetTitle("M^{jet}");
  p_m_incl->GetYaxis()->SetTitle("prob.");
  g_m_incl->GetXaxis()->SetTitle("M^{jet}");
  g_m_incl->GetYaxis()->SetTitle("prob.");
   
  p_m_lead->Scale(1/(double)p_m_lead->Integral());
  p_m_incl->Scale(1/(double)p_m_incl->Integral());
  g_m_lead->Scale(1/(double)g_m_lead->Integral());
  g_m_incl->Scale(1/(double)g_m_incl->Integral());

  p_m_incl->GetXaxis()->SetRangeUser(0,11);
  p_m_lead->GetXaxis()->SetRangeUser(0,11);
  g_m_incl->GetXaxis()->SetRangeUser(0,11);
  g_m_lead->GetXaxis()->SetRangeUser(0,11);

  TLegend *tmass = new TLegend(0.7,0.7,0.9,0.8); tmass->SetBorderSize(0);
  tmass->AddEntry((TObject*)0, "Ch+Ne jets", "");
  tmass->AddEntry(p_m_incl, "PYTHIA6 inclusive", "p");
  tmass->AddEntry(p_m_lead, "PYTHIA6 leading", "p");
  tmass->AddEntry(g_m_incl, "PYTHIA6+GEANT inclusive","p");
  tmass->AddEntry(g_m_lead, "PYTHIA6+GEANT leading","p");
  tmass->Draw();

  TCanvas * cmass = new TCanvas("cmass","cmass",800,800); cmass->cd(); cmass->SetLogy();
  p_m_incl->Draw(); p_m_lead->Draw("same"); g_m_incl->Draw("same"); g_m_lead->Draw("same"); tmass->Draw("same");
  
  cmass->SaveAs((out_path + "mass" + filetype).c_str());
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~3D HISTS~~~~~~~~~~~~~~~~~~~
  
  TH3D * p_j_PtEtaPhi_incl = (TH3D*)f->Get("p_PtEtaPhi_incl");
  TH3D * p_j_PtEtaPhi_lead = (TH3D*)f->Get("p_PtEtaPhi_lead");
  TH3D * p_c_PtEtaPhi_incl = (TH3D*)f->Get("p_cons_PtEtaPhi_incl");
  TH3D * p_c_PtEtaPhi_lead = (TH3D*)f->Get("p_cons_PtEtaPhi_lead");
  TH3D * g_j_PtEtaPhi_incl = (TH3D*)f->Get("g_PtEtaPhi_incl");
  TH3D * g_j_PtEtaPhi_lead = (TH3D*)f->Get("g_PtEtaPhi_lead");
  TH3D * g_c_PtEtaPhi_incl = (TH3D*)f->Get("g_cons_PtEtaPhi_incl");
  TH3D * g_c_PtEtaPhi_lead = (TH3D*)f->Get("g_cons_PtEtaPhi_lead");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   //marker, axis, scale - incl blue, lead red 

  //~~~~~~~~~~~~~~~~~~~JET PHI~~~~~~~~~~~~~~~~~~~~
  TH1 * p_j_phi_incl = p_j_PtEtaPhi_incl->Project3D("z");
  TH1 * p_j_phi_lead = p_j_PtEtaPhi_lead->Project3D("z");
  TH1 * g_j_phi_incl = g_j_PtEtaPhi_incl->Project3D("z");
  TH1 * g_j_phi_lead = g_j_PtEtaPhi_lead->Project3D("z");
  
  p_j_phi_incl->SetTitle(""); p_j_phi_lead->SetTitle("");
  g_j_phi_incl->SetTitle(""); g_j_phi_lead->SetTitle("");

  p_j_phi_incl->SetMarkerStyle(20);
  p_j_phi_lead->SetMarkerStyle(20);
  p_j_phi_incl->SetMarkerColor(kBlue);
  p_j_phi_lead->SetMarkerColor(kRed);
  p_j_phi_incl->SetLineColor(kBlue);
  p_j_phi_lead->SetLineColor(kRed);
  g_j_phi_incl->SetMarkerStyle(21);
  g_j_phi_lead->SetMarkerStyle(21);
  g_j_phi_incl->SetMarkerColor(kBlue);
  g_j_phi_lead->SetMarkerColor(kRed);
  g_j_phi_incl->SetLineColor(kBlue);
  g_j_phi_lead->SetLineColor(kRed);

  p_j_phi_incl->GetXaxis()->SetTitle("#phi^{jet}");
  p_j_phi_incl->GetYaxis()->SetTitle("arb.");
  g_j_phi_incl->GetXaxis()->SetTitle("#phi^{jet}");
  g_j_phi_incl->GetYaxis()->SetTitle("arb.");
   
  p_j_phi_incl->Scale(1/(double)p_j_phi_incl->Integral());
  p_j_phi_lead->Scale(1/(double)p_j_phi_lead->Integral());
  g_j_phi_incl->Scale(1/(double)g_j_phi_incl->Integral());
  g_j_phi_lead->Scale(1/(double)g_j_phi_lead->Integral());

  p_j_phi_incl->GetYaxis()->SetRangeUser(0,0.2);
  p_j_phi_lead->GetYaxis()->SetRangeUser(0,0.2);
  g_j_phi_incl->GetYaxis()->SetRangeUser(0,0.2);
  g_j_phi_lead->GetYaxis()->SetRangeUser(0,0.2);


  TLegend *tjphi = new TLegend(0.7,0.7,0.9,0.8); tjphi->SetBorderSize(0);
  tjphi->AddEntry((TObject*)0, "Ch+Ne jets", "");
  tjphi->AddEntry(p_j_phi_incl, "PYTHIA6 inclusive","p");
  tjphi->AddEntry(p_j_phi_lead, "PYTHIA6 leading","p");
  tjphi->AddEntry(g_j_phi_incl, "PYTHIA6+GEANT inclusive","p");
  tjphi->AddEntry(g_j_phi_lead, "PYTHIA6+GEANT leading","p");

  TCanvas * cjphi = new TCanvas("cjphi","cjphi",800,800); cjphi->cd();
  p_j_phi_incl->Draw(); p_j_phi_lead->Draw("same"); g_j_phi_incl->Draw("same"); g_j_phi_lead->Draw("same"); tjphi->Draw("same");

  cjphi->SaveAs((out_path + "jet_phi" + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~JET ETA~~~~~~~~~~~~~~~~~~~~
  
  TH1 * p_j_eta_incl = p_j_PtEtaPhi_incl->Project3D("y");
  TH1 * p_j_eta_lead = p_j_PtEtaPhi_lead->Project3D("y");
  TH1 * g_j_eta_incl = g_j_PtEtaPhi_incl->Project3D("y");
  TH1 * g_j_eta_lead = g_j_PtEtaPhi_lead->Project3D("y");                            

  p_j_eta_incl->SetTitle(""); p_j_eta_lead->SetTitle("");
  g_j_eta_incl->SetTitle(""); g_j_eta_lead->SetTitle("");
   
  p_j_eta_incl->SetMarkerStyle(20);
  p_j_eta_lead->SetMarkerStyle(20);
  p_j_eta_incl->SetMarkerColor(kBlue);
  p_j_eta_lead->SetMarkerColor(kRed);
  p_j_eta_incl->SetLineColor(kBlue);
  p_j_eta_lead->SetLineColor(kRed);
  g_j_eta_incl->SetMarkerStyle(21);
  g_j_eta_lead->SetMarkerStyle(21);
  g_j_eta_incl->SetMarkerColor(kBlue);
  g_j_eta_lead->SetMarkerColor(kRed);
  g_j_eta_incl->SetLineColor(kBlue);
  g_j_eta_lead->SetLineColor(kRed);

  p_j_eta_incl->GetXaxis()->SetTitle("#eta^{jet}");
  p_j_eta_incl->GetYaxis()->SetTitle("arb.");
  g_j_eta_incl->GetXaxis()->SetTitle("#eta^{jet}");
  g_j_eta_incl->GetYaxis()->SetTitle("arb.");

  p_j_eta_incl->Scale(1/(double)p_j_eta_incl->Integral());
  p_j_eta_lead->Scale(1/(double)p_j_eta_lead->Integral());
  g_j_eta_incl->Scale(1/(double)g_j_eta_incl->Integral());
  g_j_eta_lead->Scale(1/(double)g_j_eta_lead->Integral());

  TLegend *tjeta = new TLegend(0.7,0.7,0.9,0.8); tjeta->SetBorderSize(0);
  tjeta->AddEntry((TObject*)0, "Ch+Ne jets", "");
  tjeta->AddEntry(p_j_eta_incl, "PYTHIA6 inclusive","p");
  tjeta->AddEntry(p_j_eta_lead, "PYTHIA6 leading","p");
  tjeta->AddEntry(g_j_eta_incl, "PYTHIA6+GEANT inclusive","p");
  tjeta->AddEntry(g_j_eta_lead, "PYTHIA6+GEANT leading","p");

  TCanvas * cjeta = new TCanvas("cjeta","cjeta",800,800); cjeta->cd();
  p_j_eta_incl->Draw(); p_j_eta_lead->Draw("same"); g_j_eta_incl->Draw("same"); g_j_eta_lead->Draw("same"); tjeta->Draw("same");

  cjeta->SaveAs((out_path + "jet_eta" + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~JET PT~~~~~~~~~~~~~~~~~~~~                                                                                                                                                         

  TH1 * p_j_pt_incl = p_j_PtEtaPhi_incl->Project3D("x");
  TH1 * p_j_pt_lead = p_j_PtEtaPhi_lead->Project3D("x");
  TH1 * g_j_pt_incl = g_j_PtEtaPhi_incl->Project3D("x");
  TH1 * g_j_pt_lead = g_j_PtEtaPhi_lead->Project3D("x");

  p_j_pt_incl->SetTitle(""); p_j_pt_lead->SetTitle("");
  g_j_pt_incl->SetTitle(""); g_j_pt_lead->SetTitle("");

  p_j_pt_incl->SetMarkerStyle(20);
  p_j_pt_lead->SetMarkerStyle(20);
  p_j_pt_incl->SetMarkerColor(kBlue);
  p_j_pt_lead->SetMarkerColor(kRed);
  p_j_pt_incl->SetLineColor(kBlue);
  p_j_pt_lead->SetLineColor(kRed);
  g_j_pt_incl->SetMarkerStyle(21);
  g_j_pt_lead->SetMarkerStyle(21);
  g_j_pt_incl->SetMarkerColor(kBlue);
  g_j_pt_lead->SetMarkerColor(kRed);
  g_j_pt_incl->SetLineColor(kBlue);
  g_j_pt_lead->SetLineColor(kRed);

  p_j_pt_incl->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  p_j_pt_incl->GetYaxis()->SetTitle("arb.");
  g_j_pt_incl->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  g_j_pt_incl->GetYaxis()->SetTitle("arb.");
  
  p_j_pt_incl->Scale(1/(double)p_j_pt_incl->Integral());
  p_j_pt_lead->Scale(1/(double)p_j_pt_lead->Integral());
  g_j_pt_incl->Scale(1/(double)g_j_pt_incl->Integral());
  g_j_pt_lead->Scale(1/(double)g_j_pt_lead->Integral());

  TLegend *tjpt = new TLegend(0.7,0.7,0.9,0.8); tjpt->SetBorderSize(0);
  tjpt->AddEntry((TObject*)0, "Ch+Ne jets", "");
  tjpt->AddEntry(p_j_pt_incl, "PYTHIA6 inclusive","p");
  tjpt->AddEntry(p_j_pt_lead, "PYTHIA6 leading","p");
  tjpt->AddEntry(g_j_pt_incl, "PYTHIA6+GEANT inclusive","p");
  tjpt->AddEntry(g_j_pt_lead, "PYTHIA6+GEANT leading","p");

  TCanvas * cjpt = new TCanvas("cjpt","cjpt",800,800); cjpt->cd(); cjpt->SetLogy();
  p_j_pt_incl->Draw(); p_j_pt_lead->Draw("same"); g_j_pt_incl->Draw("same"); g_j_pt_lead->Draw("same"); tjpt->Draw("same");

  cjpt->SaveAs((out_path + "jet_pt" + filetype).c_str());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  //~~~~~~~~~~~~~~~~~~~CONS PHI~~~~~~~~~~~~~~~~~~~~
  TH1 * p_c_phi_incl = p_c_PtEtaPhi_incl->Project3D("z");
  TH1 * p_c_phi_lead = p_c_PtEtaPhi_lead->Project3D("z");
  TH1 * g_c_phi_incl = g_c_PtEtaPhi_incl->Project3D("z");
  TH1 * g_c_phi_lead = g_c_PtEtaPhi_lead->Project3D("z");
  
  p_c_phi_incl->SetTitle(""); p_c_phi_lead->SetTitle("");
  g_c_phi_incl->SetTitle(""); g_c_phi_lead->SetTitle("");

  p_c_phi_incl->SetMarkerStyle(20);
  p_c_phi_lead->SetMarkerStyle(20);
  p_c_phi_incl->SetMarkerColor(kBlue);
  p_c_phi_lead->SetMarkerColor(kRed);
  p_c_phi_incl->SetLineColor(kBlue);
  p_c_phi_lead->SetLineColor(kRed);
  g_c_phi_incl->SetMarkerStyle(21);
  g_c_phi_lead->SetMarkerStyle(21);
  g_c_phi_incl->SetMarkerColor(kBlue);
  g_c_phi_lead->SetMarkerColor(kRed);
  g_c_phi_incl->SetLineColor(kBlue);
  g_c_phi_lead->SetLineColor(kRed);

  p_c_phi_incl->GetXaxis()->SetTitle("#phi^{cons}");
  p_c_phi_incl->GetYaxis()->SetTitle("arb.");
  g_c_phi_incl->GetXaxis()->SetTitle("#phi^{cons}");
  g_c_phi_incl->GetYaxis()->SetTitle("arb.");
   
  p_c_phi_incl->Scale(1/(double)p_c_phi_incl->Integral());
  p_c_phi_lead->Scale(1/(double)p_c_phi_lead->Integral());
  g_c_phi_incl->Scale(1/(double)g_c_phi_incl->Integral());
  g_c_phi_lead->Scale(1/(double)g_c_phi_lead->Integral());

  p_c_phi_incl->GetYaxis()->SetRangeUser(0,0.2);
  p_c_phi_lead->GetYaxis()->SetRangeUser(0,0.2);
  g_c_phi_incl->GetYaxis()->SetRangeUser(0,0.2);
  g_c_phi_lead->GetYaxis()->SetRangeUser(0,0.2);

  TLegend *tcphi = new TLegend(0.7,0.7,0.9,0.8); tcphi->SetBorderSize(0);
  tcphi->AddEntry((TObject*)0, "Ch+Ne jet constituents", "");
  tcphi->AddEntry(p_c_phi_incl, "PYTHIA6 inclusive","p");
  tcphi->AddEntry(p_c_phi_lead, "PYTHIA6 leading","p");
  tcphi->AddEntry(g_c_phi_incl, "PYTHIA6+GEANT inclusive","p");
  tcphi->AddEntry(g_c_phi_lead, "PYTHIA6+GEANT leading","p");

  TCanvas * ccphi = new TCanvas("ccphi","ccphi",800,800); ccphi->cd();
  p_c_phi_incl->Draw(); p_c_phi_lead->Draw("same"); g_c_phi_incl->Draw("same"); g_c_phi_lead->Draw("same"); tcphi->Draw("same");

  ccphi->SaveAs((out_path + "cons_phi" + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~CONS ETA~~~~~~~~~~~~~~~~~~~~
  
  TH1 * p_c_eta_incl = p_c_PtEtaPhi_incl->Project3D("y");
  TH1 * p_c_eta_lead = p_c_PtEtaPhi_lead->Project3D("y");
  TH1 * g_c_eta_incl = g_c_PtEtaPhi_incl->Project3D("y");
  TH1 * g_c_eta_lead = g_c_PtEtaPhi_lead->Project3D("y");                            

  p_c_eta_incl->SetTitle(""); p_c_eta_lead->SetTitle("");
  g_c_eta_incl->SetTitle(""); g_c_eta_lead->SetTitle("");
   
  p_c_eta_incl->SetMarkerStyle(20);
  p_c_eta_lead->SetMarkerStyle(20);
  p_c_eta_incl->SetMarkerColor(kBlue);
  p_c_eta_lead->SetMarkerColor(kRed);
  p_c_eta_incl->SetLineColor(kBlue);
  p_c_eta_lead->SetLineColor(kRed);
  g_c_eta_incl->SetMarkerStyle(21);
  g_c_eta_lead->SetMarkerStyle(21);
  g_c_eta_incl->SetMarkerColor(kBlue);
  g_c_eta_lead->SetMarkerColor(kRed);
  g_c_eta_incl->SetLineColor(kBlue);
  g_c_eta_lead->SetLineColor(kRed);

  p_c_eta_incl->GetXaxis()->SetTitle("#eta^{cons}");
  p_c_eta_incl->GetYaxis()->SetTitle("arb.");
  g_c_eta_incl->GetXaxis()->SetTitle("#eta^{cons}");
  g_c_eta_incl->GetYaxis()->SetTitle("arb.");

  p_c_eta_incl->Scale(1/(double)p_c_eta_incl->Integral());
  p_c_eta_lead->Scale(1/(double)p_c_eta_lead->Integral());
  g_c_eta_incl->Scale(1/(double)g_c_eta_incl->Integral());
  g_c_eta_lead->Scale(1/(double)g_c_eta_lead->Integral());

  TLegend *tceta = new TLegend(0.63,0.59,0.83,0.69); tceta->SetBorderSize(0);
  tceta->AddEntry((TObject*)0, "Ch+Ne jet constituents","");
  tceta->AddEntry(p_c_eta_incl, "PYTHIA6 inclusive","p");
  tceta->AddEntry(p_c_eta_lead, "PYTHIA6 leading","p");
  tceta->AddEntry(g_c_eta_incl, "PYTHIA6+GEANT inclusive","p");
  tceta->AddEntry(g_c_eta_lead, "PYTHIA6+GEANT leading","p");

  TCanvas * cceta = new TCanvas("cceta","cceta",800,800); cceta->cd();
  p_c_eta_incl->Draw(); p_c_eta_lead->Draw("same"); g_c_eta_incl->Draw("same"); g_c_eta_lead->Draw("same"); tceta->Draw("same");

  cceta->SaveAs((out_path + "cons_eta" + filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~CONS PT~~~~~~~~~~~~~~~~~~~~                                                                                                             

  TH1 * p_c_pt_incl = p_c_PtEtaPhi_incl->Project3D("x");
  TH1 * p_c_pt_lead = p_c_PtEtaPhi_lead->Project3D("x");
  TH1 * g_c_pt_incl = g_c_PtEtaPhi_incl->Project3D("x");
  TH1 * g_c_pt_lead = g_c_PtEtaPhi_lead->Project3D("x");

  p_c_pt_incl->SetTitle(""); p_c_pt_lead->SetTitle("");
  g_c_pt_incl->SetTitle(""); g_c_pt_lead->SetTitle("");

  p_c_pt_incl->SetMarkerStyle(20);
  p_c_pt_lead->SetMarkerStyle(20);
  p_c_pt_incl->SetMarkerColor(kBlue);
  p_c_pt_lead->SetMarkerColor(kRed);
  p_c_pt_incl->SetLineColor(kBlue);
  p_c_pt_lead->SetLineColor(kRed);
  g_c_pt_incl->SetMarkerStyle(21);
  g_c_pt_lead->SetMarkerStyle(21);
  g_c_pt_incl->SetMarkerColor(kBlue);
  g_c_pt_lead->SetMarkerColor(kRed);
  g_c_pt_incl->SetLineColor(kBlue);
  g_c_pt_lead->SetLineColor(kRed);

  p_c_pt_incl->GetXaxis()->SetTitle("p_{T}^{cons} [GeV]");
  p_c_pt_incl->GetYaxis()->SetTitle("arb.");
  g_c_pt_incl->GetXaxis()->SetTitle("p_{T}^{cons} [GeV]");
  g_c_pt_incl->GetYaxis()->SetTitle("arb.");
  
  p_c_pt_incl->Scale(1/(double)p_c_pt_incl->Integral());
  p_c_pt_lead->Scale(1/(double)p_c_pt_lead->Integral());
  g_c_pt_incl->Scale(1/(double)g_c_pt_incl->Integral());
  g_c_pt_lead->Scale(1/(double)g_c_pt_lead->Integral());

  TLegend *tcpt = new TLegend(0.7,0.7,0.9,0.8); tcpt->SetBorderSize(0);
  tcpt->AddEntry((TObject*)0, "Ch+Ne jet constituents", "");
  tcpt->AddEntry(p_c_pt_incl, "PYTHIA6 inclusive","p");
  tcpt->AddEntry(p_c_pt_lead, "PYTHIA6 leading","p");
  tcpt->AddEntry(g_c_pt_incl, "PYTHIA6+GEANT inclusive","p");
  tcpt->AddEntry(g_c_pt_lead, "PYTHIA6+GEANT leading","p");

  TCanvas * ccpt = new TCanvas("ccpt","ccpt",800,800); ccpt->cd(); ccpt->SetLogy();
  p_c_pt_incl->Draw(); p_c_pt_lead->Draw("same"); g_c_pt_incl->Draw("same"); g_c_pt_lead->Draw("same"); tcpt->Draw("same");

  ccpt->SaveAs((out_path + "cons_pt" + filetype).c_str());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /*
  //REMOVE SOON!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
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
  */

  return;
}
