#include <string>
#include <iostream>

//CHECKING THAT TRIGGER AND RECOIL JET SELECTION WORKS IN BOTH DATA AND SIMULATION                                             

using namespace std;

void trig_recoil_initial_checks () {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HISTOGRAMS AND CANVASES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  string dir = "~/jetmass/";
  string simin = "out/sim/";
  string datain = "out/data/";
  string file = "full.root";
  string out = "~/jetmass/plots/trig_rec_initial_checks/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "rec";
  string flag3 = "cons";
    
  TCanvas *cpt = new TCanvas("cpt","cpt",800,800); cpt->SetLogy();
  TCanvas *ceta = new TCanvas("ceta","ceta",800,800); ceta->SetLogy();
  TCanvas *cphi = new TCanvas("cphi","cphi",800,800);
  TCanvas *cm = new TCanvas("cm","cm",800,800); cm->SetLogy();
  
  TFile* simFile = new TFile( (dir + simin + file).c_str(), "READ");
  TH3D *ptetaphi_py = (TH3D*) simFile->Get(("PtEtaPhi_" + flag1 + "_" + flag2 + "_" + flag3 + "_" + "py").c_str());
  TH3D *ptetaphi_ge = (TH3D*) simFile->Get(("PtEtaPhi_" + flag1 + "_" + flag2 + "_" + flag3 + "_" + "ge").c_str());
  TH1D *m_py = (TH1D*) simFile->Get(("m_" + flag1 + "_" + flag2 + "_" + flag3 + "_" + "py").c_str());
  TH1D *m_ge = (TH1D*) simFile->Get(("m_" + flag1 + "_" + flag2 + "_" + flag3 + "_" + "ge").c_str());
   
    
  TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  TH3D *ptetaphi_dat= (TH3D*) dataFile->Get(("PtEtaPhi_" + flag1 + "_" + flag2 + "_" + flag3).c_str());
   TH1D *m_dat = (TH1D*) dataFile->Get(("m_" + flag1 + "_" + flag2 + "_" + flag3).c_str());
  

  //string pyfull = py_m_pt->GetName();
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PROJECTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    TH1 *pt_py = ptetaphi_py->Project3D("x"); TH1 *pt_ge = ptetaphi_ge->Project3D("x"); TH1 *pt_dat = ptetaphi_dat->Project3D("x");
    pt_py->SetTitle(""); pt_ge->SetTitle(""); pt_dat->SetTitle("");
    TH1 *eta_py = ptetaphi_py->Project3D("y"); TH1 *eta_ge = ptetaphi_ge->Project3D("y"); TH1 *eta_dat = ptetaphi_dat->Project3D("y");
    eta_py->SetTitle(""); eta_ge->SetTitle(""); eta_dat->SetTitle("");
     TH1 *phi_py = ptetaphi_py->Project3D("z"); TH1 *phi_ge = ptetaphi_ge->Project3D("z"); TH1 *phi_dat = ptetaphi_dat->Project3D("z");
    phi_py->SetTitle(""); phi_ge->SetTitle(""); phi_dat->SetTitle("");
     
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TLEGENDS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    TLegend *tpt = new TLegend(0.244,0.15,0.57,0.25); tpt->SetBorderSize(0);
    tpt->AddEntry((TObject*)0,(flag1 + " " + flag2 + " " + flag3).c_str(), "");
    tpt->AddEntry(pt_py, "PYTHIA6", "p");
    tpt->AddEntry(pt_ge, "PYTHIA6+GEANT", "p");
    tpt->AddEntry(pt_dat, "STAR pp run6 - 200 GeV", "p");

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRETTIFYING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //axes, markers, ranges, scaling
    pt_py->Scale(1/(double)pt_py->Integral()); pt_ge->Scale(1/(double)pt_ge->Integral()); pt_dat->Scale(1/(double)pt_dat->Integral());
    pt_py->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]"); pt_ge->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]"); pt_dat->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
    pt_py->GetYaxis()->SetTitle("prob."); pt_ge->GetYaxis()->SetTitle("prob."); pt_dat->GetYaxis()->SetTitle("prob.");
    pt_py->GetYaxis()->SetTitleOffset(1.3); pt_ge->GetYaxis()->SetTitleOffset(1.3); pt_dat->GetYaxis()->SetTitleOffset(1.3);
    pt_py->SetMarkerColor(kBlue); pt_ge->SetMarkerColor(kRed); pt_dat->SetMarkerColor(kViolet);
    pt_py->SetLineColor(kBlue); pt_ge->SetLineColor(kRed); pt_dat->SetLineColor(kViolet);
    pt_py->SetMarkerStyle(20); pt_ge->SetMarkerStyle(21); pt_dat->SetMarkerStyle(24);
    //    pt_py->GetXaxis()->SetRangeUser(0,40); pt_ge->GetXaxis()->SetRangeUser(0,40); pt_dat->GetXaxis()->SetRangeUser(0,40);
    
    eta_py->Scale(1/(double)eta_py->Integral()); eta_ge->Scale(1/(double)eta_ge->Integral()); eta_dat->Scale(1/(double)eta_dat->Integral());
    eta_py->GetXaxis()->SetTitle("#eta^{jet}"); eta_ge->GetXaxis()->SetTitle("#eta^{jet}"); eta_dat->GetXaxis()->SetTitle("#eta^{jet}");
    eta_py->GetYaxis()->SetTitle("prob."); eta_ge->GetYaxis()->SetTitle("prob."); eta_dat->GetYaxis()->SetTitle("prob.");
    eta_py->GetYaxis()->SetTitleOffset(1.3); eta_ge->GetYaxis()->SetTitleOffset(1.3); eta_dat->GetYaxis()->SetTitleOffset(1.3);
    eta_py->SetMarkerColor(kBlue); eta_ge->SetMarkerColor(kRed); eta_dat->SetMarkerColor(kViolet);
    eta_py->SetLineColor(kBlue); eta_ge->SetLineColor(kRed); eta_dat->SetLineColor(kViolet);
    eta_py->SetMarkerStyle(20); eta_ge->SetMarkerStyle(21); eta_dat->SetMarkerStyle(24);
    //eta_py->GetXaxis()->SetRangeUser(0,40); eta_ge->GetXaxis()->SetRangeUser(0,40); eta_dat->GetXaxis()->SetRangeUser(0,40);

    phi_py->Scale(1/(double)phi_py->Integral()); phi_ge->Scale(1/(double)phi_ge->Integral()); phi_dat->Scale(1/(double)phi_dat->Integral());
    phi_py->GetXaxis()->SetTitle("#phi^{jet}"); phi_ge->GetXaxis()->SetTitle("#phi^{jet}"); phi_dat->GetXaxis()->SetTitle("#phi^{jet}");
    phi_py->GetYaxis()->SetTitle("prob."); phi_ge->GetYaxis()->SetTitle("prob."); phi_dat->GetYaxis()->SetTitle("prob.");
    phi_py->GetYaxis()->SetTitleOffset(1.3); phi_ge->GetYaxis()->SetTitleOffset(1.3); phi_dat->GetYaxis()->SetTitleOffset(1.3);
    phi_py->SetMarkerColor(kBlue); phi_ge->SetMarkerColor(kRed); phi_dat->SetMarkerColor(kViolet);
    phi_py->SetLineColor(kBlue); phi_ge->SetLineColor(kRed); phi_dat->SetLineColor(kViolet);
    phi_py->SetMarkerStyle(20); phi_ge->SetMarkerStyle(21); phi_dat->SetMarkerStyle(24);
    //phi_py->GetXaxis()->SetRangeUser(0,40); phi_ge->GetXaxis()->SetRangeUser(0,40); phi_dat->GetXaxis()->SetRangeUser(0,40);
    phi_py->GetYaxis()->SetRangeUser(0,0.2); phi_ge->GetYaxis()->SetRangeUser(0,0.2); phi_dat->GetYaxis()->SetRangeUser(0,0.2);


    
    m_py->Scale(1/(double)m_py->Integral()); m_ge->Scale(1/(double)m_ge->Integral()); m_dat->Scale(1/(double)m_dat->Integral());
    m_py->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); m_ge->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); m_dat->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    m_py->GetYaxis()->SetTitle("prob."); m_ge->GetYaxis()->SetTitle("prob."); m_dat->GetYaxis()->SetTitle("prob.");
    m_py->GetYaxis()->SetTitleOffset(1.3); m_ge->GetYaxis()->SetTitleOffset(1.3); m_dat->GetYaxis()->SetTitleOffset(1.3);
    m_py->SetMarkerColor(kBlue); m_ge->SetMarkerColor(kRed); m_dat->SetMarkerColor(kViolet);
    m_py->SetLineColor(kBlue); m_ge->SetLineColor(kRed); m_dat->SetLineColor(kViolet);
    m_py->SetMarkerStyle(20); m_ge->SetMarkerStyle(21); m_dat->SetMarkerStyle(24);
    m_py->GetXaxis()->SetRangeUser(0,16); m_ge->GetXaxis()->SetRangeUser(0,16); m_dat->GetXaxis()->SetRangeUser(0,16);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DRAWING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    cpt->cd(); pt_py->Draw("same"); pt_ge->Draw("same"); pt_dat->Draw("same"); tpt->Draw("same");
    ceta->cd(); eta_py->Draw("same"); eta_ge->Draw("same"); eta_dat->Draw("same"); tpt->Draw("same");
    cphi->cd(); phi_py->Draw("same"); phi_ge->Draw("same"); phi_dat->Draw("same"); tpt->Draw("same");
    cm->cd(); m_py->Draw("same"); m_ge->Draw("same"); m_dat->Draw("same"); tpt->Draw("same");
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SAVING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
     cpt->SaveAs((out + "pt_" + flag1 + "_" + flag2 + "_" + flag3 + filetype).c_str());
     ceta->SaveAs((out + "eta_" + flag1 + "_" + flag2 + "_" + flag3 + filetype).c_str());
     cphi->SaveAs((out + "phi_" + flag1 + "_" + flag2 + "_" + flag3 + filetype).c_str());
     cm->SaveAs((out + "m_" + flag1 + "_" + flag2 + "_" + flag3 + filetype).c_str());
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 
    
  return;
}
