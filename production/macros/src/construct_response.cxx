#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "math.h"

#include "funcs.hh"

using namespace std; using namespace Analysis;


int main() {
  string dir = "~/jetmass/production/";
  string in = "";
  string file_on = "decays_on.root";
  string file_off = "decays_off.root";
  string out = "~/jetmass/production/macros/hists/";
  string filetype = ".pdf";

  TFile* onFile = new TFile( (dir + in + file_on).c_str(), "READ");
  TFile* offFile = new TFile( (dir + in + file_off).c_str(), "READ");
  
  //  RooUnfoldResponse *m_decay_res = new RooUnfoldResponse(20,0,10,20,0,10,"m_decay_res","");  
  //std::vector<RooUnfoldResponse*> res_vec;
  //res_vec.push_back(m_decay_res);
  
  TH1D* hon1520 = new TH1D("hon1520","",10,0,10);
  TH1D* hon2025 = new TH1D("hon2025","",10,0,10);
  TH1D* hon2530 = new TH1D("hon2530","",10,0,10);
  TH1D* hon3040 = new TH1D("hon3040","",10,0,10);
  TH1D* hon4060 = new TH1D("hon4060","",10,0,10);

  TH1D* pton = new TH1D("pton","",15,5,80);
  
  //  TH1D* hon_clone = (TH1D*) hon->Clone();
  TH1D* hoff1520 = new TH1D("hoff1520","",10,0,10);
  TH1D* hoff2025 = new TH1D("hoff2025","",10,0,10);
  TH1D* hoff2530 = new TH1D("hoff2530","",10,0,10);
  TH1D* hoff3040 = new TH1D("hoff3040","",10,0,10);
  TH1D* hoff4060 = new TH1D("hoff4060","",10,0,10);

  TH1D* ptoff = new TH1D("ptoff","",15,5,80);
  
  HistFromTree(onFile, pton, hon1520,hon2025,hon2530,hon3040,hon4060/*offFile, res_vec*/);  
  HistFromTree(offFile, ptoff, hoff1520,hoff2025,hoff2530,hoff3040,hoff4060);
  
  hon1520->Scale(1/(double)hon1520->Integral()); hoff1520->Scale(1/(double)hoff1520->Integral());
  hon2025->Scale(1/(double)hon2025->Integral()); hoff2025->Scale(1/(double)hoff2025->Integral());
  hon2530->Scale(1/(double)hon2530->Integral()); hoff2530->Scale(1/(double)hoff2530->Integral());
  hon3040->Scale(1/(double)hon3040->Integral()); hoff3040->Scale(1/(double)hoff3040->Integral());
  hon4060->Scale(1/(double)hon4060->Integral()); hoff4060->Scale(1/(double)hoff4060->Integral());

  pton->Scale(1/(double)pton->Integral()); ptoff->Scale(1/(double)ptoff->Integral());

  double binwidth_pt = (pton->GetXaxis()->GetXmax() - pton->GetXaxis()->GetXmin()) / (double) pton->GetXaxis()->GetNbins();
  double binwidth_mass = (hon1520->GetXaxis()->GetXmax() - hon1520->GetXaxis()->GetXmin()) / (double) hon1520->GetXaxis()->GetNbins();
  hon1520->Scale(1/(double)binwidth_mass);
  hon2025->Scale(1/(double)binwidth_mass);
  hon2530->Scale(1/(double)binwidth_mass);
  hon3040->Scale(1/(double)binwidth_mass);
  hon4060->Scale(1/(double)binwidth_mass);
  hoff1520->Scale(1/(double)binwidth_mass);
  hoff2025->Scale(1/(double)binwidth_mass);
  hoff2530->Scale(1/(double)binwidth_mass);
  hoff3040->Scale(1/(double)binwidth_mass);
  hoff4060->Scale(1/(double)binwidth_mass);

  pton->Scale(1/(double)binwidth_pt); ptoff->Scale(1/(double)binwidth_pt);
  
  TH1D* hratio1520 = (TH1D*) hoff1520->Clone("hratio1520");
  TH1D* hratio2025 = (TH1D*) hoff2025->Clone("hratio2025");
  TH1D* hratio2530 = (TH1D*) hoff2530->Clone("hratio2530");
  TH1D* hratio3040 = (TH1D*) hoff3040->Clone("hratio3040");
  TH1D* hratio4060 = (TH1D*) hoff4060->Clone("hratio4060");
  
  TH1D* ptratio = (TH1D*) ptoff->Clone("ptratio");
  
  hratio1520->Divide(hon1520);
  hratio2025->Divide(hon2025);
  hratio2530->Divide(hon2530);
  hratio3040->Divide(hon3040);
  hratio4060->Divide(hon4060);
  
  ptratio->Divide(pton);
  
  TFile *fout = new TFile( ( out + "ratio.root" ).c_str() ,"RECREATE");
  
  //m_decay_res->Write();
  hratio1520->Write();
  hratio2025->Write();
  hratio2530->Write();
  hratio3040->Write();
  hratio4060->Write();

  hoff1520->Write(); hon1520->Write();
  hoff2025->Write(); hon2025->Write();
  hoff2530->Write(); hon2530->Write();
  hoff3040->Write(); hon3040->Write();
  hoff4060->Write(); hon4060->Write();
  
  pton->Write(); ptoff->Write(); ptratio->Write();
  
  fout->Close();
  
  return 0;
}
