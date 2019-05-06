#include "Plots.h"
#include "RooUnfold.h"
#include <string>
#include <iostream>

using namespace std;

void quicktest () {
  TFile *f = new TFile("~/jetmass/out/matching/test2.root","READ");
  
  RooUnfoldResponse* res = (RooUnfoldResponse*) f->Get("pt_res_coarse");
  
  TH2D* hres = (TH2D*) res->Hresponse();
  TH1D* gen_matched = (TH1D*) hres->ProjectionY("gen_matched");
  TH1D* truth = (TH1D*) res->Htruth();
  
  TH1D* hratio = (TH1D*) gen_matched->Clone("hratio");
  hratio->Divide(truth);
  
  Prettify1D(hratio,kBlue,kOpenCircle,3,kBlue,"Gen jet p_{T} [GeV/c]","Gen matched / Gen truth",5,80,0,1.2);
  
  hratio->Draw();
  
  return;
}
