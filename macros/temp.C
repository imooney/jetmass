#include "Plots.h"

using namespace std;

void temp () {
  TFile *fmat = new TFile ("~/jetmass/out/matching/temp_jetptmin10_JP2_evtEtCut_on.root","READ");
  TFile *finc = new TFile ("~/jetmass/out/sim/py/temp_0_on.root","READ");
  
  RooUnfoldResponse* res = (RooUnfoldResponse*) fmat->Get("pt_response");//"pt_res_coarse");
  TH2D* res2D = (TH2D*) res->Hresponse();
  TH1D* genmat = (TH1D*) res2D->ProjectionY("genmat");
  
  TTree* tinc = (TTree*) finc->Get("event"); 
  TH1D* hinc = new TH1D("hinc","",15,5,80);
  
  tinc->Draw("Pt>>hinc","weight");
  genmat->Draw("same");
  
  TH1D* hratio = (TH1D*) genmat->Clone("hratio");
  
  hratio->Divide(hinc);
  hratio->GetYaxis()->SetRangeUser(0,1.4);
  
  TCanvas *crat = new TCanvas ("crat","crat",1000,800);
  crat->cd();
 
  hratio->Draw();

  TLine * one = new TLine(5,1,80,1);
  TLine * ninety = new TLine(5,0.9,80,0.9); ninety->SetLineStyle(kDashed);
  
  one->Draw("same"); ninety->Draw("same");
  
  for (int i = 1; i <= hratio->GetNbinsX(); ++ i) {
    cout << "bin " << i << "'s factor is " << hratio->GetBinContent(i) << endl;
  }
  
  //crat->SaveAs("~/temp/jetptmin10_JP2_evtEtCut_trackcuts_NEF_on.pdf");
 
  return;
}
