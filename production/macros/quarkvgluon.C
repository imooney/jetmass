
#include "Plots.h"

using namespace std;

void quarkvgluon () {
  gStyle->SetPalette(kPastel);

  TFile *f = new TFile ("~/jetmass/production/macros/hists/hists_new.root","READ");
  
  TH2D* un_q_2D = (TH2D*) f->Get("m_v_pt_un_q");
  TH2D* un_g_2D = (TH2D*) f->Get("m_v_pt_un_g");
  TH2D* un_n_2D = (TH2D*) f->Get("m_v_pt_un_neither");
  TH2D* un_a_2D = (TH2D*) f->Get("m_v_pt_un");

  TH2D* dec_q_2D = (TH2D*) f->Get("m_v_pt_dec_q");
  TH2D* dec_g_2D = (TH2D*) f->Get("m_v_pt_dec_g");
  TH2D* dec_n_2D = (TH2D*) f->Get("m_v_pt_dec_neither");
  TH2D* dec_a_2D = (TH2D*) f->Get("m_v_pt_dec");

  const int nBins = 5;
  TAxis* yaxis = un_q_2D->GetYaxis();
  double ranges[nBins+1] = {(double) yaxis->FindBin(15),(double) yaxis->FindBin(20),(double) yaxis->FindBin(25),(double)yaxis->FindBin(30),(double)yaxis->FindBin(40),(double)yaxis->FindBin(60)};
  string pts[nBins+1] = {"15","20","25","30","40","60"}; 
  
  vector<TH1D*> un_q = Projection2D(un_q_2D,nBins,ranges,"x");
  vector<TH1D*> un_g = Projection2D(un_g_2D,nBins,ranges,"x");
  vector<TH1D*> un_n = Projection2D(un_n_2D,nBins,ranges,"x");
  vector<TH1D*> un_a = Projection2D(un_a_2D,nBins,ranges,"x");
  
  vector<TH1D*> dec_q = Projection2D(dec_q_2D,nBins,ranges,"x");
  vector<TH1D*> dec_g = Projection2D(dec_g_2D,nBins,ranges,"x");
  vector<TH1D*> dec_n = Projection2D(dec_n_2D,nBins,ranges,"x");
  vector<TH1D*> dec_a = Projection2D(dec_a_2D,nBins,ranges,"x");
  
  for (int i = 0; i < nBins; ++ i) {
    cout << "uq: " << un_q[i]->Integral() / (double) un_a[i]->Integral() << endl;
    cout << "ug: " << un_g[i]->Integral() / (double) un_a[i]->Integral() << endl;
    cout << "un: " << un_n[i]->Integral() / (double) un_a[i]->Integral() << endl;
    //cout << "a: " << un_a[i]->Integral() << endl;
    cout << "dq: " << dec_q[i]->Integral() / (double) dec_a[i]->Integral() << endl;
    cout << "dg: " << dec_g[i]->Integral() / (double) dec_a[i]->Integral() << endl;
    cout << "dn: " << dec_n[i]->Integral() / (double) dec_a[i]->Integral() << endl;
    
    Prettify1DwLineStyle(un_q[i], kDashed, 5,"M [GeV/c^{2}]","tmp",0,10,-1,-1);
    Prettify1DwLineStyle(un_g[i], kDashed, 5,"M [GeV/c^{2}]","tmp",0,10,-1,-1);
    Prettify1DwLineStyle(un_n[i], kDashed, 5,"M [GeV/c^{2}]","tmp",0,10,-1,-1);
    Prettify1DwLineStyle(un_a[i], kDashed, 5,"M [GeV/c^{2}]","tmp",0,10,-1,-1);
    
    Prettify1DwLineStyle(dec_q[i], kSolid, 5,"M [GeV/c^{2}]","tmp",0,10,-1,-1);
    Prettify1DwLineStyle(dec_g[i], kSolid, 5,"M [GeV/c^{2}]","tmp",0,10,-1,-1);
    Prettify1DwLineStyle(dec_n[i], kSolid, 5,"M [GeV/c^{2}]","tmp",0,10,-1,-1);
    Prettify1DwLineStyle(dec_a[i], kSolid, 5,"M [GeV/c^{2}]","tmp",0,10,-1,-1);
  }
  
  TLegend *t1 = new TLegend (0.6,0.4,0.85,0.7); t1->SetBorderSize(0);
  t1->AddEntry(un_q[0],"quark jets","l");
  t1->AddEntry(un_g[0],"gluon jets","l");
  t1->AddEntry(un_n[0],"neither","l");
  t1->AddEntry(un_a[0],"all jets","l");
  TLegend *t2 = new TLegend (0.6,0.4,0.85,0.7); t2->SetBorderSize(0);
  t2->AddEntry(dec_q[0],"quark jets","l");
  t2->AddEntry(dec_g[0],"gluon jets","l");
  t2->AddEntry(dec_n[0],"neither","l");
  t2->AddEntry(dec_a[0],"all jets","l");

  TH1D* hdummy = new TH1D("hdummy","tmp",1,0,10); hdummy->GetYaxis()->SetRangeUser(-1,-1);
  TH1D* hdummy2 = new TH1D("hdummy2","tmp",1,0,10); hdummy2->GetYaxis()->SetRangeUser(-1,-1);
  TLatex *p = new TLatex(); TLatex *slice = new TLatex();

  TCanvas *cun = new TCanvas("cun","cun",1400,1000); cun->cd();
  DivideCanvas(cun,"0",3,2);
  
  cun->cd(1); hdummy->Draw(); p = PanelTitle();
  for (int i = 0; i < nBins; ++ i) {
    cun->cd(i+2);
    un_q[i]->Draw("c PLC  PMC"); un_g[i]->Draw("csame PLC  PMC"); un_n[i]->Draw("csame PLC  PMC"); un_a[i]->Draw("csame PLC  PMC");
    if (i == 0) {t1->Draw("same");} slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }

  TCanvas *cdec = new TCanvas("cdec","cdec",1400,1000); cdec->cd();
  DivideCanvas(cdec,"0",3,2);
  
  cdec->cd(1); hdummy->Draw(); p = PanelTitle();
  for (int i = 0; i < nBins; ++ i) {
    cdec->cd(i+2);
    dec_q[i]->Draw("c PLC  PMC"); dec_g[i]->Draw("csame PLC  PMC"); dec_n[i]->Draw("csame PLC  PMC"); dec_a[i]->Draw("csame PLC  PMC");
    if (i == 0) {t2->Draw("same");} slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
  
  //  cun->SaveAs("~/jetmass/production/plots/qvg/qvg_py8_undecayed.pdf");
  // cdec->SaveAs("~/jetmass/production/plots/qvg/qvg_py8_decayed.pdf");
  
  return;
}
