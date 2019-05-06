
#include "Plots.h"

using namespace std;


//note: the vector hists is assumed to follow the ordering: p8 decayed, p8 undecayed, p6 (same order), h7 (same order)
void plot_from2D(string obs, string xtitle, string units, vector<TH2D*> hists, vector<double> axes, vector<double> lcoords) {
  const int nBins = 5;
  TAxis* yax = (TAxis*) hists[0]->GetYaxis();
  double bins[nBins+1] = {(double)yax->FindBin(15),(double)yax->FindBin(20),(double)yax->FindBin(25),(double)yax->FindBin(30),(double)yax->FindBin(40),(double)yax->FindBin(60)};
  string pts[nBins+1] = {"15","20","25","30","40","60"};
  
  vector<TH1D*> obs_p8on = Projection2D(hists[0],nBins,bins,"x");
  vector<TH1D*> obs_p8off = Projection2D(hists[1],nBins,bins,"x");
  vector<TH1D*> obs_p6on = Projection2D(hists[2],nBins,bins,"x");
  vector<TH1D*> obs_p6off = Projection2D(hists[3],nBins,bins,"x");
  vector<TH1D*> obs_h7on = Projection2D(hists[4],nBins,bins,"x");
  vector<TH1D*> obs_h7off = Projection2D(hists[5],nBins,bins,"x");

  for(int i = 0; i < nBins; ++i) {
    Prettify1DwLineStyle(obs_p8on[i],kBlue,kSolid,5,(xtitle + " " + units).c_str(), ("1/N_{j} dN/d"+xtitle).c_str(),axes[0],axes[1],axes[2],axes[3]);
    Prettify1DwLineStyle(obs_p8off[i],kBlue,kDashed,5,(xtitle + " " + units).c_str(), ("1/N_{j} dN/d"+xtitle).c_str(),axes[0],axes[1],axes[2],axes[3]);
    Prettify1DwLineStyle(obs_p6on[i],kRed,kSolid,5,(xtitle + " " + units).c_str(), ("1/N_{j} dN/d"+xtitle).c_str(),axes[0],axes[1],axes[2],axes[3]);
    Prettify1DwLineStyle(obs_p6off[i],kRed,kDashed,5,(xtitle + " " + units).c_str(), ("1/N_{j} dN/d"+xtitle).c_str(),axes[0],axes[1],axes[2],axes[3]);
    Prettify1DwLineStyle(obs_h7on[i],kCyan,kSolid,5,(xtitle + " " + units).c_str(), ("1/N_{j} dN/d"+xtitle).c_str(),axes[0],axes[1],axes[2],axes[3]);
    Prettify1DwLineStyle(obs_h7off[i],kCyan,kDashed,5,(xtitle + " " + units).c_str(), ("1/N_{j} dN/d"+xtitle).c_str(),axes[0],axes[1],axes[2],axes[3]);
  }
  
  TH1D* hdummy = new TH1D(("hdummy_"+obs).c_str(),"",1,axes[0],axes[1]); hdummy->GetYaxis()->SetRangeUser(axes[2],axes[3]);
  
  TLegend *t1 = new TLegend(lcoords[0],lcoords[1],lcoords[2],lcoords[3]); t1->SetBorderSize(0);
  t1->AddEntry(obs_p8on[0],"P8 decayed","l");
  t1->AddEntry(obs_p8off[0],"P8 undecayed","l");
  t1->AddEntry(obs_p6on[0],"P6 decayed","l");
  TLegend *t2 = new TLegend(lcoords[0],lcoords[1],lcoords[2],lcoords[3]); t2->SetBorderSize(0);
  t2->AddEntry(obs_p6off[0],"P6 undecayed","l");
  t2->AddEntry(obs_h7on[0],"H7 decayed","l");
  t2->AddEntry(obs_h7off[0],"H7 undecayed","l");

  TLatex *p = new TLatex(); TLatex *slice = new TLatex();
  TCanvas *c = new TCanvas (("c"+obs).c_str(),("c"+obs).c_str(),1400,1000); c->cd();
  DivideCanvas(c,"0",3,2);
  
  c->cd(1); hdummy->Draw(); p = PanelTitle();
  for (int i = 0; i < nBins; ++ i) {
    c->cd(i+2);
    obs_p8on[i]->Draw("c"); obs_p8off[i]->Draw("csame");
    obs_p6on[i]->Draw("csame"); obs_p6off[i]->Draw("csame");
    obs_h7on[i]->Draw("csame"); obs_h7off[i]->Draw("csame");
    if (i == 0) {t1->Draw("same");} slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    if (i == 1) {t2->Draw("same");}
  }
  
  c->SaveAs(("~/jetmass/plots/decays/"+obs+"_allsim.pdf").c_str());
  
  return;
}


void decayscompare() {
  
  TFile *f = new TFile("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved.root","READ");
  
  vector<TH2D*> mvpt = {(TH2D*)f->Get("mvpt_p8on"),(TH2D*)f->Get("mvpt_p8off"),(TH2D*)f->Get("mvpt_p6on"),
			(TH2D*)f->Get("mvpt_p6off"),(TH2D*)f->Get("mvpt_h7on"),(TH2D*)f->Get("mvpt_h7off")};
  vector<TH2D*> mgvpt = {(TH2D*)f->Get("mgvpt_p8on"),(TH2D*)f->Get("mgvpt_p8off"),(TH2D*)f->Get("mgvpt_p6on"),
			(TH2D*)f->Get("mgvpt_p6off"),(TH2D*)f->Get("mgvpt_h7on"),(TH2D*)f->Get("mgvpt_h7off")};
  vector<TH2D*> zgvpt = {(TH2D*)f->Get("zgvpt_p8on"),(TH2D*)f->Get("zgvpt_p8off"),(TH2D*)f->Get("zgvpt_p6on"),
			(TH2D*)f->Get("zgvpt_p6off"),(TH2D*)f->Get("zgvpt_h7on"),(TH2D*)f->Get("zgvpt_h7off")};
  vector<TH2D*> rgvpt = {(TH2D*)f->Get("rgvpt_p8on"),(TH2D*)f->Get("rgvpt_p8off"),(TH2D*)f->Get("rgvpt_p6on"),
			(TH2D*)f->Get("rgvpt_p6off"),(TH2D*)f->Get("rgvpt_h7on"),(TH2D*)f->Get("rgvpt_h7off")};
  
  vector<double> m_axes = {0,10,0,0.5};
  vector<double> zr_axes = {0,0.5,0,8};
  
  vector<double> lcoords = {0.55,0.55,0.85,0.75};
  vector<double> lrcoords = {0.15,0.55,0.45,0.75};
  
  plot_from2D("m", "M", "[GeV/c^{2}]", mvpt, m_axes, lcoords);
  plot_from2D("mg", "M_{g}", "[GeV/c^{2}]", mgvpt, m_axes, lcoords);
  plot_from2D("zg","z_{g}", "", zgvpt, zr_axes, lcoords);
  plot_from2D("rg","R_{g}", "", rgvpt, zr_axes, lrcoords);
  
  return;
}
