#include "RooUnfoldResponse.h"
#include "Plots.h"

using namespace std;

void kinematic_eff() {
  const string path = "~/jetmass/out/matching/";
  const string file_wo = "full_w_o_bin_drop_w_o_pt_cut.root";
  const string file_w = "full_w_o_bin_drop_w_pt_cut.root";
  
  TFile *fwo = new TFile((path+file_wo).c_str(),"READ");
  TFile *fw = new TFile ((path+file_w).c_str(),"READ");
  TFile *fdat = new TFile("~/jetmass/macros/hists/hists_w_o_bin_drop.root","READ");
  
  RooUnfoldResponse *pt_wo_res = (RooUnfoldResponse*) fwo->Get("pt_response");
  RooUnfoldResponse *pt_w_res = (RooUnfoldResponse*) fw->Get("pt_res_coarse");

  TH1D* misses = (TH1D*) fwo->Get("misses");
  TH1D* miss_frac = (TH1D*) misses->Clone("miss_frac");
  
  TH1D* tru_wo = (TH1D*) pt_wo_res->Htruth(); tru_wo->SetTitle("");
  TH1D* gen_proj_wo = (TH1D*) pt_wo_res->Hresponse()->ProjectionY("gen_proj_wo", pt_wo_res->Hresponse()->GetXaxis()->GetFirst(), pt_wo_res->Hresponse()->GetXaxis()->GetLast());
  miss_frac->Divide(tru_wo);
  miss_frac->Sumw2();

  TH1D* hdiff = new TH1D("hdiff","",15,5,80);
  for (int i = 1; i < hdiff->GetNbinsX(); ++ i) {
    hdiff->SetBinContent(i, tru_wo->GetBinContent(i) - gen_proj_wo->GetBinContent(i));
    hdiff->SetBinError(i, gen_proj_wo->GetBinError(i));
  }
  hdiff->Sumw2();
  
  TCanvas *cmis = new TCanvas("cmis","cmis",1400,1000); cmis->cd(); cmis->SetLogy();
  Prettify1D(tru_wo,kBlue,kOpenCircle,3,kBlue,"p_{T} [GeV/c]","counts",5,80,1e-11,1e-1);
  Prettify1D(gen_proj_wo,kRed,kOpenCircle,3,kRed,"p_{T} [GeV/c]","counts",5,80,1e-11,1e-1);
  Prettify1D(misses,kViolet,kFullCircle,3,kViolet,"p_{T} [GeV/c]","counts",5,80,1e-11,1e-1);
  Prettify1D(miss_frac,kBlue,kOpenSquare,3,kBlue,"p_{T} [GeV/c]","counts",5,80,0,1.1);
  Prettify1D(hdiff,kViolet,kOpenStar,4,kViolet,"p_{T} [GeV/c]","counts",5,80,1e-11,1e-1);

  TLegend *tm = new TLegend(0.4,0.5,0.8,0.8); tm->SetBorderSize(0);
  tm->AddEntry(tru_wo,"Truth (match+miss)","p");
  tm->AddEntry(gen_proj_wo,"Gen.-level matches","p");
  tm->AddEntry(misses,"Misses","p");
  tm->AddEntry(hdiff,"Truth - gen.-level matches","p");

  // Upper plot will be in pad1                                                                                                                              
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined                                                                                              
  pad1->SetGridx();         // Vertical grid                                                                                                                
  pad1->Draw();             // Draw the upper pad: pad1                                                                                                     
  pad1->cd();               // pad1 becomes the current pad                                                                                                 
  pad1->SetLogy();
  tru_wo->SetStats(0);          // No statistics on upper plot                                                                                             
  tru_wo->Draw();                                                                                                              
  gen_proj_wo->Draw("same");
  misses->Draw("same");
  hdiff->Draw("same");
  tm->Draw("same");

  // lower plot will be in pad                                                                                                                               
  cmis->cd();          // Go back to the main canvas before defining pad2                                                                                  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid                                                                                                                        
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad   

  miss_frac->Draw("ep");       // Draw the ratio plot                                                                                                         
  TLine *one = new TLine(5,1,80,1); one->SetLineStyle(kDashed);
  TLine *half = new TLine(5,0.5,80,0.5); half->SetLineStyle(kDashed);
  TLine *quarter = new TLine(5,0.25,80,0.25); quarter->SetLineStyle(kDashed);

  one->Draw("same"); half->Draw("same"); quarter->Draw("same");

  // Y axis p plot settings                                                                                                                                 
  tru_wo->GetYaxis()->SetTitleSize(20);
  tru_wo->GetYaxis()->SetTitleFont(43);
  tru_wo->GetYaxis()->SetTitleOffset(1.55);

  // Ratio plot settings                                                                                                                              
  miss_frac->SetTitle(""); // Remove the ratio title                                                                                                               
  // Y axis ratio plot settings
  miss_frac->GetYaxis()->SetTitle("misses / truth");
  miss_frac->GetYaxis()->SetNdivisions(505);
  miss_frac->GetYaxis()->SetTitleSize(20);
  miss_frac->GetYaxis()->SetTitleFont(43);
  miss_frac->GetYaxis()->SetTitleOffset(1.55);
  miss_frac->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)                                                                    
  miss_frac->GetYaxis()->SetLabelSize(15);

  // X axis ratio plot settings                                                                                                                       
  miss_frac->GetXaxis()->SetTitleSize(20);
  miss_frac->GetXaxis()->SetTitleFont(43);
  miss_frac->GetXaxis()->SetTitleOffset(4.);
  miss_frac->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)                                                                     
  miss_frac->GetXaxis()->SetLabelSize(15);

  cout << "check: " << misses->GetBinContent(5) << " " << tru_wo->GetBinContent(5) - gen_proj_wo->GetBinContent(5) << endl;
  cout << "check: " << misses->GetBinContent(6) << " " << tru_wo->GetBinContent(6) - gen_proj_wo->GetBinContent(6) << endl;
  cout << "check: " << misses->GetBinContent(7) << " " << tru_wo->GetBinContent(7) - gen_proj_wo->GetBinContent(7) << endl;

  for (int i = 1; i < tru_wo->GetNbinsX(); ++ i) {
    cout << "scale factor for " << 5*i << " to " << 5*(i+1) << " GeV: " << tru_wo->GetBinContent(i) / (double) gen_proj_wo->GetBinContent(i) << endl;
  }
  
  //gen_proj_wo->Draw(); tru_wo->Draw("same"); misses->Draw("same"); miss_frac->Draw("same"); hdiff->Draw("same"); tm->Draw("same");

  cmis->SaveAs("~/misses_check.pdf");
  
  TH1D *pt_data = (TH1D*) fdat->Get("pt_d_pt15");

  RooUnfoldBayes *unfold1D = new RooUnfoldBayes(pt_w_res, pt_data, 4, false, "unfold1D","");
  TH1D *reco = (TH1D*) unfold1D->Hreco((RooUnfold::ErrorTreatment) 3);

  TH2D* hpt_wo = (TH2D*) pt_wo_res->Hresponse();
  TH2D* hpt_w = (TH2D*) pt_w_res->Hresponse();
  
  TAxis* hpt_wo_x = (TAxis*) hpt_wo->GetXaxis();
  TAxis* hpt_wo_y = (TAxis*) hpt_wo->GetYaxis();
  TAxis* hpt_w_x = (TAxis*) hpt_w->GetXaxis();
  TAxis* hpt_w_y = (TAxis*) hpt_w->GetYaxis();
  
  //normalizing the rows:
  vector<double> wo_ints; vector<double> w_ints;
  
  for (int j = 1; j <= hpt_wo->GetNbinsY(); ++ j) {
    wo_ints.push_back(hpt_wo->Integral(1,hpt_wo->GetNbinsX(),hpt_wo_y->FindBin(5*j),hpt_wo_y->FindBin(5*j)));
  }
  for (int j = 1; j <= hpt_w->GetNbinsY(); ++ j) {
    w_ints.push_back(hpt_w->Integral(1,hpt_w->GetNbinsX(),hpt_w_y->FindBin(5*j),hpt_w_y->FindBin(5*j)));
  }

  for (int j = 1; j <= hpt_wo->GetNbinsY(); ++ j) {
    for (int i = 1; i <= hpt_wo->GetNbinsX(); ++ i) {
      if (wo_ints[j-1] != 0) {
	hpt_wo->SetBinContent(i,j,hpt_wo->GetBinContent(i,j)/(double)wo_ints[j-1]);
      }
    }
  }
  for (int j = 1; j <= hpt_w->GetNbinsY(); ++ j) {
    for (int i = 1; i <= hpt_w->GetNbinsX(); ++ i) {
      if (w_ints[j-1] != 0) {
	hpt_w->SetBinContent(i,j,hpt_w->GetBinContent(i,j)/(double)w_ints[j-1]);
      }
    }
  }
  
  TCanvas *cres_wo = new TCanvas("cres_wo","cres_wo",1400,1000); cres_wo->SetLogz();
  TCanvas *cres_w = new TCanvas("cres_w","cres_w",1400,1000); cres_w->SetLogz();
  
  TLine *diag_wo = new TLine(5,5,80,80);
  TLine *diag_w = new TLine(15,15,60,60);

  hpt_wo->SetTitle(""); hpt_w->SetTitle("");
  hpt_wo->GetXaxis()->SetTitle("p^{det.}_{T} [GeV/c]"); hpt_w->GetXaxis()->SetTitle("p^{det.}_{T} [GeV/c]");
  hpt_wo->GetYaxis()->SetTitle("p^{gen.}_{T} [GeV/c]"); hpt_w->GetYaxis()->SetTitle("p^{gen.}_{T} [GeV/c]");
  
  cres_wo->cd(); hpt_wo->Draw("colz"); diag_wo->Draw("same");
  cres_w->cd(); hpt_w->Draw("colz"); diag_w->Draw("same");

  cres_wo->SaveAs("~/res_no_cuts_normalized_rows.pdf"); cres_w->SaveAs("~/res_w_cuts_normalized_rows.pdf");
  
  //projecting and taking the ratio
  
  vector<TH1D*> pt_wo; vector<TH1D*> pt_w;
  
  const int nBins = 1;
  vector<int> bins_w = {1,/*2,3,4,6,*/10}; vector<int> bins_wo = {3,/*4,5,6,8,*/12};
  int pts[nBins+1] = {15,/*20,25,30,40,*/60};
  
  for (int i = 0; i < bins_w.size() - 1; ++ i) {
    pt_wo.push_back((TH1D*) hpt_wo->ProjectionY(("pt"+to_string(pts[i])+to_string(pts[i+1])+"_wo").c_str(),bins_wo[i],bins_wo[i+1]));
    pt_w.push_back((TH1D*) hpt_w->ProjectionY(("pt"+to_string(pts[i])+to_string(pts[i+1])+"_w").c_str(),bins_w[i],bins_w[i+1]));
    // cout << pt_wo[i-1]->Integral() << " " << pt_w[i-1]->Integral() << pt_wo[i-1]->Integral() / (double) pt_w[i-1]->Integral() << endl;
  }
  cout << "test! wo, w" << endl;
    cout << "test! 40 - 45 " << pt_wo[0]->Integral(pt_wo[0]->GetXaxis()->FindBin(40), pt_wo[0]->GetXaxis()->FindBin(45)) << " " << pt_w[0]->Integral(pt_w[0]->GetXaxis()->FindBin(40), pt_w[0]->GetXaxis()->FindBin(45)) << endl;
  cout << "test! 45 - 50 " << pt_wo[0]->Integral(pt_wo[0]->GetXaxis()->FindBin(45), pt_wo[0]->GetXaxis()->FindBin(50)) << " " << pt_w[0]->Integral(pt_w[0]->GetXaxis()->FindBin(45), pt_w[0]->GetXaxis()->FindBin(50)) << endl;
  cout << "test! 50 - 55 " << pt_wo[0]->Integral(pt_wo[0]->GetXaxis()->FindBin(50), pt_wo[0]->GetXaxis()->FindBin(55)) << " " << pt_w[0]->Integral(pt_w[0]->GetXaxis()->FindBin(50), pt_w[0]->GetXaxis()->FindBin(55)) << endl;
  cout << "test! 55 - 60 " << pt_wo[0]->Integral(pt_wo[0]->GetXaxis()->FindBin(55), pt_wo[0]->GetXaxis()->FindBin(60)) << " " << pt_w[0]->Integral(pt_w[0]->GetXaxis()->FindBin(55), pt_w[0]->GetXaxis()->FindBin(60)) << endl;
  cout << "test! 60 - 65 " << pt_wo[0]->Integral(pt_wo[0]->GetXaxis()->FindBin(60), pt_wo[0]->GetXaxis()->FindBin(65)) << " " << pt_w[0]->Integral(pt_w[0]->GetXaxis()->FindBin(60), pt_w[0]->GetXaxis()->FindBin(65)) << endl;
  cout << "test! 65 - 70 " << pt_wo[0]->Integral(pt_wo[0]->GetXaxis()->FindBin(65), pt_wo[0]->GetXaxis()->FindBin(70)) << " " << pt_w[0]->Integral(pt_w[0]->GetXaxis()->FindBin(65), pt_w[0]->GetXaxis()->FindBin(70)) << endl;  
  cout << "test! 70 - 75 " << pt_wo[0]->Integral(pt_wo[0]->GetXaxis()->FindBin(70), pt_wo[0]->GetXaxis()->FindBin(75)) << " " << pt_w[0]->Integral(pt_w[0]->GetXaxis()->FindBin(70), pt_w[0]->GetXaxis()->FindBin(75)) << endl;

  vector<TH1D*> ratios;
  for (int i = 0; i < pt_w.size(); ++ i) {
    ratios.push_back((TH1D*)pt_w[i]->Clone(("r"+to_string(pts[i])+to_string(pts[i+1])).c_str()));
    ratios[i]->Divide(pt_wo[i]);
    Prettify1D(ratios[i],kBlue,kOpenCircle,2,kBlue,"p^{gen.}_{T} [GeV/c]", "w/ cuts / w/o cuts",5,80,0.9,2);
    ratios[i]->SetTitle("");
    
    //    ratios[i]->GetYaxis()->SetMoreLogLabels();
    //ratios[i]->GetYaxis()->SetNdivisions(999, kFALSE);
  }
  
  TCanvas *c = new TCanvas("c","c",1400,1000);// c->SetLogy();
  //DivideCanvas(c,"y",3,2);
  
  TLatex *slice = new TLatex(); TLatex *p = new TLatex();

  TLine *t = new TLine(5,1,80,1); t->SetLineStyle(kDashed);

  TH1D* hdummy = new TH1D("hdummy","",1,5,80); hdummy->GetYaxis()->SetRangeUser(0.9,2.1);
  
  //c->cd(1); hdummy->Draw(); p = PanelTitle(); t->Draw("same");
  
  for (int i = 0; i < ratios.size(); ++ i) {
    /*c->cd(i+2);*/ ratios[i]->Draw(); t->Draw("same");/*slice->DrawLatexNDC(0.4,0.8,(to_string(pts[i])+" < p^{det.}_{T} < "+to_string(pts[i+1])+" GeV/c").c_str()); */
  } 

  c->SaveAs("~/jetmass/plots/kinematicefficiency_wholerange.pdf");
  
  TCanvas *ccomp = new TCanvas("ccomp","ccomp",1400,1000);
  ccomp->cd(); ccomp->SetLogy();
  Prettify1D(pt_data,kRed,kOpenCircle,3,kRed,"p_{T} [GeV/c]", "counts",15,60,1e-1,1e8);
  Prettify1D(reco,kBlue,kOpenCircle,3,kBlue,"p_{T} [GeV/c]", "counts",15,60,1e-1,1e8);

  TLegend *leg = new TLegend(0.45,0.6,0.8,0.8); leg->SetBorderSize(0);
  leg->AddEntry(pt_data,"raw data","p");
  leg->AddEntry(reco,"unfolded data","p");
  
  pt_data->Draw(); reco->Draw("same"); leg->Draw("same");

  ccomp->SaveAs("~/jetmass/plots/unfolding/counts_compare_1D.pdf");
  
  TAxis* daxis = pt_data->GetXaxis(); TAxis* raxis = reco->GetXaxis();
  
  cout << "Counts from 15 - 20 GeV - data, unfolded: " << reco->Integral(raxis->FindBin(15),raxis->FindBin(20)) /(double) pt_data->Integral(daxis->FindBin(15),daxis->FindBin(20)) << endl;
  cout << "Counts from 20 - 25 GeV - data, unfolded: " << reco->Integral(raxis->FindBin(20),raxis->FindBin(25)) /(double) pt_data->Integral(daxis->FindBin(20),daxis->FindBin(25)) << endl;
  cout << "Counts from 25 - 30 GeV - data, unfolded: " << reco->Integral(raxis->FindBin(25),raxis->FindBin(30)) /(double) pt_data->Integral(daxis->FindBin(25),daxis->FindBin(30)) << endl;
  cout << "Counts from 30 - 40 GeV - data, unfolded: " << reco->Integral(raxis->FindBin(30),raxis->FindBin(40)) /(double) pt_data->Integral(daxis->FindBin(30),daxis->FindBin(40)) << endl;
  cout << "Counts from 40 - 60 GeV - data, unfolded: " << reco->Integral(raxis->FindBin(40),raxis->FindBin(60)) /(double) pt_data->Integral(daxis->FindBin(40),daxis->FindBin(60)) << endl;
  return;
}
