//Isaac Mooney - 9/17/2018
//This macro plots the p6 mass detector resolution as a function of p6 mass, for various bins of p6 pT. This plot is "corrected" to a resolution of 1 for all mass. I then re-run the code with the "corrected" resolution and plot the new mu/sigma as a function of p6 mass for various bins of p6 pT. It is fit with a function, and this function is used to "correct" the p8 mass to detector-level, after which we construct a p8 response and unfold.

#include <iostream>
#include <string>
#include <math.h>
  
#include "Plots.h"

using namespace std;

TF1* fit_funcs(TProfile* mres_profy, const std::string flag, TFile* funcs, const int i, const bool groom) {
  std::string groomstring; if (groom) {groomstring = "groom_";} else {groomstring = "";}
  std::string fit_name = (groomstring+"scale_"+flag+"_"+to_string(i)).c_str();
  
  TF1 *res; 
  //ungroomed!
  if (!groom && flag == "lo" && i == 2) {
    res = new TF1(fit_name.c_str(),"[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",1,10);
    res->SetParameters(4.31806e-01,6.15115e-01,-2.26141e-01,3.62100e-02,-2.05147e-03);
  }
  else if (!groom && flag == "hi" && i == 3) { res = new TF1(fit_name.c_str(),"pol6",1,10);}
  else if (!groom) { res = new TF1(fit_name.c_str(), "pol4",1,10);}

  //groomed!
  if (groom && flag == "hi" && (i == 2 || i == 3)) {
    res = new TF1(fit_name.c_str(),"[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",0,10);
    res->SetParameters(5.47610e-01,4.33022e-01,-1.30284e-01,1.71122e-02,-7.80361e-04);
  }
  else if (groom && flag == "lo" && i == 1) {
    res = new TF1(fit_name.c_str(),"pol9",0,10);
  }
  else if (groom) { res = new TF1(fit_name.c_str(), "pol4",0,10);}
  //TF1 *res = new TF1(fit_name.c_str(),"[0]*TMath::ATan([1]*x + [2])+[3]",0,10);

  //if (flag == "lo" && i == 2) {mres_profy->Fit("chebyshev9","MEI");}
  mres_profy->Fit(fit_name.c_str(),"MEI");
  res->Draw("same");
  
  funcs->cd();
  res->Write();
  
  return res;
}


void p6res(TH3D* res_m_pt_p6, TH3D* res_gen_m_pt_p6, const std::string out, const std::string filetype, const std::string flag, TFile* funcs, const bool fitting, const bool groom) {
  TLine * one = new TLine (0,1,10,1); one->SetLineStyle(kDashed);
  TLine * fivepercent = new TLine (0,1.05,10,1.05); fivepercent->SetLineStyle(kDashed);
  TLine * minusfivepercent = new TLine (0,0.95,10,0.95); minusfivepercent->SetLineStyle(kDashed);
  // funcs->cd();
  
  TCanvas *cscale = MakeCanvas(((string) "cscale"+ (string) res_m_pt_p6->GetName()).c_str(),"0",1200,800);
  DivideCanvas(cscale, "0", 3, 2);

  TCanvas *cres;
  if (!fitting) {
    cres = MakeCanvas(((string) "cres"+ (string) res_m_pt_p6->GetName()).c_str(),"0",1200,800);
    DivideCanvas(cres, "0", 4, 2);
  }
  
  const int nBins_z = 5; const int nBins_z_gen = 7; const int nBins_y = 10;
  //reminder: ranges is full of bin NUMBERS not VALUES
  int ranges_z[nBins_z+1] = {1,2,3,4,6,10};
  int ranges_z_gen[nBins_z_gen+1] = {1,2,3,4,5,6,8,12};
  double corresp_pts[nBins_z+1] = {15,20,25,30,40,60}; 
  double corresp_pts_gen[nBins_z_gen+1] = {5,10,15,20,25,30,40,60};
  double ranges_y[nBins_y+1] = {1,2,3,4,5,6,7,8,9,10};
  double masses[nBins_y+1] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
  double mass_err[nBins_y+1] = {};
  
  std::vector<TH2D*> mres_v_m = Projection3D (res_m_pt_p6, nBins_z, ranges_z, "yx");  
  std::vector<TH2D*> mres_v_m_gen = Projection3D(res_gen_m_pt_p6, nBins_z_gen, ranges_z_gen, "yx");
  
  TLegend *tslices[nBins_z]; TLegend *tslices_gen[nBins_z];
  for (int i = 0; i < nBins_z; ++ i) {
    tslices[i] = SliceLegend(((to_string(corresp_pts[i])).substr(0,2) + " < p_{T}^{det} < " + (to_string(corresp_pts[i + 1])).substr(0,2) + " GeV/c").c_str(), 0.2,0.15,0.9,0.3);
  }
  for (int i = 0; i < nBins_z_gen; ++i) {
    tslices_gen[i] = SliceLegend(((to_string(corresp_pts_gen[i])).substr(0,2) + " < p_{T}^{gen} < " + (to_string(corresp_pts_gen[i + 1])).substr(0,2) + " GeV/c").c_str(), 0.2,0.75,0.9,0.85);
  }
  
  string eta_low = "-0.6"; string eta_high = "0.6";
  if (((string)res_m_pt_p6->GetName()).find("lo") != string::npos) {eta_high = "-0.2";}
  else if (((string)res_m_pt_p6->GetName()).find("mid") != string::npos) {eta_low = "-0.2"; eta_high = "0.2";}
  else {eta_low = "0.2";}
  cscale->cd(1); TLatex *t = new TLatex();
  t->SetTextAlign(11);
  //  t->SetTextFont(63);                                                                                                                           
  t->SetTextSizePixels(26);
  t->DrawLatex(0.1,0.9, "pp 200 GeV run12 JP2");
  t->DrawLatex(0.1,0.75, "anti-k_{T}, R = 0.4");
  t->DrawLatex(0.1,0.6, ("Ch+Ne jets, " + eta_low + " < #eta < " + eta_high).c_str());
  
  if (!fitting) {
    cres->cd(1);
    t->DrawLatex(0.1,0.9, "pp 200 GeV run12 JP2");
    t->DrawLatex(0.1,0.75, "anti-k_{T}, R = 0.4");
    t->DrawLatex(0.1,0.6, ("Ch+Ne jets, " + eta_low + " < #eta < " + eta_high).c_str());
  }
  
  for(int i = 0; i < mres_v_m.size(); ++ i) {
    TProfile* mres_profy = Profile2D(mres_v_m[i], "y");
    PrettifyTProfile(mres_profy,kBlue,kFullCircle,1,kBlue,"M^{det} [GeV/c^{2}]","<M^{det} / M^{gen}>", -1,-1,0.5,1.5);
    cscale->cd(i+2);
    mres_profy->Draw();
    tslices[i]->Draw("same");
    one->Draw("same"); fivepercent->Draw("same"); minusfivepercent->Draw("same");
    
    if (fitting) {
      fit_funcs(mres_profy,flag,funcs,i, groom);
    }
  }
  
   if (!fitting) {//i.e. we've already done the jet energy scale correction
     for(int i = 0; i < mres_v_m_gen.size(); ++ i) {
       std::vector<TH1D*> mres_projxs = Projection2D(mres_v_m_gen[i], nBins_y, ranges_y, "x");
       double sigmameans[nBins_y+1]; double errs[nBins_y+1];
       for (int j = 0; j < mres_projxs.size(); ++ j) {
	 sigmameans[j] = mres_projxs[j]->GetRMS() / (double) mres_projxs[j]->GetMean();
	 double rel_mean_err = mres_projxs[j]->GetMeanError() / (double) mres_projxs[j]->GetMean();
	 double rel_rms_err = mres_projxs[j]->GetRMSError() / (double) mres_projxs[j]->GetRMS();
	 errs[j] = sqrt((rel_mean_err*rel_mean_err) + (rel_rms_err*rel_rms_err));
       }
       
       TGraphErrors *m_v_res = new TGraphErrors(nBins_y,masses,sigmameans,mass_err,errs); m_v_res->SetTitle("");
       PrettifyTGraph(m_v_res, kBlue, kFullCircle, 1, kBlue, "M^{gen} [GeV/c^{2}]", "#sigma/#mu (M^{det} / M^{gen})", -1,-1,0,0.5);
       cres->cd(i+2);
       m_v_res->Draw("AP");
       tslices_gen[i]->Draw("same");
       
       std::string fit_name = ("res_"+flag+"_"+to_string(i)).c_str();
       TF1 *sigmafit = new TF1(fit_name.c_str(),"[0]"/*"sqrt([0]*[0]+([1]/sqrt(x))*([1]/sqrt(x))+([2]/x)*([2]/x))"*/,0,10);
       m_v_res->Fit(fit_name.c_str(),"MI");
       sigmafit->Draw("same");
     }
   }
   
   std::string filename;
   if (groom && !fitting) {filename = "corrected_JMg";} else if (!fitting) {filename = "corrected_JM";}
   if (groom && fitting) {filename = "JMg";} else if (fitting) {filename = "JM";}
   
   cscale->SaveAs((out+filename+"S"+flag+filetype).c_str());
   if (!fitting){ cres->SaveAs((out+filename+"R"+flag+filetype).c_str()); }
   return;
}



void detector_mass_res () {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILES&STYLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                                       
  string dir = "~/jetmass/";
  string in = "macros/hists/";
  string funcpath = "macros/funcs/";
  string infile = "hists.root";
  string funcs = "funcs.root";
  string out = "~/jetmass/plots/Y12_detector_mass_response_correction/";
  string filetype = ".pdf";
  string flag1 = "full";
  string flag2 = "incl";

  TFile* inFile = new TFile( (dir + in + infile).c_str(), "READ");
  TFile* funcFile = new TFile( (dir + funcpath + funcs).c_str(), "RECREATE");

  //MASS
  TH3D* res_m_pt_p6_lo = (TH3D*) inFile->Get("mass_res_v_det_m_pt_eta_lo");
  TH3D* res_m_pt_p6_mid = (TH3D*) inFile->Get("mass_res_v_det_m_pt_eta_mid");
  TH3D* res_m_pt_p6_hi = (TH3D*) inFile->Get("mass_res_v_det_m_pt_eta_hi");
  
  TH3D* res_m_pt_p6_lo_corr = (TH3D*) inFile->Get("mass_res_v_det_m_pt_eta_lo_corr");
  TH3D* res_m_pt_p6_mid_corr = (TH3D*) inFile->Get("mass_res_v_det_m_pt_eta_mid_corr");
  TH3D* res_m_pt_p6_hi_corr = (TH3D*) inFile->Get("mass_res_v_det_m_pt_eta_hi_corr");
  
  TH3D* res_gen_m_pt_p6_lo_corr = (TH3D*) inFile->Get("mass_res_v_gen_m_pt_eta_lo_corr");
  TH3D* res_gen_m_pt_p6_mid_corr = (TH3D*) inFile->Get("mass_res_v_gen_m_pt_eta_mid_corr");
  TH3D* res_gen_m_pt_p6_hi_corr = (TH3D*) inFile->Get("mass_res_v_gen_m_pt_eta_hi_corr");
  
  //GROOMED MASS
  TH3D* res_mg_pt_p6_lo = (TH3D*) inFile->Get("mg_res_v_det_mg_pt_eta_lo");
  TH3D* res_mg_pt_p6_mid = (TH3D*) inFile->Get("mg_res_v_det_mg_pt_eta_mid");
  TH3D* res_mg_pt_p6_hi = (TH3D*) inFile->Get("mg_res_v_det_mg_pt_eta_hi");
  
  TH3D* res_mg_pt_p6_lo_corr = (TH3D*) inFile->Get("mg_res_v_det_mg_pt_eta_lo_corr");
  TH3D* res_mg_pt_p6_mid_corr = (TH3D*) inFile->Get("mg_res_v_det_mg_pt_eta_mid_corr");
  TH3D* res_mg_pt_p6_hi_corr = (TH3D*) inFile->Get("mg_res_v_det_mg_pt_eta_hi_corr");
  
  TH3D* res_gen_mg_pt_p6_lo_corr = (TH3D*) inFile->Get("mg_res_v_gen_mg_pt_eta_lo_corr");
  TH3D* res_gen_mg_pt_p6_mid_corr = (TH3D*) inFile->Get("mg_res_v_gen_mg_pt_eta_mid_corr");
  TH3D* res_gen_mg_pt_p6_hi_corr = (TH3D*) inFile->Get("mg_res_v_gen_mg_pt_eta_hi_corr");
  
  
  TH3D* dummy = new TH3D("dummy","",1,0,1,1,0,1,1,0,1);
  
  
  p6res(res_m_pt_p6_lo, dummy, out, filetype, "lo",funcFile, 1,0);
  p6res(res_m_pt_p6_mid, dummy, out, filetype, "mid",funcFile, 1,0);
  p6res(res_m_pt_p6_hi, dummy, out, filetype, "hi",funcFile, 1,0);
  
  p6res(res_m_pt_p6_lo_corr, res_gen_m_pt_p6_lo_corr, out, filetype, "lo",funcFile, 0,0);
  p6res(res_m_pt_p6_mid_corr, res_gen_m_pt_p6_mid_corr, out, filetype, "mid",funcFile, 0,0);
  p6res(res_m_pt_p6_hi_corr, res_gen_m_pt_p6_hi_corr, out, filetype, "hi",funcFile, 0,0);
  
  p6res(res_mg_pt_p6_lo, dummy, out, filetype, "lo",funcFile, 1,1);
  p6res(res_mg_pt_p6_mid, dummy, out, filetype, "mid",funcFile, 1,1);
  p6res(res_mg_pt_p6_hi, dummy, out, filetype, "hi",funcFile, 1,1);
    
  p6res(res_mg_pt_p6_lo_corr, res_gen_mg_pt_p6_lo_corr, out, filetype, "lo",funcFile, 0,1);
  p6res(res_mg_pt_p6_mid_corr, res_gen_mg_pt_p6_mid_corr, out, filetype, "mid",funcFile, 0,1);
  p6res(res_mg_pt_p6_hi_corr, res_gen_mg_pt_p6_hi_corr, out, filetype, "hi",funcFile, 0,1);
  
  
  return;
}
