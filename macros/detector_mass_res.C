//Isaac Mooney - 9/17/2018
//This macro plots the p6 mass detector resolution as a function of p6 mass, for various bins of p6 pT. This plot is "corrected" to a resolution of 1 for all mass. I then re-run the code with the "corrected" resolution and plot the new mu/sigma as a function of p6 mass for various bins of p6 pT. It is fit with a function, and this function is used to "correct" the p8 mass to detector-level, after which we construct a p8 response and unfold.

#include <iostream>
#include <string>
#include <math.h>
  
#include "Plots.h"

using namespace std;

void p6res(TH3D* res_m_pt_p6, const std::string out, const std::string filetype, const std::string flag, TFile* funcs, const bool fitting) {
  TLine * one = new TLine (0,1,10,1); one->SetLineStyle(kDashed);
  funcs->cd();
  
  TCanvas *cres = MakeCanvas(((string) "cres"+ (string) res_m_pt_p6->GetName()).c_str(),"0",1200,800);
  DivideCanvas(cres, "0", 3, 2);
  
  const int nBins_z = 5; const int nBins_y = 20;
  //reminder: ranges is full of bin NUMBERS not VALUES
  int ranges_z[nBins_z+1] = {1,2,3,4,6,10}; 
  double corresp_pts[nBins_z+1] = {15,20,25,30,40,60}; 
  double ranges_y[nBins_y+1] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  double masses[nBins_y+1] = {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};
  double mass_err[nBins_y+1] = {};
  
  std::vector<TH2D*> mres_v_m = Projection3D (res_m_pt_p6, nBins_z, ranges_z, "yx");  

  TLegend *tslices[nBins_z];
  for (int i = 0; i < nBins_z; ++ i) {
    tslices[i] = SliceLegend(((to_string(corresp_pts[i])).substr(0,2) + " < p_{T}^{det} < " + (to_string(corresp_pts[i + 1])).substr(0,2) + " GeV/c").c_str(), 0.2,0.15,0.9,0.3);
  }
  
  string eta_low = "-0.6"; string eta_high = "0.6";
  if (((string)res_m_pt_p6->GetName()).find("lo") != string::npos) {eta_high = "-0.2";}
  else if (((string)res_m_pt_p6->GetName()).find("mid") != string::npos) {eta_low = "-0.2"; eta_high = "0.2";}
  else {eta_low = "0.2";}
  cres->cd(1); TLatex *t = new TLatex();
  t->SetTextAlign(11);
  //  t->SetTextFont(63);                                                                                                                           
  t->SetTextSizePixels(26);
  t->DrawLatex(0.1,0.9, "pp 200 GeV run12 MinBias");
  t->DrawLatex(0.1,0.75, "anti-k_{T}, R = 0.4");
  t->DrawLatex(0.1,0.6, ("Ch+Ne jets, " + eta_low + " < #eta < " + eta_high).c_str());
  
  for(int i = 0; i < mres_v_m.size(); ++ i) {
    TProfile* mres_profy = Profile2D(mres_v_m[i], "y");
    PrettifyTProfile(mres_profy,kBlue,kFullCircle,1,kBlue,"M^{det} [GeV/c^{2}]","M^{det} / M^{gen}", -1,-1,0,2);
    cres->cd(i+2);
    mres_profy->Draw();
    tslices[i]->Draw("same");
    one->Draw("same");
    
    if (fitting) {
      std::string fit_name = ("res_"+flag+"_"+to_string(i)).c_str();
      
      TF1 *res = new TF1(fit_name.c_str(),"[0]*TMath::ATan([1]*x + [2])+[3]",0,10);
      if (flag == "lo") {
	if      (i == 0) {res->SetParameters(3.06080e-01, 1.46152e+00,0, 7.18850e-01);
			  
			  res->FixParameter(0,1.32872e+01);res->FixParameter(1,8.35957e+00);res->FixParameter(2,2.34800e+01);res->FixParameter(3,-1.95327e+01);}
	else if (i == 1) {res->SetParameters(9.77184e-01, 3.48758e+00,0, -3.30133e-01);
	
	  res->FixParameter(0,9.79988e-01);res->FixParameter(1,1.46052e+00);res->FixParameter(2,1.52811e+00);res->FixParameter(3,-3.28350e-01);}
	else if (i == 2) {res->SetParameters(1.05,10,0.12,-0.5);
	  res->FixParameter(0,1.05);res->FixParameter(2,1.2e-01); res->FixParameter(3,-5e-01);res->FixParameter(1,10);
	  res->FixParameter(0,1.05000e+00);res->FixParameter(1,1.00000e+01);res->FixParameter(2,1.20000e-01);res->FixParameter(3,-5.00000e-01);}
	else if (i == 3) {res->SetParameters(9.77184e-01, 10,-2, 0);
	  res->FixParameter(1,10);res->FixParameter(0,0.75);
	  res->FixParameter(0,7.50000e-01);res->FixParameter(1,1.00000e+01);res->FixParameter(2,-6.33119e+00);res->FixParameter(3,-9.28106e-02);}
	else if (i == 4) {res->SetParameters(9.77184e-01, 10,-2, 0);
	  res->FixParameter(1,10);res->FixParameter(2,-2);
	  res->FixParameter(0,1.11611e+00);res->FixParameter(1,1.00000e+01);res->FixParameter(2,-2.00000e+00);res->FixParameter(3,-6.38040e-01);}
      }
      else if (flag == "mid") {
	if      (i == 0) {res->SetParameters(3.06080e-01, 1.46152e+00,0, 7.18850e-01);
	  
	  res->FixParameter(0,3.12410e-01);res->FixParameter(1,4.04520e-01);res->FixParameter(2,3.31368e-01);res->FixParameter(3,7.45250e-01);}
	else if (i == 1) {res->SetParameters(9.77184e-01, 10,-2, 0);
	  res->FixParameter(1,10);res->FixParameter(2,-2);
	  res->FixParameter(0,8.22198e-01);res->FixParameter(1,1.00000e+01);res->FixParameter(2,-2.00000e+00);res->FixParameter(3,-2.41397e-01);}
	else if (i == 2) {res->SetParameters(9.77184e-01, 10,-2, 0);
	  res->FixParameter(1,10);res->FixParameter(2,-2);
	  res->FixParameter(0,7.02658e-01);res->FixParameter(1,1.00000e+01);res->FixParameter(2,-2.00000e+00);res->FixParameter(3,-3.82504e-02);}
	else if (i == 3) {res->SetParameters(4.29, 3.28,-1, -5.62);
	  res->FixParameter(0,4.29);res->FixParameter(1,3.28);
	  res->FixParameter(0,4.29000e+00);res->FixParameter(1,3.28000e+00);res->FixParameter(2,2.77480e+00);res->FixParameter(3,-5.45843e+00);}
	else if (i == 4) {res->SetParameters(9.77184e-01, 10,-2, 0);
	  res->FixParameter(1,3);res->FixParameter(2,-2);
	  res->FixParameter(0,8.65640e-01);res->FixParameter(1,3.00000e+00);res->FixParameter(2,-2.00000e+00);res->FixParameter(3,-1.71121e-01);}
      }
      else if (flag == "hi") {
	if      (i == 0) {res->SetParameters(9.77184e-01, 10,-2, 0);
	  res->FixParameter(1,10);res->FixParameter(2,-2);
	  res->FixParameter(0,8.35238e-01);res->FixParameter(1,1.00000e+01);res->FixParameter(2,-2.00000e+00);res->FixParameter(3,-2.20800e-01);}
	else if (i == 1) {res->SetParameters(0.2, 1, 0, 0.75);
	  res->FixParameter(0,0.2); res->FixParameter(1,1); res->FixParameter(2,0); res->FixParameter(3,0.75);
	  res->FixParameter(0,2.00000e-01);res->FixParameter(1,1.00000e+00);res->FixParameter(2,0.00000e+00);res->FixParameter(3,7.50000e-01);}
	else if (i == 2) {res->SetParameters(0.4, 1.5, 0, 0.6);
	  res->FixParameter(1,1.5); res->FixParameter(2,0); res->FixParameter(3,0.6);
	  res->FixParameter(0,2.99066e-01);res->FixParameter(1,1.50000e+00);res->FixParameter(2,0.00000e+00);res->FixParameter(3,6.00000e-01);}
	else if (i == 3) {res->SetParameters(0.6, 1.5,0, 0.25);
	  res->FixParameter(0,0.6);res->FixParameter(1,2);res->FixParameter(2,0); res->FixParameter(3,0.15);
	  res->FixParameter(0,6.00000e-01);res->FixParameter(1,2.00000e+00);res->FixParameter(2,0.00000e+00);res->FixParameter(3,1.50000e-01);}
	else if (i == 4) {res->SetParameters(9.77184e-01, 10,-2, 0);
	  res->FixParameter(1,10);res->FixParameter(2,-2);
	  res->FixParameter(0,1.63914e+00);res->FixParameter(1,1.00000e+01);res->FixParameter(2,-2.00000e+00);res->FixParameter(3,-1.48141e+00);}
      }
      mres_profy->Fit(fit_name.c_str(),"MI");
      res->Draw("same");
      
      res->Write();
    }
  }
  
  cres->SaveAs((out+res_m_pt_p6->GetName()+filetype).c_str());
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

  TH3D* res_m_pt_p6_lo = (TH3D*) inFile->Get("mass_res_v_py_m_pt_eta_lo");
  TH3D* res_m_pt_p6_mid = (TH3D*) inFile->Get("mass_res_v_py_m_pt_eta_mid");
  TH3D* res_m_pt_p6_hi = (TH3D*) inFile->Get("mass_res_v_py_m_pt_eta_hi");
  
  TH3D* res_m_pt_p6_lo_corr = (TH3D*) inFile->Get("mass_res_v_py_m_pt_eta_lo_corr");
  TH3D* res_m_pt_p6_mid_corr = (TH3D*) inFile->Get("mass_res_v_py_m_pt_eta_mid_corr");
  TH3D* res_m_pt_p6_hi_corr = (TH3D*) inFile->Get("mass_res_v_py_m_pt_eta_hi_corr");

  p6res(res_m_pt_p6_lo, out, filetype, "lo",funcFile, 1);
  p6res(res_m_pt_p6_mid, out, filetype, "mid",funcFile, 1);
  p6res(res_m_pt_p6_hi, out, filetype, "hi",funcFile, 1);
    
  p6res(res_m_pt_p6_lo_corr, out, filetype, "lo",funcFile, 0);
  p6res(res_m_pt_p6_mid_corr, out, filetype, "mid",funcFile, 0);
  p6res(res_m_pt_p6_hi_corr, out, filetype, "hi",funcFile, 0);
  
  return;
}
