//  Plots.h
//  Created by Isaac Mooney on 7/23/18.
//  A bunch of functions to make plotting more automated and less repetitive

#ifndef Plots_h
#define Plots_h

#include "TROOT.h"
#include <string>
#include <iostream>
#include <vector>

//constructs canvas
//second argument is: 0 = no log, 1 = logx, 2 = logy, 3 = logz (or combinations thereof)
TCanvas * MakeCanvas(const std::string can_name, const std::string log_scale, const double xdim, const double ydim) {
  TCanvas * can = new TCanvas ((can_name).c_str(),(can_name).c_str(),800,800);
    //set desired axes' log scales
    
  if (log_scale.find("1") != std::string::npos || log_scale.find("x") != std::string::npos || log_scale.find("X") != std::string::npos) {can->SetLogx();}
  if (log_scale.find("2") != std::string::npos || log_scale.find("y") != std::string::npos || log_scale.find("Y") != std::string::npos) {can->SetLogy();}
  if (log_scale.find("3") != std::string::npos || log_scale.find("z") != std::string::npos || log_scale.find("Z") != std::string::npos) {can->SetLogz();}
  if (log_scale.find("0") != std::string::npos) {can->SetLogx(0); can->SetLogy(0); can->SetLogz(0);}
    
    return can;
}

//divides existing canvas into multiple (rows*cols) panels
void DivideCanvas(TCanvas *can, const std::string log_scale, const unsigned rows, const unsigned cols) {
    can->Divide(rows, cols, 0, 0);
    for (unsigned i = 1; i <= cols*rows; ++ i) {
        can->cd(i);
        if (log_scale.find("1") != std::string::npos || log_scale.find("x") != std::string::npos || log_scale.find("X") != std::string::npos) {gPad->SetLogx();}
        if (log_scale.find("2") != std::string::npos || log_scale.find("y") != std::string::npos || log_scale.find("Y") != std::string::npos) {gPad->SetLogy();}
        if (log_scale.find("3") != std::string::npos || log_scale.find("z") != std::string::npos || log_scale.find("Z") != std::string::npos) {gPad->SetLogz();}
        if (log_scale.find("0") != std::string::npos) {gPad->SetLogx(0); gPad->SetLogy(0); gPad->SetLogz(0);}
    }
    
    return;
}

//projects a 2D histogram in desired ranges and returns an array of the (1D) projections on desired axis
std::vector<TH1D*> Projection2D (TH2D * hist2D, const int nBins, double * ranges, const std::string axis) {
  std::vector<TH1D*> proj1Ds;
    for (int i = 0; i < nBins; ++ i) {
      std::string low = std::to_string(ranges[i]);
      std::string high = std::to_string(ranges[i+1]);
      std::string low_rough = low.substr(0,2);
      std::string high_rough = high.substr(0,2);
      if (low_rough.substr(1,2) == ".") {low_rough = low_rough.substr(0,1);}
      if (high_rough.substr(1,2) == ".") {high_rough = high_rough.substr(0,1);}
      if (axis == "x" || axis == "X" || axis == "1") {
	proj1Ds.push_back(hist2D->ProjectionX((hist2D->GetName() + axis + low_rough + high_rough).c_str(),ranges[i],ranges[i+1]));
      }
      else if (axis == "y" || axis == "Y" || axis == "2") {
	proj1Ds.push_back(hist2D->ProjectionY((hist2D->GetName() + axis + low_rough + high_rough).c_str(),ranges[i],ranges[i+1]));
      }
      else {
	std::cerr << "Improper axis given for projections. Exiting." << std::endl; exit(1);
      }
      proj1Ds[i]->SetTitle("");
    }
    return proj1Ds;
}

//prettify 1D histogram
void Prettify1D (TH1D * hist, const Color_t markColor, const Style_t markStyle, const double markSize, const Color_t lineColor,
		 const std::string xTitle, const std::string yTitle, const double lowx, const double highx, const double lowy, const double highy) {
  if (yTitle.find("arb") == std::string::npos) { //|| yTitle.find("prob") != std::string::npos) {
    hist->Scale(1/(double)hist->Integral());
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else {
    cout << "Not scaling by histogram " << hist->GetName() << "'s integral & bin width! If this is a cross section measurement, this will be a problem! Fix by passing histogram y-axis title containing anything but 'arb'" << endl;
  }
  hist->SetMarkerColor(markColor); hist->SetMarkerStyle(markStyle); hist->SetMarkerSize(markSize); hist->SetLineColor(lineColor);
  hist->GetXaxis()->SetTitle((xTitle).c_str()); hist->GetYaxis()->SetTitle((yTitle).c_str());
  if (highx != -1) {
    hist->GetXaxis()->SetRangeUser(lowx, highx);
  }
  if (highy != -1) {
    hist->GetYaxis()->SetRangeUser(lowy, highy);
  }
  return;
}

//prettify 1D histogram to be drawn as line with "C" option
void Prettify1DwLineStyle(TH1D * hist, const Color_t lineColor, const Style_t lineStyle, const double lineWidth,
			  const std::string xTitle, const std::string yTitle, const double lowx, const double highx, const double lowy, const double highy) {
  if (yTitle.find("arb") == std::string::npos) {
    hist->Scale(1/(double)hist->Integral());
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else {
    std::cout << "Not scaling by histogram " << hist->GetName() << "'s integral! If this is not a cross section measurement, this will be a problem! Fix by passing histogram y-axis title containing 'prob' or 'arb'" << std::endl;
  }
  hist->SetLineColor(lineColor); hist->SetLineStyle(lineStyle); hist->SetLineWidth(lineWidth);
  hist->GetXaxis()->SetTitle((xTitle).c_str()); hist->GetYaxis()->SetTitle((yTitle).c_str());
  if (highx != -1) {
    hist->GetXaxis()->SetRangeUser(lowx, highx);
  }
  if (highy != -1) {
    hist->GetYaxis()->SetRangeUser(lowy, highy);
  }
  hist->Sumw2(0);
  return;
}

//prettify TGraphErrors
void PrettifyTGraph (TGraphErrors * gr, const Color_t markColor, const Style_t markStyle, const double markSize, const Color_t lineColor,
		     const std::string xTitle, const std::string yTitle, const double lowx, const double highx, const double lowy, const double highy) {
  gr->SetMarkerColor(markColor); gr->SetMarkerStyle(markStyle); gr->SetMarkerSize(markSize); gr->SetLineColor(lineColor);
    gr->SetTitle("");
    gr->GetXaxis()->SetTitle((xTitle).c_str()); gr->GetYaxis()->SetTitle((yTitle).c_str());
    gr->GetXaxis()->SetRangeUser(lowx, highx);
    if (highx != -1) {
        gr->GetXaxis()->SetRangeUser(lowx, highx);
    }
    if (highy != -1) {
        gr->GetYaxis()->SetRangeUser(lowy, highy);
    }
    return;
}

//prettify 2D histogram
void Prettify2D (TH2D * hist, const std::string xTitle, const std::string yTitle,
                 const double lowx, const double highx, const double lowy, const double highy, const double lowz, const double highz) {
  hist->GetXaxis()->SetTitle((xTitle).c_str()); hist->GetYaxis()->SetTitle((yTitle).c_str());
    if (highx != -1) {
        hist->GetXaxis()->SetRangeUser(lowx, highx);
    }
    if (highy != -1) {
        hist->GetYaxis()->SetRangeUser(lowy, highy);
    }
    if (highz != -1) {
        hist->GetZaxis()->SetRangeUser(lowz, highz);
    }
    return;
}

//construct legend with data information
TLegend * TitleLegend(const double xlow, const double ylow, const double xhigh, const double yhigh) {
    TLegend * leg = new TLegend(xlow,ylow,xhigh,yhigh);
    leg->SetBorderSize(0);
    leg->AddEntry((TObject*)0,"pp 200 GeV run12 JP2", "");
    leg->AddEntry((TObject*)0,"anti-k_{T}, R = 0.4", "");
    leg->AddEntry((TObject*)0, "Ch+Ne jets, |#eta| < 0.6", "");
    return leg;
}

TLatex * PanelTitle() {
  TLatex *t = new TLatex();
  t->SetTextAlign(11);
  //  t->SetTextFont(63);
  t->SetTextSizePixels(26);
  t->DrawLatex(0.1,0.9, "pp 200 GeV run12 JP2");
  t->DrawLatex(0.1,0.75, "anti-k_{T}, R = 0.4");
  t->DrawLatex(0.1,0.6, "Ch+Ne jets, |#eta| < 0.6");
  return t;
}

//constructs a slim legend useful for one panel of a projection denoting the range projected over
TLegend * SliceLegend(const std::string title, const double xlow, const double ylow, const double xhigh, const double yhigh) {
    TLegend * leg = new TLegend(xlow,ylow,xhigh,yhigh);
    leg->SetBorderSize(0);
    leg->AddEntry((TObject*)0,(title).c_str(),"");
    return leg;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~COMMON PLOTS: e.g. resolution, observable in bins of pT, ...~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


void Resolution (TH2D* DeltaObsvPt, const std::string out, const std::string filetype, const std::string yTitle, const bool groom) {
  std::string kind = "";
  if (groom == 1) {kind = "SD ";}
  std::string cName = (std::string) "c" + (std::string) DeltaObsvPt->GetName();
  TCanvas * cres = MakeCanvas(cName, "z", 600,800);
  DivideCanvas(cres,"z",1,3);

  const int nDeltas = 11;
  double ranges_delta[nDeltas+1] = {0,1,2,3,4,5,6,7,8,9,10,11};
  vector<TH1D*> delta_projs = Projection2D(DeltaObsvPt, nDeltas, ranges_delta, "Y");
  double means[nDeltas]; double rms[nDeltas]; double meanerr[nDeltas]; double rmserr[nDeltas];
  for (int i = 0; i < nDeltas; ++ i) {
    means[i] = delta_projs[i]->GetMean();
    rms[i] = delta_projs[i]->GetRMS();
    meanerr[i] = delta_projs[i]->GetMeanError();
    rmserr[i] = delta_projs[i]->GetRMSError();
  }
  /*
  double x[nDeltas] = {7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5};
  double y[nDeltas] = {means[0],means[1],means[2],means[3],means[4],means[5],means[6],means[7],means[8],means[9],means[10]};
  double xerr[nDeltas] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  double yerr[nDeltas] = {rms[0],rms[1],rms[2],rms[3],rms[4],rms[5],rms[6],rms[7],rms[8],rms[9],rms[10]};
  */
  double x[nDeltas] = {7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5};
  double ymean[nDeltas] = {means[0],means[1],means[2],means[3],means[4],means[5],means[6],means[7],means[8],means[9],means[10]};
  double xerr[nDeltas] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  double yerrmean[nDeltas] = {meanerr[0],meanerr[1],meanerr[2],meanerr[3],meanerr[4],meanerr[5],meanerr[6],meanerr[7],meanerr[8],meanerr[9],meanerr[10]};

  double yrms[nDeltas] = {rms[0],rms[1],rms[2],rms[3],rms[4],rms[5],rms[6],rms[7],rms[8],rms[9],rms[10]};
  double yerrrms[nDeltas] = {rmserr[0],rmserr[1],rmserr[2],rmserr[3],rmserr[4],rmserr[5],rmserr[6],rmserr[7],rmserr[8],rmserr[9],rmserr[10]};

  TGraphErrors *mean = new TGraphErrors(nDeltas,x,ymean,xerr,yerrmean);
  TGraphErrors *sigma = new TGraphErrors(nDeltas,x,yrms,xerr,yerrrms);
  
  PrettifyTGraph(mean, kRed, kFullCircle, 1, kBlack, "Gen. p^{jet}_{T} [GeV/c]", "#mu", 5, 60, -2, 2);
  PrettifyTGraph(sigma, kRed, kFullCircle, 1, kBlack, "Gen. p^jet}_{T} [GeV/c]","#sigma",5,60,0,0.5);
  
  Prettify2D(DeltaObsvPt, "Gen. p^{jet}_{T} [GeV/c]", yTitle, 5,60,-1,-1,1e-11,1e-3);

  cres->cd(1); DeltaObsvPt->Draw("colz");
  cres->cd(2); mean->Draw("APZ");
  cres->cd(3); sigma->Draw("APZ");
  
  cres->SaveAs((out + DeltaObsvPt->GetName() + filetype).c_str()); 
  
  return;
}


void ObservablePtSlices(TH2D* ObsvPt_d, TH2D* ObsvPt_py, TH2D* ObsvPt_ge, TH2D* ObsvPt_p8, const std::string out, const std::string filetype, std::string xTitle, const bool matched, const bool groom, double* pt_bins) {  
  //NOTE: this function is also used for the matched jet resolution for various observables, so "ObsvPt_ge" does not have the same meaning in this context. The other two 2Ds are unused.
  ObsvPt_d->SetTitle(""); ObsvPt_py->SetTitle(""); ObsvPt_ge->SetTitle(""); ObsvPt_p8->SetTitle("");
  std::string kind = "";
  /*  if (matched == 1) {if (groom == 1) {kind = "gen-SD ";} else {*/kind = "gen-";//}}
  std::string cName = (std::string) "c" + (std::string) ObsvPt_ge->GetName();
  TCanvas * c = MakeCanvas(cName,"0",800,600);
  DivideCanvas(c,"0",3,2);
  
  const int nHists = 5;
  //double gen_pt_bins[nHists+1] = {pt_bins[0] + 2,pt_bins[1] + 2,pt_bins[2] + 2,pt_bins[3] + 2, pt_bins[4] + 2,pt_bins[5] + 2}; //pT is usually 11 5 GeV bins from 5 to 60 GeV, so e.g. bin 0 = 5, bin 1 = 10, etc. 
  double corresp_pts[nHists+1] ={15,20,25,30,40,60};
  double gen_pts[nHists+1] = {2,3,4,5,7,11};
  
  vector<TH1D*> py_projs = Projection2D (ObsvPt_py, nHists, gen_pts, "X");
  //  vector<TH1D*> ge_projs = Projection2D (ObsvPt_ge, nHists, pt_bins, "X");
  //  vector<TH1D*> d_projs = Projection2D (ObsvPt_d, nHists, pt_bins, "X");
  vector<TH1D*> p8_projs = Projection2D (ObsvPt_p8, nHists, gen_pts, "X");
 
  double lox, hix, loy, hiy;
  
  if ((xTitle.find("z_") == std::string::npos && xTitle.find("r_") == std::string::npos && xTitle.find("Z_") == std::string::npos && xTitle.find("R_") == std::string::npos) || xTitle.find(" / ") != std::string::npos) {lox = -1; hix = -1;} else {lox = 0; hix = 0.5;}
  
  if (xTitle.find("M_{c") != std::string::npos) {lox = -0.5; hix = 5;}

  if (xTitle.find(" / ") != std::string::npos) { hiy = 0.18; } else { hiy = 0.5; } 
  
  if (xTitle.find("p_{T,g} / ") != std::string::npos) {hiy = 20;}
  
  std::string yTitle;

  if (matched == 1) {hiy /= (double) 0.039; /*std::string sub = xTitle.substr(xTitle.find(' ')); yTitle = ("1/N_{pair} dN_{pair}/d" + sub.substr(1,sub.length())).c_str();*/ yTitle = "1/N_{pair} dN_{pair}/d(ratio)";}
  else {
    yTitle = ("1/N_{pair} dN_{pair}/d" + xTitle.substr(xTitle.find(' '))).c_str();
  }

  if (xTitle.find("z_") != std::string::npos || xTitle.find("R_") != std::string::npos) {xTitle = (xTitle.substr(0,5)).c_str(); hiy = 7;}
  
  for(int i = 0; i < nHists; ++ i) {
    Prettify1DwLineStyle(py_projs[i], kGreen, kSolid, 5, xTitle, yTitle,lox, hix, loy, hiy);
    //Prettify1D(ge_projs[i], kBlue, kOpenCircle, 1, kBlue, xTitle, yTitle,lox, hix, 0, hiy);
    //Prettify1D(d_projs[i], kBlack, kFullStar, 2, kBlack, xTitle, yTitle,lox, hix, 0, hiy);
    Prettify1DwLineStyle(p8_projs[i], kGreen, kDashed, 5, xTitle, yTitle, lox, hix, loy, hiy);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TLEGEND~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                                   
  //TLegend *t = TitleLegend(0.1,0.15,0.8,0.8);
  //TLatex *t = PanelTitle();
  TLegend *pan_leg = new TLegend(0.1,0.15,0.8,0.45); pan_leg->SetBorderSize(0);
  if (matched == 0) {
    pan_leg->AddEntry(py_projs[0], "PYTHIA6", "l");
    //pan_leg->AddEntry(ge_projs[0], "PYTHIA6+GEANT", "p");
    //pan_leg->AddEntry(d_projs[0], "Raw data", "p");
    pan_leg->AddEntry(p8_projs[0], "PYTHIA8", "l");
  }

  TLegend *tslices[nHists];
  for (int i = 0; i < nHists; ++ i) {
    tslices[i] = SliceLegend(((to_string(corresp_pts[i])).substr(0,2) + " < p_{T}^{" + kind + "jet} < " + (to_string(corresp_pts[i + 1])).substr(0,2) + " GeV/c").c_str(), 0.13,0.8,0.9,0.95);
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                                 
  
  c->cd(1); TLatex *t = PanelTitle(); if(matched==0) { pan_leg->Draw(); }
  for (int i = 0; i < nHists; ++ i) {
    c->cd(i+2);
    if (matched == 0) {
      py_projs[i]->Draw("C"); p8_projs[i]->Draw("C,same"); //d_projs[i]->Draw("same");
    }
    /*ge_projs[i]->Draw("same");*/ tslices[i]->Draw("same");
  }
  cout << "vi" << endl;
  //  c->SaveAs((out + ObsvPt_ge->GetName()  + "_slices" + filetype).c_str());
  return;
}

void Response(TFile *matchFile, const std::string resName, TH1D* ObsPy, TH1D* ObsGe, const std::string xTitle, const std::string out, const std::string filetype, const std::string flag) {
  ObsPy->SetTitle(""); ObsGe->SetTitle("");
  std::string cresName = (std::string) "cres" + (std::string) resName;
  std::string cmesName = (std::string) "cmes" + (std::string) ObsGe->GetName();
  TCanvas * cres = MakeCanvas(cresName.c_str(),"z",800,800);
  string logyn;
  if (xTitle.find("p_{T") != std::string::npos) { logyn = "y";} else{ logyn = "0";}
  TCanvas * cmes = MakeCanvas(cmesName.c_str(),logyn,800,800);
  
  double lox = -2; double hix = -2; double loy = -2; double hiy = -2;
  if (flag != "pt") {loy = 0; hiy = 7;}
  if (flag == "mass") {lox = -1; hix = -1; loy = 0; hiy = 2;} else if (flag == "pt") {loy = -1; hiy = -1;} else {lox = 0; hix = 0.5;}
  
  //if (xTitle.find("z_") == std::string::npos && xTitle.find("r_") == std::string::npos && xTitle.find("Z_") == std::string::npos && xTitle.find("R_") == std::string::npos) { lox = -1; hix = -1; } else {lox = 0; hix = 0.5;}

  std::string yTitle = "";
  if (flag == "pt") {
    yTitle = "arb.";
    ObsGe->Scale(1/(double)ObsGe->Integral());
    double rat_scale = ObsGe->GetBinContent(1) / (double) ObsPy->GetBinContent(1);
    ObsPy->Scale(rat_scale);
  }
  else {yTitle = ("1/N_{j} dN_{j}/d" + xTitle.substr(0,xTitle.find(' '))).c_str();}

  Prettify1D(ObsPy, kGreen, kOpenCircle, 2, kBlack, xTitle, yTitle, lox, hix,loy, hiy);
  Prettify1D(ObsGe, kBlue, kFullCircle, 2, kBlack, xTitle, yTitle, lox, hix,loy, hiy);
  RooUnfoldResponse *res = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  TH2D* hres = (TH2D*) res->Hresponse(); hres->SetTitle("");

  double lor = -1; double hir = -1;
  if (flag == "pt") {lor = 15; hir = 60;}
  
  Prettify2D(hres, ("Det. " + xTitle).c_str(), ("Gen. " + xTitle).c_str(), lor,hir,-1,-1,1e-11,1e-4);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TLEGEND~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                  
  TLegend *tmes;
  if (xTitle.find("R_") != std::string::npos || xTitle.find("r_") != std::string::npos) {
    tmes = TitleLegend(0.15, 0.57, 0.45, 0.88);
  }
  else { tmes = TitleLegend(0.5,0.57,0.8,0.88); }
  tmes->AddEntry(ObsPy, "PYTHIA6", "p");
  tmes->AddEntry(ObsGe, "PYTHIA6+GEANT", "p");
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                             
  cres->cd(); hres->Draw("colz");
  cmes->cd(); ObsPy->Draw(); ObsGe->Draw("same"); tmes->Draw("same");

  cres->SaveAs((out + resName + filetype).c_str());
  cmes->SaveAs((out + ObsGe->GetName() + "_and_py" + filetype).c_str());
  return;
}

#endif /* Plots_h */
