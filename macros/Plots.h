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

//make variable bin-size histogram from tree
TH1D* HistFromTree(const std::string hist_name, const int nBins, double * edges, TTree * tree, double & branch, double & weight) {
  TH1D * hist = new TH1D((hist_name).c_str(), "", nBins, edges);
    
    for (int i = 0; i < tree->GetEntries(); ++ i) {
        tree->GetEntry(i);
        hist->Fill(branch, weight);
    }
    return hist;
}

//make variable bin-size histogram from tree
TH1D* HistFromTree(const std::string hist_name, const int nBins, const double lo, const double hi, TTree * tree, double & branch, double & weight) {
  TH1D * hist = new TH1D((hist_name).c_str(), "", nBins, lo, hi);
    
    for (int i = 0; i < tree->GetEntries(); ++ i) {
        tree->GetEntry(i);
        hist->Fill(branch, weight);
    }
    return hist;
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
  if (yTitle.find("arb") != std::string::npos || yTitle.find("prob") != std::string::npos) {
    hist->Scale(1/(double)hist->Integral());
  }
  else {
    cout << "Not scaling by histogram " << hist->GetName() << "'s integral! If this is not a cross section measurement, this will be a problem! Fix by passing histogram y-axis title containing 'prob' or 'arb'" << endl;
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
  if (yTitle.find("arb") != std::string::npos || yTitle.find("prob") != std::string::npos) {
    hist->Scale(1/(double)hist->Integral());
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
  DivideCanvas(cres,"z",1,2);

  const int nDeltas = 11;
  double ranges_delta[nDeltas+1] = {0,1,2,3,4,5,6,7,8,9,10,11};
  vector<TH1D*> delta_projs = Projection2D(DeltaObsvPt, nDeltas, ranges_delta, "Y");
  double means[nDeltas]; double rms[nDeltas];
  for (int i = 0; i < nDeltas; ++ i) {
    means[i] = delta_projs[i]->GetMean();
    rms[i] = delta_projs[i]->GetRMS();
  }

  double x[nDeltas] = {7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5};
  double y[nDeltas] = {means[0],means[1],means[2],means[3],means[4],means[5],means[6],means[7],means[8],means[9],means[10]};
  double xerr[nDeltas] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  double yerr[nDeltas] = {rms[0],rms[1],rms[2],rms[3],rms[4],rms[5],rms[6],rms[7],rms[8],rms[9],rms[10]};

  TGraphErrors *gr = new TGraphErrors(nDeltas,x,y,xerr,yerr);

  PrettifyTGraph(gr, kRed, kFullCircle, 1, kBlack, ("Gen. p^{" + kind + "jet}_{T} [GeV/c]").c_str(), "#mu #pm #sigma", 5, 60, -2, 2);

  Prettify2D(DeltaObsvPt, ("Gen. p^{" + kind + "jet}_{T} [GeV/c]").c_str(), yTitle, 5,60,-1,-1,1e-11,1e-3);

  cres->cd(1); DeltaObsvPt->Draw("colz");
  cres->cd(2); gr->Draw("APZ");

  cres->SaveAs((out + DeltaObsvPt->GetName() + filetype).c_str()); 
  
  return;
}


void ObservablePtSlices(TH2D* ObsvPt_d, TH2D* ObsvPt_py, TH2D* ObsvPt_ge, const std::string out, const std::string filetype, const std::string xTitle, const bool matched, const bool groom) {

  //NOTE: this function is also used for the matched jet resolution for various observables, so "ObsvPt_ge" does not have the same meaning in this context. The other two 2Ds are unused.
  ObsvPt_d->SetTitle(""); ObsvPt_py->SetTitle(""); ObsvPt_ge->SetTitle("");
  std::string kind = "";
  if (matched == 1) {if (groom == 1) {kind = "gen-SD ";} else {kind = "gen-";}}
  std::string cName = (std::string) "c" + (std::string) ObsvPt_ge->GetName();
  TCanvas * c = MakeCanvas(cName,"0",800,600);
  DivideCanvas(c,"0",4,2);

  const int nHists = 7;
  double pt_bins[nHists+1] = {0,1,2,3,4,5,7,11}; //pT is usually 11 5 GeV bins from 5 to 60 GeV, so e.g. bin 0 = 5, bin 1 = 10, etc. 
  double corresp_pts[nHists+1] ={5,10,15,20,25,30,40,60};
  vector<TH1D*> py_projs = Projection2D (ObsvPt_py, nHists, pt_bins, "X");
  vector<TH1D*> ge_projs = Projection2D (ObsvPt_ge, nHists, pt_bins, "X");
  vector<TH1D*> d_projs = Projection2D (ObsvPt_d, nHists, pt_bins, "X");

  for(int i = 0; i < nHists; ++ i) {
    Prettify1DwLineStyle(py_projs[i], kGreen, kDashed, 5, xTitle, "prob.",-1,-1, 0, 0.5);
    Prettify1D(ge_projs[i], kBlue, kOpenCircle, 1, kBlue, xTitle, "prob.",-1,-1, 0, 0.5);
    Prettify1D(d_projs[i], kBlack, kFullStar, 2, kBlack, xTitle, "prob.",-1,-1, 0, 0.5);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TLEGEND~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                                   
  //TLegend *t = TitleLegend(0.1,0.15,0.8,0.8);
  //TLatex *t = PanelTitle();
  TLegend *pan_leg = new TLegend(0.1,0.15,0.8,0.45); pan_leg->SetBorderSize(0);
  if (matched == 0) {
    pan_leg->AddEntry(py_projs[0], "PYTHIA6", "l");
    pan_leg->AddEntry(ge_projs[0], "PYTHIA6+GEANT", "p");
    pan_leg->AddEntry(d_projs[0], "Raw data", "p");
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
      py_projs[i]->Draw("C"); d_projs[i]->Draw("same");
    }
    ge_projs[i]->Draw("same"); tslices[i]->Draw("same");
  }

  c->SaveAs((out + ObsvPt_ge->GetName()  + "_slices" + filetype).c_str());
  return;
}

void Response(TFile *matchFile, const std::string resName, TH1D* ObsPy, TH1D* ObsGe, const std::string xTitle, const std::string out, const std::string filetype) {
  ObsPy->SetTitle(""); ObsGe->SetTitle("");
  std::string cresName = (std::string) "cres" + (std::string) resName;
  std::string cmesName = (std::string) "cmes" + (std::string) ObsGe->GetName();
  TCanvas * cres = MakeCanvas(cresName.c_str(),"z",800,800);
  TCanvas * cmes = MakeCanvas(cmesName.c_str(),"0",800,800);

  Prettify1D(ObsPy, kGreen, kOpenCircle, 2, kBlack, xTitle, "prob.", -1,-1,-1,-1);
  Prettify1D(ObsGe, kBlue, kFullCircle, 2, kBlack, xTitle, "prob.", -1,-1,-1,-1);

  RooUnfoldResponse *res = (RooUnfoldResponse*) matchFile->Get(resName.c_str());
  TH2D* hres = (TH2D*) res->Hresponse(); hres->SetTitle("");

  Prettify2D(hres, ("Det. " + xTitle).c_str(), ("Gen. " + xTitle).c_str(), -1,-1,-1,-1,1e-11,1e-3);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TLEGEND~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                  
  TLegend *tmes = TitleLegend(0.5,0.57,0.8,0.88);
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