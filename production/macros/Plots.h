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
  TCanvas * can = new TCanvas ((can_name).c_str(),(can_name).c_str(),xdim,ydim);
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

//projects a 3D histogram in desired ranges and returns an array of the (2D) projections on desired axis
//For 3D plots, have to use SetRange instead of SetRangeUser for some stupid reason. => pass bin numbers instead of bin values to the function 
std::vector<TH2D*> Projection3D (TH3D * hist3D, const int nBins, int * ranges, const std::string axis) {
  std::vector<TH2D*> proj2Ds;
  for (int i = 0; i < nBins; ++ i) {
    std::string low = std::to_string(ranges[i]);
    std::string high = std::to_string(ranges[i+1]);
    std::string low_rough = low.substr(0,2);
    std::string high_rough = high.substr(0,2);
    if (low_rough.substr(1,2) == ".") {low_rough = low_rough.substr(0,1);}
    if (high_rough.substr(1,2) == ".") {high_rough = high_rough.substr(0,1);}
    if (axis == "xy") {
      hist3D->GetZaxis()->SetRange(ranges[i],ranges[i+1]);
      proj2Ds.push_back((TH2D*) hist3D->Project3D((hist3D->GetName() + axis +"e"+ low_rough + high_rough).c_str()));
    }
    else if (axis == "yx") {
      hist3D->GetZaxis()->SetRange(ranges[i],ranges[i+1]);
      proj2Ds.push_back((TH2D*) hist3D->Project3D((hist3D->GetName() + axis +"e"+ low_rough + high_rough).c_str()));
    }
    else if (axis == "yz") {
      hist3D->GetXaxis()->SetRange(ranges[i],ranges[i+1]);
      proj2Ds.push_back((TH2D*) hist3D->Project3D((hist3D->GetName() + axis +"e"+ low_rough + high_rough).c_str()));
    }
    else if (axis == "zy") {
      hist3D->GetXaxis()->SetRange(ranges[i],ranges[i+1]);
      proj2Ds.push_back((TH2D*) hist3D->Project3D((hist3D->GetName() + axis +"e"+ low_rough + high_rough).c_str()));
    }
    else if (axis == "xz") {
      hist3D->GetYaxis()->SetRange(ranges[i],ranges[i+1]);
      proj2Ds.push_back((TH2D*) hist3D->Project3D((hist3D->GetName() + axis +"e"+ low_rough + high_rough).c_str()));
    }
    else if (axis == "zx") {
      hist3D->GetYaxis()->SetRange(ranges[i],ranges[i+1]);
      proj2Ds.push_back((TH2D*) hist3D->Project3D((hist3D->GetName() + axis +"e"+ low_rough + high_rough).c_str()));
    }
    else {
      std::cerr << "Improper axis given for projections. Exiting." << std::endl; exit(1);
    }
    proj2Ds[i]->SetTitle("");
  }
  return proj2Ds;
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
  
  if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") == std::string::npos) { //|| yTitle.find("prob") != std::string::npos) {
    hist->Scale(1/(double)hist->Integral());
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") != std::string::npos) {
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
  
  if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") == std::string::npos) { //|| yTitle.find("prob") != std::string::npos) {
    hist->Scale(1/(double)hist->Integral());
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") != std::string::npos) {
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else {
    cout << "Not scaling by histogram " << hist->GetName() << "'s integral & bin width! If this is a cross section measurement, this will be a problem! Fix by passing histogram y-axis title containing anything but 'arb'" << endl;
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

//prettify 1D histogram to be drawn as line with "C" option //OVERLOADED
void Prettify1DwLineStyle(TH1D * hist, const Style_t lineStyle, const double lineWidth,
			  const std::string xTitle, const std::string yTitle, const double lowx, const double highx, const double lowy, const double highy) {
  
  if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") == std::string::npos) { //|| yTitle.find("prob") != std::string::npos) {
    hist->Scale(1/(double)hist->Integral());
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") != std::string::npos) {
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else {
    cout << "Not scaling by histogram " << hist->GetName() << "'s integral & bin width! If this is a cross section measurement, this will be a problem! Fix by passing histogram y-axis title containing anything but 'arb'" << endl;
  }
  
  hist->SetLineStyle(lineStyle); hist->SetLineWidth(lineWidth);
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
    leg->AddEntry((TObject*)0,"pp 200 GeV" /*#pi_{0} trigger > 5.4 GeV/c*/,"");
    leg->AddEntry((TObject*)0,"anti-k_{T}, R = 0.4"/*2"*/,"");
    leg->AddEntry((TObject*)0,"Ch+Ne jets, |#eta| < 0.6"/*"Ch rec. jets, |#eta| < 0.8, p^{cons.}_{T} < 30 GeV/c"*/,"");
    leg->AddEntry((TObject*)0,"Constituents assigned m_{PDG}"/*"#phi_{rec.} #in [3#pi/4, 5#pi/4] w/r/t #phi_{trig.}"*/,"");
    
    return leg;
}

TLatex * PanelTitle() {
  TLatex *t = new TLatex();
  t->SetTextAlign(11);
  //  t->SetTextFont(63);
  t->SetTextSize(0.07);
  t->DrawLatexNDC(0.2,0.65, "Pythia8 - pp 200 GeV"/*#pi_{0} trigger > 5.4 GeV/c"*/);
  t->DrawLatexNDC(0.2,0.55, "anti-k_{T}, R = 0.4"/*2"*/);
  t->DrawLatexNDC(0.2,0.45, "Ch+Ne jets, |#eta| < 0.6"/*"Ch rec. jets, |#eta| < 0.8"*/);
  t->DrawLatexNDC(0.2,0.35,"Constituents assigned m_{PDG}"/*"#phi_{rec.} within 1/4 of #phi_{trig.} + #pi"*/);

  return t;
}

//constructs a slim legend useful for one panel of a projection denoting the range projected over
TLegend * SliceLegend(const std::string title, const double xlow, const double ylow, const double xhigh, const double yhigh) {
    TLegend * leg = new TLegend(xlow,ylow,xhigh,yhigh);
    leg->SetBorderSize(0);
    leg->AddEntry((TObject*)0,(title).c_str(),"");
    return leg;
}

#endif /* Plots_h */
