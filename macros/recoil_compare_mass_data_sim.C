#include <string>
#include <iostream>

//COMPARING INCLUSIVE JET MASS FOR PYTHIA, PYTHIA+GEANT, & DATA                                                                                                                  

using namespace std;

void recoil_compare_mass_data_sim () {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HISTOGRAMS AND CANVASES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  string sim_dir = "~/jetmass/";
  string data_dir = "~/jetmass/";
  string sim_in = "out/sim/";
  string data_in = "out/data/";
  string sim_file = "full.root";
  string data_file = "full.root";
  string out = "~/jetmass/plots/massplots/";
  string filetype = ".pdf";
    
  TCanvas *cfull = new TCanvas("cfull","cfull",800,800);    //1D full range
  TCanvas *cpympt = new TCanvas("cpympt","cpympt",800,800); //2D pythia
  TCanvas *cgempt = new TCanvas("cgempt","cgempt",800,800); //2D geant
  TCanvas *cdmpt = new TCanvas("cdmpt","cdmpt",800,800);    //2D data
  TCanvas *cslices = new TCanvas("cslices","cslices",1200,800); //canvas to be partitioned for pt slices
  TCanvas *cratios = new TCanvas("cratios","cratios",1200,800); //to be partitioned for ratios of the slices
  TCanvas *cmean = new TCanvas("cmean","cmean",800,800);    //means for each pt slice
  TCanvas *crms = new TCanvas("crms","crms",800,800);       //rms for each pt slice
  
  cfull->SetLogy(); cpympt->SetLogz(); cgempt->SetLogz(); cdmpt->SetLogz();
  
  cslices->Divide(3,2,0,0); cratios->Divide(3,2,0,0);
  
  TFile* simFile = new TFile( (sim_dir + sim_in + sim_file).c_str(), "READ");
    //MASS v PT
  TH2D *py_m_pt = (TH2D*) simFile->Get("m_v_pt_full_trig_jet_py");
  TH2D *ge_m_pt = (TH2D*) simFile->Get("m_v_pt_full_rec_jet_ge");
  
  TFile* dataFile = new TFile( (data_dir + data_in + data_file).c_str(), "READ");
    //MASS v PT
  TH2D *d_m_pt = (TH2D*) dataFile->Get("m_v_pt_full_rec_jet");
  
  string pyfull = py_m_pt->GetName(); string gefull = ge_m_pt->GetName(); string dfull = d_m_pt->GetName();
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MASS FOR PT SLICES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  //bins go from e.g. 9 to 14 but correspond to values 10 to 15 for pT.
  TH1D* py_10_15 = py_m_pt->ProjectionX((pyfull + "1015").c_str(),9,14);
  TH1D* py_15_20 = py_m_pt->ProjectionX((pyfull + "1520").c_str(),14,19);
  TH1D* py_20_25 = py_m_pt->ProjectionX((pyfull + "2025").c_str(),19,24);
  TH1D* py_25_30 = py_m_pt->ProjectionX((pyfull + "2530").c_str(),24,29);
  TH1D* py_30_plus = py_m_pt->ProjectionX((pyfull + "30plus").c_str(),29,1000);
    TH1D* py_full = py_m_pt->ProjectionX((pyfull + "all").c_str(),0,1000); //MASS FOR FULL PT RANGE
    
  TH1D* ge_10_15 = ge_m_pt->ProjectionX((gefull + "1015").c_str(),9,14);
  TH1D* ge_15_20 = ge_m_pt->ProjectionX((gefull + "1520").c_str(),14,19);
  TH1D* ge_20_25 = ge_m_pt->ProjectionX((gefull + "2025").c_str(),19,24);
  TH1D* ge_25_30 = ge_m_pt->ProjectionX((gefull + "2530").c_str(),24,29);
  TH1D* ge_30_plus = ge_m_pt->ProjectionX((gefull + "30plus").c_str(),29,1000);
    TH1D* ge_full = ge_m_pt->ProjectionX((gefull + "all").c_str(),0,1000); //MASS FOR FULL PT RANGE
    
  TH1D* d_10_15 = d_m_pt->ProjectionX((dfull + "1015").c_str(),9,14);
  TH1D* d_15_20 = d_m_pt->ProjectionX((dfull + "1520").c_str(),14,19);
  TH1D* d_20_25 = d_m_pt->ProjectionX((dfull + "2025").c_str(),19,24);
  TH1D* d_25_30 = d_m_pt->ProjectionX((dfull + "2530").c_str(),24,29);
  TH1D* d_30_plus = d_m_pt->ProjectionX((dfull + "30plus").c_str(),29,1000);
    TH1D* d_full = d_m_pt->ProjectionX((dfull + "all").c_str(),0,1000); //MASS FOR FULL PT RANGE

    ge_10_15->Scale(1/(double)ge_10_15->Integral()); d_10_15->Scale(1/(double)d_10_15->Integral());
    ge_15_20->Scale(1/(double)ge_15_20->Integral()); d_15_20->Scale(1/(double)d_15_20->Integral());
    ge_20_25->Scale(1/(double)ge_20_25->Integral()); d_20_25->Scale(1/(double)d_20_25->Integral());
    ge_25_30->Scale(1/(double)ge_25_30->Integral()); d_25_30->Scale(1/(double)d_25_30->Integral());
    ge_30_plus->Scale(1/(double)ge_30_plus->Integral()); d_30_plus->Scale(1/(double)d_30_plus->Integral());
    
    TH1D* ratio_10_15 = (TH1D*)d_10_15->Clone("ratio"); ratio_10_15->Divide(ge_10_15);
    TH1D* ratio_15_20 = (TH1D*)d_15_20->Clone("ratio"); ratio_15_20->Divide(ge_15_20);
    TH1D* ratio_20_25 = (TH1D*)d_20_25->Clone("ratio"); ratio_20_25->Divide(ge_20_25);
    TH1D* ratio_25_30 = (TH1D*)d_25_30->Clone("ratio"); ratio_25_30->Divide(ge_25_30);
    TH1D* ratio_30_plus = (TH1D*)d_30_plus->Clone("ratio"); ratio_30_plus->Divide(ge_30_plus);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MASS MEAN & RMS FOR PT SLICES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    double py_mean_10_15 = py_10_15->GetMean();  double ge_mean_10_15 = ge_10_15->GetMean(); double d_mean_10_15 = d_10_15->GetMean();
    double py_mean_15_20 = py_15_20->GetMean();  double ge_mean_15_20 = ge_15_20->GetMean(); double d_mean_15_20 = d_15_20->GetMean();
    double py_mean_20_25 = py_20_25->GetMean();  double ge_mean_20_25 = ge_20_25->GetMean(); double d_mean_20_25 = d_20_25->GetMean();
    double py_mean_25_30 = py_25_30->GetMean();  double ge_mean_25_30 = ge_25_30->GetMean(); double d_mean_25_30 = d_25_30->GetMean();
    double py_mean_30_plus = py_30_plus->GetMean();  double ge_mean_30_plus = ge_30_plus->GetMean(); double d_mean_30_plus = d_30_plus->GetMean();
    
    double py_rms_10_15 = py_10_15->GetRMS();  double ge_rms_10_15 = ge_10_15->GetRMS(); double d_rms_10_15 = d_10_15->GetRMS();
    double py_rms_15_20 = py_15_20->GetRMS();  double ge_rms_15_20 = ge_15_20->GetRMS(); double d_rms_15_20 = d_15_20->GetRMS();
    double py_rms_20_25 = py_20_25->GetRMS();  double ge_rms_20_25 = ge_20_25->GetRMS(); double d_rms_20_25 = d_20_25->GetRMS();
    double py_rms_25_30 = py_25_30->GetRMS();  double ge_rms_25_30 = ge_25_30->GetRMS(); double d_rms_25_30 = d_25_30->GetRMS();
    double py_rms_30_plus = py_30_plus->GetRMS();  double ge_rms_30_plus = ge_30_plus->GetRMS(); double d_rms_30_plus = d_30_plus->GetRMS();
    
    Double_t pt_array[5] = {12.5,17.5,22.5,27.5,32.5};//FIX LATER! THE LAST POINT AT 32.5 SHOULD BE ADJUSTED SOMEHOW!!! HINT: what's the mean pt in the pt hist zoomed in on 30+?
    Double_t py_m_mean[5] = {py_mean_10_15, py_mean_15_20, py_mean_20_25, py_mean_25_30, py_mean_30_plus};
    Double_t ge_m_mean[5] = {ge_mean_10_15, ge_mean_15_20, ge_mean_20_25, ge_mean_25_30, ge_mean_30_plus};
    Double_t d_m_mean[5] = {d_mean_10_15, d_mean_15_20, d_mean_20_25, d_mean_25_30, d_mean_30_plus};
    Double_t pt_errs[5] = {2.5,2.5,2.5,2.5,2.5}; //LAST POINT SHOULD ALSO BE ADJUSTED HERE ACCORDINGLY
    Double_t py_m_mean_errs[5] = {py_10_15->GetMeanError(), py_15_20->GetMeanError(), py_20_25->GetMeanError(), py_25_30->GetMeanError(), py_30_plus->GetMeanError()};
    Double_t ge_m_mean_errs[5] = {ge_10_15->GetMeanError(), ge_15_20->GetMeanError(), ge_20_25->GetMeanError(), ge_25_30->GetMeanError(), ge_30_plus->GetMeanError()};
    Double_t d_m_mean_errs[5] = {d_10_15->GetMeanError(), d_15_20->GetMeanError(), d_20_25->GetMeanError(), d_25_30->GetMeanError(), d_30_plus->GetMeanError()};
    Double_t py_m_rms_errs[5] = {py_10_15->GetRMSError(), py_15_20->GetRMSError(), py_20_25->GetRMSError(), py_25_30->GetRMSError(), py_30_plus->GetRMSError()};
    Double_t ge_m_rms_errs[5] = {ge_10_15->GetRMSError(), ge_15_20->GetRMSError(), ge_20_25->GetRMSError(), ge_25_30->GetRMSError(), ge_30_plus->GetRMSError()};
    Double_t d_m_rms_errs[5] = {d_10_15->GetRMSError(), d_15_20->GetRMSError(), d_20_25->GetRMSError(), d_25_30->GetRMSError(), d_30_plus->GetRMSError()};
    Double_t py_m_rms[5] = {py_rms_10_15, py_rms_15_20, py_rms_20_25, py_rms_25_30, py_rms_30_plus};
    Double_t ge_m_rms[5] = {ge_rms_10_15, ge_rms_15_20, ge_rms_20_25, ge_rms_25_30, ge_rms_30_plus};
    Double_t d_m_rms[5] = {d_rms_10_15, d_rms_15_20, d_rms_20_25, d_rms_25_30, d_rms_30_plus};
    
    py_mean_m_v_pt_gr = new TGraphErrors(5,pt_array,py_m_mean,pt_errs,py_m_mean_errs);
    ge_mean_m_v_pt_gr = new TGraphErrors(5,pt_array,ge_m_mean,pt_errs,ge_m_mean_errs);
    d_mean_m_v_pt_gr = new TGraphErrors(5,pt_array,d_m_mean,pt_errs,d_m_mean_errs);
    py_rms_m_v_pt_gr = new TGraphErrors(5,pt_array,py_m_rms,pt_errs,py_m_rms_errs);
    ge_rms_m_v_pt_gr = new TGraphErrors(5,pt_array,ge_m_rms,pt_errs,ge_m_rms_errs);
    d_rms_m_v_pt_gr = new TGraphErrors(5,pt_array,d_m_rms,pt_errs,d_m_rms_errs);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TLEGENDS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    TLegend *tfull = new TLegend(0.52,0.7,0.84,0.8); tfull->SetBorderSize(0);
    tfull->AddEntry(py_full, "PYTHIA6", "p");
    tfull->AddEntry(ge_full, "PYTHIA6+GEANT", "p");
    tfull->AddEntry(d_full, "STAR pp run6 - 200 GeV", "p");
    
    TLegend *tmean = new TLegend(0.2,0.6,0.5,0.8); tmean->SetBorderSize(0);
    tmean->AddEntry(py_mean_m_v_pt_gr, "PYTHIA6", "p");
    tmean->AddEntry(ge_mean_m_v_pt_gr, "PYTHIA6+GEANT", "p");
    tmean->AddEntry(d_mean_m_v_pt_gr, "STAR pp run6 - 200 GeV", "p");
    
    TLegend *trms = new TLegend(0.2,0.6,0.5,0.8); trms->SetBorderSize(0);
    trms->AddEntry(py_rms_m_v_pt_gr, "PYTHIA6", "p");
    trms->AddEntry(ge_rms_m_v_pt_gr, "PYTHIA6+GEANT", "p");
    trms->AddEntry(d_rms_m_v_pt_gr, "STAR pp run6 - 200 GeV", "p");
    
    TLegend *tslices = new TLegend(0.2,0.2,0.8,0.7); tslices->SetBorderSize(0);
    tslices->AddEntry(py_10_15, "PYTHIA6", "p");
    tslices->AddEntry(ge_10_15, "PYTHIA6+GEANT", "p");
    tslices->AddEntry(d_10_15, "STAR pp run6 - 200 GeV", "p");
    tslices->AddEntry((TObject*)0, "Ch+Ne trigger jets, |#eta| < 0.6","");
    tslices->AddEntry((TObject*)0, "anti-kt, R = 0.4", "");
    
    TLegend *t1015 = new TLegend(0.49,0.75,0.9,0.9); t1015->SetBorderSize(0);
    t1015->AddEntry((TObject*)0, "10 < p_{T} < 15 GeV", "");
    TLegend *t1520 = new TLegend(0.49,0.75,0.9,0.9); t1520->SetBorderSize(0);
    t1520->AddEntry((TObject*)0, "15 < p_{T} < 20 GeV", "");
    TLegend *t2025 = new TLegend(0.49,0.75,0.9,0.9); t2025->SetBorderSize(0);
    t2025->AddEntry((TObject*)0, "20 < p_{T} < 25 GeV", "");
    TLegend *t2530 = new TLegend(0.49,0.75,0.9,0.9); t2530->SetBorderSize(0);
    t2530->AddEntry((TObject*)0, "25 < p_{T} < 30 GeV", "");
    TLegend *t30plus = new TLegend(0.6,0.76,0.9,0.9); t30plus->SetBorderSize(0);
    t30plus->AddEntry((TObject*)0, "p_{T} > 30 GeV", "");
    
    TLegend *tratios = new TLegend(0.005,0.35,0.92,0.66); tratios->SetBorderSize(0);
    tratios->AddEntry((TObject*)0, "Ch+Ne trigger jets, |#eta| < 0.6","");
    tratios->AddEntry((TObject*)0, "anti-kt, R = 0.4", "");
    tratios->AddEntry((TObject*)0, "STAR pp run6 200 GeV / PYTHIA6+GEANT", "");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRETTIFYING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //full mass
        //axes
    py_full->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); py_full->GetYaxis()->SetTitle("prob.");
    ge_full->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ge_full->GetYaxis()->SetTitle("prob.");
    d_full->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); d_full->GetYaxis()->SetTitle("prob.");
    py_full->GetYaxis()->SetTitleOffset(1.3); ge_full->GetYaxis()->SetTitleOffset(1.3); d_full->GetYaxis()->SetTitleOffset(1.3);
    py_full->GetXaxis()->SetRangeUser(0,15); ge_full->GetXaxis()->SetRangeUser(0,15); d_full->GetXaxis()->SetRangeUser(0,15);
        //markers
    py_full->SetMarkerStyle(20); ge_full->SetMarkerStyle(21); d_full->SetMarkerStyle(22);
    py_full->SetMarkerColor(kBlue); ge_full->SetMarkerColor(kRed); d_full->SetMarkerColor(kViolet);
    py_full->SetLineColor(kBlue); ge_full->SetLineColor(kRed); d_full->SetLineColor(kViolet);
        //scaling
    py_full->Scale(1/(double)py_full->Integral()); ge_full->Scale(1/(double)ge_full->Integral()); d_full->Scale(1/(double)d_full->Integral());
    
    //2D m v pt
        //axes
    py_m_pt->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); py_m_pt->GetYaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
    ge_m_pt->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ge_m_pt->GetYaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
    d_m_pt->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); d_m_pt->GetYaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
    py_m_pt->GetXaxis()->SetRangeUser(0,10); ge_m_pt->GetXaxis()->SetRangeUser(0,10); d_m_pt->GetXaxis()->SetRangeUser(0,10);
    py_m_pt->GetYaxis()->SetRangeUser(0,60); ge_m_pt->GetYaxis()->SetRangeUser(0,60); d_m_pt->GetYaxis()->SetRangeUser(0,60);
    py_m_pt->GetYaxis()->SetTitleOffset(1.3); ge_m_pt->GetYaxis()->SetTitleOffset(1.3); d_m_pt->GetYaxis()->SetTitleOffset(1.3);    

    //mean mass graph
        //axes
    py_mean_m_v_pt_gr->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    ge_mean_m_v_pt_gr->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    d_mean_m_v_pt_gr->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    py_mean_m_v_pt_gr->GetYaxis()->SetTitle("<M^{jet}> [GeV/c^{2}]");
    ge_mean_m_v_pt_gr->GetYaxis()->SetTitle("<M^{jet}> [GeV/c^{2}]");
    d_mean_m_v_pt_gr->GetYaxis()->SetTitle("<M^{jet}> [GeV/c^{2}]");
    py_mean_m_v_pt_gr->GetYaxis()->SetTitleOffset(1.3);
    ge_mean_m_v_pt_gr->GetYaxis()->SetTitleOffset(1.3);
    d_mean_m_v_pt_gr->GetYaxis()->SetTitleOffset(1.3);
    py_mean_m_v_pt_gr->SetTitle("");
    ge_mean_m_v_pt_gr->SetTitle("");
    d_mean_m_v_pt_gr->SetTitle("");
        //markers
    py_mean_m_v_pt_gr->SetMarkerStyle(20);
    ge_mean_m_v_pt_gr->SetMarkerStyle(20);
    d_mean_m_v_pt_gr->SetMarkerStyle(20);
    py_mean_m_v_pt_gr->SetMarkerColor(kBlue);
    ge_mean_m_v_pt_gr->SetMarkerColor(kRed);
    d_mean_m_v_pt_gr->SetMarkerColor(kViolet);
    
    //mean rms graph
        //axes
    py_rms_m_v_pt_gr->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    ge_rms_m_v_pt_gr->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    d_rms_m_v_pt_gr->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    py_rms_m_v_pt_gr->GetYaxis()->SetTitle("#sigma(M^{jet}) [GeV/c^{2}]");
    ge_rms_m_v_pt_gr->GetYaxis()->SetTitle("#sigma(M^{jet}) [GeV/c^{2}]");
    d_rms_m_v_pt_gr->GetYaxis()->SetTitle("#sigma(M^{jet}) [GeV/c^{2}]");
    py_rms_m_v_pt_gr->GetYaxis()->SetTitleOffset(1.3);
    ge_rms_m_v_pt_gr->GetYaxis()->SetTitleOffset(1.3);
    d_rms_m_v_pt_gr->GetYaxis()->SetTitleOffset(1.3);
    py_rms_m_v_pt_gr->SetTitle("");
    ge_rms_m_v_pt_gr->SetTitle("");
    d_rms_m_v_pt_gr->SetTitle("");
        //markers
    py_rms_m_v_pt_gr->SetMarkerStyle(20);
    ge_rms_m_v_pt_gr->SetMarkerStyle(20);
    d_rms_m_v_pt_gr->SetMarkerStyle(20);
    py_rms_m_v_pt_gr->SetMarkerColor(kBlue);
    ge_rms_m_v_pt_gr->SetMarkerColor(kRed);
    d_rms_m_v_pt_gr->SetMarkerColor(kViolet);
    
    //slices
        //axes
    py_10_15->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ge_10_15->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); d_10_15->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    py_15_20->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ge_15_20->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); d_15_20->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    py_20_25->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ge_20_25->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); d_20_25->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    py_25_30->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ge_25_30->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); d_25_30->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    py_30_plus->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ge_30_plus->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); d_30_plus->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    py_10_15->GetYaxis()->SetTitle("prob."); ge_10_15->GetYaxis()->SetTitle("prob."); d_10_15->GetYaxis()->SetTitle("prob.");
    py_15_20->GetYaxis()->SetTitle("prob."); ge_15_20->GetYaxis()->SetTitle("prob."); d_15_20->GetYaxis()->SetTitle("prob.");
    py_20_25->GetYaxis()->SetTitle("prob."); ge_20_25->GetYaxis()->SetTitle("prob."); d_20_25->GetYaxis()->SetTitle("prob.");
    py_25_30->GetYaxis()->SetTitle("prob."); ge_25_30->GetYaxis()->SetTitle("prob."); d_25_30->GetYaxis()->SetTitle("prob.");
    py_30_plus->GetYaxis()->SetTitle("prob."); ge_30_plus->GetYaxis()->SetTitle("prob."); d_30_plus->GetYaxis()->SetTitle("prob.");
    py_10_15->GetXaxis()->SetRangeUser(0,15); ge_10_15->GetXaxis()->SetRangeUser(0,15); d_10_15->GetXaxis()->SetRangeUser(0,15);
    py_15_20->GetXaxis()->SetRangeUser(0,15); ge_15_20->GetXaxis()->SetRangeUser(0,15); d_15_20->GetXaxis()->SetRangeUser(0,15);
    py_20_25->GetXaxis()->SetRangeUser(0,15); ge_20_25->GetXaxis()->SetRangeUser(0,15); d_20_25->GetXaxis()->SetRangeUser(0,15);
    py_25_30->GetXaxis()->SetRangeUser(0,15); ge_25_30->GetXaxis()->SetRangeUser(0,15); d_25_30->GetXaxis()->SetRangeUser(0,15);
    py_30_plus->GetXaxis()->SetRangeUser(0,15); ge_30_plus->GetXaxis()->SetRangeUser(0,15); d_30_plus->GetXaxis()->SetRangeUser(0,15);
        //markers
    py_10_15->SetMarkerStyle(20); ge_10_15->SetMarkerStyle(21); d_10_15->SetMarkerStyle(22);
    py_15_20->SetMarkerStyle(20); ge_15_20->SetMarkerStyle(21); d_15_20->SetMarkerStyle(22);
    py_20_25->SetMarkerStyle(20); ge_20_25->SetMarkerStyle(21); d_20_25->SetMarkerStyle(22);
    py_25_30->SetMarkerStyle(20); ge_25_30->SetMarkerStyle(21); d_25_30->SetMarkerStyle(22);
    py_30_plus->SetMarkerStyle(20); ge_30_plus->SetMarkerStyle(21); d_30_plus->SetMarkerStyle(22);
    py_10_15->SetMarkerColor(kBlue); ge_10_15->SetMarkerColor(kRed); d_10_15->SetMarkerColor(kViolet);
    py_15_20->SetMarkerColor(kBlue); ge_15_20->SetMarkerColor(kRed); d_15_20->SetMarkerColor(kViolet);
    py_20_25->SetMarkerColor(kBlue); ge_20_25->SetMarkerColor(kRed); d_20_25->SetMarkerColor(kViolet);
    py_25_30->SetMarkerColor(kBlue); ge_25_30->SetMarkerColor(kRed); d_25_30->SetMarkerColor(kViolet);
    py_30_plus->SetMarkerColor(kBlue); ge_30_plus->SetMarkerColor(kRed); d_30_plus->SetMarkerColor(kViolet);
    py_10_15->SetLineColor(kBlue); ge_10_15->SetLineColor(kRed); d_10_15->SetLineColor(kViolet);
    py_15_20->SetLineColor(kBlue); ge_15_20->SetLineColor(kRed); d_15_20->SetLineColor(kViolet);
    py_20_25->SetLineColor(kBlue); ge_20_25->SetLineColor(kRed); d_20_25->SetLineColor(kViolet);
    py_25_30->SetLineColor(kBlue); ge_25_30->SetLineColor(kRed); d_25_30->SetLineColor(kViolet);
    py_30_plus->SetLineColor(kBlue); ge_30_plus->SetLineColor(kRed); d_30_plus->SetLineColor(kViolet);
        //scaling
    py_10_15->Scale(1/(double)py_10_15->Integral());
    py_15_20->Scale(1/(double)py_15_20->Integral());
    py_20_25->Scale(1/(double)py_20_25->Integral());
    py_25_30->Scale(1/(double)py_25_30->Integral());
    py_30_plus->Scale(1/(double)py_30_plus->Integral());

    //ratios
    //axes
    ratio_10_15->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ratio_15_20->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    ratio_20_25->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]"); ratio_25_30->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    ratio_30_plus->GetXaxis()->SetTitle("M^{jet} [GeV/c^{2}]");
    ratio_10_15->GetYaxis()->SetTitle("Data / GEANT"); ratio_15_20->GetYaxis()->SetTitle("Data / GEANT");
    ratio_20_25->GetYaxis()->SetTitle("Data / GEANT"); ratio_25_30->GetYaxis()->SetTitle("Data / GEANT");
    ratio_30_plus->GetYaxis()->SetTitle("Data / GEANT");
    ratio_10_15->GetYaxis()->SetRangeUser(-0.5,2.5); ratio_15_20->GetYaxis()->SetRangeUser(-0.5,2.5);
    ratio_20_25->GetYaxis()->SetRangeUser(-0.5,2.5); ratio_25_30->GetYaxis()->SetRangeUser(-0.5,2.5);
    ratio_30_plus->GetYaxis()->SetRangeUser(-0.5,2.5);
    ratio_10_15->GetXaxis()->SetRangeUser(0,15); ratio_15_20->GetXaxis()->SetRangeUser(0,15);
    ratio_20_25->GetXaxis()->SetRangeUser(0,15); ratio_25_30->GetXaxis()->SetRangeUser(0,15);
    ratio_30_plus->GetXaxis()->SetRangeUser(0,15);
    ratio_10_15->GetYaxis()->SetTitleOffset(1.3); ratio_15_20->GetYaxis()->SetTitleOffset(1.3);
    ratio_20_25->GetYaxis()->SetTitleOffset(1.3); ratio_25_30->GetYaxis()->SetTitleOffset(1.3);
    ratio_30_plus->GetYaxis()->SetTitleOffset(1.3);
    //markers
    ratio_10_15->SetMarkerStyle(20); ratio_10_15->SetMarkerColor(kBlue); ratio_10_15->SetLineColor(kBlue);
    ratio_15_20->SetMarkerStyle(20); ratio_15_20->SetMarkerColor(kBlue); ratio_15_20->SetLineColor(kBlue);
    ratio_20_25->SetMarkerStyle(20); ratio_20_25->SetMarkerColor(kBlue); ratio_20_25->SetLineColor(kBlue);
    ratio_25_30->SetMarkerStyle(20); ratio_25_30->SetMarkerColor(kBlue); ratio_25_30->SetLineColor(kBlue);
    ratio_30_plus->SetMarkerStyle(20); ratio_30_plus->SetMarkerColor(kBlue); ratio_30_plus->SetLineColor(kBlue);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DRAWING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //full mass
    cfull->cd(); py_full->Draw(); ge_full->Draw("same"); d_full->Draw("same"); tfull->Draw("same");
    //2Ds
    cpympt->cd(); py_m_pt->Draw("colz"); cgempt->cd(); ge_m_pt->Draw("colz"); cdmpt->cd(); d_m_pt->Draw("colz");
    //mean
    cmean->cd(); py_mean_m_v_pt_gr->Draw("APZ"); ge_mean_m_v_pt_gr->Draw("PZsame"); d_mean_m_v_pt_gr->Draw("PZsame"); tmean->Draw("same");
    //rms
    crms->cd(); py_rms_m_v_pt_gr->Draw("APZ"); ge_rms_m_v_pt_gr->Draw("PZsame"); d_rms_m_v_pt_gr->Draw("PZsame"); trms->Draw("same");
    //slices
    cslices->cd(); cslices->cd(1); tslices->Draw(); gPad->SetTickx(2); gPad->SetTicky(2);
    cslices->cd(2); gPad->SetLogy(); py_10_15->Draw(); ge_10_15->Draw("same"); d_10_15->Draw("same"); t1015->Draw("same");
    cslices->cd(3); gPad->SetLogy(); py_15_20->Draw(); ge_15_20->Draw("same"); d_15_20->Draw("same"); t1520->Draw("same");
    cslices->cd(4); gPad->SetLogy(); py_20_25->Draw(); ge_20_25->Draw("same"); d_20_25->Draw("same"); t2025->Draw("same");
    cslices->cd(5); gPad->SetLogy(); py_25_30->Draw(); ge_25_30->Draw("same"); d_25_30->Draw("same"); t2530->Draw("same");
    cslices->cd(6); gPad->SetLogy(); py_30_plus->Draw(); ge_30_plus->Draw("same"); d_30_plus->Draw("same"); t30plus->Draw("same");
    //ratios
    cratios->cd(); cratios->cd(1); tratios->Draw();
    cratios->cd(2); ratio_10_15->Draw(); t1015->Draw("same");
    cratios->cd(3); ratio_15_20->Draw(); t1520->Draw("same");
    cratios->cd(4); ratio_20_25->Draw(); t2025->Draw("same");
    cratios->cd(5); ratio_25_30->Draw(); t2530->Draw("same");
    cratios->cd(6); ratio_30_plus->Draw(); t30plus->Draw("same");
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SAVING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    /*    
    cfull->SaveAs((out + "jet_mass_full_trig" + filetype).c_str());
    cpympt->SaveAs((out + "py_mass_v_pt_trig" + filetype).c_str());
    cgempt->SaveAs((out + "ge_mass_v_pt_trig" + filetype).c_str());
    cdmpt->SaveAs((out + "data_mass_v_pt_trig" + filetype).c_str());
    cmean->SaveAs((out + "mass_means_trig" + filetype).c_str());
    crms->SaveAs((out + "mass_rms_trig" + filetype).c_str());
    cslices->SaveAs((out + "mass_slices_trig" + filetype).c_str());
    cratios->SaveAs((out + "mass_ratios_geant_data_trig" + filetype).c_str());
    */
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 
    
  return;
}
