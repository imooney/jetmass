//! test unfolding macro to check with Isaac's code 

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>


using namespace std;

const int MaxnBayes2D = 1;

void test_unfolding(int RADIUS = 4){

  std::string trainname = Form("GeantToPythia_Match_ppRun12_200GeV_NoEff_NoBg_JP2_v3_largeBins_R04.root", RADIUS);

  TFile* ftrain = TFile::Open( Form("%s", trainname.c_str()));

  TFile * fout = new TFile("test_unfolding.root","RECREATE");

  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovToy;
  
  TH1D* IncJetpTMeasMCClosure1D=0;
  TH1D* IncJetpTMeasTestMCClosure1D=0;
  TH1D* IncJetpTTruthMCClosure1D=0;
  IncJetpTMeasMCClosure1D = (TH1D*)ftrain->Get("IncJetpTMeasMCClosure1D");
  IncJetpTMeasTestMCClosure1D = (TH1D*)ftrain->Get("IncJetpTMeasTestMCClosure1D");
  IncJetpTTruthMCClosure1D = (TH1D*)ftrain->Get("IncJetpTTruthMCClosure1D");

  RooUnfoldResponse* IncPtResponse_MCClosure = (RooUnfoldResponse*) ftrain->Get("IncPtResponse_MCClosure");

  TH1D** IncBayesUnfolded_JetPt_MCClosure_OppSide = new TH1D*[MaxnBayes2D];
  TH1D** IncBayesUnfolded_JetPt_MCClosure_SameSide = new TH1D*[MaxnBayes2D];
  
  TH1D** hMCClosure_UnfoldedJetPt_OppSide = new TH1D*[MaxnBayes2D];
  TH1D** hMCClosure_UnfoldedJetPt_SameSide = new TH1D*[MaxnBayes2D];

  for ( int nBayes2D=0; nBayes2D<1; ++nBayes2D ){

    RooUnfoldBayes IncBayesUnfold_Pt_MCClosure_OppSide ( IncPtResponse_MCClosure, IncJetpTMeasMCClosure1D, 4);
    IncBayesUnfold_Pt_MCClosure_OppSide.SetVerbose(1);
    IncBayesUnfolded_JetPt_MCClosure_OppSide[nBayes2D] = (TH1D*) IncBayesUnfold_Pt_MCClosure_OppSide.Hreco( errorTreatment );
    IncBayesUnfolded_JetPt_MCClosure_OppSide[nBayes2D]->SetName(Form("IncBayesUnfolded_JetPt_MCClosure_OppSide_iter%d", nBayes2D));
    //! Do the MC Closure Tests - take ratio with input Gen level spectra
    hMCClosure_UnfoldedJetPt_OppSide[nBayes2D] = (TH1D*)IncBayesUnfolded_JetPt_MCClosure_OppSide[nBayes2D]->Clone(Form("hMCClosure_UnfoldedJetPt_OppSide_iter%d", nBayes2D));
    hMCClosure_UnfoldedJetPt_OppSide[nBayes2D]->Divide(IncJetpTTruthMCClosure1D);

    RooUnfoldBayes IncBayesUnfold_Pt_MCClosure_SameSide ( IncPtResponse_MCClosure, IncJetpTMeasTestMCClosure1D, 4);
    IncBayesUnfold_Pt_MCClosure_SameSide.SetVerbose(1);
    IncBayesUnfolded_JetPt_MCClosure_SameSide[nBayes2D] = (TH1D*) IncBayesUnfold_Pt_MCClosure_SameSide.Hreco( errorTreatment );
    IncBayesUnfolded_JetPt_MCClosure_SameSide[nBayes2D]->SetName(Form("IncBayesUnfolded_JetPt_MCClosure_SameSide_iter%d", nBayes2D));
    //! Do the MC Closure Tests - take ratio with input Gen level spectra
    hMCClosure_UnfoldedJetPt_SameSide[nBayes2D] = (TH1D*)IncBayesUnfolded_JetPt_MCClosure_SameSide[nBayes2D]->Clone(Form("hMCClosure_UnfoldedJetPt_SameSide_iter%d", nBayes2D));
    hMCClosure_UnfoldedJetPt_SameSide[nBayes2D]->Divide(IncJetpTTruthMCClosure1D);
    
  }

  fout->Write();
  fout->Close();
  
}
