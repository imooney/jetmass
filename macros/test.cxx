#include <ctime>
#include <iostream>
#include <math.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TRandom.h>

using namespace std;

void TreetoHist (TFile *f, string system, string trig, vector<TH1D*> evts, vector<TH1D*> trks, vector<TH1D*> tows, vector<TH2D*> evts2D, vector<TH2D*> trks2D, vector<TH2D*> tows2D, vector<TH3D*> trks3D, vector<TH3D*> tows3D, double *meanEt, double *meanEtg2GeV) {  
  
  //these things will be used for the average Et/towerID determination
  const int nTows = 4800;
  int counts[nTows] = {0}; //initially none of the towers have fired
  int countsg2GeV[nTows] = {0};

  double evt_vtx = -9999; double vpdvz = -9999; double vzdiff = -9999; double bbc_coinc = -9999; double runID = -9999;
  double n_trks = -9999; double n_tows = -9999; double n_vertices = -9999; double n_globals = -9999;
  vector<double> * trackPt = 0; vector<double> * trackEta = 0; vector<double> * trackPhi = 0; vector<double> * trackDCA = 0;
  vector<double> * trackNhits = 0; vector<double> * trackNhitsposs = 0; //these are different sizes than other track vectors!
  vector<double> * towerEt = 0; vector<double> * towerEta = 0; vector<double> * towerPhi = 0; vector<double> * towerId = 0;
  
  TTree *t = (TTree*) f->Get(("QA"+trig).c_str()); //string = BBCMB, VPDMB, or JP2 for now  
  t->SetBranchAddress("n_trks", &n_trks);
  t->SetBranchAddress("n_tows", &n_tows);
  t->SetBranchAddress("bbc_coinc", &bbc_coinc);
  t->SetBranchAddress("evt_vtx", &evt_vtx); //tpcvz
  if (system != "pp") {
    t->SetBranchAddress("vpdvz", &vpdvz); //vpdvz
    t->SetBranchAddress("vzdiff", &vzdiff); // abs(tpcvz - vpdvz)
  }
  t->SetBranchAddress("runID", &runID);
  t->SetBranchAddress("n_vertices", &n_vertices);
  t->SetBranchAddress("n_globals", &n_globals);
  t->SetBranchAddress("trackPt", &trackPt);
  t->SetBranchAddress("trackEta", &trackEta);
  t->SetBranchAddress("trackPhi", &trackPhi);
  t->SetBranchAddress("trackDCA", &trackDCA);
  t->SetBranchAddress("trackNhits", &trackNhits);
  t->SetBranchAddress("trackNhitsposs", &trackNhitsposs);
  t->SetBranchAddress("towerEta", &towerEta);
  t->SetBranchAddress("towerPhi", &towerPhi);
  t->SetBranchAddress("towerEt", &towerEt);
  t->SetBranchAddress("towerId", &towerId);

  cout << ("RUNNING OVER "+system+trig+"! Entries: ").c_str() << t->GetEntries() << endl;
  const clock_t begin_time = clock();
  for (int i = 0; i < t->GetEntries(); ++ i) {
    if (i % 1000000 == 0) { cout << "still chuggin. " << i << endl; std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;}
    t->GetEntry(i);

    evts[0]->Fill(n_trks);
    evts[1]->Fill(n_tows);
    evts[2]->Fill(bbc_coinc);
    // cout<<"evt_vtx = "<<evt_vtx<<endl;
    evts[3]->Fill(evt_vtx);
    if (system != "pp") {
      evts[4]->Fill(vpdvz);
      evts[5]->Fill(vzdiff);
    }
    evts[6]->Fill(runID);
    evts[7]->Fill(n_vertices);
    evts[8]->Fill(n_globals);
    evts2D[0]->Fill(bbc_coinc, evt_vtx);
    evts2D[1]->Fill(bbc_coinc, n_trks);
    evts2D[2]->Fill(bbc_coinc, n_tows);
    evts2D[3]->Fill(bbc_coinc, n_vertices);
    evts2D[4]->Fill(bbc_coinc, n_globals);
    evts2D[5]->Fill(runID, evt_vtx);
    evts2D[6]->Fill(runID, bbc_coinc);

    for (int j = 0; j < trackNhits->size(); ++ j) {//these have different size than other track vectors
      trks[4]->Fill(trackNhits->at(j));
      trks[5]->Fill(trackNhitsposs->at(j));
      trks[6]->Fill(trackNhits->at(j) /(double) trackNhitsposs->at(j));
    }

    for (int j = 0; j < trackPt->size(); ++ j) { //all vectors of doubles in the branches should have the same size     
      trks[0]->Fill(trackPt->at(j));
      trks[1]->Fill(trackEta->at(j));
      trks[2]->Fill(trackPhi->at(j));
      trks[3]->Fill(trackDCA->at(j));
      if ((trackEta->at(j)<0 && trackEta->at(j)>-1) && (trackPhi->at(j)<0 && trackPhi->at(j)>-1)) {trks[7]->Fill(trackPt->at(j));trks[9]->Fill(trackDCA->at(j));}
      else {trks[8]->Fill(trackPt->at(j));trks[10]->Fill(trackDCA->at(j));}
      trks2D[0]->Fill(trackEta->at(j),trackPhi->at(j));
      trks2D[1]->Fill(trackDCA->at(j),n_trks);
      trks2D[2]->Fill(bbc_coinc, trackDCA->at(j));
      trks2D[3]->Fill(bbc_coinc, trackPt->at(j));
      trks2D[4]->Fill(runID, trackPt->at(j));
      trks2D[5]->Fill(runID, trackDCA->at(j));
      trks3D[0]->Fill(runID, trackEta->at(j), trackPhi->at(j));
      trks3D[1]->Fill(trackPt->at(j), trackEta->at(j), trackPhi->at(j));
    }//! track loop 
    for (int j = 0; j < towerEt->size(); ++ j) {
      tows[0]->Fill(towerEt->at(j));
      tows[1]->Fill(towerEta->at(j));
      tows[2]->Fill(towerPhi->at(j));
      tows[3]->Fill(towerId->at(j));
      tows2D[0]->Fill(towerEta->at(j),towerPhi->at(j));
      tows2D[1]->Fill(bbc_coinc, towerEt->at(j));
      tows2D[2]->Fill(runID, towerEt->at(j));
      tows2D[3]->Fill(runID, towerId->at(j));
      tows2D[4]->Fill(runID, towerEt->at(j));
      tows3D[0]->Fill(runID, towerEta->at(j), towerPhi->at(j));
      tows3D[1]->Fill(towerEt->at(j), towerEta->at(j), towerPhi->at(j));
      
      counts[((int) towerId->at(j))-1] ++;
      meanEt[((int) towerId->at(j))-1] += towerEt->at(j);
      if (towerEt->at(j) > 2) {
	countsg2GeV[((int) towerId->at(j))-1] ++;
	meanEtg2GeV[((int) towerId->at(j))-1] += towerEt->at(j);
      }
    }//! tower size loop
    
  }//! event loop
  
  for (int i = 0; i < nTows; ++ i) {
    if (counts[i] != 0) {
      meanEt[i] /= (double) counts[i];
    }
    else {
      meanEt[i] = 0;
    }
    if (countsg2GeV[i] != 0) {
      meanEtg2GeV[i] /= (double) countsg2GeV[i];
    }
    else {
      meanEtg2GeV[i] = 0;
    }
  }
  
  //! needs to be outside the event loop 
  t->ResetBranchAddresses();
  
  return;
}

int main () {
  
  TFile *fpAu = new TFile ("~/jetmass/out/QA/full_pAu_vpdvzcut_fixed.root","READ");
  TFile *fpp = new TFile ("~/jetmass/out/QA/full_pp.root","READ");
  TFile *fpAuBBCMB = new TFile("~/jetmass/out/QA/full_pAuBBCMB_vpdvzcut_fixed.root","READ");
  
  //MB hists

  TH1D* hn_trks_pAuBBCMB = new TH1D("hn_trks_pAuBBCMB","hn_trks_pAuBBCMB",150,0,150);
  TH1D* hn_tows_pAuBBCMB = new TH1D("hn_tows_pAuBBCMB","hn_tows_pAuBBCMB",250,0,500);
  TH1D* hbbc_coinc_pAuBBCMB = new TH1D("hbbc_coinc_pAuBBCMB","hbbc_coinc_pAuBBCMB",200,0,2600000);
  TH1D* hevt_vtx_pAuBBCMB = new TH1D("hevt_vtx_pAuBBCMB","hevt_vtx_pAuBBCMB",100,-35,35);
  TH1D* hvpdvz_pAuBBCMB = new TH1D("hvpdvz_pAuBBCMB","hvpdvz_pAuBBCMB",100,-35,35);
  TH1D* hvzdiff_pAuBBCMB = new TH1D("hvzdiff_pAuBBCMB","hvzdiff_pAuBBCMB",100,0,5);
  TH1D* hrunId_pAuBBCMB = new TH1D("hrunId_pAuBBCMB","hrunId_pAuBBCMB",1000,16125000,16159025);
  TH1D* hn_vertices_pAuBBCMB = new TH1D("hn_vertices_pAuBBCMB","hn_vertices_pAuBBCMB",40,0,40);
  TH1D* hn_globals_pAuBBCMB = new TH1D("hn_globals_pAuBBCMB","hn_globals_pAuBBCMB",200,0,4000);
  

  //raw track info (before QA cuts):
  TH1D* htrackNhits_pAuBBCMB = new TH1D("htrackNhits_pAuBBCMB","htrackNhits_pAuBBCMB",50,0,50);
  TH1D* htrackNhitsposs_pAuBBCMB = new TH1D("htrackNhitsposs_pAuBBCMB","htrackNhitsposs_pAuBBCMB",50,0,50);
  TH1D* htrackNhitsratio_pAuBBCMB = new TH1D("htrackNhitsratio_pAuBBCMB","htrackNhitsratio_pAuBBCMB",100,0,1);
  //
  
  TH1D* htrackPt_pAuBBCMB = new TH1D("htrackPt_pAuBBCMB","htrackPt_pAuBBCMB",200,0,40);
  TH1D* htrackEta_pAuBBCMB = new TH1D("htrackEta_pAuBBCMB","htrackEta_pAuBBCMB",40,-1,1);
  TH1D* htrackPhi_pAuBBCMB = new TH1D("htrackPhi_pAuBBCMB","htrackPhi_pAuBBCMB",120,-M_PI,M_PI);
  TH1D* htrackDCA_pAuBBCMB = new TH1D("htrackDCA_pAuBBCMB","htrackDCA_pAuBBCMB",200,0,1);
  

  //in/out of dead zone in TPC
  TH1D* htrackPt_indeadzone_pAuBBCMB = new TH1D("htrackPt_indeadzone_pAuBBCMB","htrackPt_indeadzone_pAuBBCMB",200,0,40);
  TH1D* htrackPt_outdeadzone_pAuBBCMB = new TH1D("htrackPt_outdeadzone_pAuBBCMB","htrackPt_outdeadzone_pAuBBCMB",200,0,40);
  TH1D* htrackDCA_indeadzone_pAuBBCMB = new TH1D("htrackDCA_indeadzone_pAuBBCMB","htrackDCA_indeadzone_pAuBBCMB",200,0,1);
  TH1D* htrackDCA_outdeadzone_pAuBBCMB = new TH1D("htrackDCA_outdeadzone_pAuBBCMB","htrackDCA_outdeadzone_pAuBBCMB",200,0,1); 
  //TH1D* htowerEt_indeadzone_pAuBBCMB = new TH1D("htowerEt_indeadzone_pAuBBCMB","htowerEt_indeadzone_pAuBBCMB",200,0,40);
  //TH1D* htowerEt_outdeadzone_pAuBBCMB = new TH1D("htowerEt_outdeadzone_pAuBBCMB","htowerEt_outdeadzone_pAuBBCMB",200,0,40);
  //TH1D* htrackNhits_indeadzone_pAuBBCMB = new TH1D("htrackNhits_indeadzone_pAuBBCMB","htrackNhits_indeadzone_pAuBBCMB",50,0,50);
  //TH1D* htrackNhits_outdeadzone_pAuBBCMB = new TH1D("htrackNhits_outdeadzone_pAuBBCMB","htrackNhits_outdeadzone_pAuBBCMB",50,0,50);
 

  TH1D* htowerEt_pAuBBCMB = new TH1D("htowerEt_pAuBBCMB","htowerEt_pAuBBCMB",200,0,40);
  TH1D* htowerEta_pAuBBCMB = new TH1D("htowerEta_pAuBBCMB","htowerEta_pAuBBCMB",40,-1,1);
  TH1D* htowerPhi_pAuBBCMB = new TH1D("htowerPhi_pAuBBCMB","htowerPhi_pAuBBCMB",120,-M_PI,M_PI);
  TH1D* htowerId_pAuBBCMB = new TH1D("htowerId_pAuBBCMB","htowerId_pAuBBCMB",4800,0.5,4800.5);

  TH2D* htrackEta_Phi_pAuBBCMB = new TH2D("htrackEta_Phi_pAuBBCMB","pAu BBCMB;#eta_{trk};#phi_{trk}",40,-1,1,120,-M_PI,M_PI);
  TH2D* htrackDCA_n_trks_pAuBBCMB = new TH2D("htrackDCA_n_trks_pAuBBCMB","pAu BBCMB;track DCA;N_{trk}",200,0,1,150,0,150);
  
  TH2D* htowerEta_Phi_pAuBBCMB = new TH2D("htowerEta_Phi_pAuBBCMB","pAu BBCMB;#eta_{tow};#phi_{tow}",40,-1,1,120,-M_PI,M_PI);

  TH3D* htrackPt_Eta_Phi_pAuBBCMB = new TH3D("htrackPt_Eta_Phi_pAuBBCMB","pAu BBCMB;#p^{trk}_{T};#eta_{trk};#phi_{trk}",200,0,40,40,-1,1,120,-M_PI,M_PI);
  TH3D* htowerEt_Eta_Phi_pAuBBCMB = new TH3D("htowerEt_Eta_Phi_pAuBBCMB","pAu BBCMB;#E^{tow}_{T};#eta_{tow};#phi_{tow}",200,0,40,40,-1,1,120,-M_PI,M_PI);

  TH2D* hbbc_coinc_evt_vtx_pAuBBCMB = new TH2D("hbbc_coinc_evt_vtx_pAuBBCMB","pAu BBCMB;BBC coincidence;v_{z} [cm]",200,0,2600000,100,-35,35);
  TH2D* hbbc_coinc_n_trks_pAuBBCMB = new TH2D("hbbc_coinc_n_trks_pAuBBCMB","pAu BBCMB;BBC coincidence;N_{trk}",200,0,2600000,150,0,150);
  TH2D* hbbc_coinc_n_tows_pAuBBCMB = new TH2D("hbbc_coinc_n_tows_pAuBBCMB","pAu BBCMB;BBC coincidence;N_{tow}",200,0,2600000,250,0,500);
  TH2D* hbbc_coinc_n_vertices_pAuBBCMB = new TH2D("hbbc_coinc_n_vertices_pAuBBCMB","pAu BBCMB;BBC coincidence;N_{vertices}",200,0,2600000,40,0,40);
  TH2D* hbbc_coinc_n_globals_pAuBBCMB = new TH2D("hbbc_coinc_n_globals_pAuBBCMB","pAu BBCMB;BBC coincidence;N_{globals}",200,0,2600000,200,0,4000);
  TH2D* hbbc_coinc_trackDCA_pAuBBCMB = new TH2D("hbbc_coinc_trackDCA_pAuBBCMB","pAu BBCMB;BBC coincidence;track DCA",200,0,2600000,200,0,1);
  TH2D* hbbc_coinc_trackPt_pAuBBCMB = new TH2D("hbbc_coinc_trackPt_pAuBBCMB","pAu BBCMB;BBC coincidence;p^{trk}_{T} [GeV/c]",200,0,2600000,200,0,40);
  TH2D* hbbc_coinc_towerEt_pAuBBCMB = new TH2D("hbbc_coinc_towerEt_pAuBBCMB","pAu BBCMB;BBC coincidence;E^{tow}_{T} [GeV]",200,0,2600000,200,0,40);

  TH2D* hrunId_trackPt_pAuBBCMB = new TH2D("hrunId_trackPt_pAuBBCMB","pAu BBCMB;run ID; p^{trk}_{T} [GeV/c]",1000,16125000,16159025,200,0,40);
  TH2D* hrunId_towerEt_pAuBBCMB = new TH2D("hrunId_towerEt_pAuBBCMB","pAu BBCMB;run ID; E^{tow}_{T} [GeV]",1000,16125000,16159025,200,0,40);
  //temp:
  TH2D* hrunId_towerEt_pAuBBCMB_dummy = new TH2D("hrunId_towerEt_pAuBBCMB_dummy","pAu BBCMB;run ID; E^{tow}_{T} [GeV]",1000,16125000,16159025,200,0,40);

  TH2D* hrunId_towerId_pAuBBCMB = new TH2D("hrunId_towerId_pAuBBCMB","pAu BBCMB;run ID; tower ID",1000,16125000,16159025,4800,0.5,4800.5);

  TH2D* hrunId_evt_vtx_pAuBBCMB = new TH2D("hrunId_evt_vtx_pAuBBCMB","pAu BBCMB; runID; v_{z} [cm]",1000,16125000,16159025,100,-35,35);
  TH2D* hrunId_bbc_coinc_pAuBBCMB = new TH2D("hrunId_bbc_coinc_pAuBBCMB","pAu BBCMB; runID; BBC coincidence [kHz]",1000,16125000,16159025,200,0,2600000);
  TH2D* hrunId_trackDCA_pAuBBCMB = new TH2D("hrunId_trackDCA_pAuBBCMB","pAu BBCMB; runID; track DCA",1000,16125000,16159025,200,0,1);
  

  TH3D* hrunId_trackEta_Phi_pAuBBCMB = new TH3D("hrunId_trackEta_Phi_pAuBBCMB","pAu BBCMB;runID;#eta_{trk};#phi_{trk}",1000,16125000,16159025,40,-1,1,120,-M_PI,M_PI);
  TH3D* hrunId_towerEta_Phi_pAuBBCMB = new TH3D("hrunId_towerEta_Phi_pAuBBCMB","pAu BBCMB;runID;#eta_{tow};#phi_{tow}",1000,16125000,16159025,40,-1,1,120,-M_PI,M_PI);

  vector<TH1D*> evts_pAuBBCMB = {hn_trks_pAuBBCMB, hn_tows_pAuBBCMB, hbbc_coinc_pAuBBCMB, hevt_vtx_pAuBBCMB, hvpdvz_pAuBBCMB, hvzdiff_pAuBBCMB, hrunId_pAuBBCMB, hn_vertices_pAuBBCMB, hn_globals_pAuBBCMB};
  vector<TH1D*> trks_pAuBBCMB = {htrackPt_pAuBBCMB, htrackEta_pAuBBCMB, htrackPhi_pAuBBCMB, htrackDCA_pAuBBCMB,htrackNhits_pAuBBCMB,htrackNhitsposs_pAuBBCMB,htrackNhitsratio_pAuBBCMB,htrackPt_indeadzone_pAuBBCMB,htrackPt_outdeadzone_pAuBBCMB,htrackDCA_indeadzone_pAuBBCMB,htrackDCA_outdeadzone_pAuBBCMB/*,htrackNhits_indeadzone_pAuBBCMB,htrackNhits_outdeadzone_pAuBBCMB*/};
  vector<TH1D*> tows_pAuBBCMB = {htowerEt_pAuBBCMB, htowerEta_pAuBBCMB, htowerPhi_pAuBBCMB, htowerId_pAuBBCMB/*,htowerEt_indeadzone_pAuBBCMB,htowerEt_outdeadzone_pAuBBCMB*/};

  vector<TH2D*> evts2D_pAuBBCMB = {hbbc_coinc_evt_vtx_pAuBBCMB,hbbc_coinc_n_trks_pAuBBCMB,hbbc_coinc_n_tows_pAuBBCMB,hbbc_coinc_n_vertices_pAuBBCMB,hbbc_coinc_n_globals_pAuBBCMB,hrunId_evt_vtx_pAuBBCMB,hrunId_bbc_coinc_pAuBBCMB};

  vector<TH2D*> trks2D_pAuBBCMB = {htrackEta_Phi_pAuBBCMB, htrackDCA_n_trks_pAuBBCMB,hbbc_coinc_trackDCA_pAuBBCMB,hbbc_coinc_trackPt_pAuBBCMB,hrunId_trackPt_pAuBBCMB,hrunId_trackDCA_pAuBBCMB};
  vector<TH2D*> tows2D_pAuBBCMB = {htowerEta_Phi_pAuBBCMB,hbbc_coinc_towerEt_pAuBBCMB,hrunId_towerEt_pAuBBCMB,hrunId_towerId_pAuBBCMB,hrunId_towerEt_pAuBBCMB_dummy};
    
  vector<TH3D*> trks3D_pAuBBCMB = {hrunId_trackEta_Phi_pAuBBCMB, htrackPt_Eta_Phi_pAuBBCMB};

  vector<TH3D*> tows3D_pAuBBCMB = {hrunId_towerEta_Phi_pAuBBCMB, htowerEt_Eta_Phi_pAuBBCMB};
  

  //JP2 hists

  TH1D* hn_trks_pAuJP2 = new TH1D("hn_trks_pAuJP2","hn_trks_pAuJP2",150,0,150);
  TH1D* hn_tows_pAuJP2 = new TH1D("hn_tows_pAuJP2","hn_tows_pAuJP2",250,0,500);
  TH1D* hbbc_coinc_pAuJP2 = new TH1D("hbbc_coinc_pAuJP2","hbbc_coinc_pAuJP2",200,0,2600000);
  TH1D* hevt_vtx_pAuJP2 = new TH1D("hevt_vtx_pAuJP2","hevt_vtx_pAuJP2",100,-35,35);
  TH1D* hvpdvz_pAuJP2 = new TH1D("hvpdvz_pAuJP2","hvpdvz_pAuJP2",100,-35,35);
  TH1D* hvzdiff_pAuJP2 = new TH1D("hvzdiff_pAuJP2","hvzdiff_pAuJP2",100,0,5);
  TH1D* hrunId_pAuJP2 = new TH1D("hrunId_pAuJP2","hrunId_pAuJP2",1000,16125000,16159025);
  TH1D* hn_vertices_pAuJP2 = new TH1D("hn_vertices_pAuJP2","hn_vertices_pAuJP2",40,0,40);
  TH1D* hn_globals_pAuJP2 = new TH1D("hn_globals_pAuJP2","hn_globals_pAuJP2",200,0,4000);

  TH1D* hn_trks_ppJP2 = new TH1D("hn_trks_ppJP2","hn_trks_ppJP2",150,0,150);
  TH1D* hn_tows_ppJP2 = new TH1D("hn_tows_ppJP2","hn_tows_ppJP2",250,0,500);
  TH1D* hbbc_coinc_ppJP2 = new TH1D("hbbc_coinc_ppJP2","hbbc_coinc_ppJP2",200,0,2600000);
  TH1D* hevt_vtx_ppJP2 = new TH1D("hevt_vtx_ppJP2","hevt_vtx_ppJP2",100,-35,35);
  TH1D* hvpdvz_ppJP2 = new TH1D("hvpdvz_ppJP2","hvpdvz_ppJP2",100,-35,35);
  TH1D* hvzdiff_ppJP2 = new TH1D("hvzdiff_ppJP2","hvzdiff_ppJP2",100,0,5);
  TH1D* hrunId_ppJP2 = new TH1D("hrunId_ppJP2","hrunId_ppJP2",1000,16125000,16159025);//this is actually a dummy for now
  TH1D* hn_vertices_ppJP2 = new TH1D("hn_vertices_ppJP2","hn_vertices_ppJP2",40,0,40);
  TH1D* hn_globals_ppJP2 = new TH1D("hn_globals_ppJP2","hn_globals_ppJP2",200,0,4000);
 
  //raw track info (before QA cuts):
  TH1D* htrackNhits_pAuJP2 = new TH1D("htrackNhits_pAuJP2","htrackNhits_pAuJP2",50,0,50);
  TH1D* htrackNhitsposs_pAuJP2 = new TH1D("htrackNhitsposs_pAuJP2","htrackNhitsposs_pAuJP2",50,0,50);
  TH1D* htrackNhitsratio_pAuJP2 = new TH1D("htrackNhitsratio_pAuJP2","htrackNhitsratio_pAuJP2",100,0,1);

  TH1D* htrackNhits_ppJP2 = new TH1D("htrackNhits_ppJP2","htrackNhits_ppJP2",50,0,50);
  TH1D* htrackNhitsposs_ppJP2 = new TH1D("htrackNhitsposs_ppJP2","htrackNhitsposs_ppJP2",50,0,50);
  TH1D* htrackNhitsratio_ppJP2 = new TH1D("htrackNhitsratio_ppJP2","htrackNhitsratio_ppJP2",100,0,1);


  TH1D* htrackPt_pAuJP2 = new TH1D("htrackPt_pAuJP2","htrackPt_pAuJP2",200,0,40);
  TH1D* htrackEta_pAuJP2 = new TH1D("htrackEta_pAuJP2","htrackEta_pAuJP2",40,-1,1);
  TH1D* htrackPhi_pAuJP2 = new TH1D("htrackPhi_pAuJP2","htrackPhi_pAuJP2",120,-M_PI,M_PI);
  TH1D* htrackDCA_pAuJP2 = new TH1D("htrackDCA_pAuJP2","htrackDCA_pAuJP2",200,0,1);
 
    //in/out of dead zone in TPC
  TH1D* htrackPt_indeadzone_pAuJP2 = new TH1D("htrackPt_indeadzone_pAuJP2","htrackPt_indeadzone_pAuJP2",200,0,40);
  TH1D* htrackPt_outdeadzone_pAuJP2 = new TH1D("htrackPt_outdeadzone_pAuJP2","htrackPt_outdeadzone_pAuJP2",200,0,40);
  TH1D* htrackDCA_indeadzone_pAuJP2 = new TH1D("htrackDCA_indeadzone_pAuJP2","htrackDCA_indeadzone_pAuJP2",200,0,1);
  TH1D* htrackDCA_outdeadzone_pAuJP2 = new TH1D("htrackDCA_outdeadzone_pAuJP2","htrackDCA_outdeadzone_pAuJP2",200,0,1); 
  //TH1D* htowerEt_indeadzone_pAuJP2 = new TH1D("htowerEt_indeadzone_pAuJP2","htowerEt_indeadzone_pAuJP2",200,0,40);
  //TH1D* htowerEt_outdeadzone_pAuJP2 = new TH1D("htowerEt_outdeadzone_pAuJP2","htowerEt_outdeadzone_pAuJP2",200,0,40);
  //TH1D* htrackNhits_indeadzone_pAuJP2 = new TH1D("htrackNhits_indeadzone_pAuJP2","htrackNhits_indeadzone_pAuJP2",50,0,50);
  //TH1D* htrackNhits_outdeadzone_pAuJP2 = new TH1D("htrackNhits_outdeadzone_pAuJP2","htrackNhits_outdeadzone_pAuJP2",50,0,50);
 

   //in/out of dead zone in TPC
  TH1D* htrackPt_indeadzone_ppJP2 = new TH1D("htrackPt_indeadzone_ppJP2","htrackPt_indeadzone_ppJP2",200,0,40);
  TH1D* htrackPt_outdeadzone_ppJP2 = new TH1D("htrackPt_outdeadzone_ppJP2","htrackPt_outdeadzone_ppJP2",200,0,40);
  TH1D* htrackDCA_indeadzone_ppJP2 = new TH1D("htrackDCA_indeadzone_ppJP2","htrackDCA_indeadzone_ppJP2",200,0,1);
  TH1D* htrackDCA_outdeadzone_ppJP2 = new TH1D("htrackDCA_outdeadzone_ppJP2","htrackDCA_outdeadzone_ppJP2",200,0,1); 
  //TH1D* htowerEt_indeadzone_ppJP2 = new TH1D("htowerEt_indeadzone_ppJP2","htowerEt_indeadzone_ppJP2",200,0,40);
  //TH1D* htowerEt_outdeadzone_ppJP2 = new TH1D("htowerEt_outdeadzone_ppJP2","htowerEt_outdeadzone_ppJP2",200,0,40);
  //TH1D* htrackNhits_indeadzone_ppJP2 = new TH1D("htrackNhits_indeadzone_ppJP2","htrackNhits_indeadzone_ppJP2",50,0,50);
  //TH1D* htrackNhits_outdeadzone_ppJP2 = new TH1D("htrackNhits_outdeadzone_ppJP2","htrackNhits_outdeadzone_ppJP2",50,0,50);
 
 
  TH1D* htrackPt_ppJP2 = new TH1D("htrackPt_ppJP2","htrackPt_ppJP2",200,0,40);
  TH1D* htrackEta_ppJP2 = new TH1D("htrackEta_ppJP2","htrackEta_ppJP2",40,-1,1);
  TH1D* htrackPhi_ppJP2 = new TH1D("htrackPhi_ppJP2","htrackPhi_ppJP2",120,-M_PI,M_PI);
  TH1D* htrackDCA_ppJP2 = new TH1D("htrackDCA_ppJP2","htrackDCA_ppJP2",200,0,1);

  TH1D* htowerEt_pAuJP2 = new TH1D("htowerEt_pAuJP2","htowerEt_pAuJP2",200,0,40);
  TH1D* htowerEta_pAuJP2 = new TH1D("htowerEta_pAuJP2","htowerEta_pAuJP2",40,-1,1);
  TH1D* htowerPhi_pAuJP2 = new TH1D("htowerPhi_pAuJP2","htowerPhi_pAuJP2",120,-M_PI,M_PI);
  TH1D* htowerId_pAuJP2 = new TH1D("htowerId_pAuJP2","htowerId_pAuJP2",4800,0.5,4800.5);

  TH1D* htowerEt_ppJP2 = new TH1D("htowerEt_ppJP2","htowerEt_ppJP2",200,0,40);
  TH1D* htowerEta_ppJP2 = new TH1D("htowerEta_ppJP2","htowerEta_ppJP2",40,-1,1);
  TH1D* htowerPhi_ppJP2 = new TH1D("htowerPhi_ppJP2","htowerPhi_ppJP2",120,-M_PI,M_PI);
  TH1D* htowerId_ppJP2 = new TH1D("htowerId_ppJP2","htowerId_ppJP2",4800,0.5,4800.5);

  TH2D* htrackEta_Phi_pAuJP2 = new TH2D("htrackEta_Phi_pAuJP2","pAu JP2;#eta_{trk};#phi_{trk}",40,-1,1,120,-M_PI,M_PI);
  TH2D* htrackDCA_n_trks_pAuJP2 = new TH2D("htrackDCA_n_trks_pAuJP2","pAu JP2;track DCA;N_{trk}",200,0,1,150,0,150);
  
  TH2D* htrackEta_Phi_ppJP2 = new TH2D("htrackEta_Phi_ppJP2","pp JP2;#eta_{trk};#phi_{trk}",40,-1,1,120,-M_PI,M_PI);
  TH2D* htrackDCA_n_trks_ppJP2 = new TH2D("htrackDCA_n_trks_ppJP2","pp JP2;track DCA;N_{trk}",200,0,1,150,0,150);
  
  TH2D* htowerEta_Phi_pAuJP2 = new TH2D("htowerEta_Phi_pAuJP2","pAu JP2;#eta_{tow};#phi_{tow}",40,-1,1,120,-M_PI,M_PI);
  
  TH2D* htowerEta_Phi_ppJP2 = new TH2D("htowerEta_Phi_ppJP2","pp JP2;#eta_{tow};#phi_{tow}",40,-1,1,120,-M_PI,M_PI);

  TH3D* htrackPt_Eta_Phi_pAuJP2 = new TH3D("htrackPt_Eta_Phi_pAuJP2","pAu JP2;#p^{trk}_{T};#eta_{trk};#phi_{trk}",200,0,40,40,-1,1,120,-M_PI,M_PI);
  TH3D* htowerEt_Eta_Phi_pAuJP2 = new TH3D("htowerEt_Eta_Phi_pAuJP2","pAu JP2;#E^{tow}_{T};#eta_{tow};#phi_{tow}",200,0,40,40,-1,1,120,-M_PI,M_PI);

  TH3D* htrackPt_Eta_Phi_ppJP2 = new TH3D("htrackPt_Eta_Phi_ppJP2","pp JP2;#p^{trk}_{T};#eta_{trk};#phi_{trk}",200,0,40,40,-1,1,120,-M_PI,M_PI);
  TH3D* htowerEt_Eta_Phi_ppJP2 = new TH3D("htowerEt_Eta_Phi_ppJP2","pp JP2;#E^{tow}_{T};#eta_{tow};#phi_{tow}",200,0,40,40,-1,1,120,-M_PI,M_PI);

  TH2D* hbbc_coinc_evt_vtx_pAuJP2 = new TH2D("hbbc_coinc_evt_vtx_pAuJP2","pAu JP2;BBC coincidence;v_{z} [cm]",200,0,2600000,100,-35,35);
  TH2D* hbbc_coinc_n_trks_pAuJP2 = new TH2D("hbbc_coinc_n_trks_pAuJP2","pAu JP2;BBC coincidence;N_{trk}",200,0,2600000,150,0,150);
  TH2D* hbbc_coinc_n_tows_pAuJP2 = new TH2D("hbbc_coinc_n_tows_pAuJP2","pAu JP2;BBC coincidence;N_{tow}",200,0,2600000,250,0,500);
  TH2D* hbbc_coinc_n_vertices_pAuJP2 = new TH2D("hbbc_coinc_n_vertices_pAuJP2","pAu JP2;BBC coincidence;N_{vertices}",200,0,2600000,40,0,40);//!~!
  TH2D* hbbc_coinc_n_globals_pAuJP2 = new TH2D("hbbc_coinc_n_globals_pAuJP2","pAu JP2;BBC coincidence;N_{globals}",200,0,2600000,200,0,4000);
  TH2D* hbbc_coinc_trackDCA_pAuJP2 = new TH2D("hbbc_coinc_trackDCA_pAuJP2","pAu JP2;BBC coincidence;track DCA",200,0,2600000,200,0,1);
  TH2D* hbbc_coinc_trackPt_pAuJP2 = new TH2D("hbbc_coinc_trackPt_pAuJP2","pAu JP2;BBC coincidence;p^{trk}_{T} [GeV/c]",200,0,2600000,200,0,40);
  TH2D* hbbc_coinc_towerEt_pAuJP2 = new TH2D("hbbc_coinc_towerEt_pAuJP2","pAu JP2;BBC coincidence;E^{tow}_{T} [GeV]",200,0,2600000,200,0,40);

  TH2D* hbbc_coinc_evt_vtx_ppJP2 = new TH2D("hbbc_coinc_evt_vtx_ppJP2","pp JP2;BBC coincidence;v_{z} [cm]",200,0,2600000,100,-35,35);
  TH2D* hbbc_coinc_n_trks_ppJP2 = new TH2D("hbbc_coinc_n_trks_ppJP2","pp JP2;BBC coincidence;N_{trk}",200,0,2600000,150,0,150);
  TH2D* hbbc_coinc_n_tows_ppJP2 = new TH2D("hbbc_coinc_n_tows_ppJP2","pp JP2;BBC coincidence;N_{tow}",200,0,2600000,250,0,500);
  TH2D* hbbc_coinc_n_vertices_ppJP2 = new TH2D("hbbc_coinc_n_vertices_ppJP2","pp JP2;BBC coincidence;N_{vertices}",200,0,2600000,40,0,40);//!~!
  TH2D* hbbc_coinc_n_globals_ppJP2 = new TH2D("hbbc_coinc_n_globals_ppJP2","pp JP2;BBC coincidence;N_{globals}",200,0,2600000,200,0,4000);
  TH2D* hbbc_coinc_trackDCA_ppJP2 = new TH2D("hbbc_coinc_trackDCA_ppJP2","pp JP2;BBC coincidence;track DCA",200,0,2600000,200,0,1);
  TH2D* hbbc_coinc_trackPt_ppJP2 = new TH2D("hbbc_coinc_trackPt_ppJP2","pp JP2;BBC coincidence;p^{trk}_{T} [GeV/c]",200,0,2600000,200,0,40);
  TH2D* hbbc_coinc_towerEt_ppJP2 = new TH2D("hbbc_coinc_towerEt_ppJP2","pp JP2;BBC coincidence;E^{tow}_{T} [GeV]",200,0,2600000,200,0,40);


  TH2D* hrunId_trackPt_pAuJP2 = new TH2D("hrunId_trackPt_pAuJP2","pAu JP2;run ID; p^{trk}_{T} [GeV/c]",1000,16125000,16159025,200,0,40);
  TH2D* hrunId_towerEt_pAuJP2 = new TH2D("hrunId_towerEt_pAuJP2","pAu JP2;run ID; E^{tow}_{T} [GeV]",1000,16125000,16159025,200,0,40);
  //temp:
  TH2D* hrunId_towerEt_pAuJP2_badzone = new TH2D("hrunId_towerEt_pAuJP2_badzone","pAu JP2;run ID; E^{tow}_{T} [GeV]",100,16135000,16135100,200,0,40);

  TH2D* hrunId_towerId_pAuJP2 = new TH2D("hrunId_towerId_pAuJP2","pAu JP2;run ID; tower ID",1000,16125000,16159025,4800,0.5,4800.5);
  
  TH2D* hrunId_trackPt_ppJP2 = new TH2D("hrunId_trackPt_ppJP2","pp JP2;run ID; p^{trk}_{T} [GeV/c]",1000,16125000,16159025,200,0,40);
  TH2D* hrunId_towerEt_ppJP2 = new TH2D("hrunId_towerEt_ppJP2","pp JP2;run ID; E^{tow}_{T} [GeV]",1000,16125000,16159025,200,0,40);
  //temp:
  TH2D* hrunId_towerEt_ppJP2_dummy = new TH2D("hrunId_towerEt_ppJP2_dummy","pp JP2;run ID; E^{tow}_{T} [GeV]",1000,16125000,16159025,200,0,40);
  
  TH2D* hrunId_towerId_ppJP2 = new TH2D("hrunId_towerId_ppJP2","pp JP2;run ID; tower ID",1000,16125000,16159025,4800,0.5,4800.5);
  
  TH3D* hrunId_trackEta_Phi_pAuJP2 = new TH3D("hrunId_trackEta_Phi_pAuJP2","pAu JP2;runID;#eta_{trk};#phi_{trk}",1000,16125000,16159025,40,-1,1,120,-M_PI,M_PI);
  TH3D* hrunId_towerEta_Phi_pAuJP2 = new TH3D("hrunId_towerEta_Phi_pAuJP2","pAu JP2;runID;#eta_{tow};#phi_{tow}",1000,16125000,16159025,40,-1,1,120,-M_PI,M_PI);

  TH3D* hrunId_trackEta_Phi_ppJP2 = new TH3D("hrunId_trackEta_Phi_ppJP2","pp JP2;runID;#eta_{trk};#phi_{trk}",1000,16125000,16159025,40,-1,1,120,-M_PI,M_PI);
  TH3D* hrunId_towerEta_Phi_ppJP2 = new TH3D("hrunId_towerEta_Phi_ppJP2","pp JP2;runID;#eta_{tow};#phi_{tow}",1000,16125000,16159025,40,-1,1,120,-M_PI,M_PI);

  TH2D* hrunId_evt_vtx_pAuJP2 = new TH2D("hrunId_evt_vtx_pAuJP2","pAu JP2; runID; v_{z} [cm]",1000,16125000,16159025,100,-35,35);
  TH2D* hrunId_bbc_coinc_pAuJP2 = new TH2D("hrunId_bbc_coinc_pAuJP2","pAu JP2; runID; BBC coincidence [kHz]",1000,16125000,16159025,200,0,2600000);
  TH2D* hrunId_trackDCA_pAuJP2 = new TH2D("hrunId_trackDCA_pAuJP2","pAu JP2; runID; track DCA",1000,16125000,16159025,200,0,1);
  
  //reminder:pp hists for run ID dependence are dummies
  TH2D* hrunId_evt_vtx_ppJP2 = new TH2D("hrunId_evt_vtx_ppJP2","pp JP2; runID; v_{z} [cm]",1000,16125000,16159025,100,-35,35);
  TH2D* hrunId_bbc_coinc_ppJP2 = new TH2D("hrunId_bbc_coinc_ppJP2","pp JP2; runID; BBC coincidence [kHz]",1000,16125000,16159025,200,0,2600000);
  TH2D* hrunId_trackDCA_ppJP2 = new TH2D("hrunId_trackDCA_ppJP2","pp JP2; runID; track DCA",1000,16125000,16159025,200,0,1);
  

  vector<TH1D*> evts_pAuJP2 = {hn_trks_pAuJP2, hn_tows_pAuJP2, hbbc_coinc_pAuJP2, hevt_vtx_pAuJP2, hvpdvz_pAuJP2, hvzdiff_pAuJP2, hrunId_pAuJP2, hn_vertices_pAuJP2, hn_globals_pAuJP2};
  vector<TH1D*> trks_pAuJP2 = {htrackPt_pAuJP2, htrackEta_pAuJP2, htrackPhi_pAuJP2, htrackDCA_pAuJP2,htrackNhits_pAuJP2,htrackNhitsposs_pAuJP2,htrackNhitsratio_pAuJP2,htrackPt_indeadzone_pAuJP2,htrackPt_outdeadzone_pAuJP2,htrackDCA_indeadzone_pAuJP2,htrackDCA_outdeadzone_pAuJP2/*,htrackNhits_indeadzone_pAuJP2,htrackNhits_outdeadzone_pAuJP2*/};
  vector<TH1D*> tows_pAuJP2 = {htowerEt_pAuJP2, htowerEta_pAuJP2, htowerPhi_pAuJP2, htowerId_pAuJP2/*,htowerEt_indeadzone_pAuJP2,htowerEt_outdeadzone_pAuJP2*/};

  vector<TH1D*> evts_ppJP2 = {hn_trks_ppJP2, hn_tows_ppJP2, hbbc_coinc_ppJP2, hevt_vtx_ppJP2, hvpdvz_ppJP2, hvzdiff_ppJP2, hrunId_ppJP2, hn_vertices_ppJP2, hn_globals_ppJP2};   
  vector<TH1D*> trks_ppJP2 = {htrackPt_ppJP2, htrackEta_ppJP2, htrackPhi_ppJP2, htrackDCA_ppJP2,htrackNhits_ppJP2,htrackNhitsposs_ppJP2,htrackNhitsratio_ppJP2,htrackPt_indeadzone_ppJP2,htrackPt_outdeadzone_ppJP2,htrackDCA_indeadzone_ppJP2,htrackDCA_outdeadzone_ppJP2/*,htrackNhits_indeadzone_ppJP2,htrackNhits_outdeadzone_ppJP2*/};
  vector<TH1D*> tows_ppJP2 = {htowerEt_ppJP2, htowerEta_ppJP2, htowerPhi_ppJP2, htowerId_ppJP2/*,htowerEt_indeadzone_ppJP2,htowerEt_outdeadzone_ppJP2*/};

  vector<TH2D*> evts2D_pAuJP2 = {hbbc_coinc_evt_vtx_pAuJP2,hbbc_coinc_n_trks_pAuJP2,hbbc_coinc_n_tows_pAuJP2,hbbc_coinc_n_vertices_pAuJP2,hbbc_coinc_n_globals_pAuJP2,hrunId_evt_vtx_pAuJP2,hrunId_bbc_coinc_pAuJP2};
  vector<TH2D*> evts2D_ppJP2 = {hbbc_coinc_evt_vtx_ppJP2,hbbc_coinc_n_trks_ppJP2,hbbc_coinc_n_tows_ppJP2,hbbc_coinc_n_vertices_ppJP2,hbbc_coinc_n_globals_ppJP2,hrunId_evt_vtx_ppJP2,hrunId_bbc_coinc_ppJP2};

  vector<TH2D*> trks2D_pAuJP2 = {htrackEta_Phi_pAuJP2, htrackDCA_n_trks_pAuJP2,hbbc_coinc_trackDCA_pAuJP2,hbbc_coinc_trackPt_pAuJP2,hrunId_trackPt_pAuJP2,hrunId_trackDCA_pAuJP2};
  vector<TH2D*> tows2D_pAuJP2 = {htowerEta_Phi_pAuJP2,hbbc_coinc_towerEt_pAuJP2,hrunId_towerEt_pAuJP2,hrunId_towerId_pAuJP2,hrunId_towerEt_pAuJP2_badzone};
  
  vector<TH2D*> trks2D_ppJP2 = {htrackEta_Phi_ppJP2, htrackDCA_n_trks_ppJP2,hbbc_coinc_trackDCA_ppJP2,hbbc_coinc_trackPt_ppJP2,hrunId_trackPt_ppJP2,hrunId_trackDCA_ppJP2};
  vector<TH2D*> tows2D_ppJP2 = {htowerEta_Phi_ppJP2,hbbc_coinc_towerEt_ppJP2,hrunId_towerEt_ppJP2,hrunId_towerId_ppJP2,hrunId_towerEt_ppJP2_dummy};
  
  vector<TH3D*> trks3D_pAuJP2 = {hrunId_trackEta_Phi_pAuJP2, htrackPt_Eta_Phi_pAuJP2};
  vector<TH3D*> trks3D_ppJP2 = {hrunId_trackEta_Phi_ppJP2, htrackPt_Eta_Phi_ppJP2};

  vector<TH3D*> tows3D_pAuJP2 = {hrunId_towerEta_Phi_pAuJP2, htowerEt_Eta_Phi_pAuJP2};
  vector<TH3D*> tows3D_ppJP2 = {hrunId_towerEta_Phi_ppJP2, htowerEt_Eta_Phi_ppJP2};
  
  const int n_tows = 4800;
  double ids_pp[n_tows], ids_pAu[n_tows], ids_pAuBBCMB[n_tows];
  double meanEt_pp[n_tows] = {0}; double meanEt_pAu[n_tows] = {0}; double meanEt_pAuBBCMB[n_tows];
  double meanEtg2GeV_pp[n_tows] = {0}; double meanEtg2GeV_pAu[n_tows] = {0}; double meanEtg2GeV_pAuBBCMB[n_tows] = {0};
  
  for (int i = 0; i < n_tows; ++ i) {
    ids_pp[i] = i+1; ids_pAu[i] = i+1; ids_pAuBBCMB[i] = i+1;
  }
  //calling functions!
  TreetoHist (fpAuBBCMB, "pAu", "BBCMB", evts_pAuBBCMB, trks_pAuBBCMB, tows_pAuBBCMB, evts2D_pAuBBCMB, trks2D_pAuBBCMB, tows2D_pAuBBCMB, trks3D_pAuBBCMB,  tows3D_pAuBBCMB, meanEt_pAuBBCMB, meanEtg2GeV_pAuBBCMB);
  
  //temporarily changing the tree string to BBCMB because I messed up in naming it. Change back later after running again!
  TreetoHist (fpAu, "pAu", /*"JP2"*/"BBCMB", evts_pAuJP2, trks_pAuJP2, tows_pAuJP2, evts2D_pAuJP2, trks2D_pAuJP2, tows2D_pAuJP2, trks3D_pAuJP2,  tows3D_pAuJP2, meanEt_pAu, meanEtg2GeV_pAu);

  TreetoHist (fpp, "pp", "JP2", evts_ppJP2, trks_ppJP2, tows_ppJP2, evts2D_ppJP2, trks2D_ppJP2, tows2D_ppJP2, trks3D_ppJP2, tows3D_ppJP2, meanEt_pp, meanEtg2GeV_pp);
  
  TGraph* htowerId_meanEt_pAuJP2 = new TGraph(n_tows,ids_pAu,meanEt_pAu); htowerId_meanEt_pAuJP2->SetName("htowerId_meanEt_pAuJP2");
  TGraph* htowerId_meanEt_pAuBBCMB = new TGraph(n_tows,ids_pAuBBCMB,meanEt_pAuBBCMB); htowerId_meanEt_pAuBBCMB->SetName("htowerId_meanEt_pAuBBCMB");
  TGraph* htowerId_meanEt_ppJP2 = new TGraph(n_tows,ids_pp,meanEt_pp); htowerId_meanEt_ppJP2->SetName("htowerId_meanEt_ppJP2");
  TGraph* htowerId_meanEtg2GeV_pAuJP2 = new TGraph(n_tows,ids_pAu,meanEtg2GeV_pAu); htowerId_meanEtg2GeV_pAuJP2->SetName("htowerId_meanEtg2GeV_pAuJP2");
  TGraph* htowerId_meanEtg2GeV_pAuBBCMB = new TGraph(n_tows,ids_pAuBBCMB,meanEtg2GeV_pAuBBCMB); htowerId_meanEtg2GeV_pAuBBCMB->SetName("htowerId_meanEtg2GeV_pAuBBCMB");
  TGraph* htowerId_meanEtg2GeV_ppJP2 = new TGraph(n_tows,ids_pp,meanEtg2GeV_pp); htowerId_meanEtg2GeV_ppJP2->SetName("htowerId_meanEtg2GeV_ppJP2");


  TFile *fout = new TFile ("~/jetmass/macros/QA_pp_pAu_w_MB_and_vpdvzcut_3.root","RECREATE");
  
  htowerId_meanEt_pAuJP2->Write();
  htowerId_meanEt_pAuBBCMB->Write();
  htowerId_meanEt_ppJP2->Write();
  htowerId_meanEtg2GeV_pAuJP2->Write();
  htowerId_meanEtg2GeV_pAuBBCMB->Write();
  htowerId_meanEtg2GeV_ppJP2->Write();

  for (int i = 0; i < evts_pAuJP2.size(); ++ i) {
    evts_pAuJP2[i]->Write(); evts_pAuBBCMB[i]->Write();
    evts_ppJP2[i]->Write();
  }
  
  for (int i = 0; i < trks_pAuJP2.size(); ++ i) {
    trks_pAuJP2[i]->Write(); trks_pAuBBCMB[i]->Write();
    trks_ppJP2[i]->Write();
  }
  
  for (int i = 0; i < tows_pAuJP2.size(); ++ i) {
    tows_pAuJP2[i]->Write(); tows_pAuBBCMB[i]->Write();
    tows_ppJP2[i]->Write();
  }

  for (int i = 0; i < evts2D_pAuJP2.size(); ++ i) {
    evts2D_pAuJP2[i]->Write(); evts2D_pAuBBCMB[i]->Write();
    evts2D_ppJP2[i]->Write();
  }
  
  for (int i = 0; i < trks2D_pAuJP2.size(); ++ i) {
    trks2D_pAuJP2[i]->Write(); trks2D_pAuBBCMB[i]->Write();
    trks2D_ppJP2[i]->Write();
  }

  for (int i = 0; i < tows2D_pAuJP2.size(); ++ i) {
    tows2D_pAuJP2[i]->Write(); tows2D_pAuBBCMB[i]->Write();
    tows2D_ppJP2[i]->Write();
  }

  for (int i = 0; i < trks3D_pAuJP2.size(); ++ i) {
    trks3D_pAuJP2[i]->Write(); trks3D_pAuBBCMB[i]->Write();
    trks3D_ppJP2[i]->Write();
  }

  for (int i = 0; i < tows3D_pAuJP2.size(); ++ i) {
    tows3D_pAuJP2[i]->Write(); tows3D_pAuBBCMB[i]->Write();
    tows3D_ppJP2[i]->Write();
  }
  
  fout->Close();
  
  return 0;
}
