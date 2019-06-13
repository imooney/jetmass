#include <ctime>

using namespace std;

void TreetoHist (TFile *f, string system, string trig, vector<TH1D*> evts, vector<TH1D*> trks, vector<TH1D*> tows, vector<TH2D*> trks2D, vector<TH2D*> tows2D, double *meanEt, double *meanEtg2GeV) {
  //these things will be used for the average Et/towerID determination

  const int nTows = 4800;
  int counts[nTows] = {0}; //initially none of the towers have fired
  int countsg2GeV[nTows] = {0};

  double evt_vtx = -9999; double bbc_coinc = -9999;
  double n_trks = -9999; double n_tows = -9999;
  vector<double> * trackPt = 0; vector<double> * trackEta = 0; vector<double> * trackPhi = 0; vector<double> * trackDCA = 0;
  vector<double> * towerEt = 0; vector<double> * towerEta = 0; vector<double> * towerPhi = 0; vector<double> * towerId = 0;
  
  TTree *t = (TTree*) f->Get(("QA"+trig).c_str()); //string = BBCMB, VPDMB, or JP2 for now  
  t->SetBranchAddress("n_trks", &n_trks);
  t->SetBranchAddress("n_tows", &n_tows);
  t->SetBranchAddress("bbc_coinc", &bbc_coinc);
  t->SetBranchAddress("evt_vtx", &evt_vtx); //vpdvz                                                                                                
  t->SetBranchAddress("trackPt", &trackPt);
  t->SetBranchAddress("trackEta", &trackEta);
  t->SetBranchAddress("trackPhi", &trackPhi);
  t->SetBranchAddress("trackDCA", &trackDCA);
  t->SetBranchAddress("towerEta", &towerEta);
  t->SetBranchAddress("towerPhi", &towerPhi);
  t->SetBranchAddress("towerEt", &towerEt);
  t->SetBranchAddress("towerId", &towerId);

  cout << ("RUNNING OVER "+system+trig+"! Entries: ").c_str() << t->GetEntries() << endl;
  const clock_t begin_time = clock();
  for (int i = 0; i < t->GetEntries(); ++ i) {
    if (i % 100000 == 0) { cout << "still chuggin. " << i << endl; std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;}
    t->GetEntry(i);

    evts[0]->Fill(n_trks);
    evts[1]->Fill(n_tows);
    evts[2]->Fill(bbc_coinc);
    // cout<<"evt_vtx = "<<evt_vtx<<endl;
    evts[3]->Fill(evt_vtx);
        
    for (int j = 0; j < trackPt->size(); ++ j) { //all vectors of doubles in the branches should have the same size     
      trks[0]->Fill(trackPt->at(j));
      trks[1]->Fill(trackEta->at(j));
      trks[2]->Fill(trackPhi->at(j));
      trks[3]->Fill(trackDCA->at(j));
      trks2D[0]->Fill(trackEta->at(j),trackPhi->at(j));
      trks2D[1]->Fill(trackDCA->at(j),n_trks);
    }//! track loop 
    for (int j = 0; j < towerEt->size(); ++ j) {
      tows[0]->Fill(towerEt->at(j));
      tows[1]->Fill(towerEta->at(j));
      tows[2]->Fill(towerPhi->at(j));
      tows[3]->Fill(towerId->at(j));
      tows2D[0]->Fill(towerEta->at(j),towerPhi->at(j));
      
      counts[((int) towerId->at(j))-1] ++;
      meanEt[((int) towerId->at(j))-1] += towerEt->at(j);
      if (towerEt->at(j) > 2) {
	countsg2GeV[((int) towerId->at(j))-1] ++;
	meanEtg2GeV[((int) towerId->at(j))-1] += towerEt->at(j);
      }
    }//! tower size loop
    
  }//! event loop
  
  for (int i = 0; i < nTows; ++ i) {
    meanEt[i] /= (double) counts[i];
    meanEtg2GeV[i] /= (double) countsg2GeV[i];
  }
  
  //! needs to be outside the event loop 
  t->ResetBranchAddresses();
  
  return;
}

void tohist () {
  TFile *fpAu = new TFile ("~/jetmass/out/QA/full_pAu.root","READ");
  TFile *fpp = new TFile ("~/jetmass/out/QA/full_pp.root","READ");
  
  TH1D* hn_trks_pAuJP2 = new TH1D("hn_trks_pAuJP2","hn_trks_pAuJP2",250,0,500);
  TH1D* hn_tows_pAuJP2 = new TH1D("hn_tows_pAuJP2","hn_tows_pAuJP2",150,0,150);
  TH1D* hbbc_coinc_pAuJP2 = new TH1D("hbbc_coinc_pAuJP2","hbbc_coinc_pAuJP2",200,0,2600000);
  TH1D* hevt_vtx_pAuJP2 = new TH1D("hevt_vtx_pAuJP2","hevt_vtx_pAuJP2",100,-35,35);
  
  TH1D* hn_trks_ppJP2 = new TH1D("hn_trks_ppJP2","hn_trks_ppJP2",250,0,500);
  TH1D* hn_tows_ppJP2 = new TH1D("hn_tows_ppJP2","hn_tows_ppJP2",150,0,150);
  TH1D* hbbc_coinc_ppJP2 = new TH1D("hbbc_coinc_ppJP2","hbbc_coinc_ppJP2",200,0,2600000);
  TH1D* hevt_vtx_ppJP2 = new TH1D("hevt_vtx_ppJP2","hevt_vtx_ppJP2",100,-35,35);
 
  TH1D* htrackPt_pAuJP2 = new TH1D("htrackPt_pAuJP2","htrackPt_pAuJP2",200,0,40);
  TH1D* htrackEta_pAuJP2 = new TH1D("htrackEta_pAuJP2","htrackEta_pAuJP2",200,-2,2);
  TH1D* htrackPhi_pAuJP2 = new TH1D("htrackPhi_pAuJP2","htrackPhi_pAuJP2",200,-3.15,3.15);
  TH1D* htrackDCA_pAuJP2 = new TH1D("htrackDCA_pAuJP2","htrackDCA_pAuJP2",200,0,1);
  
  TH1D* htrackPt_ppJP2 = new TH1D("htrackPt_ppJP2","htrackPt_ppJP2",200,0,40);
  TH1D* htrackEta_ppJP2 = new TH1D("htrackEta_ppJP2","htrackEta_ppJP2",200,-2,2);
  TH1D* htrackPhi_ppJP2 = new TH1D("htrackPhi_ppJP2","htrackPhi_ppJP2",200,-3.15,3.15);
  TH1D* htrackDCA_ppJP2 = new TH1D("htrackDCA_ppJP2","htrackDCA_ppJP2",200,0,1);

  TH1D* htowerEt_pAuJP2 = new TH1D("htowerEt_pAuJP2","htowerEt_pAuJP2",200,0,40);
  TH1D* htowerEta_pAuJP2 = new TH1D("htowerEta_pAuJP2","htowerEta_pAuJP2",200,-2,2);
  TH1D* htowerPhi_pAuJP2 = new TH1D("htowerPhi_pAuJP2","htowerPhi_pAuJP2",200,-3.15,3.15);
  TH1D* htowerId_pAuJP2 = new TH1D("htowerId_pAuJP2","htowerId_pAuJP2",4801,0,4801);

  TH1D* htowerEt_ppJP2 = new TH1D("htowerEt_ppJP2","htowerEt_ppJP2",200,0,40);
  TH1D* htowerEta_ppJP2 = new TH1D("htowerEta_ppJP2","htowerEta_ppJP2",200,-2,2);
  TH1D* htowerPhi_ppJP2 = new TH1D("htowerPhi_ppJP2","htowerPhi_ppJP2",200,-3.15,3.15);
  TH1D* htowerId_ppJP2 = new TH1D("htowerId_ppJP2","htowerId_ppJP2",4801,0,4801);

  TH2D* htrackEta_Phi_pAuJP2 = new TH2D("htrackEta_Phi_pAuJP2","pAu JP2;#eta_{trk};#phi_{trk}",100,-1,1,200,-3.15,3.15);
  TH2D* htrackDCA_n_trks_pAuJP2 = new TH2D("htrackDCA_n_trks_pAuJP2","pAu JP2;track DCA;N_{trk}",200,0,1,250,0,500);
  
  TH2D* htrackEta_Phi_ppJP2 = new TH2D("htrackEta_Phi_ppJP2","pp JP2;#eta_{trk};#phi_{trk}",100,-1,1,200,-3.15,3.15);
  TH2D* htrackDCA_n_trks_ppJP2 = new TH2D("htrackDCA_n_trks_ppJP2","pp JP2;track DCA;N_{trk}",200,0,1,250,0,500);
  
  TH2D* htowerEta_Phi_pAuJP2 = new TH2D("htowerEta_Phi_pAuJP2","pAu JP2;#eta_{tow};#phi_{tow}",100,-1,1,200,-3.15,3.15);
  
  TH2D* htowerEta_Phi_ppJP2 = new TH2D("htowerEta_Phi_ppJP2","pp JP2;#eta_{tow};#phi_{tow}",100,-1,1,200,-3.15,3.15);

  vector<TH1D*> evts_pAuJP2 = {hn_trks_pAuJP2, hn_tows_pAuJP2, hbbc_coinc_pAuJP2, hevt_vtx_pAuJP2};
  vector<TH1D*> trks_pAuJP2 = {htrackPt_pAuJP2, htrackEta_pAuJP2, htrackPhi_pAuJP2, htrackDCA_pAuJP2};
  vector<TH1D*> tows_pAuJP2 = {htowerEt_pAuJP2, htowerEta_pAuJP2, htowerPhi_pAuJP2, htowerId_pAuJP2};

  vector<TH1D*> evts_ppJP2 = {hn_trks_ppJP2, hn_tows_ppJP2, hbbc_coinc_ppJP2, hevt_vtx_ppJP2};
  vector<TH1D*> trks_ppJP2 = {htrackPt_ppJP2, htrackEta_ppJP2, htrackPhi_ppJP2, htrackDCA_ppJP2};
  vector<TH1D*> tows_ppJP2 = {htowerEt_ppJP2, htowerEta_ppJP2, htowerPhi_ppJP2, htowerId_ppJP2};

  vector<TH2D*> trks2D_pAuJP2 = {htrackEta_Phi_pAuJP2, htrackDCA_n_trks_pAuJP2};
  vector<TH2D*> tows2D_pAuJP2 = {htowerEta_Phi_pAuJP2};
  
  vector<TH2D*> trks2D_ppJP2 = {htrackEta_Phi_ppJP2, htrackDCA_n_trks_ppJP2};
  vector<TH2D*> tows2D_ppJP2 = {htowerEta_Phi_ppJP2};
  
  const int n_tows = 4800;
  double ids_pp[n_tows], ids_pAu[n_tows];
  double meanEt_pp[n_tows] = {0}; double meanEt_pAu[n_tows] = {0};
  double meanEtg2GeV_pp[n_tows] = {0}; double meanEtg2GeV_pAu[n_tows] = {0};
  
  for (int i = 0; i < n_tows; ++ i) {
    ids_pp[i] = i+1; ids_pAu[i] = i+1;
  }
  
  TreetoHist (fpAu, "pAu", "JP2", evts_pAuJP2, trks_pAuJP2, tows_pAuJP2, trks2D_pAuJP2, tows2D_pAuJP2, meanEt_pAu, meanEtg2GeV_pAu);
  //  TreetoHist (fpAu, "pAu", "VPDMB", evts_pAuVPDMB, trks_pAuVPDMB, tows_pAuVPDMB);
  TreetoHist (fpp, "pp", "JP2", evts_ppJP2, trks_ppJP2, tows_ppJP2, trks2D_ppJP2, tows2D_ppJP2, meanEt_pp, meanEtg2GeV_pp);
  
  TGraph* htowerId_meanEt_pAuJP2 = new TGraph(n_tows,ids_pAu,meanEt_pAu);
  TGraph* htowerId_meanEt_ppJP2 = new TGraph(n_tows,ids_pp,meanEt_pp);
  TGraph* htowerId_meanEtg2GeV_pAuJP2 = new TGraph(n_tows,ids_pAu,meanEtg2GeV_pAu);
  TGraph* htowerId_meanEtg2GeV_ppJP2 = new TGraph(n_tows,ids_pp,meanEtg2GeV_pp);

  TFile *fout = new TFile ("~/jetmass/macros/QA_pp_pAu_test.root","RECREATE");
  
  htowerId_meanEt_pAuJP2->Write();
  htowerId_meanEt_ppJP2->Write();
  htowerId_meanEtg2GeV_pAuJP2->Write();
  htowerId_meanEtg2GeV_ppJP2->Write();

  for (int i = 0; i < evts_pAuJP2.size(); ++ i) {
    evts_pAuJP2[i]->Write();
    evts_ppJP2[i]->Write();
  }
  
  for (int i = 0; i < trks_pAuJP2.size(); ++ i) {
    trks_pAuJP2[i]->Write();
    trks_ppJP2[i]->Write();
  }
  
  for (int i = 0; i < tows_pAuJP2.size(); ++ i) {
    tows_pAuJP2[i]->Write();
    tows_ppJP2[i]->Write();
  }

  for (int i = 0; i < trks2D_pAuJP2.size(); ++ i) {
    trks2D_pAuJP2[i]->Write();
    trks2D_ppJP2[i]->Write();
  }

  for (int i = 0; i < tows2D_pAuJP2.size(); ++ i) {
    tows2D_pAuJP2[i]->Write();
    tows2D_ppJP2[i]->Write();
  }
  
  fout->Close();
  
  return;
  
}
