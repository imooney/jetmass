
using namespace std;

void tohist () {
  TFile *f = new TFile ("~/jetmass/out/MB/full.root","READ");

  TH1D* hn_trks_BBCMB = new TH1D("hn_trks_BBCMB","hn_trks_BBCMB",250,0,250);
  TH1D* hn_tows_BBCMB = new TH1D("hn_tows_BBCMB","hn_tows_BBCMB",150,0,150);
  TH1D* hbbc_coinc_BBCMB = new TH1D("hbbc_coinc_BBCMB","hbbc_coinc_BBCMB",1000,0,3000);
  TH1D* hevt_vtx_BBCMB = new TH1D("hevt_vtx_BBCMB","hevt_vtx_BBCMB",100,-35,35);
  
  TH1D* hn_trks_VPDMB = new TH1D("hn_trks_VPDMB","hn_trks_VPDMB",250,0,250);
  TH1D* hn_tows_VPDMB = new TH1D("hn_tows_VPDMB","hn_tows_VPDMB",150,0,150);
  TH1D* hbbc_coinc_VPDMB = new TH1D("hbbc_coinc_VPDMB","hbbc_coinc_VPDMB",1000,0,3000);
  TH1D* hevt_vtx_VPDMB = new TH1D("hevt_vtx_VPDMB","hevt_vtx_VPDMB",100,-35,35);
 
  TH1D* htrackPt_BBCMB = new TH1D("htrackPt_BBCMB","htrackPt_BBCMB",200,0,40);
  TH1D* htrackEta_BBCMB = new TH1D("htrackEta_BBCMB","htrackEta_BBCMB",200,-2,2);
  TH1D* htrackPhi_BBCMB = new TH1D("htrackPhi_BBCMB","htrackPhi_BBCMB",200,-3.15,3.15);
  TH1D* htrackDCA_BBCMB = new TH1D("htrackDCA_BBCMB","htrackDCA_BBCMB",50,0,5);

  TH1D* htrackPt_VPDMB = new TH1D("htrackPt_VPDMB","htrackPt_VPDMB",200,0,40);
  TH1D* htrackEta_VPDMB = new TH1D("htrackEta_VPDMB","htrackEta_VPDMB",200,-2,2);
  TH1D* htrackPhi_VPDMB = new TH1D("htrackPhi_VPDMB","htrackPhi_VPDMB",200,-3.15,3.15);
  TH1D* htrackDCA_VPDMB = new TH1D("htrackDCA_VPDMB","htrackDCA_VPDMB",50,0,5);

  TH1D* htowerEt_BBCMB = new TH1D("htowerEt_BBCMB","htowerEt_BBCMB",200,0,40);
  TH1D* htowerEta_BBCMB = new TH1D("htowerEta_BBCMB","htowerEta_BBCMB",200,-2,2);
  TH1D* htowerPhi_BBCMB = new TH1D("htowerPhi_BBCMB","htowerPhi_BBCMB",200,-3.15,3.15);
  TH1D* htowerId_BBCMB = new TH1D("htowerId_BBCMB","htowerId_BBCMB",4801,0,4801);

  TH1D* htowerEt_VPDMB = new TH1D("htowerEt_VPDMB","htowerEt_VPDMB",200,0,40);
  TH1D* htowerEta_VPDMB = new TH1D("htowerEta_VPDMB","htowerEta_VPDMB",200,-2,2);
  TH1D* htowerPhi_VPDMB = new TH1D("htowerPhi_VPDMB","htowerPhi_VPDMB",200,-3.15,3.15);
  TH1D* htowerId_VPDMB = new TH1D("htowerId_VPDMB","htowerId_VPDMB",4801,0,4801);
  
  
  double evt_vtx_BBCMB = -9999; double bbc_coinc_BBCMB = -9999;
  double n_trks_BBCMB = -9999; double n_tows_BBCMB = -9999;
  vector<double> * trackPt_BBCMB = 0; vector<double> * trackEta_BBCMB = 0; vector<double> * trackPhi_BBCMB = 0; vector<double> * trackDCA_BBCMB = 0;
  vector<double> * towerEt_BBCMB = 0; vector<double> * towerEta_BBCMB = 0; vector<double> * towerPhi_BBCMB = 0; vector<double> * towerId_BBCMB = 0;

  double evt_vtx_VPDMB = -9999; double bbc_coinc_VPDMB = -9999;
  double n_trks_VPDMB = -9999; double n_tows_VPDMB = -9999;
  vector<double> * trackPt_VPDMB = 0; vector<double> * trackEta_VPDMB = 0; vector<double> * trackPhi_VPDMB = 0; vector<double> * trackDCA_VPDMB = 0;
  vector<double> * towerEt_VPDMB = 0; vector<double> * towerEta_VPDMB = 0; vector<double> * towerPhi_VPDMB = 0; vector<double> * towerId_VPDMB = 0;

  TTree *tVPD = (TTree*) f->Get("QAVPDMB");  
  tVPD->SetBranchAddress("n_trks", &n_trks_VPDMB);
  tVPD->SetBranchAddress("n_tows", &n_tows_VPDMB);
  tVPD->SetBranchAddress("bbc_coinc", &bbc_coinc_VPDMB);
  tVPD->SetBranchAddress("evt_vtx", &evt_vtx_VPDMB); //vpdvz                                                                                                

  tVPD->SetBranchAddress("trackPt", &trackPt_VPDMB);
  tVPD->SetBranchAddress("trackEta", &trackEta_VPDMB);
  tVPD->SetBranchAddress("trackPhi", &trackPhi_VPDMB);
  tVPD->SetBranchAddress("trackDCA", &trackDCA_VPDMB);
  tVPD->SetBranchAddress("towerEta", &towerEta_VPDMB);
  tVPD->SetBranchAddress("towerPhi", &towerPhi_VPDMB);
  tVPD->SetBranchAddress("towerEt", &towerEt_VPDMB);
  tVPD->SetBranchAddress("towerId", &towerId_VPDMB);


  TTree *tBBC = (TTree*) f->Get("QABBCMB");
  tBBC->SetBranchAddress("n_trks", &n_trks_BBCMB);
  tBBC->SetBranchAddress("n_tows", &n_tows_BBCMB);
  tBBC->SetBranchAddress("bbc_coinc", &bbc_coinc_BBCMB);
  tBBC->SetBranchAddress("evt_vtx", &evt_vtx_BBCMB); //vpdvz                                                                                                      
  tBBC->SetBranchAddress("trackPt", &trackPt_BBCMB);
  tBBC->SetBranchAddress("trackEta", &trackEta_BBCMB);
  tBBC->SetBranchAddress("trackPhi", &trackPhi_BBCMB);
  tBBC->SetBranchAddress("trackDCA", &trackDCA_BBCMB);
  tBBC->SetBranchAddress("towerEta", &towerEta_BBCMB);
  tBBC->SetBranchAddress("towerPhi", &towerPhi_BBCMB);
  tBBC->SetBranchAddress("towerEt", &towerEt_BBCMB);
  tBBC->SetBranchAddress("towerId", &towerId_BBCMB);

  cout << "RUNNING OVER BBCMB! Entries: " << tBBC->GetEntries() << endl;
  for (int i = 0; i < tBBC->GetEntries(); ++ i) {
    if (i % 1000000 == 0) { cout << "still chuggin. " << i << endl;}
    tBBC->GetEntry(i);

    hn_trks_BBCMB->Fill(n_trks_BBCMB);
    hn_tows_BBCMB->Fill(n_tows_BBCMB);
    hbbc_coinc_BBCMB->Fill(bbc_coinc_BBCMB);
    // cout<<"evt_vtx_BBCMB = "<<evt_vtx_BBCMB<<endl;
    hevt_vtx_BBCMB->Fill(evt_vtx_BBCMB);
    
    for (int j = 0; j < trackPt_BBCMB->size(); ++ j) { //all vectors of doubles in the branches should have the same size     
      htrackPt_BBCMB->Fill(trackPt_BBCMB->at(j));
      htrackEta_BBCMB->Fill(trackEta_BBCMB->at(j));
      htrackPhi_BBCMB->Fill(trackPhi_BBCMB->at(j));
      htrackDCA_BBCMB->Fill(trackDCA_BBCMB->at(j));
    }//! track loop 
    for (int j = 0; j < towerEt_BBCMB->size(); ++ j) {
      htowerEt_BBCMB->Fill(towerEt_BBCMB->at(j));
      htowerEta_BBCMB->Fill(towerEta_BBCMB->at(j));
      htowerPhi_BBCMB->Fill(towerPhi_BBCMB->at(j));
      htowerId_BBCMB->Fill(towerId_BBCMB->at(j));
    }//! tower size loop 
  }//! event loop

  //! needs to be outside the event loop 
  //tBBC->ResetBranchAddresses();
  
  cout << "RUNNING OVER VPDMB! Entries: " << tVPD->GetEntries() << endl;
  for (int i = 0; i < tVPD->GetEntries(); ++ i) {
    if (i % 1000000 == 0) { cout << "still chuggin. " << i << endl;}
    tVPD->GetEntry(i);
    
    hn_trks_VPDMB->Fill(n_trks_VPDMB);
    hn_tows_VPDMB->Fill(n_tows_VPDMB);
    hbbc_coinc_VPDMB->Fill(bbc_coinc_VPDMB);
    hevt_vtx_VPDMB->Fill(evt_vtx_VPDMB);
    
    for (int j = 0; j < trackPt_VPDMB->size(); ++ j) { //all vectors of doubles in the branches should have the same size     
      htrackPt_VPDMB->Fill(trackPt_VPDMB->at(j));
      htrackEta_VPDMB->Fill(trackEta_VPDMB->at(j));
      htrackPhi_VPDMB->Fill(trackPhi_VPDMB->at(j));
      htrackDCA_VPDMB->Fill(trackDCA_VPDMB->at(j));
    }
    for (int j = 0; j < towerEt_VPDMB->size(); ++ j) {
      htowerEt_VPDMB->Fill(towerEt_VPDMB->at(j));
      htowerEta_VPDMB->Fill(towerEta_VPDMB->at(j));
      htowerPhi_VPDMB->Fill(towerPhi_VPDMB->at(j));
      htowerId_VPDMB->Fill(towerId_VPDMB->at(j));
    }
  }
  
  tVPD->ResetBranchAddresses();
 
  TFile *fout = new TFile ("~/jetmass/macros/pAu_hists_test.root","RECREATE");
  
  hn_trks_BBCMB->Write();
  hn_tows_BBCMB->Write();
  hbbc_coinc_BBCMB->Write();
  hevt_vtx_BBCMB->Write();
  
  hn_trks_VPDMB->Write();
  hn_tows_VPDMB->Write();
  hbbc_coinc_VPDMB->Write();
  hevt_vtx_VPDMB->Write();
 
  htrackPt_BBCMB->Write();
  htrackEta_BBCMB->Write();
  htrackPhi_BBCMB->Write();
  htrackDCA_BBCMB->Write();

  htrackPt_VPDMB->Write();
  htrackEta_VPDMB->Write();
  htrackPhi_VPDMB->Write();
  htrackDCA_VPDMB->Write();

  htowerEt_BBCMB->Write();
  htowerEta_BBCMB->Write();
  htowerPhi_BBCMB->Write();
  htowerId_BBCMB->Write();

  htowerEt_VPDMB->Write();
  htowerEta_VPDMB->Write();
  htowerPhi_VPDMB->Write();
  htowerId_VPDMB->Write();

  fout->Close();
  
  return;
  
}
