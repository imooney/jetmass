
#include "Plots.h"

using namespace std;

void p6p8compare () {
  
  const string p6path = "~/jetmass/macros/hists/";
  const string p6decpath = "~/jetmass/production/macros/p6_dec_hists/";
  const string p8path = "~/jetmass/production/macros/hists/";
  const string herpath = "~/jetmass/production/macros/herwig_hists/";
  const string file = "hists.root";
  const string out = "~/jetmass/production/plots/p6p8compare/";
  const string filetype = ".pdf";
  
  TFile * p8 = new TFile((p8path + file).c_str(),"READ");
  TFile * p6 = new TFile((p6path + file).c_str(),"READ");
  TFile * p6dec = new TFile((p6decpath + file).c_str(),"READ");
  TFile * her = new TFile((herpath+file).c_str(),"READ");
  
  TH2D* m_pt_un = (TH2D*) p8->Get("m_v_pt_un");
  TH2D* m_pt_dec = (TH2D*) p8->Get("m_v_pt_dec");
  TH2D* m_pt_6 = (TH2D*) p6->Get("m_v_pt_p");
  
  TH2D* conspT8 = (TH2D*) p8->Get("conspT_v_ptp8");
  TH2D* conspT8dec = (TH2D*) p8->Get("conspT_v_ptp8dec");
  TH2D* conspT6 = (TH2D*) p6->Get("conspT_v_pt");
  TH2D* consGirth8 = (TH2D*) p8->Get("consGirth_v_ptp8");
  TH2D* consGirth8dec = (TH2D*) p8->Get("consGirth_v_ptp8dec");
  TH2D* consGirth6 = (TH2D*) p6->Get("consGirth_v_pt");
  TH2D* consDist8 = (TH2D*) p8->Get("consDist_v_ptp8");
  TH2D* consDist8dec = (TH2D*) p8->Get("consDist_v_ptp8dec");
  TH2D* consDist6 = (TH2D*) p6->Get("consDist_v_pt");
  TH2D* jetMult8 = (TH2D*) p8->Get("jetMult_v_ptp8");
  TH2D* jetMult6 = (TH2D*) p6->Get("jetMult_v_pt");
  TH2D* jetMult6dec = (TH2D*) p6dec->Get("jetMult_v_ptp6dec");
  TH2D* jetMult8dec = (TH2D*) p8->Get("jetMult_v_ptp8dec");
  TH2D* jetGirth8 = (TH2D*) p8->Get("jetGirth_v_ptp8");
  TH2D* jetGirth8dec = (TH2D*) p8->Get("jetGirth_v_ptp8dec");
  TH2D* jetGirth6 = (TH2D*) p6->Get("jetGirth_v_pt");
  
  
  TH2D* jetEta8dec = (TH2D*) p8->Get("eta_dec_v_pt");
  TH2D* jetEta8 = (TH2D*) p8->Get("eta_un_v_pt");
  TH2D* jetEta6 = (TH2D*) p6->Get("eta_v_pt_p");
  TH2D* jetEta6Ge = (TH2D*) p6->Get("eta_v_pt_g");
  TH2D* jetEtaDat = (TH2D*) p6->Get("eta_v_pt_d");
  TH2D* jetEtaHer = (TH2D*) her->Get("eta_her_v_pt");
  
  TH2D* jetPhi8dec = (TH2D*) p8->Get("phi_dec_v_pt");
  TH2D* jetPhi8 = (TH2D*) p8->Get("phi_un_v_pt");
  TH2D* jetPhi6 = (TH2D*) p6->Get("phi_v_pt_p");
  TH2D* jetPhi6Ge = (TH2D*) p6->Get("phi_v_pt_g");
  TH2D* jetPhiDat = (TH2D*) p6->Get("phi_v_pt_d");
  TH2D* jetPhiHer = (TH2D*) her->Get("phi_her_v_pt");

  TH1D* jetPt6 = (TH1D*) p6->Get("pt_p_var_bin");
  TH1D* jetPt8 = (TH1D*) p8->Get("pt_coarse_un");
  TH1D* jetPt8dec = (TH1D*) p8->Get("pt_coarse_dec");
  TH1D* jetPt6Ge = (TH1D*) p6->Get("pt_g_var_bin");
  TH1D* jetPtDat = (TH1D*) p6->Get("pt_d_var_bin");
  TH1D* jetPtHer = (TH1D*) her->Get("pt_coarse_her");
  
  vector<TH1D*> unvec; vector<TH1D*> decvec; vector<TH1D*> p6vec;
  
  vector<TH1D*> p8conspT; vector<TH1D*> p8decconspT; vector<TH1D*> p6conspT;
  vector<TH1D*> p8consGirth; vector<TH1D*> p8decconsGirth; vector<TH1D*> p6consGirth;
  vector<TH1D*> p8consDist; vector<TH1D*> p8decconsDist; vector<TH1D*> p6consDist;
  vector<TH1D*> p8jetMult; vector<TH1D*> p6jetMult; vector<TH1D*> p6decjetMult; vector<TH1D*> p8decjetMult;
  vector<TH1D*> p8jetGirth; vector<TH1D*> p8decjetGirth; vector<TH1D*> p6jetGirth;
  vector<TH1D*> p8jetEta; vector<TH1D*> p6jetEta; vector<TH1D*> p8decjetEta; vector<TH1D*> p6GejetEta; vector<TH1D*> datjetEta; vector<TH1D*> herjetEta;
  vector<TH1D*> p8jetPhi; vector<TH1D*> p6jetPhi; vector<TH1D*> p8decjetPhi; vector<TH1D*> p6GejetPhi; vector<TH1D*> datjetPhi; vector<TH1D*> herjetPhi;
  
  const int nBins = 5;
  double ranges[nBins+1] = {3,4,5,6,8,12};
  double ranges_d[nBins+1] = {1,2,3,4,6,10};
  std::string pts[nBins+1] = {"15","20","25","30","40","60"};
  
  unvec = Projection2D (m_pt_un, nBins, ranges, "x");
  decvec = Projection2D (m_pt_dec, nBins, ranges, "x");
  p6vec = Projection2D (m_pt_6, nBins, ranges, "x");
  
  p8conspT = Projection2D(conspT8, nBins, ranges, "x");
  p8decconspT = Projection2D(conspT8dec, nBins, ranges, "x");
  p6conspT = Projection2D(conspT6, nBins, ranges, "x");
  p8consGirth = Projection2D(consGirth8, nBins, ranges, "x");
  p8decconsGirth = Projection2D(consGirth8dec, nBins, ranges, "x");
  p6consGirth = Projection2D(consGirth6, nBins, ranges, "x");
  p8consDist = Projection2D(consDist8, nBins, ranges, "x");
  p8decconsDist = Projection2D(consDist8dec, nBins, ranges, "x");
  p6consDist = Projection2D(consDist6, nBins, ranges, "x");
  p8jetMult = Projection2D(jetMult8, nBins, ranges, "x");
  p6jetMult = Projection2D(jetMult6, nBins, ranges, "x");
  p6decjetMult = Projection2D(jetMult6dec, nBins, ranges, "x");
  p8decjetMult = Projection2D(jetMult8dec, nBins, ranges, "x");
  p8jetGirth = Projection2D(jetGirth8, nBins, ranges, "x");
  p8decjetGirth = Projection2D(jetGirth8dec, nBins, ranges, "x");
  p6jetGirth = Projection2D(jetGirth6, nBins, ranges, "x");
  
  p8jetPhi = Projection2D(jetPhi8, nBins, ranges, "x");
  p6jetPhi = Projection2D(jetPhi6, nBins, ranges, "x");
  p8decjetPhi = Projection2D(jetPhi8dec, nBins, ranges, "x");
  herjetPhi = Projection2D(jetPhiHer, nBins, ranges, "x");
  p6GejetPhi = Projection2D(jetPhi6Ge, nBins, ranges_d, "x");
  datjetPhi = Projection2D(jetPhiDat, nBins, ranges_d, "x");
  
  p8jetEta = Projection2D(jetEta8, nBins, ranges, "x");
  p6jetEta = Projection2D(jetEta6, nBins, ranges, "x");
  p8decjetEta = Projection2D(jetEta8dec, nBins, ranges, "x");
  herjetEta = Projection2D(jetEtaHer, nBins, ranges, "x");
  p6GejetEta = Projection2D(jetEta6Ge, nBins, ranges_d, "x");
  datjetEta = Projection2D(jetEtaDat, nBins, ranges_d, "x");
  
  
  vector<double> njetsp6; vector<double> njetsp8; vector<double> njetsp8dec;
  for (int i = 0; i < p8jetMult.size(); ++ i) {
    njetsp8.push_back(p8jetMult[i]->Integral());
    njetsp8dec.push_back(p8decjetMult[i]->Integral());
    njetsp6.push_back(p6jetMult[i]->Integral());
  }
  
  TH1D* jetPt8clone = (TH1D*) jetPt8->Clone("jetPt8clone");
  TH1D* jetPt8decclone = (TH1D*) jetPt8dec->Clone("jetPt8decclone");
  
  TH1D* jetPt8nostyle = (TH1D*) jetPt8clone->Clone("jetPt8nostyle");
  TH1D* jetPt8decnostyle = (TH1D*) jetPt8decclone->Clone("jetPt8decnostyle");
  TH1D* jetPt6nostyle = (TH1D*) jetPt6->Clone("jetPt6nostyle");
  TH1D* jetPtHernostyle = (TH1D*) jetPtHer->Clone("jetPtHernostyle");
  TH1D* jetPt6Genostyle = (TH1D*) jetPt6Ge->Clone("jetPt6Genostyle");
  TH1D* jetPtDatnostyle = (TH1D*) jetPtDat->Clone("jetPtDatnostyle");
  
  jetPt8nostyle->Scale(1/(double)jetPt8nostyle->Integral());
  jetPt8decnostyle->Scale(1/(double)jetPt8decnostyle->Integral());
  jetPt6nostyle->Scale(1/(double)jetPt6nostyle->Integral());
  jetPt6Genostyle->Scale(1/(double)jetPt6Genostyle->Integral());
  jetPtHernostyle->Scale(1/(double)jetPtHernostyle->Integral());
  jetPtDatnostyle->Scale(1/(double)jetPtDatnostyle->Integral());
  

  Prettify1DwLineStyle(jetPt6, kRed, kDashed, 5, "p_{T,jet} [GeV/c]","arb.",5,60,-1,-1);
  Prettify1DwLineStyle(jetPt8clone, kBlue, kDashed, 5, "p_{T,jet} [GeV/c]","arb.",5,60,-1,-1);
  Prettify1DwLineStyle(jetPt8decclone, kBlue, kSolid, 5, "p_{T,jet} [GeV/c]","arb.",5,60,-1,-1);  
  Prettify1DwLineStyle(jetPtHer, kCyan, kSolid, 5, "p_{T,jet} [GeV/c]","arb.",5,60,-1,-1);  
  Prettify1D(jetPt6Ge, kRed, kOpenCircle,2,kRed, "p_{T,jet} [GeV/c]","arb.",5,60,-1,-1);
  Prettify1D(jetPtDat, kBlack, kOpenStar, 2.5,kBlack, "p_{T,jet} [GeV/c]","arb.",5,60,-1,-1);
  

  double rat_scale = jetPt6->GetBinContent(6) / (double) jetPt8->GetBinContent(6);
  double rat_scale_dec = jetPt6->GetBinContent(6) / (double) jetPt8dec->GetBinContent(6);
  jetPt8->Scale(rat_scale);
  jetPt8dec->Scale(rat_scale_dec);
  
  for (int i = 0; i < nBins; ++ i) {
    Prettify1DwLineStyle (unvec[i], kBlue, kDashed, 5, "M [GeV/c^{2}]", "1/N_{jets} dN/dM", -1, -1, 0,0.4);
    Prettify1DwLineStyle (decvec[i], kBlue, kSolid, 5, "M [GeV/c^{2}]", "1/N_{jets} dN/dM", -1, -1, 0,0.4);
    Prettify1DwLineStyle (p6vec[i], kRed, kDashed, 5, "M [GeV/c^{2}]", "1/N_{jets} dN/dM", -1, -1, 0,0.4);

    p8consGirth[i]->Scale(1/(double)p8consGirth[i]->Integral());//njetsp8[i]);
    p8decconsGirth[i]->Scale(1/(double)p8decconsGirth[i]->Integral());//njetsp8dec[i]);
    p6consGirth[i]->Scale(1/(double)p6consGirth[i]->Integral());//njetsp6[i]);
    p8consDist[i]->Scale(1/(double)p8consDist[i]->Integral());//(double)njetsp8[i]);
    p8decconsDist[i]->Scale(1/(double)p8decconsDist[i]->Integral());//(double)njetsp8dec[i]);
    p6consDist[i]->Scale(1/(double)p6consDist[i]->Integral());//(double)njetsp6[i]);
    p8conspT[i]->Scale(1/(double)p8conspT[i]->Integral());//njetsp8[i]);
    p8decconspT[i]->Scale(1/(double)p8decconspT[i]->Integral());//njetsp8dec[i]);
    p6conspT[i]->Scale(1/(double)p6conspT[i]->Integral());//njetsp6[i]);

    Prettify1DwLineStyle (p8conspT[i], kBlue, kDashed, 5, "p_{T,cons}","1/N_{cons} dN_{cons}/dp_{T}",-1,-1,1e-5,5e-1);
    Prettify1DwLineStyle (p8decconspT[i], kBlue, kSolid, 5, "p_{T,cons}","1/N_{cons} dN_{cons}/dp_{T}",-1,-1,1e-5,5e-1);
    Prettify1DwLineStyle (p6conspT[i], kRed, kDashed, 5, "p_{T,cons}","1/N_{cons} dN_{cons}/dp_{T}",-1,-1,1e-5,5e-1);
    Prettify1DwLineStyle (p8consGirth[i], kBlue, kDashed, 5, "|r_{i}|p_{T,i}/p_{T,jet}","1/N_{cons} dN_{cons}/dg_{i}",0,0.06,0,65);
    Prettify1DwLineStyle (p8decconsGirth[i], kBlue, kSolid, 5, "|r_{i}|p_{T,i}/p_{T,jet}","1/N_{cons} dN_{cons}/dg_{i}",0,0.06,0,65);
    Prettify1DwLineStyle (p6consGirth[i], kRed, kDashed, 5, "|r_{i}|p_{T,i}/p_{T,jet}","1/N_{cons} dN_{cons}/dg_{i}",0,0.06,0,65);
    Prettify1DwLineStyle (p8consDist[i], kBlue, kDashed, 5, "|r_{i}|","1/N_{cons} dN_{cons}/dr_{i}",-1,-1,0,6);
    Prettify1DwLineStyle (p8decconsDist[i], kBlue, kSolid, 5, "|r_{i}|","1/N_{cons} dN_{cons}/dr_{i}",-1,-1,0,6);
    Prettify1DwLineStyle (p6consDist[i], kRed, kDashed, 5, "|r_{i}|","1/N_{cons} dN_{cons}/dr_{i}",-1,-1,0,6);
    Prettify1DwLineStyle (p8jetMult[i], kBlue, kDashed, 5, "N_{c}","1/N_{j} dN_{j}/dN_{c}",-1,-1,0,0.25);
    Prettify1DwLineStyle (p6jetMult[i], kRed, kDashed, 5, "N_{c}","1/N_{j} dN_{j}/dN_{c}",-1,-1,0,0.25);
    Prettify1DwLineStyle (p6decjetMult[i], kRed, kSolid, 5, "N_{c}","1/N_{j} dN_{j}/dN_{c}",-1,-1,0,0.25);
    Prettify1DwLineStyle (p8decjetMult[i], kBlue, kSolid, 5, "N_{c}","1/N_{j} dN_{j}/dN_{c}",-1,-1,0,0.25);
    Prettify1DwLineStyle (p8jetGirth[i], kBlue, kDashed, 5, "#sum |r_{i}|p_{T,i}/p_{T,jet}","1/N dN/dg",-1,-1,0,16);
    Prettify1DwLineStyle (p8decjetGirth[i], kBlue, kSolid, 5, "#sum |r_{i}|p_{T,i}/p_{T,jet}","1/N dN/dg",-1,-1,0,16);
    Prettify1DwLineStyle (p6jetGirth[i], kRed, kDashed, 5, "#sum |r_{i}|p_{T,i}/p_{T,jet}","1/N dN/dg",-1,-1,0,16);
    Prettify1DwLineStyle (p8jetPhi[i], kBlue, kDashed, 5, "#phi_{jet}","1/N dN/d#phi",-1,-1,0.1,0.25);
    Prettify1DwLineStyle (p6jetPhi[i], kRed, kDashed, 5, "#phi_{jet}","1/N dN/d#phi",-1,-1,0.1,0.25);
    Prettify1DwLineStyle (p8decjetPhi[i], kBlue, kSolid, 5, "#phi_{jet}","1/N dN/d#phi",-1,-1,0.1,0.25);
    Prettify1DwLineStyle (herjetPhi[i], kCyan, kSolid, 5, "#phi_{jet}","1/N dN/d#phi",-1,-1,0.1,0.25);
    Prettify1D (p6GejetPhi[i], kRed, kOpenCircle, 2, kRed, "#phi_{jet}","1/N dN/d#phi",-1,-1,0.1,0.25);
    Prettify1D (datjetPhi[i], kBlack, kOpenStar, 2.5, kBlack, "#phi_{jet}","1/N dN/d#phi",-1,-1,0.1,0.25);
    Prettify1DwLineStyle (p8jetEta[i], kBlue, kDashed, 5, "#eta_{jet}","1/N dN/d#eta",-1,-1,0.6,1.2);
    Prettify1DwLineStyle (p6jetEta[i], kRed, kDashed, 5, "#eta_{jet}","1/N dN/d#eta",-1,-1,0.6,1.2);
    Prettify1DwLineStyle (p8decjetEta[i], kBlue, kSolid, 5, "#eta_{jet}","1/N dN/d#eta",-1,-1,0.6,1.2);
    Prettify1DwLineStyle (herjetEta[i], kCyan, kSolid, 5, "#eta_{jet}","1/N dN/d#eta",-1,-1,0.6,1.2);
    Prettify1D (p6GejetEta[i], kRed, kOpenCircle, 2, kRed, "#eta_{jet}","1/N dN/d#eta",-1,-1,0.6,1.2);
    Prettify1D (datjetEta[i], kBlack, kOpenStar, 2.5, kBlack, "#eta_{jet}","1/N dN/d#eta",-1,-1,0.6,1.2);
  }
   
  TLegend *title = new TLegend(0.55,0.5,0.75,0.7); title->SetBorderSize(0);//TitleLegend(0.2,0.2,0.8,0.8);
  title->AddEntry(decvec[0],"P8 decayed","l");
  title->AddEntry(unvec[0],"P8 undecayed","l");
  title->AddEntry(p6vec[0],"P6 undecayed","l");
  
  TLatex * tslices = new TLatex(); tslices->SetTextAlign(11);
  TLatex *p;
  
  
  /*
  TCanvas *ccomp = MakeCanvas("ccomp","0",1400,1000);
  DivideCanvas(ccomp,"0",3,2);
  
  TH1D* hdummy = new TH1D("hdummy",";M [GeV/c^{2}];1/N_{jets} dN/dM",1,0,10);
  
  hdummy->GetYaxis()->SetRangeUser(0,0.4);
  
  ccomp->cd(1); hdummy->Draw(); p=PanelTitle();
  for (int i = 0; i < nBins; ++ i) {
    ccomp->cd(i+2); unvec[i]->Draw("C"); decvec[i]->Draw("Csame"); p6vec[i]->Draw("Csame"); tslices->DrawLatexNDC(0.55,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if (i==0){title->Draw("same");}
  }
    
  ccomp->SaveAs((out+"mass_compare_p6p8"+filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~constituent girth~~~~~~~~~~~~~~~~~~~~~~//
  
  TCanvas *ccg = MakeCanvas("ccg","0",1400,1000);
  DivideCanvas(ccg,"0",3,2);
  
  TH1D* hdummy_1 = new TH1D("hdummy_1",";|r_{i}|p_{T,i}/p_{T,jet};1/N_{cons} dN_{cons}/dg_{i}",1,0,0.06);

  hdummy_1->GetYaxis()->SetRangeUser(0,65);
  
  ccg->cd(1); hdummy_1->Draw(); p=PanelTitle();
  for (int i = 0; i < nBins; ++i ) {
    ccg->cd(i+2); p8consGirth[i]->Draw("C"); p6consGirth[i]->Draw("Csame"); p8decconsGirth[i]->Draw("Csame"); tslices->DrawLatexNDC(0.55,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if (i==0){title->Draw("same");}
  }
  
  ccg->SaveAs((out+"consGirth_compare_p6p8"+filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~constituent distance~~~~~~~~~~~~~~~~~~~~//
  TCanvas *ccd = MakeCanvas("ccd","0",1400,1000);
  DivideCanvas(ccd,"0",3,2);
  
  TH1D* hdummy_2 = new TH1D("hdummy_2",";|r_{i}|;1/N_{cons} dN_{cons}/dr_{i}",1,0,0.4);

  hdummy_2->GetYaxis()->SetRangeUser(0,6);
  
  ccd->cd(1); hdummy_2->Draw(); p=PanelTitle();
  for (int i = 0; i < nBins; ++i ) {
    ccd->cd(i+2); p8consDist[i]->Draw("C"); p6consDist[i]->Draw("Csame"); p8decconsDist[i]->Draw("Csame"); tslices->DrawLatexNDC(0.55,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if (i==0){title->Draw("same");}
  }
  
  ccd->SaveAs((out+"consDist_compare_p6p8"+filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~constituent pT~~~~~~~~~~~~~~~~~~~~//
  TCanvas *ccpt = MakeCanvas("ccpt","0",1400,1000);
  DivideCanvas(ccpt,"0",3,2);
 
  TH1D* hdummy_3 = new TH1D("hdummy_3",";p_{T,cons};1/N_{cons} dN_{cons}/dp_{T}",1,0.2,30.2);
  hdummy_3->GetYaxis()->SetRangeUser(1e-5,5e-1);
  
  ccpt->cd(1); gPad->SetLogy(); hdummy_3->Draw(); p=PanelTitle();
  for (int i = 0; i < nBins; ++i ) {
    ccpt->cd(i+2); gPad->SetLogy(); p8conspT[i]->Draw("C"); p6conspT[i]->Draw("Csame"); p8decconspT[i]->Draw("Csame");  tslices->DrawLatexNDC(0.55,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if(i==0){title->Draw("same");}
  }
  
  ccpt->SaveAs((out+"conspT_compare_p6p8"+filetype).c_str());
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~jet multiplicity~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  TCanvas *cjm = MakeCanvas("cjm","0",1400,1000);
  DivideCanvas(cjm,"0",3,2);

  TH1D* hdummy_4 = new TH1D("hdummy_4",";N_{c};1/N_{j} dN_{j}/dN_{c}",1,0,25);
  
  hdummy_4->GetYaxis()->SetRangeUser(0,0.25);
  
  cjm->cd(1); hdummy_4->Draw(); p=PanelTitle();
  for (int i = 0; i < nBins; ++i ) {
    cjm->cd(i+2); p8jetMult[i]->Draw("C"); p8decjetMult[i]->Draw("Csame"); p6jetMult[i]->Draw("Csame");  tslices->DrawLatexNDC(0.55,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if(i==0){title->Draw("same");}
    //p6decjetMult[i]->Draw("Csame");
  }
  
  cjm->SaveAs((out+"jetMult_compare_p6p8"+filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~jet girth~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  TCanvas *cjg = MakeCanvas("cjg","0",1400,1000);
  DivideCanvas(cjg,"0",3,2);
  
  TH1D* hdummy_5 = new TH1D("hdummy_5",";#sum |r_{i}|p_{T,i}/p_{T,jet};1/N dN/dg",1,0,0.3);
  hdummy_5->GetYaxis()->SetRangeUser(0,16);
 
  cjg->cd(1); hdummy_5->Draw(); p=PanelTitle();
  for (int i = 0; i < nBins; ++i ) {
    cjg->cd(i+2); p8jetGirth[i]->Draw("C"); p6jetGirth[i]->Draw("Csame"); p8decjetGirth[i]->Draw("Csame");  tslices->DrawLatexNDC(0.55,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if(i==0){title->Draw("same");}
  }
  
  cjg->SaveAs((out+"jetGirth_compare_p6p8"+filetype).c_str());
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~jet eta~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TLegend *stitle = new TLegend(0.05,0.75,0.3,0.95); stitle->SetBorderSize(0);
  stitle->AddEntry(decvec[0],"P8 decayed","l");
  stitle->AddEntry(unvec[0],"P8 undecayed","l");
  stitle->AddEntry(p6vec[0],"P6 undecayed","l");
  stitle->AddEntry(herjetEta[0],"H7 decayed","l");
  
  TLegend *dtitle = new TLegend(0.05,0.75,0.3,0.9); dtitle->SetBorderSize(0);
  dtitle->AddEntry(p6GejetEta[0],"P6+Ge","p");
  dtitle->AddEntry(datjetEta[0],"Raw data","p");
  
  TCanvas *cje = MakeCanvas("cje","0",1400,1000);
  DivideCanvas(cje,"0",3,2);
  
  TH1D* hdummy_6 = new TH1D("hdummy_6",";#eta_{jet};1/N dN/d#eta",1,-1,1);
  hdummy_6->GetYaxis()->SetRangeUser(0.6,1.2);
  
  cje->cd(1); hdummy_6->Draw(); p=PanelTitle();
  for (int i = 0; i < nBins; ++i ) {
    cje->cd(i+2); p8jetEta[i]->Draw("C"); p6jetEta[i]->Draw("Csame"); p8decjetEta[i]->Draw("Csame"); herjetEta[i]->Draw("Csame"); p6GejetEta[i]->Draw("same"); datjetEta[0]->Draw("same"); tslices->DrawLatexNDC(0.55,0.85,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if(i==0){dtitle->Draw("same");} if(i==4){stitle->Draw("same");}
  }
  
  cje->SaveAs((out+"jetEta_compare_p6p8"+filetype).c_str());
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~jet phi~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
  TCanvas *cjp = MakeCanvas("cjp","0",1400,1000);
  DivideCanvas(cjp,"0",3,2);
  
  TH1D* hdummy_7 = new TH1D("hdummy_7",";#phi_{jet};1/N dN/d#phi",1,0,2*M_PI);
  hdummy_7->GetYaxis()->SetRangeUser(0.1,0.25);
  
  cjp->cd(1); hdummy_7->Draw(); p=PanelTitle();
  for (int i = 0; i < nBins; ++i ) {
    cjp->cd(i+2); p8jetPhi[i]->Draw("C"); p6jetPhi[i]->Draw("Csame"); p8decjetPhi[i]->Draw("Csame"); herjetPhi[i]->Draw("Csame"); p6GejetPhi[i]->Draw("same"); datjetPhi[0]->Draw("same"); tslices->DrawLatexNDC(0.55,0.85,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if(i==0){dtitle->Draw("same");} if(i==4){stitle->Draw("same");}
  }
  
  cjp->SaveAs((out+"jetPhi_compare_p6p8"+filetype).c_str());
  */
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~jet pt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  TCanvas *c = MakeCanvas("c","0",1000,1200);

  jetPt6->Scale(1/(double)jetPt6->Integral(),"width");
  jetPt8clone->Scale(1/(double)jetPt8clone->Integral(),"width");
  jetPt8decclone->Scale(1/(double)jetPt8decclone->Integral(),"width");
  jetPt6Ge->Scale(1/(double)jetPt6Ge->Integral(),"width");
  jetPtDat->Scale(1/(double)jetPtDat->Integral(),"width");
  jetPtHer->Scale(1/(double)jetPtHer->Integral(),"width");

  jetPt6->GetYaxis()->SetRangeUser(1e-7,3e-1);
  jetPt8clone->GetYaxis()->SetRangeUser(1e-7,3e-1);
  jetPt8decclone->GetYaxis()->SetRangeUser(1e-7,3e-1);
  jetPtHer->GetYaxis()->SetRangeUser(1e-7,3e-1);
  jetPtDat->GetYaxis()->SetRangeUser(1e-7,3e-1);
  jetPt6Ge->GetYaxis()->SetRangeUser(1e-7,3e-1);
  
  double rat1 = jetPt6->GetBinContent(4) / (double) jetPt8clone->GetBinContent(4);
  double rat2 = jetPt6->GetBinContent(4) / (double) jetPt8decclone->GetBinContent(4);
  double rat3 = jetPt6->GetBinContent(4) / (double) jetPt6Ge->GetBinContent(4);
  double rat4 = jetPt6->GetBinContent(4) / (double) jetPtDat->GetBinContent(4);
  double rat5 = jetPt6->GetBinContent(4) / (double) jetPtHer->GetBinContent(4);
  
  jetPt8clone->Scale(rat1);
  jetPt8decclone->Scale(rat2);
  jetPt6Ge->Scale(rat3);
  jetPtDat->Scale(rat4);
  jetPtHer->Scale(rat5);

  jetPt8nostyle->Scale(rat1);
  jetPt8decnostyle->Scale(rat2);
  jetPt6Genostyle->Scale(rat3);
  jetPtDatnostyle->Scale(rat4);
  jetPtHernostyle->Scale(rat5);
  
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetLogy();
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  jetPt6->SetStats(0);          // No statistics on upper plot

  jetPt6->Draw("Csame"); jetPt8clone->Draw("Csame"); jetPt8decclone->Draw("Csame"); jetPtHer->Draw("Csame"); jetPt6Ge->Draw("same"); jetPtDat->Draw("same"); 

  TLatex *t = new TLatex();
  t->SetTextAlign(11);                                                                                                         
  t->DrawLatexNDC(0.15,0.3, "pp 200 GeV");
  t->DrawLatexNDC(0.15,0.2, "anti-k_{T}, R = 0.4");
  t->DrawLatexNDC(0.15,0.1, "Ch+Ne jets, |#eta| < 0.6");

  // lower plot will be in pad
  c->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  // Define the ratio plot
  TH1D *p6p8 = (TH1D*)jetPt8nostyle->Clone("p6p8");
  p6p8->SetStats(0);      // No statistics on lower plot
  p6p8->Divide(jetPt6nostyle);
  Prettify1D(p6p8,kGreen+3, kOpenSquare,2,kGreen+3,"p_{jet,T} [GeV/c]","ratios",5,60,0.8,1.2);
  p6p8->SetTitle(""); // Remove the ratio title
  
  TH1D *herp8 = (TH1D*)jetPt8nostyle->Clone("herp8");
  herp8->SetStats(0);      // No statistics on lower plot
  herp8->Divide(jetPtHernostyle);
  Prettify1D(herp8,kMagenta, kOpenSquare,2,kMagenta,"p_{jet,T} [GeV/c]","ratios",5,60,0.8,1.2);
  herp8->SetTitle(""); // Remove the ratio title

  TH1D *dg = (TH1D*)jetPtDatnostyle->Clone("dg");
  dg->SetStats(0);      // No statistics on lower plot
  dg->Divide(jetPt6Genostyle);
  Prettify1D(dg,kOrange, kOpenSquare,2,kOrange,"p_{jet,T} [GeV/c]","ratios",5,60,0.8,1.2);
  dg->SetTitle(""); // Remove the ratio title

  TH1D *decun = (TH1D*)jetPt8decnostyle->Clone("decun");
  decun->SetStats(0);      // No statistics on lower plot
  decun->Divide(jetPt8nostyle);
  Prettify1D(decun,kBlue, kOpenSquare,2,kBlue,"p_{jet,T} [GeV/c]","ratios",5,60,0.8,1.2);
  decun->SetTitle(""); // Remove the ratio title
  
  TLegend *rats = new TLegend(0.3,0.6,0.55,0.8); rats->SetBorderSize(0);
  rats->AddEntry(decun,"P8 decayed / undecayed","p");
  rats->AddEntry(dg,"Raw data / P6+Ge","p");
  rats->AddEntry(herp8,"P8 / H7","p");
  rats->AddEntry(p6p8, "P8 undecayed / P6 undecayed","p");
  
  TLegend *specs = new TLegend(0.6,0.2,0.75,0.59); specs->SetBorderSize(0);
  specs->AddEntry(jetPt6, "P6 undecayed","l"); specs->AddEntry(jetPt8clone,"P8 undecayed","l");
  specs->AddEntry(jetPt8decclone, "P6 decayed","l"); specs->AddEntry(jetPtHer,"H7 decayed","l");
  specs->AddEntry(jetPt6Ge,"P6+Ge","p"); specs->AddEntry(jetPtDat,"Raw data","p");
  
  pad1->cd();
  gStyle->SetLegendTextSize(0.04);
  rats->Draw("same"); specs->Draw("same");
  pad2->cd();
  
  p6p8->Draw("p");       // Draw the ratio plot
  herp8->Draw("psame");
  dg->Draw("psame");
  decun->Draw("same");

  TLine *fifteen_up = new TLine(15,1e-7,15,3e-1);
  TLine *fifteen_down = new TLine(15,0.8,15,1.2);
  TLine *one = new TLine(5,1,60,1);
  
  fifteen_up->SetLineStyle(kDashed); fifteen_down->SetLineStyle(kDashed);
  one->SetLineStyle(kDashed);
  
  fifteen_down->Draw("same"); one->Draw("same");
  pad1->cd(); fifteen_up->Draw("same");
  pad2->cd();
  
  // Y axis ratio plot settings
  p6p8->GetYaxis()->SetNdivisions(405,1);
  p6p8->GetYaxis()->SetTitleSize(20);
  p6p8->GetYaxis()->SetTitleFont(43);
  p6p8->GetYaxis()->SetTitleOffset(1.4);
  p6p8->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  p6p8->GetYaxis()->SetLabelSize(15);

  // X axis ratio plot settings
  p6p8->GetXaxis()->SetTitleSize(20);
  p6p8->GetXaxis()->SetTitleFont(43);
  p6p8->GetXaxis()->SetTitleOffset(4.);
  p6p8->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  p6p8->GetXaxis()->SetLabelSize(15);

  c->SaveAs((out+"jetPt_compare_all_curves"+filetype).c_str()); 
  
  return;
}
