
using namespace std;
/*
//loops over the unweighted response object's response histogram and finds any bins where there are fewer than 20 counts. These statistics-limited bins are dropped from the corresponding weighted response.                                                                     
void DropLowStatsBins(RooUnfoldResponse* weighted, RooUnfoldResponse* unweighted) {                                              
  TH2D* unweighted_res = (TH2D*) unweighted->Hresponse();                                                                                                                                                                                                           
  for (int i = 0; i <= unweighted_res->GetXaxis()->GetNbins(); ++ i) {                                                            
    for (int j = 0; j <= unweighted_res->GetYaxis()->GetNbins(); ++ j) {                                                 
      if (unweighted_res->GetBinContent(i,j) < 20) {                                                                             
	weighted->Hresponse()->SetBinContent(i,j,0);                                                                           
	weighted->Hresponse()->SetBinError(i,j,0);                                                                               
      }                                                                                                                      
    }                                                                                                                       
  }                                                                                                                                                       
  return;                                                                                                                                                 
}  
*/

//loops over the 2D spectra to find bins with bin content 0 (low-stats bins should already have been dropped in the spectra before this function is called). These bins are passed to the response, to be dropped there as well.
//PS this assumes mass (or whatever observable) is on the x-axis!
vector<int> DropLowStatsBins(RooUnfoldResponse* weighted, TH2D* det_unweighted) {//, TH2D* gen) {
  //TH2D* unweighted_res = (TH2D*) unweighted->Hresponse();

  vector<int> bins_to_drop;

  //this method requires the binning for the TH2 to be the same as for the RooUnfoldResponse
  for (int i = 1; i <= det_unweighted->GetXaxis()->GetNbins(); ++ i) {//mass bins in det-level spectrum
    for (int j = 1; j <= det_unweighted->GetYaxis()->GetNbins(); ++ j) {//pT bins in det-level spectrum
      if (det_unweighted->GetBinContent(i,j) < 20) {
	cout << "For mass bin " << i << " and pT bin " << j << ": " << endl;
	cout << "Removing column: " << (det_unweighted->GetXaxis()->GetNbins()*(j-1))+i << endl;
	bins_to_drop.push_back((det_unweighted->GetXaxis()->GetNbins()*(j-1))+i);
	for (int k = 1; k <= weighted->Hresponse()->GetYaxis()->GetNbins(); ++ k) {//gen bins in response
	  weighted->Hresponse()->SetBinContent((det_unweighted->GetXaxis()->GetNbins()*(j-1))+i,k,0);
	  weighted->Hresponse()->SetBinError((det_unweighted->GetXaxis()->GetNbins()*(j-1))+i,k,0);                                      
	}                                                                                                                      
      }                                                                                                                       
    }
  }
  return bins_to_drop;
}


//loops over the unweighted 1D histogram and finds any bins where there are fewer than 20 counts. These statistics-limited bins are dropped from the corresponding weighted histogram.
void DropLowStatsBins(TH1D* weighted, TH1D* unweighted) {                                                                         
  for (int i = 1; i <= unweighted->GetXaxis()->GetNbins(); ++ i) {                                                         
    if (unweighted->GetBinContent(i) < 20) {                                                                                  
      weighted->SetBinContent(i,0);                                                                                          
      weighted->SetBinError(i,0);                                                                                             
    }                                                                                                                            
  }                                                                                                                               

  return;                                                                                                                                                 
}                                                                                                                                                         
                                                                                                                                                            
//loops over the unweighted 2D histogram and finds any bins where there are fewer than 20 counts. These statistics-limited bins are dropped from the corresponding weighted histogram.
void DropLowStatsBins(TH2D* weighted, TH2D* unweighted) {                                                                
  for (int i = 1; i <= unweighted->GetXaxis()->GetNbins(); ++ i) {                                                        
    for (int j = 1; j <= unweighted->GetYaxis()->GetNbins(); ++ j) {                                                            
      if (unweighted->GetBinContent(i,j) < 20) {                                                                             
	weighted->SetBinContent(i,j,0);                                                                                      
	weighted->SetBinError(i,j,0);                                                                                        
      }                                                                                                                  
    }                                                                                                                            
  }                                                                                                                                                       
  return;                                                                                                                                                 
} 

void DropBins(RooUnfoldResponse *response, vector<int> to_drop) {
  for (int i = 0; i < to_drop.size(); ++ i) { //make sure you index from 0 to size - 1 here since you're using the vector.
    for (int j = 1; j <= response->Hresponse()->GetYaxis()->GetNbins(); ++ j) { //indexing response here so 1 to size.
      response->Hresponse()->SetBinContent(to_drop[i], j, 0);
      response->Hresponse()->SetBinError(to_drop[i], j, 0);
    }
  }
  
  return;
}

void bin_drop() {
  const string syst_path = "~/jetmass/out/systematics/";
  const string match_path = "~/jetmass/out/matching/";
  const string clos_path = "~/jetmass/out/closure/";
  const string data_path = "~/jetmass/macros/hists/";
  const string file_in = "full_w_o_bin_drop_R04.root";
  const string data_file_in = "hists_w_o_bin_drop_R04.root";
  const string file_out = "full_R04.root";
  const string data_file_out = "hists_R04.root";
  
  TFile *match_in = new TFile((match_path+file_in).c_str(),"READ");
  TFile *syst_in = new TFile((syst_path+file_in).c_str(),"READ");
  TFile *clos_in = new TFile((clos_path+file_in).c_str(),"READ");
  TFile *data_in = new TFile((data_path+data_file_in).c_str(),"READ");

  //~~~closure - must be done first
  RooUnfoldResponse *pt_m_response_even = (RooUnfoldResponse*) clos_in->Get("pt_m_response");
  RooUnfoldResponse *pt_mg_response_even = (RooUnfoldResponse*) clos_in->Get("pt_mg_response");
  RooUnfoldResponse *pt_m_response_odd = (RooUnfoldResponse*) clos_in->Get("pt_m_response_odd");
  RooUnfoldResponse *pt_mg_response_odd = (RooUnfoldResponse*) clos_in->Get("pt_mg_response_odd");
  RooUnfoldResponse *pt_m_response_even_counts = (RooUnfoldResponse*) clos_in->Get("pt_m_response_counts");
  RooUnfoldResponse *pt_mg_response_even_counts = (RooUnfoldResponse*) clos_in->Get("pt_mg_response_counts");
  RooUnfoldResponse *pt_m_response_odd_counts = (RooUnfoldResponse*) clos_in->Get("pt_m_response_odd_counts");
  RooUnfoldResponse *pt_mg_response_odd_counts = (RooUnfoldResponse*) clos_in->Get("pt_mg_response_odd_counts");
  RooUnfoldResponse *pt_even = (RooUnfoldResponse*) clos_in->Get("pt_even");
  RooUnfoldResponse *pt_odd = (RooUnfoldResponse*) clos_in->Get("pt_odd");
  
  TH1D *pt_gen_odd = (TH1D*) clos_in->Get("pt_gen_odd");
  TH1D *pt_det_odd = (TH1D*) clos_in->Get("pt_det_odd");
  TH1D *pt_gen_even = (TH1D*) clos_in->Get("pt_gen_even");
  TH1D *pt_det_even = (TH1D*) clos_in->Get("pt_det_even");
  TH1D *pt_gen_odd_counts = (TH1D*) clos_in->Get("pt_gen_odd_counts");
  TH1D *pt_det_odd_counts = (TH1D*) clos_in->Get("pt_det_odd_counts");
  TH1D *pt_gen_even_counts = (TH1D*) clos_in->Get("pt_gen_even_counts");
  TH1D *pt_det_even_counts = (TH1D*) clos_in->Get("pt_det_even_counts");

  TH2D *pt_m_gen_odd = (TH2D*) clos_in->Get("pt_m_gen_odd");
  TH2D *pt_m_det_odd = (TH2D*) clos_in->Get("pt_m_det_odd");
  TH2D *pt_m_gen_even = (TH2D*) clos_in->Get("pt_m_gen_even");
  TH2D *pt_m_det_even = (TH2D*) clos_in->Get("pt_m_det_even");
  TH2D *pt_m_gen_odd_counts = (TH2D*) clos_in->Get("pt_m_gen_odd_counts");
  TH2D *pt_m_det_odd_counts = (TH2D*) clos_in->Get("pt_m_det_odd_counts");
  TH2D *pt_m_gen_even_counts = (TH2D*) clos_in->Get("pt_m_gen_even_counts");
  TH2D *pt_m_det_even_counts = (TH2D*) clos_in->Get("pt_m_det_even_counts");
  TH2D *pt_mg_gen_odd = (TH2D*) clos_in->Get("pt_mg_gen_odd");
  TH2D *pt_mg_det_odd = (TH2D*) clos_in->Get("pt_mg_det_odd");
  TH2D *pt_mg_gen_even = (TH2D*) clos_in->Get("pt_mg_gen_even");
  TH2D *pt_mg_det_even = (TH2D*) clos_in->Get("pt_mg_det_even");
  TH2D *pt_mg_gen_odd_counts = (TH2D*) clos_in->Get("pt_mg_gen_odd_counts");
  TH2D *pt_mg_det_odd_counts = (TH2D*) clos_in->Get("pt_mg_det_odd_counts");
  TH2D *pt_mg_gen_even_counts = (TH2D*) clos_in->Get("pt_mg_gen_even_counts");
  TH2D *pt_mg_det_even_counts = (TH2D*) clos_in->Get("pt_mg_det_even_counts");

                                                          
  DropLowStatsBins(pt_gen_odd, pt_gen_odd_counts);                                                                                
  DropLowStatsBins(pt_det_odd, pt_det_odd_counts);                                                                                
  DropLowStatsBins(pt_gen_even, pt_gen_even_counts);                                                                              
  DropLowStatsBins(pt_det_even, pt_det_even_counts);        

  DropLowStatsBins(pt_m_gen_odd, pt_m_gen_odd_counts);                                                                            
  DropLowStatsBins(pt_m_det_odd, pt_m_det_odd_counts);                                                                            
  DropLowStatsBins(pt_m_gen_even, pt_m_gen_even_counts);                                                                          
  DropLowStatsBins(pt_m_det_even, pt_m_det_even_counts); 
  
  vector<int> to_drop_m; vector<int> to_drop_mg;
  
  //remember to call this block of after spectra have been pruned already
  to_drop_m = DropLowStatsBins(pt_m_response_even, pt_m_det_even_counts);//, pt_m_gen_even);
  to_drop_mg = DropLowStatsBins(pt_mg_response_even, pt_mg_det_even_counts);//, pt_mg_gen_even);
  cout << "A" << endl;
  DropBins(pt_m_response_odd, to_drop_m);//, pt_m_det_odd_counts);//, pt_m_gen_odd);
  DropBins(pt_mg_response_odd, to_drop_mg);//, pt_mg_det_odd_counts);//, pt_mg_gen_odd); 
  //~~~
  
  //~~~matching
  RooUnfoldResponse *pt_m_response = (RooUnfoldResponse*) match_in->Get("pt_m_response");
  RooUnfoldResponse *pt_mg_response = (RooUnfoldResponse*) match_in->Get("pt_mg_response");
  RooUnfoldResponse *pt_m_response_counts = (RooUnfoldResponse*) match_in->Get("pt_m_response_counts");
  RooUnfoldResponse *pt_mg_response_counts = (RooUnfoldResponse*) match_in->Get("pt_mg_response_counts");
  cout << "B" << endl;
  DropBins(pt_m_response, to_drop_m);//, pt_m_response_counts);   
  DropBins(pt_mg_response, to_drop_mg);//, pt_mg_response_counts);   
  
  //~~~


  //~~~systematics
  RooUnfoldResponse *pt_m_res_nom = (RooUnfoldResponse*) syst_in->Get("pt_m_res_nom");
  RooUnfoldResponse *pt_m_res_TS = (RooUnfoldResponse*) syst_in->Get("pt_m_res_TS");
  RooUnfoldResponse *pt_m_res_TU = (RooUnfoldResponse*) syst_in->Get("pt_m_res_TU");
  RooUnfoldResponse *pt_m_res_HC50 = (RooUnfoldResponse*) syst_in->Get("pt_m_res_HC50");
  RooUnfoldResponse *pt_m_res_HC0 = (RooUnfoldResponse*) syst_in->Get("pt_m_res_HC0");
  RooUnfoldResponse *pt_m_res_DS = (RooUnfoldResponse*) syst_in->Get("pt_m_res_DS");
  RooUnfoldResponse *pt_m_res_GS = (RooUnfoldResponse*) syst_in->Get("pt_m_res_GS");
  RooUnfoldResponse *pt_m_res_MS = (RooUnfoldResponse*) syst_in->Get("pt_m_res_MS");
  RooUnfoldResponse *pt_mg_res_nom = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_nom");
  RooUnfoldResponse *pt_mg_res_TS = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_TS");
  RooUnfoldResponse *pt_mg_res_TU = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_TU");
  RooUnfoldResponse *pt_mg_res_HC50 = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_HC50");
  RooUnfoldResponse *pt_mg_res_HC0 = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_HC0");
  RooUnfoldResponse *pt_mg_res_DS = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_DS");
  RooUnfoldResponse *pt_mg_res_GS = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_GS");
  RooUnfoldResponse *pt_mg_res_MS = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_MS");
  RooUnfoldResponse *pt_m_res_nom_counts = (RooUnfoldResponse*) syst_in->Get("pt_m_res_nom_counts");
  RooUnfoldResponse *pt_m_res_TS_counts = (RooUnfoldResponse*) syst_in->Get("pt_m_res_TS_counts");
  RooUnfoldResponse *pt_m_res_TU_counts = (RooUnfoldResponse*) syst_in->Get("pt_m_res_TU_counts");
  RooUnfoldResponse *pt_m_res_HC50_counts = (RooUnfoldResponse*) syst_in->Get("pt_m_res_HC50_counts");
  RooUnfoldResponse *pt_m_res_HC0_counts = (RooUnfoldResponse*) syst_in->Get("pt_m_res_HC0_counts");
  RooUnfoldResponse *pt_m_res_DS_counts = (RooUnfoldResponse*) syst_in->Get("pt_m_res_DS_counts");
  RooUnfoldResponse *pt_m_res_GS_counts = (RooUnfoldResponse*) syst_in->Get("pt_m_res_GS_counts");
  RooUnfoldResponse *pt_m_res_MS_counts = (RooUnfoldResponse*) syst_in->Get("pt_m_res_MS_counts");
  RooUnfoldResponse *pt_mg_res_nom_counts = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_nom_counts");
  RooUnfoldResponse *pt_mg_res_TS_counts = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_TS_counts");
  RooUnfoldResponse *pt_mg_res_TU_counts = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_TU_counts");
  RooUnfoldResponse *pt_mg_res_HC50_counts = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_HC50_counts");
  RooUnfoldResponse *pt_mg_res_HC0_counts = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_HC0_counts");
  RooUnfoldResponse *pt_mg_res_DS_counts = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_DS_counts");
  RooUnfoldResponse *pt_mg_res_GS_counts = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_GS_counts");
  RooUnfoldResponse *pt_mg_res_MS_counts = (RooUnfoldResponse*) syst_in->Get("pt_mg_res_MS_counts");
  cout << "C" << endl;
  DropBins(pt_m_res_nom,to_drop_m);// pt_m_res_nom_counts);
  cout<< "i" <<endl;
  DropBins(pt_mg_res_nom,to_drop_mg);// pt_mg_res_nom_counts);
  DropBins(pt_m_res_TS,to_drop_m);// pt_m_res_TS_counts);
  DropBins(pt_mg_res_TS,to_drop_mg);// pt_mg_res_TS_counts);
  DropBins(pt_m_res_TU,to_drop_m);// pt_m_res_TU_counts);
  DropBins(pt_mg_res_TU,to_drop_mg);// pt_mg_res_TU_counts); 
  DropBins(pt_m_res_HC50,to_drop_m);// pt_m_res_HC50_counts);
  DropBins(pt_mg_res_HC50,to_drop_mg);// pt_mg_res_HC50_counts);
  DropBins(pt_m_res_HC0,to_drop_m);// pt_m_res_HC0_counts);
  DropBins(pt_mg_res_HC0,to_drop_mg);// pt_mg_res_HC0_counts);
  DropBins(pt_m_res_DS,to_drop_m);// pt_m_res_DS_counts);
  DropBins(pt_mg_res_DS,to_drop_mg);// pt_mg_res_DS_counts);
  cout << "ii" <<endl;
  DropBins(pt_m_res_GS,to_drop_m);// pt_m_res_GS_counts);
  cout << "iia" << endl;
  DropBins(pt_mg_res_GS,to_drop_mg);// pt_mg_res_GS_counts);   
  cout << "iii" << endl;
  DropBins(pt_m_res_MS,to_drop_m);// pt_m_res_MS_counts);
  cout << "iv" << endl;
  DropBins(pt_mg_res_MS,to_drop_mg);// pt_mg_res_MS_counts);   

  cout <<"D"<<endl;
  //~~~

  //~~~data
  TH2D *m_v_pt_d = (TH2D*) data_in->Get("m_v_pt_d");
  TH2D *m_v_pt_g = (TH2D*) data_in->Get("m_v_pt_g");
  TH2D *m_v_pt_p = (TH2D*) data_in->Get("m_v_pt_p");
  TH2D *m_v_pt_d_counts = (TH2D*) data_in->Get("m_v_pt_d_counts");
  TH2D *m_v_pt_g_counts = (TH2D*) data_in->Get("m_v_pt_g_counts");
  TH2D *m_v_pt_p_counts = (TH2D*) data_in->Get("m_v_pt_p_counts");
  TH2D *mg_v_pt_d = (TH2D*) data_in->Get("mg_v_pt_d");
  TH2D *mg_v_pt_g = (TH2D*) data_in->Get("mg_v_pt_g");
  TH2D *mg_v_pt_p = (TH2D*) data_in->Get("mg_v_pt_p");
  TH2D *mg_v_pt_d_counts = (TH2D*) data_in->Get("mg_v_pt_d_counts");
  TH2D *mg_v_pt_g_counts = (TH2D*) data_in->Get("mg_v_pt_g_counts");
  TH2D *mg_v_pt_p_counts = (TH2D*) data_in->Get("mg_v_pt_p_counts");
  
  DropLowStatsBins(m_v_pt_d, m_v_pt_d_counts);                                                                             
  DropLowStatsBins(m_v_pt_g, m_v_pt_g_counts);                                                                             
  DropLowStatsBins(m_v_pt_p, m_v_pt_p_counts);
  DropLowStatsBins(mg_v_pt_d, mg_v_pt_d_counts);                                                                             
  DropLowStatsBins(mg_v_pt_g, mg_v_pt_g_counts);                                                                             
  DropLowStatsBins(mg_v_pt_p, mg_v_pt_p_counts);

  //~~~

  TFile *match_out = new TFile((match_path+file_out).c_str(),"RECREATE");
  TFile *syst_out = new TFile((syst_path+file_out).c_str(),"RECREATE");
  TFile *clos_out = new TFile((clos_path+file_out).c_str(),"RECREATE");
  TFile *data_out = new TFile((data_path+data_file_out).c_str(),"RECREATE");

  match_out->cd();
  pt_m_response->Write(); pt_mg_response->Write();
  
  clos_out->cd();
  pt_even->Write(); pt_odd->Write();
  pt_m_response_even->Write(); pt_mg_response_even->Write();
  pt_m_response_odd->Write(); pt_mg_response_odd->Write();
  pt_gen_odd->Write(); pt_gen_even->Write(); pt_det_odd->Write(); pt_det_even->Write();
  pt_m_gen_odd->Write(); pt_m_gen_even->Write(); pt_m_det_odd->Write(); pt_m_det_even->Write();
  pt_mg_gen_odd->Write(); pt_mg_gen_even->Write(); pt_mg_det_odd->Write(); pt_mg_det_even->Write();
  
  syst_out->cd();
  pt_m_res_nom->Write(); pt_mg_res_nom->Write();
  pt_m_res_TS->Write(); pt_mg_res_TS->Write();
  pt_m_res_TU->Write(); pt_mg_res_TU->Write();
  pt_m_res_HC50->Write(); pt_mg_res_HC50->Write();
  pt_m_res_HC0->Write(); pt_mg_res_HC0->Write();
  pt_m_res_DS->Write(); pt_mg_res_DS->Write();
  pt_m_res_GS->Write(); pt_mg_res_GS->Write();
  pt_m_res_MS->Write(); pt_mg_res_MS->Write();

  data_out->cd();
  m_v_pt_d->Write(); m_v_pt_g->Write(); m_v_pt_p->Write();
  mg_v_pt_d->Write(); mg_v_pt_g->Write(); mg_v_pt_p->Write();

  cout << "Dropped mass bins: " << endl;
  for (int i = 0; i < to_drop_m.size() - 1; ++ i) {
    cout << to_drop_m[i] << ", ";
  }
  cout << to_drop_m[to_drop_m.size() - 1] << endl;
  
  cout << "Dropped Mg bins: " << endl;
  for (int i = 0; i < to_drop_mg.size() - 1; ++ i) {
    cout << to_drop_mg[i] << ", ";
  }
  cout << to_drop_mg[to_drop_mg.size() - 1] << endl;
  
  
  return;
}
