//  functions.cxx
//  Veronica Verkest May 13, 2018
//  Adapted by Isaac Mooney June, 2018

#include "params.hh"
#include "funcs.hh"
//#include "ktTrackEff.hh"

typedef fastjet::contrib::SoftDrop SD;

namespace Analysis {

  // -------------------------                                                                                                                                                                            
  // IO/OS Manip functionality                                                                                                                                                                            
  // -------------------------                                                                                                                                                                             
  // Used to understand which format of input file is being used                                                                                                                                          
  // ( .root file, .txt, .list, etc )                                                                                                                                                                     
  // ---------------------------------------------------------------------                                                                                                                                 
  bool HasEnding (std::string const &full_string, std::string const &ending) {
    if (full_string.length() >= ending.length()) {
      return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
    } else {
      return false;
    }
  }

  // parse a CSV file to a set of unique entries.
  // All comments must start on their own line, and be proceeded   
  // by a pound sign (#)                                                                                                                                                              
  template <typename T>
  std::set<T> ParseCSV(std::string csv) {
    // return set                                                                                                        
    std::set<T> ret;
    std::ifstream fs(csv);
    std::string line;
    // first, split by line
    while (std::getline(fs, line)) {
      if (line.size() == 0) // reject empty lines   
	continue;
      if (line[0] == '#') // reject comments                                                                                                                                          
	continue;
      // split the string by commas                                                                                                                                                   
      std::istringstream ss(line);
      while (ss) {
	std::string str_value;
	std::getline(ss, str_value, ',');
	if (CanCast<T>(str_value)) {
	  ret.insert(CastTo<T>(str_value));
	}
      }
    }
    return ret;
  }

  template<typename T>
  bool CanCast(std::string s) {
    std::istringstream iss(s);
    T dummy;
    iss >> std::skipws >> dummy;
    return iss && iss.eof();
  }

  template<typename T>
  T CastTo(std::string s) {
    std::istringstream iss(s);
    T dummy;
    iss >> std::skipws >> dummy;
    return dummy;
  }


  void FillMatchedTree (const fastjet::PseudoJet pyjet, const fastjet::PseudoJet gejet, TTree * Tree, double &pyPt, double & gePt, double & pyM, double &geM, double &deltaPt, double &deltaM, double &Mratio, double &weight, const double mc_weight) {
    pyPt = pyjet.pt(); gePt = gejet.pt(); pyM = pyjet.m(); geM = gejet.m(); 
    deltaM = gejet.m() - pyjet.m(); deltaPt = gejet.pt() - pyjet.pt(); Mratio = gejet.m() / (double) pyjet.m(); weight = mc_weight;
    Tree->Fill();
    return;
  }

  void FillTrees ( std::vector<fastjet::PseudoJet> jets, TTree* Tree, double &jPt, double &jEta, double &jPhi, double &jM, double &jE, int &jncons, double &wt, const double weight) {
    for ( int j = 0; j< jets.size(); ++ j) {   // FILL JET INFO
      int nGhosts = 0;
      if (jets[j].pt() < 0.2) continue;
      jPt = jets[j].pt();    jEta = jets[j].eta();    jPhi = jets[j].phi();
      jE = jets[j].e();    jM = jets[j].m(); wt = weight;
      std::vector<fastjet::PseudoJet> Cons = jets[j].constituents(); //jncons = Cons.size();
      //cons without ghosts:
      for (int c = 0; c < Cons.size(); ++ c) {
	if (Cons[c].pt() < 0.2) ++nGhosts;
      }
      jncons = Cons.size() - nGhosts;
      Tree->Fill();
    }
    return;
  }

  void AnalysisSummary( int events, int pJets, int eJets, int gJets, int pgMatchedJets, int epMatchedJets, int egMatchedJets, std::string outName ) {
    std::cout << std::endl << std::endl << " Of " << events << " events: "<< std::endl;
    std::cout << pJets << " jets have been found for the Pythia6 data" << std::endl;
    std::cout << eJets << " jets have been found for the Pythia6 + Efficiency data" << std::endl;
    std::cout << gJets << " jets have been found for the Pythia6 + GEANT data" << std::endl << std::endl;
    std::cout << pgMatchedJets << " GEANT jets have been matched to Pythia6 jets" << std::endl;
    std::cout << epMatchedJets << " Pythia6+Efficiency jets have been matched to Pythia6 jets" << std::endl;
    std::cout << egMatchedJets << " GEANT jets have been matched to Pythia6+Efficiency jets" << std::endl;
    std::cout <<std::endl << "Writing to:  " << outName << std::endl << std::endl;
  }
  //remind self later: why am I passing TStarJetVector *sv? should just be declaring it before the for-loop
  void GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container, TStarJetVector *sv, std::vector<fastjet::PseudoJet> & Particles, const bool full, const bool py){
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );

      // current.set_user_index( sv->GetCharge() );
      //      if (py == 1 && full == 1) {std::cout << "starting with charge: " << sv->GetCharge() << " for particle " << i << std::endl;}
      
      if (sv->GetCharge() != 0 /*&& py == 0*/) {
	current.reset_PtYPhiM(sqrt(current.perp2()),current.rap(),current.phi(), chPionMass); //assigning pion mass to charged particles
      }/*
      if (py != 0) {
	current.reset_PtYPhiM(sqrt(current.perp2()), current.rap(), current.phi(), current.m());
      }
       */
      if ((sv->GetCharge() == 0) && (full == 0)) { continue; } // if we don't want full jets, skip neutrals
      current.set_user_index( sv->GetCharge() );
      // if (py == 1 && full == 1) {std::cout << "ending with charge: " << current.user_index() << " for particle " << i << std::endl;}
      Particles.push_back(current);
    }
    return;
  }

  void ApplyNEFSelection(const std::vector<fastjet::PseudoJet> init, std::vector<fastjet::PseudoJet> & result) {
    //Implementing a neutral energy fraction cut of 90% on inclusive jets                                                                                                                           
    for (int i = 0; i < init.size(); ++ i) {
      double towersum = 0; double ptsum = 0;
      for (int j = 0; j < init[i].constituents().size(); ++ j) {
	if (init[i].constituents()[j].user_index() == 0) {
	  towersum += init[i].constituents()[j].pt();
	}
	ptsum += init[i].constituents()[j].pt();
      }
      if (towersum / (double) ptsum < NEF_max) {
	result.push_back(init[i]);
      }
    }
    return;
  }
  
  
  double LookupRun6Xsec(TString currentfile ) {

    static const Double_t Xsec[12] = {
      1.0,        // Placeholder for 2-3
      1.30E+09, // 3-4
      3.15E+08, // 4-5
      1.37E+08, // 5-7
      2.30E+07, // 7-9
      5.53E+06, // 9-11
      2.22E+06, // 11-15
      3.90E+05, // 15-25
      1.02E+04, // 25-35
      5.01E+02, // 35-45
      2.86E+01, // 45-55
      1.46E+00 // 55-65
    };

    static const Double_t Nmc[12] = {
      1, // 2-3
      672518, // 3-4
      672447, // 4-5
      393498, // 5-7
      417659, // 7-9
      412652, // 9-11
      419030, // 11-15
      396744, // 15-25
      399919, // 25-35
      119995, // 35-45
      117999, // 45-55
      119999 // 55-65
    };

    Double_t w[12];
    for ( int i=0; i<12 ; ++i ){
      w[i] = Xsec[i] / Nmc[i];
    }

    if ( currentfile.Contains("3_4") ) return w[1];
    if ( currentfile.Contains("4_5") ) return w[2];
    if ( currentfile.Contains("5_7") ) return w[3];
    if ( currentfile.Contains("7_9") ) return w[4];
    if ( currentfile.Contains("9_11") ) return w[5];
    if ( currentfile.Contains("11_15") ) return w[6];
    if ( currentfile.Contains("15_25") ) return w[7];
    if ( currentfile.Contains("25_35") ) return w[8];
    if ( currentfile.Contains("35_45") ) return w[9];
    if ( currentfile.Contains("45_55") ) return w[10];
    if ( currentfile.Contains("55_65") ) return w[11];
    return 1;
  }

  //----------------------------------------------------------------------
  double LookupRun12Xsec( TString filename ){
  
    const int NUMBEROFPT = 11;
    // const char *PTBINS[NUMBEROFPT]={"2_3","3_4","4_5","5_7","7_9","9_11","11_15","15_20","20_25","25_35","35_-1"};
    const static float XSEC[NUMBEROFPT] = {9.00581646, 1.461908221, 0.3544350863, 0.1513760388, 0.02488645725, 0.005845846143, 0.002304880181, 0.000342661835, 4.562988397e-05, 9.738041626e-06, 5.019978175e-07};
    const static float NUMBEROFEVENT[NUMBEROFPT] = {2100295, 600300, 600300, 300289, 300289, 300289, 160295, 100302, 80293, 76303, 23307};

    const static std::vector<std::string> vptbins={"pp12Pico_pt2_3","pp12Pico_pt3_4","pp12Pico_pt4_5","pp12Pico_pt5_7","pp12Pico_pt7_9","pp12Pico_pt9_11","pp12Pico_pt11_15","pp12Pico_pt15_20","pp12Pico_pt20_25","pp12Pico_pt25_35","35_-1"};
    for ( int i=0; i<vptbins.size(); ++i ){
      if ( filename.Contains(vptbins.at(i).data())) return XSEC[i] / NUMBEROFEVENT[i];
    }

    throw std::runtime_error("Not a valid filename");
    return -1;
  
  }


  //finds potential trigger jets out of the two highest pT jets
  bool GetTriggerJet(std::vector<fastjet::PseudoJet> & triggers, const std::vector<fastjet::PseudoJet> jets) {
    triggers.clear();
    bool placeholder = 0; //to keep track of which jet was the trigger in the case of only one trigger 
    for (int i = 0; i < jets.size(); ++ i) {
      if (i == 2) {return placeholder;} //only want to look at the two highest pT jets
      for (int j = 0; j < jets[i].constituents().size(); ++ j) {
	if (jets[i].constituents()[j].pt() > det_evEtMin) {//has a trigger
	  placeholder = i;
	  triggers.push_back(jets[i]);
	  break;
	}
      }
      if (triggers.size() == 2) {return placeholder;} //we only need at most two objects in the triggers vector - one to be trigger, one to be recoil
    }
    return placeholder;
  }
  /*
  void MatchJets(const std::vector<fastjet::PseudoJet> candidates, const fastjet::PseudoJet toMatch, int & match_position) {                                                         match_position = -1;
    if (candidates.size() == 0) {
      return;
    }
    // match the leading jet                                                                                                                                                    
    fastjet::Selector selectMatchedLead = fastjet::SelectorCircle( R );
    selectMatchedLead.set_reference( toMatch );
    std::vector<fastjet::PseudoJet> matchedToLead = sorted_by_pt( selectMatchedLead( candidates ));
    if (matchedToLead.size() == 0) {return;}
    //If here, found match(es)
    for (int i = 0; i < candidates.size(); ++ i) { //finding which one was the highest pT match
      if (matchedToLead[0].delta_R(candidates[i]) < 0.0001) { //is probably the same jet
	match_position = i;
	return;
      }
    }
    return;
  }
  */
  //PLAN: get rid of "match_position". Pass two vectors instead, and have a third & fourth which are written to with the matches (tomatch & candidate). If size != 0, we have matches.
  //another consideration: need to be able to remove jets from the "candidates" vector after they've been matched. Make a copy prior to the function call? Feels risky. Do it in the function. Also make a copy candidates vector for each iteration on toMatch since it gets selected.
  //In finding which jets were the matches, we know the tomatch jet match will be the 'i'th jet since we are iterating. The candidate_copy jet should be the highest pT match, so the first one in the candidate_copy list. Geometrically match the candidate_copy jet to the nearest candidate jet.
  std::vector<int> MatchJets(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & c_matches, std::vector<fastjet::PseudoJet> & t_matches) {
    std::vector<int> match_indices;
    if (candidates_safe.size() == 0 || toMatch.size() == 0) {
      return match_indices;
    }
    std::vector<fastjet::PseudoJet> candidates = candidates_safe; 
    for (int i = 0; i < toMatch.size(); ++ i) {
      std::vector<fastjet::PseudoJet> candidates_copy = candidates;
      fastjet::Selector selectMatchedLead = fastjet::SelectorCircle( R );
      selectMatchedLead.set_reference( toMatch[i] );
      std::vector<fastjet::PseudoJet> matchedToLead = sorted_by_pt( selectMatchedLead( candidates_copy ));
      if (candidates_copy.size() == 0) { continue; } //means no match to this jet. Remove none from candidates. Continuing on to the next one.
      else { //found at least one match. Need to remove the highest pT one from candidates and add the respective jets to the match vectors.
	match_indices.push_back(i); //for using the groomed jet corresponding to this jet, later
	t_matches.push_back(toMatch[i]);
	c_matches.push_back(candidates_copy[0]);
	for (int j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
	  if (candidates_copy[0].delta_R(candidates[j]) < 0.0001) { //is probably the same jet
	    candidates.erase(candidates.begin() + j); //removing the jet from the overall list of candidates so it can't be considered next time
	    match_indices.push_back(j);
	    break; //should exit only the c_matches loop.
	  }
	}
      }
    } 
    return match_indices;
  }

  //this function is similar to the "MatchJets" function, but it instead takes a list of jets which have already been matched ("candidates_safe") and finds the jets to which they correspond in the list of unmatched jets ("toMatch") by geometrical matching. I know, it is confusing to match matched jets to unmatched jets. But that's what we're doing. When we have jets which don't have a basically perfect geometrical match to a matched jet, we know that the jet in our hands is unmatched (to a detector/particle-level complement rather than to itself in the other list, basically). 
  std::vector<int> FakesandMisses(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & unmatched) {
    std::vector<int> miss_fake_index;
    std::vector<fastjet::PseudoJet> candidates = candidates_safe; 
    for (int i = 0; i < toMatch.size(); ++ i) {
      std::vector<fastjet::PseudoJet> candidates_copy = candidates;
      fastjet::Selector selectMatchedLead = fastjet::SelectorCircle( 0.0001 ); //a "match" is now if we found the same exact jet.
      selectMatchedLead.set_reference( toMatch[i] );
      std::vector<fastjet::PseudoJet> matchedToLead = sorted_by_pt( selectMatchedLead( candidates_copy ));
      if (candidates_copy.size() == 0) { //means no match to this jet. Remove none from candidates. Add it to unmatched & continue to next one.
	miss_fake_index.push_back(i);
	unmatched.push_back(toMatch[i]);
	continue;
      } 
      else { //found at least one match. Need to remove the highest pT one from candidates
	for (int j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
	  if (candidates_copy[0].delta_R(candidates[j]) < 0.0001) { //is probably the same jet
	    candidates.erase(candidates.begin() + j); //removing the jet from the overall list of candidates so it can't be considered next time
	    break; //should exit only the candidates loop.
	  }
	}
      }
    } 
    return miss_fake_index;
  }
  

  /*
  void ConstructResponse(RooUnfoldResponse & res, const std::vector<fastjet::PseudoJet> g_Jets, const std::vector<fastjet::PseudoJet> p_Jets, const double mc_weight) {
    if (p_Jets.size() != 0) {
      int position = -1;
      MatchJets(g_Jets, p_Jets[0], position);
      if (position == -1) { //didn't find a match                                                                                                   
	res.Miss(p_Jets[0].pt(), mc_weight);
	std::cout << "MISS " << p_Jets[0].pt() << std::endl;
      }
      else { //found a match                                                                                                                        
	res.Fill(g_Jets[position].pt(), p_Jets[0].pt(), mc_weight);
	std::cout << "MATCH " << p_Jets[0].pt() << " " << g_Jets[position].pt() << std::endl;
      }
    }

    //fake rate                                                                                                                                     
    if (g_Jets.size() != 0) {
      int position = -1;
      MatchJets(p_Jets, g_Jets[0], position);
      if (position == -1) { //didn't find a match                                                                                                   
	res.Fake(g_Jets[0].pt(), mc_weight);
	std::cout << "FAKE " << g_Jets[0].pt() << " " << std::endl;
      }
    }
    return;
  }
  */

  //new construction of responses based on ~inclusive~ matching. We now need to have two vectors to be filled with matched jets. If they aren't, when looping over pythia jets, we have misses. Fill all pythia jets into the misses. Same goes when looping over geant jets with fakes. And for matches, we just fill with however many entries there are in the matched vectors.
  //MatchJets now returns a vector of pairs of indices (i,j). The first entry is the candidate's position, the second its match's position, the third the next candidate's position, the fourth its match's position, etc.
  //FakesandMisses now returns a vector of indices (i) corresponding to the indices of misses or fakes from the original candidate vector.
  void ConstructResponses(std::vector<RooUnfoldResponse*> res, const std::vector<fastjet::PseudoJet> g_Jets, const std::vector<fastjet::PseudoJet> p_Jets, const std::vector<fastjet::PseudoJet> g_Groomed, const std::vector<fastjet::PseudoJet> p_Groomed, std::vector<std::vector<double> > & g_tree, std::vector<std::vector<double> > & p_tree, const double mc_weight/*, double & nEntries, double & nFakes, double & nMisses, double & nMatches*/) {
    std::vector<fastjet::PseudoJet> g_matches; std::vector<fastjet::PseudoJet> p_matches; 
    std::vector<fastjet::PseudoJet> fakes; std::vector<fastjet::PseudoJet> misses;
    std::vector<fastjet::PseudoJet> g_sd_matches; std::vector<fastjet::PseudoJet> p_sd_matches; 
    std::vector<fastjet::PseudoJet> sd_fakes; std::vector<fastjet::PseudoJet> sd_misses;
    std::vector<int> sd_match_indices;
    std::vector<int> sd_miss_indices; std::vector<int> sd_fake_indices;
    
    if (p_Jets.size() != 0) {
      g_matches.clear(); p_matches.clear(); g_sd_matches.clear(); p_sd_matches.clear(); misses.clear(); sd_misses.clear();
      sd_match_indices.clear(); sd_miss_indices.clear();
      sd_match_indices = MatchJets(g_Jets, p_Jets, g_matches, p_matches);
      if (sd_match_indices.size() != 0) {//means there is >= 1 match. Will index the g_vector with these match indices (and the p_vector with index 0). 
	for (int i = 0; i < sd_match_indices.size(); i += 2) {
	  p_sd_matches.push_back(p_Groomed[sd_match_indices[i]]);
	}
	for (int i = 1; i < sd_match_indices.size(); i += 2) {
	  g_sd_matches.push_back(g_Groomed[sd_match_indices[i]]);
	}
      }
      if (g_matches.size() != p_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. This should never happen!" <<std::endl; exit(1);}
      if (g_sd_matches.size() != g_matches.size()) {std::cerr << "There should be as many matched groomed jets as matched ungroomed jets!" <<std::endl; exit(1);}
      if (g_matches.size() < p_Jets.size()) { //then we have misses
	sd_miss_indices = FakesandMisses(p_matches, p_Jets, misses);
	for (int i = 0; i < misses.size(); ++ i) {
	  res[0]->Miss(misses[i].pt(), mc_weight);
	  res[1]->Miss(misses[i].m(), mc_weight);
	  res[2]->Miss(misses[i].m(), misses[i].pt(), mc_weight);
	  //	  nMisses ++;
	  //	  std::cout << "MISS " << misses[i].pt() << std::endl;
	}
	for (int i = 0; i < sd_miss_indices.size(); ++ i) {
	  res[3]->Miss(p_Groomed[sd_miss_indices[i]].structure_of<SD>().symmetry(), mc_weight);
	  res[4]->Miss(p_Groomed[sd_miss_indices[i]].structure_of<SD>().delta_R(), mc_weight);
	  res[5]->Miss(p_Groomed[sd_miss_indices[i]].pt(), mc_weight);
	  res[6]->Miss(p_Groomed[sd_miss_indices[i]].m(), mc_weight);
	  res[7]->Miss(p_Groomed[sd_miss_indices[i]].structure_of<SD>().symmetry(), p_Jets[sd_miss_indices[i]].pt(), mc_weight);
	  res[8]->Miss(p_Groomed[sd_miss_indices[i]].structure_of<SD>().delta_R(), p_Jets[sd_miss_indices[i]].pt(), mc_weight);
	  res[9]->Miss(p_Groomed[sd_miss_indices[i]].pt(), p_Jets[sd_miss_indices[i]].pt(), mc_weight);
	  res[10]->Miss(p_Groomed[sd_miss_indices[i]].m(), p_Jets[sd_miss_indices[i]].pt(), mc_weight);
	}
      }
      if (g_matches.size() != 0) { //found match(es)
	for (int i = 0; i < g_matches.size(); ++ i) {
	  res[0]->Fill(g_matches[i].pt(), p_matches[i].pt(), mc_weight); //matches should be at same index in respective vectors
	  res[1]->Fill(g_matches[i].m(), p_matches[i].m(), mc_weight);
	  res[2]->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), mc_weight);
	  //  nMatches ++; nEntries ++;
	  g_tree[0].push_back(g_matches[i].pt());
	  g_tree[1].push_back(g_matches[i].m());
	  g_tree[6].push_back(g_matches[i].eta());
	  
	  p_tree[0].push_back(p_matches[i].pt());
	  p_tree[1].push_back(p_matches[i].m()); 
	  p_tree[6].push_back(p_matches[i].eta());
	  //std::cout << "MATCH p:" << p_matches[i].pt() << " g:" << g_matches[i].pt() << std::endl;
	}
	for (int i = 0; i < g_sd_matches.size(); ++ i) {
	  res[3]->Fill(g_sd_matches[i].structure_of<SD>().symmetry(), p_sd_matches[i].structure_of<SD>().symmetry(),  mc_weight);
	  res[4]->Fill(g_sd_matches[i].structure_of<SD>().delta_R(), p_sd_matches[i].structure_of<SD>().delta_R(),  mc_weight);	  
	  res[5]->Fill(g_sd_matches[i].pt(), p_sd_matches[i].pt(), mc_weight);
	  res[6]->Fill(g_sd_matches[i].m(), p_sd_matches[i].m(), mc_weight);
	  res[7]->Fill(g_sd_matches[i].structure_of<SD>().symmetry(), g_matches[i].pt(), p_sd_matches[i].structure_of<SD>().symmetry(), p_matches[i].pt(), mc_weight);
	  res[8]->Fill(g_sd_matches[i].structure_of<SD>().delta_R(), g_matches[i].pt(), p_sd_matches[i].structure_of<SD>().delta_R(), p_matches[i].pt(), mc_weight);
	  res[9]->Fill(g_sd_matches[i].pt(), g_matches[i].pt(), p_sd_matches[i].pt(), p_matches[i].pt(), mc_weight);
	  res[10]->Fill(g_sd_matches[i].m(), g_matches[i].pt(), p_sd_matches[i].m(), p_matches[i].pt(), mc_weight);
	
	  g_tree[2].push_back(g_sd_matches[i].structure_of<SD>().symmetry());
	  g_tree[3].push_back(g_sd_matches[i].structure_of<SD>().delta_R());
	  g_tree[4].push_back(g_sd_matches[i].pt());
	  g_tree[5].push_back(g_sd_matches[i].m());
	  
	  p_tree[2].push_back(p_sd_matches[i].structure_of<SD>().symmetry());
	  p_tree[3].push_back(p_sd_matches[i].structure_of<SD>().delta_R());
	  p_tree[4].push_back(p_sd_matches[i].pt());
	  p_tree[5].push_back(p_sd_matches[i].m());
	}
	
      }
    }

    //fake rate                                                
    if (g_Jets.size() != 0) {
      g_matches.clear(); p_matches.clear(); g_sd_matches.clear(); p_sd_matches.clear(); fakes.clear(); sd_fakes.clear();
      sd_match_indices.clear(); sd_fake_indices.clear();
      MatchJets(p_Jets, g_Jets, p_matches, g_matches);
      if (g_matches.size() != p_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. This should never happen!" <<std::endl; exit(1);}
      if (p_matches.size() < g_Jets.size()) { //then we have fakes
	sd_fake_indices = FakesandMisses(g_matches, g_Jets, fakes);
	for (int i = 0; i < fakes.size(); ++ i) {
	  res[0]->Fake(fakes[i].pt(), mc_weight);
	  res[1]->Fake(fakes[i].m(), mc_weight);
	  res[2]->Fake(fakes[i].m(), fakes[i].pt(), mc_weight);
	  //nFakes ++; nEntries ++;
	  //std::cout << "FAKE " << fakes[i].pt() << " " << std::endl;
	}
	for (int i = 0; i < sd_fake_indices.size(); ++ i) {
	  res[3]->Fake(g_Groomed[sd_fake_indices[i]].structure_of<SD>().symmetry(), mc_weight);
	  res[4]->Fake(g_Groomed[sd_fake_indices[i]].structure_of<SD>().delta_R(), mc_weight);
	  res[5]->Fake(g_Groomed[sd_fake_indices[i]].pt(), mc_weight);
	  res[6]->Fake(g_Groomed[sd_fake_indices[i]].m(), mc_weight);
	  res[7]->Fake(g_Groomed[sd_fake_indices[i]].structure_of<SD>().symmetry(), g_Jets[sd_fake_indices[i]].pt(), mc_weight);
	  res[8]->Fake(g_Groomed[sd_fake_indices[i]].structure_of<SD>().delta_R(), g_Jets[sd_fake_indices[i]].pt(), mc_weight);
	  res[9]->Fake(g_Groomed[sd_fake_indices[i]].pt(), g_Jets[sd_fake_indices[i]].pt(), mc_weight);
	  res[10]->Fake(g_Groomed[sd_fake_indices[i]].m(), g_Jets[sd_fake_indices[i]].pt(), mc_weight);
	}
      }
    }
    return;
  }
  
  
  void FillResponses(RooUnfoldResponse & pt_res_coarse, RooUnfoldResponse & pt_response, RooUnfoldResponse & m_response, RooUnfoldResponse & pt_m_response, const fastjet::PseudoJet g_jet, const fastjet::PseudoJet p_jet, const double mc_weight) {
    pt_res_coarse.Fill(g_jet.pt(), p_jet.pt(), mc_weight);
    pt_response.Fill(g_jet.pt(), p_jet.pt(), mc_weight);
    m_response.Fill(g_jet.m(), p_jet.m(), mc_weight);
    pt_m_response.Fill(g_jet.m(), g_jet.pt(), p_jet.m(), p_jet.pt(), mc_weight);
    return;
  }
  
  void FillResponsesSD(RooUnfoldResponse & zg_response, RooUnfoldResponse & rg_response, RooUnfoldResponse & ptg_response, RooUnfoldResponse & mg_response, RooUnfoldResponse & pt_zg_response, RooUnfoldResponse & pt_rg_response, RooUnfoldResponse & pt_ptg_response, RooUnfoldResponse & pt_mg_response, const fastjet::PseudoJet g_jet, const fastjet::PseudoJet p_jet, const fastjet::PseudoJet g_ungroom, const fastjet::PseudoJet p_ungroom, const double mc_weight) {
    zg_response.Fill(g_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), p_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), mc_weight);
    rg_response.Fill(g_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), p_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), mc_weight);
    ptg_response.Fill(g_jet.pt(), p_jet.pt(), mc_weight);
    mg_response.Fill(g_jet.m(), p_jet.m(), mc_weight);
    pt_zg_response.Fill(g_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), g_ungroom.pt(), p_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), p_ungroom.pt(), mc_weight);
    pt_rg_response.Fill(g_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), g_ungroom.pt(), p_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), p_ungroom.pt(), mc_weight);
    pt_ptg_response.Fill(g_jet.pt(), g_ungroom.pt(), p_jet.pt(), p_ungroom.pt(), mc_weight);
    pt_mg_response.Fill(g_jet.m(), g_ungroom.pt(), p_jet.m(), p_ungroom.pt(), mc_weight);
    
    return;
  }

  void Misses(RooUnfoldResponse & pt_res_coarse, RooUnfoldResponse & pt_response, RooUnfoldResponse & m_response, RooUnfoldResponse & pt_m_response, const fastjet::PseudoJet p_jet, const double mc_weight) {
    pt_res_coarse.Miss(p_jet.pt(), mc_weight);
    pt_response.Miss(p_jet.pt(), mc_weight);
    m_response.Miss(p_jet.m(), mc_weight);
    pt_m_response.Miss(p_jet.m(), p_jet.pt(), mc_weight);
    return;
  }
  
  void MissesSD(RooUnfoldResponse & zg_response, RooUnfoldResponse & rg_response, RooUnfoldResponse & ptg_response, RooUnfoldResponse & mg_response, RooUnfoldResponse & pt_zg_response, RooUnfoldResponse & pt_rg_response, RooUnfoldResponse & pt_ptg_response, RooUnfoldResponse & pt_mg_response, const fastjet::PseudoJet p_jet, const fastjet::PseudoJet p_ungroom, const double mc_weight) {
    zg_response.Miss(p_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), mc_weight);
    rg_response.Miss(p_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), mc_weight);
    ptg_response.Miss(p_jet.pt(), mc_weight);
    mg_response.Miss(p_jet.m(), mc_weight);
    pt_zg_response.Miss(p_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), p_ungroom.pt(), mc_weight);
    pt_rg_response.Miss(p_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), p_ungroom.pt(), mc_weight);
    pt_ptg_response.Miss(p_jet.pt(), p_ungroom.pt(), mc_weight);
    pt_mg_response.Miss(p_jet.m(), p_ungroom.pt(), mc_weight);
    
    return;
  }
  
  void Fakes(RooUnfoldResponse & pt_res_coarse, RooUnfoldResponse & pt_response, RooUnfoldResponse & m_response, RooUnfoldResponse & pt_m_response, const fastjet::PseudoJet g_jet, const double mc_weight) {
    pt_res_coarse.Miss(g_jet.pt(), mc_weight);
    pt_response.Miss(g_jet.pt(), mc_weight);
    m_response.Miss(g_jet.m(), mc_weight);
    pt_m_response.Miss(g_jet.m(), g_jet.pt(), mc_weight);
    return;
  }
  
  void FakesSD(RooUnfoldResponse & zg_response, RooUnfoldResponse & rg_response, RooUnfoldResponse & ptg_response, RooUnfoldResponse & mg_response, RooUnfoldResponse & pt_zg_response, RooUnfoldResponse & pt_rg_response, RooUnfoldResponse & pt_ptg_response, RooUnfoldResponse & pt_mg_response, const fastjet::PseudoJet g_jet, const fastjet::PseudoJet g_ungroom, const double mc_weight) {
    zg_response.Fake(g_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), mc_weight);
    rg_response.Fake(g_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), mc_weight);
    ptg_response.Fake(g_jet.pt(), mc_weight);
    mg_response.Fake(g_jet.m(), mc_weight);
    pt_zg_response.Fake(g_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), g_ungroom.pt(), mc_weight);
    pt_rg_response.Fake(g_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), g_ungroom.pt(), mc_weight);
    pt_ptg_response.Fake(g_jet.pt(), g_ungroom.pt(), mc_weight);
    pt_mg_response.Fake(g_jet.m(), g_ungroom.pt(), mc_weight);
    
    return;
  }
  
   
  
  //  INITIATE READER
  void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents, const std::string trig, const double vZ, const double vZDiff, const double Pt, const double Et, const double Etmin,  const double DCA, const double NFit, const double NFitRatio, const double maxEtTow, const std::string badTows, const std::string bad_run_list) {
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( true );
    reader.SetFractionHadronicCorrection( 0.9999 );
    reader.SetRejectTowerElectrons( kFALSE );

    // if bad run list is specified, add to reader
    
    if (!bad_run_list.empty()) {
      std::set<int> bad_runs = ParseCSV<int>(bad_run_list);
      for (auto run : bad_runs) {
	reader.AddMaskedRun(run);
      }
    }
        
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetTriggerSelection( trig.c_str() ); //All, MB, HT, pp, ppHT, ppJP
    evCuts->SetVertexZCut ( vZ );
    evCuts->SetVertexZDiffCut( vZDiff );
    evCuts->SetRefMultCut( refMultCut );
    evCuts->SetMaxEventPtCut( Pt );
    evCuts->SetMaxEventEtCut( Et );
    evCuts->SetMinEventEtCut( Etmin );

    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( DCA );
    trackCuts->SetMinNFitPointsCut( NFit );
    trackCuts->SetFitOverMaxPointsCut( NFitRatio );    
    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( maxEtTow );
    towerCuts->AddBadTowers( badTows );
    std::cout << badTows << std::endl;

    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1
  }
  
  //discards events on the grounds of them having jets of pT > double the high end of the pT-hat bin from which they came. Both the Py & Py+Ge event will be thrown out.
  bool DiscardEvent(const TString Filename, const std::vector<fastjet::PseudoJet> p_Jets, const std::vector<fastjet::PseudoJet> g_Jets) {
    bool bad_event = 0;
    std::string tail = ((std::string) Filename).substr(((std::string) Filename).size() - 10);
    std::string upstring = tail.substr(0,2);
    std::string upstring_copy = upstring;
    if (upstring.find("_") != std::string::npos || upstring.find("-") != std::string::npos) { if (upstring.substr(1,1) != "_") {upstring = upstring.substr(1,1);} else {upstring = upstring.substr(0,1);}}
    int upbin = std::stoi(upstring);
    if (p_Jets.size() != 0) {
      if ((p_Jets[0].pt() > 2*upbin) && upstring_copy != "-1") {
	std::cout << "due to a pythia event from " << Filename << " removing event " /*<< p_EventID << " with weight " << mc_weight*/ << " and bad jet with pt, eta, phi, and m: " << p_Jets[0].pt() << " " << p_Jets[0].eta() << " " << p_Jets[0].phi() << " " << p_Jets[0].m() << std::endl;
	bad_event = 1;
      }
    }
    if (g_Jets.size() != 0) {
      if ((g_Jets[0].pt() > 2*upbin) && upstring_copy != "-1") {
	std::cout << "due to a geant event from " << Filename << " removing event " /*<< g_EventID << " with weight " << mc_weight*/ << " and bad jet with pt, eta, phi, and m: " << g_Jets[0].pt() << " " << g_Jets[0].eta() << " " << g_Jets[0].phi() << " " << g_Jets[0].m() << std::endl;
	bad_event = 1;
      }
    }
    
    return bad_event;
  }



  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILL HISTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  void FillHistsHelper(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> &c2D, Collection<std::string, TH3D> & c3D, const std::string flag, const fastjet::PseudoJet jet, const double weight) {
    c1D.fill(("m_" + flag + "_" + "jet").c_str(), jet.m(), weight);
    c2D.fill(("m_v_pt_" + flag + "_" + "jet").c_str(), jet.m(), jet.pt(), weight);
    c2D.fill(("m_v_pt_rebin_" + flag + "_" + "jet").c_str(), jet.m(), jet.pt(), weight);
    c3D.fill(("PtEtaPhi_" + flag + "_" + "jet").c_str(), jet.pt(), jet.eta(), jet.phi(), weight);
    for (int cons = 0; cons < jet.constituents().size(); ++ cons) {
      if (jet.constituents()[cons].pt() < partMinPt) {continue;} //ignores contributions from ghosts                       
      c3D.fill(("PtEtaPhi_" + flag + "_" + "cons").c_str(), jet.constituents()[cons].pt(), jet.constituents()[cons].eta(), jet.constituents()[cons].phi(), weight);
    }
    return;
  }
  
  void FillHists(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> & c2D, Collection<std::string, TH3D> & c3D, const std::vector<fastjet::PseudoJet> jets, const double weight) {
    //leading
    FillHistsHelper(c1D, c2D, c3D, "lead", jets[0], weight);
    //subleading
    if (jets.size() > 1) {
      FillHistsHelper(c1D, c2D, c3D, "sublead", jets[1], weight);
    }
    //inclusive
    for (int i = 0; i < jets.size(); ++ i) {
      FillHistsHelper(c1D, c2D, c3D, "incl", jets[i], weight);
    }
    /*
    //trigger & recoil
    std::vector<fastjet::PseudoJet> candidates;
    bool which_one = GetTriggerJet(candidates, jets);
    if (candidates.size() == 0 || jets.size() < 2) { // means there isn't a trigger or there isn't a recoil
      return;
    }
    if (candidates.size() == 1 && jets.size() > 1) { //potential trigger
      if (fabs(fabs(jets[which_one].delta_phi_to(jets[(which_one + 1) % 2])) - Pi) < R) { //found a recoil
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "trig", jets[which_one], weight); //filling hists for trigger
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "rec", jets[(which_one + 1) % 2], weight); //filling hists for recoil
	return;
      }
    }
    if (candidates.size() == 2) {
      if (fabs(fabs(candidates[0].delta_phi_to(candidates[1])) - Pi) < R) { //trigger & recoil found!
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "trig", candidates[0], weight); //filling hists for trigger
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "rec", candidates[1], weight); //filling hists for recoil
	return;
      }
    }
    */
    return;
  }

  void FillSDHistsHelper(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> &c2D, Collection<std::string, TH3D> & c3D, const std::string flag, const fastjet::PseudoJet jet, const double weight) {
    c1D.fill(("m_" + flag + "_" + "sd").c_str(), jet.m(), weight);
    c1D.fill(("zg_" + flag + "_" + "sd").c_str(), jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), weight);
    c1D.fill(("thetag_" + flag + "_" + "sd").c_str(), jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), weight);
    c2D.fill(("m_v_pt_" + flag + "_" + "sd").c_str(), jet.m(), jet.pt(), weight);
    c3D.fill(("PtEtaPhi_" + flag + "_" + "sd").c_str(), jet.pt(), jet.eta(), jet.phi(), weight);
    return;
  }
   

  void FillSDHists(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> & c2D, Collection<std::string, TH3D> & c3D, const std::vector<fastjet::PseudoJet> jets, const double weight) {
    //leading
    FillSDHistsHelper(c1D, c2D, c3D, "lead", jets[0], weight);
    //subleading
    if (jets.size() > 1) {
      FillSDHistsHelper(c1D, c2D, c3D, "sublead", jets[1], weight);
    }
    //inclusive
    for (int i = 0; i < jets.size(); ++ i) {
      FillSDHistsHelper(c1D, c2D, c3D, "incl", jets[i], weight);
    }

  }

}
