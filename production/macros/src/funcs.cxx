
#include "funcs.hh"

namespace Analysis {
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

  //new construction of responses based on ~inclusive~ matching. We now need to have two vectors to be filled with matched jets. If they aren't, when looping over pythia jets, we have misses. Fill all pythia jets into the misses. Same goes when looping over geant jets with fakes. And for matches, we just fill with however many entries there are in the matched vectors.
  //MatchJets now returns a vector of pairs of indices (i,j). The first entry is the candidate's position, the second its match's position, the third the next candidate's position, the fourth its match's position, etc.
  //FakesandMisses now returns a vector of indices (i) corresponding to the indices of misses or fakes from the original candidate vector.
  void ConstructResponses(std::vector<RooUnfoldResponse*> res, const std::vector<fastjet::PseudoJet> g_Jets, const std::vector<fastjet::PseudoJet> p_Jets, const double mc_weight) {
    std::vector<fastjet::PseudoJet> g_matches; std::vector<fastjet::PseudoJet> p_matches; 
    std::vector<fastjet::PseudoJet> fakes; std::vector<fastjet::PseudoJet> misses;
    
    if (p_Jets.size() != 0) {
      g_matches.clear(); p_matches.clear(); misses.clear();
      MatchJets(g_Jets, p_Jets, g_matches, p_matches);
      if (g_matches.size() != p_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. This should never happen!" <<std::endl; exit(1);}
      if (g_matches.size() < p_Jets.size()) { //then we have misses
	FakesandMisses(p_matches, p_Jets, misses);
	for (int i = 0; i < misses.size(); ++ i) {
	  res[0]->Miss(misses[i].m(), mc_weight);
	  //	  nMisses ++;
	  //	  std::cout << "MISS " << misses[i].pt() << std::endl;
	}
      }
      if (g_matches.size() != 0) { //found match(es)
	for (int i = 0; i < g_matches.size(); ++ i) {
	  res[0]->Fill(g_matches[i].m(), p_matches[i].m(), mc_weight); //matches should be at same index in respective vectors
	  //  nMatches ++; nEntries ++;
	  //std::cout << "MATCH p:" << p_matches[i].pt() << " g:" << g_matches[i].pt() << std::endl;
	}
      }
    }

    //fake rate                                                
    if (g_Jets.size() != 0) {
      g_matches.clear(); p_matches.clear(); fakes.clear();
      MatchJets(p_Jets, g_Jets, p_matches, g_matches);
      if (g_matches.size() != p_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. This should never happen!" <<std::endl; exit(1);}
      if (p_matches.size() < g_Jets.size()) { //then we have fakes
	FakesandMisses(g_matches, g_Jets, fakes);
	for (int i = 0; i < fakes.size(); ++ i) {
	  res[0]->Fake(fakes[i].m(), mc_weight);
	  //nFakes ++; nEntries ++;
	  //std::cout << "FAKE " << fakes[i].pt() << " " << std::endl;
	}
      }
    }
    return;
  }  

  void HistFromTree(TFile *onfile, TH1D* pt, TH1D* hist1520,TH1D* hist2025, TH1D* hist2530, TH1D* hist3040, TH1D* hist4060/*, TFile *offfile, std::vector<RooUnfoldResponse*> res_vec*/) {
    std::vector<double> *on_Pt = 0; std::vector<double> *on_M = 0;
    std::vector<double> *on_Eta = 0; std::vector<double> *on_Phi = 0;
    /*std::vector<double> *off_Pt = 0; std::vector<double> *off_M = 0;
    std::vector<double> *off_Eta = 0; std::vector<double> *off_Phi = 0;
    */
    double on_weight = 1;/* double off_weight = 1;*/ double on_eventID = 0; /*double off_eventID = 0;*/
    
    TTree *on = (TTree*) onfile->Get("ResultTree");
    //TTree *off = (TTree*) offfile->Get("ResultTree");
    
    on->SetBranchAddress("jetpT",&on_Pt); on->SetBranchAddress("jetM",&on_M);
    on->SetBranchAddress("jeteta",&on_Eta); on->SetBranchAddress("jetphi",&on_Phi);
    on->SetBranchAddress("eventID",&on_eventID); on->SetBranchAddress("mcweight",&on_weight);
    /*off->SetBranchAddress("jetpT",&off_Pt); off->SetBranchAddress("jetM",&off_M);
    off->SetBranchAddress("jeteta",&off_Eta); off->SetBranchAddress("jetphi",&off_Phi);
    off->SetBranchAddress("eventID",&off_eventID); off->SetBranchAddress("mcweight",&off_weight);
    */
    std::cout << on->GetEntries() << " entries in tree" << std::endl;
    //std::cout << off->GetEntries() << " entries in decays = off tree" << std::endl;
    
    //note to self: PY8 with decays off is the "detector-level". We loop over this to get the fakes. Loop over "decays on" to get matches and misses.                                         
    //another note: we are going to assume that the trees have the same number of entries and that each entry is the same event (we will check this last assumption explicitly).              
    
    for (int i = 0; i < on->GetEntries(); ++ i) {
      //if (i % 1000000 == 0) { cout << "still chuggin. " << i << endl;}                                                                                                                      
      on->GetEntry(i);
      
      for (int j = 0; j < on_M->size(); ++ j) {
	pt->Fill(on_Pt->at(j),on_weight);
	if (on_Pt->at(j) > 15 && on_Pt->at(j) < 20) {
	  hist1520->Fill(on_M->at(j),on_weight);
	}
	if (on_Pt->at(j) > 20 && on_Pt->at(j) < 25) {
	  hist2025->Fill(on_M->at(j),on_weight);
	}
	if (on_Pt->at(j) > 25 && on_Pt->at(j) < 30) {
	  hist2530->Fill(on_M->at(j),on_weight);
	}
	if (on_Pt->at(j) > 30 && on_Pt->at(j) < 40) {
	  hist3040->Fill(on_M->at(j),on_weight);
	}
	if (on_Pt->at(j) > 40 && on_Pt->at(j) < 60) {
	  hist4060->Fill(on_M->at(j),on_weight);
	}
      }
      /*off->GetEntry(i);
      
      std::cout << "Event IDs: " << on_eventID << " " << off_eventID << std::endl;
      if (on_eventID != off_eventID) {std::cerr << "We've mismatched events. Exiting!" << std::endl; exit(1);}
      //now we have the same event in each. Match the jets.                                               
      std::vector<fastjet::PseudoJet> ons; std::vector<fastjet::PseudoJet> offs;
      //fill the vectors of pseudojets.                                                
      for (int j = 0; j < on_Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size
	fastjet::PseudoJet on;
	on.reset_PtYPhiM(on_Pt->at(j), on_Eta->at(j), on_Phi->at(j), on_M->at(j));
	ons.push_back(on);
      }
      for (int j = 0; j < off_Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size
	fastjet::PseudoJet off;
	off.reset_PtYPhiM(off_Pt->at(j), off_Eta->at(j), off_Phi->at(j), off_M->at(j));
	offs.push_back(off);
      }
      if (off_weight != on_weight) {std::cout << "ON V OFF: " << on_weight << " " << off_weight << std::endl; std::cerr << "Weights differ! Should never happen for the same event! Exiting!" << std::endl; exit(1);}
      
      ConstructResponses(res_vec, offs, ons, off_weight);
*/
      
					     
    }
    on->ResetBranchAddresses(); //off->ResetBranchAddresses();
    return;
  }
  
}
