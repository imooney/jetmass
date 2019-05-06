// main32.cc is a part of the PYTHIA event generator.
// Copyright (C) 2016 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a sample program showing Alpgen- or Madgraph-style MLM matching
// for Madgraph LHEF or native Alpgen format event files.
//
// Please see the 'Jet Matching Style' manual page for a description of the
// parameters and user options.

// Includes and namespace
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/CombineMatchingInput.h"
using namespace Pythia8;

//==========================================================================

int main() {
 // Generator and read in commands.
  Pythia pythia;
  pythia.readFile("parton_matching.cmnd");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  cout << "Running over " << nEvent << " events" << endl;
  int nAbort = pythia.mode("Main:timesAllowErrors");
  int nSkip  = pythia.mode("Main:spareMode1");
  
  std::cout << "stat: " << std::endl; pythia.stat();
  
  // Create UserHooks pointer. Stop if it failed. Pass pointer to Pythia.
  CombineMatchingInput combined;
  UserHooks* matching = combined.getHook(pythia);
  if (!matching) return 1;
  pythia.setUserHooksPtr(matching);
  
  
  // Initialise Pythia.
  if (!pythia.init()) {
    std::cout << "Error: could not initialise Pythia" << std::endl;
    return 1;
  };

  // Optionally skip ahead in LHEF.
  pythia.LHAeventSkip( nSkip );

  
  // Begin event loop. Optionally quit it before end of file.
  int iAbort = 0;
  for (int iEvent = 0; ;  ++iEvent) {
    std::cout << "On event " << iEvent << " of " << nEvent << std::endl;
    if (nEvent > 0 && iEvent >= nEvent) break;
    pythia.stat();
    std::cout << "So we haven't finished yet" << std::endl;
    
    // Generate events. Quit if at end of file or many failures.
    if (!pythia.next()) {
      std::cout << "no next event" << std::endl;
      if (pythia.info.atEndOfFile()) {
        std::cout << "Info: end of input file reached" << std::endl;
        break;
      }
      if (++iAbort < nAbort) continue;
      std::cout << "Abort: too many errors in generation" << std::endl;
      break;
    }
    else {std:: cout << " should be fine " << std::endl;}
    
    // Event analysis goes here.

  // End of event loop.
  }

  // Final statistics and done.
  pythia.stat();
  delete matching;
  
  return 0;
}
