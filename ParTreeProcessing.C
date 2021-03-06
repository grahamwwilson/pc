#ifndef PARTREE
#define PARTREE

#include "TH1D.h"
#include "TH2D.h"
#include "ROOT/TThreadedObject.hxx"
#include "Math/GenVector/LorentzVector.h"

#include "ROOT/TTreeProcessorMT.hxx"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <vector>
#include <string>
#include "TFile.h"
#include "histset.C"
#include "myselector.C"

#include <iostream>

int main(int argc, char *argv[])
{

// Need to make it digest both the original trees 
// and the slimmed trees at least for now.
// Code assumes that all files that are included at any one time 
// have the same tree name.

// Original name
   std::string treestring = "MyNtupleMaking/PhotonConversionsTree";
   std::string filestring;
   std::string searchstring = "_R1.root";

   //set nthreads from first input argument
   int nthreads = std::atoi (argv[1]);
   // First enable implicit multi-threading globally, so that the implicit parallelisation is on.
   // The parameter of the call specifies the number of threads to use.
   ROOT::EnableImplicitMT(nthreads);

   //load up the list of files to be processed
   std::vector<std::string_view> ifilelist{};
   for(int i=2; i<argc; i++){
      filestring = argv[i];
      std::cout << "File " << argv[i] << endl;
      if(filestring.find(searchstring) != std::string::npos){
// searchstring found. Update the name for the _R1 slimmed versions.
         treestring = "PhotonConversionsTree";
      }
	  ifilelist.push_back(std::string_view(argv[i]));
   }
	
   //Create our ThreadedHistograms and Analysis Class
   histset h;
   ROOT::TTreeProcessorMT tp(ifilelist,treestring);

   // Define the function that will process a subrange of the tree.
   // The function must receive only one parameter, a TTreeReader,
   // and it must be thread safe. To enforce the latter requirement,
   // TThreadedObject histograms will be used.
   //
   auto myFunction = [&](TTreeReader &myReader) {

	  //create our values to be read
	  myselector s;

      //copy the tree onto the selector
	  s.Init(myReader.GetTree());

      while (myReader.Next()) {
	  //synchronize the threaded readed with private selector set of variables
	  s.fReader.SetEntry(myReader.GetCurrentEntry());
	  //analyze the current set of variables with hist set  class
	  h.AnalyzeEntry(s); 
      }
   };

   // Launch the parallel processing of the tree
   tp.Process(myFunction);
	
  //automatically do all merging and writing
  h.WriteHist();
	
   return 0;
}
#endif
