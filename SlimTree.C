// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <istream>
#include <stdio.h>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>
#include <TLeafElement.h>
#include <TLorentzVector.h>

using namespace std;
#include "PhotonConversionsTree.h"

void produceReducedTree(PhotonConversionsTree& selector, std::string ofilename){
    //copy branches to output file	
    auto ofile =  new TFile(ofilename.c_str(), "RECREATE");
    auto reducedTree = selector.fChain->CloneTree();
    reducedTree->Print();
    reducedTree->Write();
    ofile->Write();
    ofile->Close();
    cout << "output file written" << endl;
}

int main(int argc, char *argv[]){

   // Display each command-line argument.
   cout << "Command-line arguments: " << endl;
   for( int i = 0; i < argc; i++ )
        cout << "  argv[" << i << "]   "
                << argv[i] << endl;

   std::string filename = argv[1];
   std::string suffix = ".root";
   std::string version = "_R1";

   std::string ifilename = filename+suffix;
   std::string ofilename = filename+version+suffix;

   cout << "input file: " << ifilename << endl;
   cout << "output file: " << ofilename << endl;

   TChain* chain;
 
   chain = (TChain*) new TChain("MyNtupleMaking/PhotonConversionsTree");
   chain->Add(ifilename.c_str());
   chain->Print();
 
   PhotonConversionsTree *tree = new PhotonConversionsTree(chain);

   produceReducedTree(*tree,ofilename);
  
}
