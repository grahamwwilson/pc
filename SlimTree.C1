// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <istream>
#include <stdio.h>
//#include <dirent.h>
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

	reducedTree->Write();
	ofile->Write();
	ofile->Close();
}

int main(){
   
   PhotonConversionsTree *tree = new PhotonConversionsTree();
 
   tree->fChain->Print();   

   produceReducedTree(*tree,"OutputTree.root");
  
}
