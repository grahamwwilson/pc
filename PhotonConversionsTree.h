//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec  5 13:42:40 2019 by ROOT version 6.18/00
// from TTree PhotonConversionsTree/PhotonConversionsTree
// found on file: Run2018_1210_Copy.root
//////////////////////////////////////////////////////////

#ifndef PhotonConversionsTree_h
#define PhotonConversionsTree_h

#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using namespace std;

class PhotonConversionsTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Bool_t          isRealData;
   UInt_t          eventNumber;
   UInt_t          runNumber;
   UInt_t          lumiSection;
   UInt_t          numberOfPV;
   vector<double>  *PV_x;
   vector<double>  *PV_y;
   vector<double>  *PV_z;
   vector<double>  *PV_xError;
   vector<double>  *PV_yError;
   vector<double>  *PV_zError;
   vector<bool>    *PV_isFake;
   UInt_t          numberOfMC_PUInfo;
   vector<unsigned int> *MC_PUInfo_bunchCrossing;
   vector<unsigned int> *MC_PUInfo_numberOfInteractions;
   Double_t        BS_x;
   Double_t        BS_y;
   Double_t        BS_z;
   Double_t        BS_zSigma;
   Double_t        BS_dxdy;
   Double_t        BS_dydz;
   Double_t        BS_xWidth;
   Double_t        BS_yWidth;
   UInt_t          numberOfMC_TrkV;
   vector<bool>    *MC_TrkV_isNuclearInteraction;
   vector<bool>    *MC_TrkV_isKaonDecay;
   vector<bool>    *MC_TrkV_isConversion;
   vector<double>  *MC_TrkV_x;
   vector<double>  *MC_TrkV_y;
   vector<double>  *MC_TrkV_z;
   vector<double>  *MC_TrkV_momentumInc_pt;
   vector<double>  *MC_TrkV_Inc_charge;
   vector<int>     *MC_TrkV_Inc_pdgId;
   vector<double>  *MC_TrkV_momentumInc_phi;
   vector<double>  *MC_TrkV_momentumInc_theta;
   vector<double>  *MC_TrkV_momentumOut_pt;
   vector<double>  *MC_TrkV_momentumOut_phi;
   vector<double>  *MC_TrkV_momentumOut_theta;
   vector<double>  *MC_TrkV_momentumOut_mass;
   vector<unsigned int> *MC_TrkV_numberOfChargedParticles_0p2;
   vector<unsigned int> *MC_TrkV_numberOfChargedParticles_0p5;
   vector<unsigned int> *MC_TrkV_numberOfChargedParticles_1p0;
   vector<unsigned int> *MC_TrkV_numberOfChargedParticles_Out0p2;
   vector<unsigned int> *MC_TrkV_numberOfChargedParticles_Out0p5;
   vector<unsigned int> *MC_TrkV_numberOfChargedParticles_Out1p0;
   vector<bool>    *MC_TrkV_isAssociatedPF;
   vector<unsigned int> *MC_TrkV_associationPCIdx;
   vector<double>  *MC_TrkV_associationPC_deltaR2d;
   vector<double>  *MC_TrkV_associationPC_deltaR3d;
   vector<double>  *MC_TrkV_associationPC_deltaR3dPerpendicular;
   vector<double>  *MC_TrkV_associationPC_deltaR3dParallel;
   UInt_t          numberOfPC;
   vector<double>  *PC_x;
   vector<double>  *PC_y;
   vector<double>  *PC_z;
   vector<double>  *PC_momentumOut_pt;
   vector<double>  *PC_momentumOut_phi;
   vector<double>  *PC_momentumOut_theta;
   vector<unsigned int> *PC_momentumOut_numberOfTracks;
   vector<double>  *PC_fitmomentumOut_pt;
   vector<double>  *PC_fitmomentumOut_phi;
   vector<double>  *PC_fitmomentumOut_theta;
   vector<double>  *PC_fitmomentumOut_mass;
   vector<double>  *PC_pairInvariantMass;
   vector<double>  *PC_pairCotThetaSeparation;
   vector<double>  *PC_distOfMinimumApproach;
   vector<double>  *PC_dPhiTracksAtVtx;
   vector<double>  *PC_vtx_chi2;
   vector<double>  *PC_vtx_ndof;
   vector<double>  *PC_vtx_normalizedChi2;
   vector<double>  *PC_vtx_sigmaxx;
   vector<double>  *PC_vtx_sigmayy;
   vector<double>  *PC_vtx_sigmazz;
   vector<double>  *PC_vtx_sigmaxy;
   vector<double>  *PC_vtx_sigmaxz;
   vector<double>  *PC_vtx_sigmayz;
   vector<vector<int> > *PC_vTrack_algo;
   vector<vector<int> > *PC_vTrack_charge;
   vector<vector<double> > *PC_vTrack_pt;
   vector<vector<double> > *PC_vTrack_eta;
   vector<vector<double> > *PC_vTrack_phi;
   vector<vector<double> > *PC_vTrack_chi2;
   vector<vector<double> > *PC_vTrack_normalizedChi2;
   vector<vector<double> > *PC_vTrack_rho;
   vector<vector<unsigned int> > *PC_vTrack_numberOfValidHits;
   vector<vector<unsigned int> > *PC_vTrack_numberOfExpectedOuterHits;
   vector<vector<unsigned int> > *PC_vTrack_closestDxyPVIdx;
   vector<vector<double> > *PC_vTrack_closestDxyPVIdx_dxy;
   vector<vector<double> > *PC_vTrack_closestDxyPVIdx_dz;
   vector<vector<unsigned int> > *PC_vTrack_closestDzPVIdx;
   vector<vector<double> > *PC_vTrack_closestDzPVIdx_dxy;
   vector<vector<double> > *PC_vTrack_closestDzPVIdx_dz;
   vector<vector<double> > *PC_fTrack_pt;
   vector<vector<double> > *PC_fTrack_eta;
   vector<vector<double> > *PC_fTrack_phi;

   // List of branches
   TBranch        *b_isRealData;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_numberOfPV;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_xError;   //!
   TBranch        *b_PV_yError;   //!
   TBranch        *b_PV_zError;   //!
   TBranch        *b_PV_isFake;   //!
   TBranch        *b_numberOfMC_PUInfo;   //!
   TBranch        *b_MC_PUInfo_bunchCrossing;   //!
   TBranch        *b_MC_PUInfo_numberOfInteractions;   //!
   TBranch        *b_BS_x;   //!
   TBranch        *b_BS_y;   //!
   TBranch        *b_BS_z;   //!
   TBranch        *b_BS_zSigma;   //!
   TBranch        *b_BS_dxdy;   //!
   TBranch        *b_BS_dydz;   //!
   TBranch        *b_BS_xWidth;   //!
   TBranch        *b_BS_yWidth;   //!
   TBranch        *b_numberOfMC_TrkV;   //!
   TBranch        *b_MC_TrkV_isNuclearInteraction;   //!
   TBranch        *b_MC_TrkV_isKaonDecay;   //!
   TBranch        *b_MC_TrkV_isConversion;   //!
   TBranch        *b_MC_TrkV_x;   //!
   TBranch        *b_MC_TrkV_y;   //!
   TBranch        *b_MC_TrkV_z;   //!
   TBranch        *b_MC_TrkV_momentumInc_pt;   //!
   TBranch        *b_MC_TrkV_Inc_charge;   //!
   TBranch        *b_MC_TrkV_Inc_pdgId;   //!
   TBranch        *b_MC_TrkV_momentumInc_phi;   //!
   TBranch        *b_MC_TrkV_momentumInc_theta;   //!
   TBranch        *b_MC_TrkV_momentumOut_pt;   //!
   TBranch        *b_MC_TrkV_momentumOut_phi;   //!
   TBranch        *b_MC_TrkV_momentumOut_theta;   //!
   TBranch        *b_MC_TrkV_momentumOut_mass;   //!
   TBranch        *b_MC_TrkV_numberOfChargedParticles_0p2;   //!
   TBranch        *b_MC_TrkV_numberOfChargedParticles_0p5;   //!
   TBranch        *b_MC_TrkV_numberOfChargedParticles_1p0;   //!
   TBranch        *b_MC_TrkV_numberOfChargedParticles_Out0p2;   //!
   TBranch        *b_MC_TrkV_numberOfChargedParticles_Out0p5;   //!
   TBranch        *b_MC_TrkV_numberOfChargedParticles_Out1p0;   //!
   TBranch        *b_MC_TrkV_isAssociatedPF;   //!
   TBranch        *b_MC_TrkV_associationPCIdx;   //!
   TBranch        *b_MC_TrkV_associationPC_deltaR2d;   //!
   TBranch        *b_MC_TrkV_associationPC_deltaR3d;   //!
   TBranch        *b_MC_TrkV_associationPC_deltaR3dPerpendicular;   //!
   TBranch        *b_MC_TrkV_associationPC_deltaR3dParallel;   //!
   TBranch        *b_numberOfPC;   //!
   TBranch        *b_PC_x;   //!
   TBranch        *b_PC_y;   //!
   TBranch        *b_PC_z;   //!
   TBranch        *b_PC_momentumOut_pt;   //!
   TBranch        *b_PC_momentumOut_phi;   //!
   TBranch        *b_PC_momentumOut_theta;   //!
   TBranch        *b_PC_momentumOut_numberOfTracks;   //!
   TBranch        *b_PC_fitmomentumOut_pt;   //!
   TBranch        *b_PC_fitmomentumOut_phi;   //!
   TBranch        *b_PC_fitmomentumOut_theta;   //!
   TBranch        *b_PC_fitmomentumOut_mass;   //!
   TBranch        *b_PC_pairInvariantMass;   //!
   TBranch        *b_PC_pairCotThetaSeparation;   //!
   TBranch        *b_PC_distOfMinimumApproach;   //!
   TBranch        *b_PC_dPhiTracksAtVtx;   //!
   TBranch        *b_PC_vtx_chi2;   //!
   TBranch        *b_PC_vtx_ndof;   //!
   TBranch        *b_PC_vtx_normalizedChi2;   //!
   TBranch        *b_PC_vtx_sigmaxx;   //!
   TBranch        *b_PC_vtx_sigmayy;   //!
   TBranch        *b_PC_vtx_sigmazz;   //!
   TBranch        *b_PC_vtx_sigmaxy;   //!
   TBranch        *b_PC_vtx_sigmaxz;   //!
   TBranch        *b_PC_vtx_sigmayz;   //!
   TBranch        *b_PC_vTrack_algo;   //!
   TBranch        *b_PC_vTrack_charge;   //!
   TBranch        *b_PC_vTrack_pt;   //!
   TBranch        *b_PC_vTrack_eta;   //!
   TBranch        *b_PC_vTrack_phi;   //!
   TBranch        *b_PC_vTrack_chi2;   //!
   TBranch        *b_PC_vTrack_normalizedChi2;   //!
   TBranch        *b_PC_vTrack_rho;   //!
   TBranch        *b_PC_vTrack_numberOfValidHits;   //!
   TBranch        *b_PC_vTrack_numberOfExpectedOuterHits;   //!
   TBranch        *b_PC_vTrack_closestDxyPVIdx;   //!
   TBranch        *b_PC_vTrack_closestDxyPVIdx_dxy;   //!
   TBranch        *b_PC_vTrack_closestDxyPVIdx_dz;   //!
   TBranch        *b_PC_vTrack_closestDzPVIdx;   //!
   TBranch        *b_PC_vTrack_closestDzPVIdx_dxy;   //!
   TBranch        *b_PC_vTrack_closestDzPVIdx_dz;   //!
   TBranch        *b_PC_fTrack_pt;   //!
   TBranch        *b_PC_fTrack_eta;   //!
   TBranch        *b_PC_fTrack_phi;   //!

   PhotonConversionsTree(TTree *tree=0);
   virtual ~PhotonConversionsTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PhotonConversionsTree_cxx
PhotonConversionsTree::PhotonConversionsTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Run2018_1210_Copy.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Run2018_1210_Copy.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Run2018_1210_Copy.root:/MyNtupleMaking");
      dir->GetObject("PhotonConversionsTree",tree);

   }
   Init(tree);
}

PhotonConversionsTree::~PhotonConversionsTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PhotonConversionsTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PhotonConversionsTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PhotonConversionsTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer

   std::cout << "Hello from Init " << std::endl;

   PV_x = 0;
   PV_y = 0;
   PV_z = 0;
   PV_xError = 0;
   PV_yError = 0;
   PV_zError = 0;
   PV_isFake = 0;
   MC_PUInfo_bunchCrossing = 0;
   MC_PUInfo_numberOfInteractions = 0;
   MC_TrkV_isNuclearInteraction = 0;
   MC_TrkV_isKaonDecay = 0;
   MC_TrkV_isConversion = 0;
   MC_TrkV_x = 0;
   MC_TrkV_y = 0;
   MC_TrkV_z = 0;
   MC_TrkV_momentumInc_pt = 0;
   MC_TrkV_Inc_charge = 0;
   MC_TrkV_Inc_pdgId = 0;
   MC_TrkV_momentumInc_phi = 0;
   MC_TrkV_momentumInc_theta = 0;
   MC_TrkV_momentumOut_pt = 0;
   MC_TrkV_momentumOut_phi = 0;
   MC_TrkV_momentumOut_theta = 0;
   MC_TrkV_momentumOut_mass = 0;
   MC_TrkV_numberOfChargedParticles_0p2 = 0;
   MC_TrkV_numberOfChargedParticles_0p5 = 0;
   MC_TrkV_numberOfChargedParticles_1p0 = 0;
   MC_TrkV_numberOfChargedParticles_Out0p2 = 0;
   MC_TrkV_numberOfChargedParticles_Out0p5 = 0;
   MC_TrkV_numberOfChargedParticles_Out1p0 = 0;
   MC_TrkV_isAssociatedPF = 0;
   MC_TrkV_associationPCIdx = 0;
   MC_TrkV_associationPC_deltaR2d = 0;
   MC_TrkV_associationPC_deltaR3d = 0;
   MC_TrkV_associationPC_deltaR3dPerpendicular = 0;
   MC_TrkV_associationPC_deltaR3dParallel = 0;
   PC_x = 0;
   PC_y = 0;
   PC_z = 0;
   PC_momentumOut_pt = 0;
   PC_momentumOut_phi = 0;
   PC_momentumOut_theta = 0;
   PC_momentumOut_numberOfTracks = 0;
   PC_fitmomentumOut_pt = 0;
   PC_fitmomentumOut_phi = 0;
   PC_fitmomentumOut_theta = 0;
   PC_fitmomentumOut_mass = 0;
   PC_pairInvariantMass = 0;
   PC_pairCotThetaSeparation = 0;
   PC_distOfMinimumApproach = 0;
   PC_dPhiTracksAtVtx = 0;
   PC_vtx_chi2 = 0;
   PC_vtx_ndof = 0;
   PC_vtx_normalizedChi2 = 0;
   PC_vtx_sigmaxx = 0;
   PC_vtx_sigmayy = 0;
   PC_vtx_sigmazz = 0;
   PC_vtx_sigmaxy = 0;
   PC_vtx_sigmaxz = 0;
   PC_vtx_sigmayz = 0;
   PC_vTrack_algo = 0;
   PC_vTrack_charge = 0;
   PC_vTrack_pt = 0;
   PC_vTrack_eta = 0;
   PC_vTrack_phi = 0;
   PC_vTrack_chi2 = 0;
   PC_vTrack_normalizedChi2 = 0;
   PC_vTrack_rho = 0;
   PC_vTrack_numberOfValidHits = 0;
   PC_vTrack_numberOfExpectedOuterHits = 0;
   PC_vTrack_closestDxyPVIdx = 0;
   PC_vTrack_closestDxyPVIdx_dxy = 0;
   PC_vTrack_closestDxyPVIdx_dz = 0;
   PC_vTrack_closestDzPVIdx = 0;
   PC_vTrack_closestDzPVIdx_dxy = 0;
   PC_vTrack_closestDzPVIdx_dz = 0;
   PC_fTrack_pt = 0;
   PC_fTrack_eta = 0;
   PC_fTrack_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetBranchStatus("*",0);   // Turn off all branches here
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("numberOfPV", &numberOfPV, &b_numberOfPV);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_xError", &PV_xError, &b_PV_xError);
   fChain->SetBranchAddress("PV_yError", &PV_yError, &b_PV_yError);
   fChain->SetBranchAddress("PV_zError", &PV_zError, &b_PV_zError);
   fChain->SetBranchAddress("PV_isFake", &PV_isFake, &b_PV_isFake);
   fChain->SetBranchAddress("numberOfMC_PUInfo", &numberOfMC_PUInfo, &b_numberOfMC_PUInfo);
   fChain->SetBranchAddress("MC_PUInfo_bunchCrossing", &MC_PUInfo_bunchCrossing, &b_MC_PUInfo_bunchCrossing);
   fChain->SetBranchAddress("MC_PUInfo_numberOfInteractions", &MC_PUInfo_numberOfInteractions, &b_MC_PUInfo_numberOfInteractions);
   fChain->SetBranchAddress("BS_x", &BS_x, &b_BS_x);
   fChain->SetBranchAddress("BS_y", &BS_y, &b_BS_y);
   fChain->SetBranchAddress("BS_z", &BS_z, &b_BS_z);
   fChain->SetBranchAddress("BS_zSigma", &BS_zSigma, &b_BS_zSigma);
   fChain->SetBranchAddress("BS_dxdy", &BS_dxdy, &b_BS_dxdy);
   fChain->SetBranchAddress("BS_dydz", &BS_dydz, &b_BS_dydz);
   fChain->SetBranchAddress("BS_xWidth", &BS_xWidth, &b_BS_xWidth);
   fChain->SetBranchAddress("BS_yWidth", &BS_yWidth, &b_BS_yWidth);
   fChain->SetBranchAddress("numberOfMC_TrkV", &numberOfMC_TrkV, &b_numberOfMC_TrkV);
   fChain->SetBranchAddress("MC_TrkV_isNuclearInteraction", &MC_TrkV_isNuclearInteraction, &b_MC_TrkV_isNuclearInteraction);
   fChain->SetBranchAddress("MC_TrkV_isKaonDecay", &MC_TrkV_isKaonDecay, &b_MC_TrkV_isKaonDecay);
   fChain->SetBranchAddress("MC_TrkV_isConversion", &MC_TrkV_isConversion, &b_MC_TrkV_isConversion);
   fChain->SetBranchAddress("MC_TrkV_x", &MC_TrkV_x, &b_MC_TrkV_x);
   fChain->SetBranchAddress("MC_TrkV_y", &MC_TrkV_y, &b_MC_TrkV_y);
   fChain->SetBranchAddress("MC_TrkV_z", &MC_TrkV_z, &b_MC_TrkV_z);
   fChain->SetBranchAddress("MC_TrkV_momentumInc_pt", &MC_TrkV_momentumInc_pt, &b_MC_TrkV_momentumInc_pt);
   fChain->SetBranchAddress("MC_TrkV_Inc_charge", &MC_TrkV_Inc_charge, &b_MC_TrkV_Inc_charge);
   fChain->SetBranchAddress("MC_TrkV_Inc_pdgId", &MC_TrkV_Inc_pdgId, &b_MC_TrkV_Inc_pdgId);
   fChain->SetBranchAddress("MC_TrkV_momentumInc_phi", &MC_TrkV_momentumInc_phi, &b_MC_TrkV_momentumInc_phi);
   fChain->SetBranchAddress("MC_TrkV_momentumInc_theta", &MC_TrkV_momentumInc_theta, &b_MC_TrkV_momentumInc_theta);
   fChain->SetBranchAddress("MC_TrkV_momentumOut_pt", &MC_TrkV_momentumOut_pt, &b_MC_TrkV_momentumOut_pt);
   fChain->SetBranchAddress("MC_TrkV_momentumOut_phi", &MC_TrkV_momentumOut_phi, &b_MC_TrkV_momentumOut_phi);
   fChain->SetBranchAddress("MC_TrkV_momentumOut_theta", &MC_TrkV_momentumOut_theta, &b_MC_TrkV_momentumOut_theta);
   fChain->SetBranchAddress("MC_TrkV_momentumOut_mass", &MC_TrkV_momentumOut_mass, &b_MC_TrkV_momentumOut_mass);
   fChain->SetBranchAddress("MC_TrkV_numberOfChargedParticles_0p2", &MC_TrkV_numberOfChargedParticles_0p2, &b_MC_TrkV_numberOfChargedParticles_0p2);
   fChain->SetBranchAddress("MC_TrkV_numberOfChargedParticles_0p5", &MC_TrkV_numberOfChargedParticles_0p5, &b_MC_TrkV_numberOfChargedParticles_0p5);
   fChain->SetBranchAddress("MC_TrkV_numberOfChargedParticles_1p0", &MC_TrkV_numberOfChargedParticles_1p0, &b_MC_TrkV_numberOfChargedParticles_1p0);
   fChain->SetBranchAddress("MC_TrkV_numberOfChargedParticles_Out0p2", &MC_TrkV_numberOfChargedParticles_Out0p2, &b_MC_TrkV_numberOfChargedParticles_Out0p2);
   fChain->SetBranchAddress("MC_TrkV_numberOfChargedParticles_Out0p5", &MC_TrkV_numberOfChargedParticles_Out0p5, &b_MC_TrkV_numberOfChargedParticles_Out0p5);
   fChain->SetBranchAddress("MC_TrkV_numberOfChargedParticles_Out1p0", &MC_TrkV_numberOfChargedParticles_Out1p0, &b_MC_TrkV_numberOfChargedParticles_Out1p0);
   fChain->SetBranchAddress("MC_TrkV_isAssociatedPF", &MC_TrkV_isAssociatedPF, &b_MC_TrkV_isAssociatedPF);
   fChain->SetBranchAddress("MC_TrkV_associationPCIdx", &MC_TrkV_associationPCIdx, &b_MC_TrkV_associationPCIdx);
   fChain->SetBranchAddress("MC_TrkV_associationPC_deltaR2d", &MC_TrkV_associationPC_deltaR2d, &b_MC_TrkV_associationPC_deltaR2d);
   fChain->SetBranchAddress("MC_TrkV_associationPC_deltaR3d", &MC_TrkV_associationPC_deltaR3d, &b_MC_TrkV_associationPC_deltaR3d);
   fChain->SetBranchAddress("MC_TrkV_associationPC_deltaR3dPerpendicular", &MC_TrkV_associationPC_deltaR3dPerpendicular, &b_MC_TrkV_associationPC_deltaR3dPerpendicular);
   fChain->SetBranchAddress("MC_TrkV_associationPC_deltaR3dParallel", &MC_TrkV_associationPC_deltaR3dParallel, &b_MC_TrkV_associationPC_deltaR3dParallel);
   fChain->SetBranchAddress("numberOfPC", &numberOfPC, &b_numberOfPC);
   fChain->SetBranchAddress("PC_x", &PC_x, &b_PC_x);
   fChain->SetBranchAddress("PC_y", &PC_y, &b_PC_y);
   fChain->SetBranchAddress("PC_z", &PC_z, &b_PC_z);
   fChain->SetBranchAddress("PC_momentumOut_pt", &PC_momentumOut_pt, &b_PC_momentumOut_pt);
   fChain->SetBranchAddress("PC_momentumOut_phi", &PC_momentumOut_phi, &b_PC_momentumOut_phi);
   fChain->SetBranchAddress("PC_momentumOut_theta", &PC_momentumOut_theta, &b_PC_momentumOut_theta);
   fChain->SetBranchAddress("PC_momentumOut_numberOfTracks", &PC_momentumOut_numberOfTracks, &b_PC_momentumOut_numberOfTracks);
   fChain->SetBranchAddress("PC_fitmomentumOut_pt", &PC_fitmomentumOut_pt, &b_PC_fitmomentumOut_pt);
   fChain->SetBranchAddress("PC_fitmomentumOut_phi", &PC_fitmomentumOut_phi, &b_PC_fitmomentumOut_phi);
   fChain->SetBranchAddress("PC_fitmomentumOut_theta", &PC_fitmomentumOut_theta, &b_PC_fitmomentumOut_theta);
   fChain->SetBranchAddress("PC_fitmomentumOut_mass", &PC_fitmomentumOut_mass, &b_PC_fitmomentumOut_mass);
   fChain->SetBranchAddress("PC_pairInvariantMass", &PC_pairInvariantMass, &b_PC_pairInvariantMass);
   fChain->SetBranchAddress("PC_pairCotThetaSeparation", &PC_pairCotThetaSeparation, &b_PC_pairCotThetaSeparation);
   fChain->SetBranchAddress("PC_distOfMinimumApproach", &PC_distOfMinimumApproach, &b_PC_distOfMinimumApproach);
   fChain->SetBranchAddress("PC_dPhiTracksAtVtx", &PC_dPhiTracksAtVtx, &b_PC_dPhiTracksAtVtx);
   fChain->SetBranchAddress("PC_vtx_chi2", &PC_vtx_chi2, &b_PC_vtx_chi2);
   fChain->SetBranchAddress("PC_vtx_ndof", &PC_vtx_ndof, &b_PC_vtx_ndof);
   fChain->SetBranchAddress("PC_vtx_normalizedChi2", &PC_vtx_normalizedChi2, &b_PC_vtx_normalizedChi2);
   fChain->SetBranchAddress("PC_vtx_sigmaxx", &PC_vtx_sigmaxx, &b_PC_vtx_sigmaxx);
   fChain->SetBranchAddress("PC_vtx_sigmayy", &PC_vtx_sigmayy, &b_PC_vtx_sigmayy);
   fChain->SetBranchAddress("PC_vtx_sigmazz", &PC_vtx_sigmazz, &b_PC_vtx_sigmazz);
   fChain->SetBranchAddress("PC_vtx_sigmaxy", &PC_vtx_sigmaxy, &b_PC_vtx_sigmaxy);
   fChain->SetBranchAddress("PC_vtx_sigmaxz", &PC_vtx_sigmaxz, &b_PC_vtx_sigmaxz);
   fChain->SetBranchAddress("PC_vtx_sigmayz", &PC_vtx_sigmayz, &b_PC_vtx_sigmayz);
   fChain->SetBranchAddress("PC_vTrack_algo", &PC_vTrack_algo, &b_PC_vTrack_algo);
   fChain->SetBranchAddress("PC_vTrack_charge", &PC_vTrack_charge, &b_PC_vTrack_charge);
   fChain->SetBranchAddress("PC_vTrack_pt", &PC_vTrack_pt, &b_PC_vTrack_pt);
   fChain->SetBranchAddress("PC_vTrack_eta", &PC_vTrack_eta, &b_PC_vTrack_eta);
   fChain->SetBranchAddress("PC_vTrack_phi", &PC_vTrack_phi, &b_PC_vTrack_phi);
   fChain->SetBranchAddress("PC_vTrack_chi2", &PC_vTrack_chi2, &b_PC_vTrack_chi2);
   fChain->SetBranchAddress("PC_vTrack_normalizedChi2", &PC_vTrack_normalizedChi2, &b_PC_vTrack_normalizedChi2);
   fChain->SetBranchAddress("PC_vTrack_rho", &PC_vTrack_rho, &b_PC_vTrack_rho);
   fChain->SetBranchAddress("PC_vTrack_numberOfValidHits", &PC_vTrack_numberOfValidHits, &b_PC_vTrack_numberOfValidHits);
   fChain->SetBranchAddress("PC_vTrack_numberOfExpectedOuterHits", &PC_vTrack_numberOfExpectedOuterHits, &b_PC_vTrack_numberOfExpectedOuterHits);
   fChain->SetBranchAddress("PC_vTrack_closestDxyPVIdx", &PC_vTrack_closestDxyPVIdx, &b_PC_vTrack_closestDxyPVIdx);
   fChain->SetBranchAddress("PC_vTrack_closestDxyPVIdx_dxy", &PC_vTrack_closestDxyPVIdx_dxy, &b_PC_vTrack_closestDxyPVIdx_dxy);
   fChain->SetBranchAddress("PC_vTrack_closestDxyPVIdx_dz", &PC_vTrack_closestDxyPVIdx_dz, &b_PC_vTrack_closestDxyPVIdx_dz);
   fChain->SetBranchAddress("PC_vTrack_closestDzPVIdx", &PC_vTrack_closestDzPVIdx, &b_PC_vTrack_closestDzPVIdx);
   fChain->SetBranchAddress("PC_vTrack_closestDzPVIdx_dxy", &PC_vTrack_closestDzPVIdx_dxy, &b_PC_vTrack_closestDzPVIdx_dxy);
   fChain->SetBranchAddress("PC_vTrack_closestDzPVIdx_dz", &PC_vTrack_closestDzPVIdx_dz, &b_PC_vTrack_closestDzPVIdx_dz);
   fChain->SetBranchAddress("PC_fTrack_pt", &PC_fTrack_pt, &b_PC_fTrack_pt);
   fChain->SetBranchAddress("PC_fTrack_eta", &PC_fTrack_eta, &b_PC_fTrack_eta);
   fChain->SetBranchAddress("PC_fTrack_phi", &PC_fTrack_phi, &b_PC_fTrack_phi);
   Notify();
}

Bool_t PhotonConversionsTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PhotonConversionsTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PhotonConversionsTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PhotonConversionsTree_cxx
