//    #include "mylocaltree.h"     //All the variable incantations needed
   	
	//always make a local copy, if its a value dereference. 
    //if you dont do this scope/dereferencing will get really weird, clunky, and unmanageable
	//have to auto& or myreader will try to register copy of the readerarray pointer

// Here we retain the old variable names as far as possible to minimize required changes

	auto& PC_x = s.Conv_vtx_X;
	auto& PC_y = s.Conv_vtx_Y;
	auto& PC_z = s.Conv_vtx_Z;
	auto& PC_vtx_chi2 = s.Conv_vtx_chi2;
	auto& PC_vtx_ndof = s.Conv_vtx_ndof;

	auto& PC_vtx_sigmaxx = s.Conv_vtx_cov_00;
	auto& PC_vtx_sigmaxy = s.Conv_vtx_cov_01;
	auto& PC_vtx_sigmayy = s.Conv_vtx_cov_11;
	auto& PC_vtx_sigmazz = s.Conv_vtx_cov_22;

	auto numberOfPC = *(s.nConv);
	auto numberOfPV = *(s.nPV);

	auto& PC_vTrack0_pt = s.Conv_Tk0_pt;
	auto& PC_vTrack0_phi = s.Conv_Tk0_phi;
	auto& PC_vTrack0_eta = s.Conv_Tk0_eta;
	auto& PC_vTrack1_pt = s.Conv_Tk1_pt;
	auto& PC_vTrack1_phi = s.Conv_Tk1_phi;
	auto& PC_vTrack1_eta = s.Conv_Tk1_eta;

// New variables
    auto& PC_mpair = s.Conv_pairInvariantMass;
    auto& PC_dcottheta = s.Conv_pairCotThetaSeparation;
    auto& PC_dmin = s.Conv_distOfMinimumApproach;
    auto& PC_dphi = s.Conv_dPhiTracksAtVtx;

//    auto isRealData = *(s.isRealData);
    bool isRealData = true;

    auto lumiSection = *(s.luminosityBlock);
//    auto& mcpu = s.MC_PUInfo_numberOfInteractions; 
    auto runNumber = *(s.run);
    auto eventNumber = *(s.event);
//    auto nMCPU = *(s.numberOfMC_PUInfo);
    int nMCPU = 0;

    auto& PC_E = s.Conv_refittedPair4Momentum_E;
    auto& PC_M = s.Conv_refittedPair4Momentum_M;
    auto& PC_Px = s.Conv_refittedPair4Momentum_Px;
    auto& PC_Py = s.Conv_refittedPair4Momentum_Py;
    auto& PC_Pz = s.Conv_refittedPair4Momentum_Pz;

//Also needed for asymmetry
    auto& PC_vTrack0_charge = s.Conv_Tk0_charge;
    auto& PC_vTrack1_charge = s.Conv_Tk1_charge;

// New variables - April 2020.
    auto& PC_vTrack0_nBefore = s.Conv_nHitsBeforeVtx_Tk0;
    auto& PC_vTrack1_nBefore = s.Conv_nHitsBeforeVtx_Tk1;
    auto& PC_nSharedHits = s.Conv_nSharedHits;

