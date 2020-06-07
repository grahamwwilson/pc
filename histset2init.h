// Now in "histset2init.h"

void histset2::init(){
//init TH1D
    TH1Manager.at(id_wtHist) = new MyTH1D("wtHist", "Weight Distribution; Weight; Events per bin", 1000, 0.0, 10.0);
    TH1Manager.at(id_wwtHist) = new MyTH1D("wwtHist", "Weighted weight Distribution; Weight; Weighted events per bin", 1000, 0.0, 10.0);
    TH1Manager.at(id_ptHist) = new MyTH1D("ptHist", "p_{T} Distribution;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 5.0);
    TH1Manager.at(id_pzHist) = new MyTH1D("pzHist", "p_{Z} Distribution;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);
    TH1Manager.at(id_numpcHist) = new MyTH1D("numpcHist", "Number of PC;;Entries per bin", 100,-0.5, 99.5);
    TH1Manager.at(id_numpvHist) = new MyTH1D("numpvHist", "Number of PV;;Entries per bin", 100,-0.5, 99.5);
    TH1Manager.at(id_numpvWHist) = new MyTH1D("numpvWHist", "Number of PV;;Entries per bin", 200,-0.5, 199.5);
    TH1Manager.at(id_numpvUWHist) = new MyTH1D("numpvUWHist", "Number of PV;;Entries per bin", 200,-0.5, 199.5);
    TH1Manager.at(id_numnopcHist) = new MyTH1D("numnopcHist", "Number of no PC;;Entries per bin", 100,-0.5, 99.5);
    TH1Manager.at(id_numpvnopcHist) = new MyTH1D("numpvnopcHist", "Number of PV (no PC);;Entries per bin", 100,-0.5, 99.5);
    TH1Manager.at(id_rerrHist) = new MyTH1D("rerrHist", "Conversion Radial Error; #Delta R (cm); Entries per 0.05 bin", 40, 0.0, 2.0);
    TH1Manager.at(id_phierrHist) = new MyTH1D("phierrHist", "Conversion Azimuthal Error;#Delta #phi; Entries per 0.002 bin",50, 0.0, 0.1);
    TH1Manager.at(id_zerrHist) = new MyTH1D("zerrHist","Conversion Z Error;#Delta z (cm); Entries per 0.1 bin", 50, 0.0, 5.0);
    TH1Manager.at(id_pfitHist) = new MyTH1D("pfitHist","Photon Conversions;Fit probability; ", 100, 0.0, 1.0);
    TH1Manager.at(id_zHist) = new MyTH1D("zHist","Photon Conversions; z (cm); ", 200, -25.0, 25.0);
    TH1Manager.at(id_costhetaHist) = new MyTH1D("costhetaHist","Photon Conversions; cos(theta); ", 200, -1.0, 1.0);
    TH1Manager.at(id_pfitHist) = new MyTH1D("pfitHist","Photon Conversions;Fit probability; ", 100, 0.0, 1.0);
    TH1Manager.at(id_r1dHist) = new MyTH1D("r1dHist","Conversion Radius No Cuts;R (cm);Entries per 0.1 bin",100, 0.0, 10.0);
    TH1Manager.at(id_r1dHist2) = new MyTH1D("r1dHist2","Conversion Radius;R (cm);Entries per 0.1 bin",250, 0.0, 25.0);
    TH1Manager.at(id_r1dHist3) = new MyTH1D("r1dHist3","Conversion Radius;R (cm);Entries per 0.1 bin",250, 0.0, 25.0);
    TH1Manager.at(id_r1dHist4) = new MyTH1D("r1dHist4","Conversion Radius;R (cm);Entries per 0.1 bin",250, 0.0, 25.0);

    TH1Manager.at(id_r1dwideHist) = new MyTH1D("r1dwideHist","Conversion Radius No Cuts;R (cm);Entries per 0.1 bin",250, 0.0, 25.0);
    TH1Manager.at(id_zcutHist) = new MyTH1D("zcutHist","Photon Conversions; z (cm); ", 250, -25.0, 25.0);
    TH1Manager.at(id_zcutHist2) = new MyTH1D("zcutHist2","Photon Conversions; z (cm); ", 500, -50.0, 50.0);
    TH1Manager.at(id_r1dcutHist) = new MyTH1D("r1dcutHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 100, 0.0, 10.0);
    TH1Manager.at(id_rendcapHist) = new MyTH1D("rendcapHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
    TH1Manager.at(id_r1dwidecutHist) = new MyTH1D("r1dwidecutHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
    TH1Manager.at(id_r1dwidecutDHist) = new MyTH1D("r1dwidecutDHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
    TH1Manager.at(id_r1dwidecutDDHist) = new MyTH1D("r1dwidecutDDHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
    TH1Manager.at(id_r1dwidecutWHist) = new MyTH1D("r1dwidecutWHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
    TH1Manager.at(id_r1dwidecutPSHist) = new MyTH1D("r1dwidecutPSHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
    TH1Manager.at(id_r1dwidecutNomHist) = new MyTH1D("r1dwidecutNomHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
    TH1Manager.at(id_r1dlowPUHist) = new MyTH1D("r1dlowPUHist","Conversion Radius: No Quality Cuts, PV #leq 16;R (cm);Entries per 0.1 bin",100,0.,10.);
    TH1Manager.at(id_r1dmedPUHist) = new MyTH1D("r1dmedPUHist","Conversion Radius: No Quality Cuts, PV #gt 16 and #lt 36;R (cm);Entries per 0.1 bin",100,0.,10.);
    TH1Manager.at(id_r1dhiPUHist) = new MyTH1D("r1dhiPUHist","Conversion Radius: No Quality Cuts, PV #geq 36;R (cm);Entries per 0.1 bin",100,0.,10.);
    TH1Manager.at(id_r1dwidelowPUHist) = new MyTH1D("r1dwidelowPUHist","Conversion Radius: No Quality Cuts, PV #leq 16;R (cm);Entries per 0.1 bin",250,0.,25.);
    TH1Manager.at(id_r1dwidemedPUHist) = new MyTH1D("r1dwidemedPUHist","Conversion Radius: No Quality Cuts, PV #gt 16 and #lt 36;R (cm);Entries per 0.1 bin",250,0.,25.);
    TH1Manager.at(id_r1dwidehiPUHist) = new MyTH1D("r1dwidehiPUHist","Conversion Radius: No Quality Cuts, PV #geq 36;R (cm);Entries per 0.1 bin",250,0.,25.);
    TH1Manager.at(id_r1dlowPUcutHist) = new MyTH1D("r1dlowPUcutHist","Conversion Radius: Quality Cuts, PV #leq 16;R (cm);Entries per 0.1 bin",100,0.,10.);
    TH1Manager.at(id_r1dmedPUcutHist) = new MyTH1D("r1dmedPUcutHist","Conversion Radius: Quality Cuts, PV #gt 16 and #lt 36;R (cm);Entries per 0.1 bin",100,0.,10.);
    TH1Manager.at(id_r1dhiPUcutHist) = new MyTH1D("r1dhiPUcutHist","Conversion Radius: Quality Cuts, PV #geq 36;R (cm);Entries per 0.1 bin",100,0.,10.);
    TH1Manager.at(id_r1dwidelowPUcutHist) = new MyTH1D("r1dwidelowPUcutHist","Conversion Radius: Quality Cuts, PV #leq 16;R (cm);Entries per 0.1 bin",250,0.,25.);
    TH1Manager.at(id_r1dwidemedPUcutHist) = new MyTH1D("r1dwidemedPUcutHist","Conversion Radius: Quality Cuts, PV #gt 16 and #lt 36;R (cm);Entries per 0.1 bin",250,0.,25.);
    TH1Manager.at(id_r1dwidehiPUcutHist) = new MyTH1D("r1dwidehiPUcutHist","Conversion Radius: Quality Cuts, PV #geq 36;R (cm);Entries per 0.1 bin",250,0.,25.);
    TH1Manager.at(id_rhobpHist) = new MyTH1D("rhobpHist","Conversion Radius w.r.t Beam Pipe Center and Quality Cuts; R (cm); Entries per 0.05 bin",100,0.,5.);
    TH1Manager.at(id_rbpHist) = new MyTH1D("rbpHist","Conversion Radius w.r.t BPIX Center and Quality Cuts; R (cm); Entries per 0.05 bin",100,0.,5.);
    TH1Manager.at(id_rnomHist) = new MyTH1D("rnomHist","Conversion Radius w.r.t Nominal Center and Quality Cuts; R (cm); Entries per 0.05 bin",100,0.,5.);
    TH1Manager.at(id_mggHist) = new MyTH1D("mggHist","Di-#gamma Mass;Mass GeV; Entries per 2.5 MeV bin", 400, 0.0, 1.0 );
    TH1Manager.at(id_mggallHist) = new MyTH1D("mggallHist","Di-#gamma Mass;Mass GeV; Entries per 2.5 MeV bin", 400, 0.0, 1.0 );
    TH1Manager.at(id_mggCutHist) = new MyTH1D("mggCutHist","Di-#gamma Mass;Mass GeV; Entries per 2.5 MeV bin", 400, 0.0,1.0 );
    TH1Manager.at(id_pTHist) = new MyTH1D("pTHist","Photon pT;pT (GeV); Entries per 0.1 GeV bin", 1000, 0.0, 100.0);
    TH1Manager.at(id_EHist) = new MyTH1D("EHist","Photon Energy;Energy (GeV); Entries per 0.1 GeV bin", 1000, 0.0, 100.0 );
    TH1Manager.at(id_phiHist) = new MyTH1D("phiHist","Photon Phi;Phi (rad); Entries per bin", 40, -PI, PI );
    TH1Manager.at(id_pTHist2) = new MyTH1D("pTHist2","Photon pT;pT (GeV); Entries per 0.1 GeV bin", 1000, 0.0, 100.0);
    TH1Manager.at(id_EHist2) = new MyTH1D("EHist2","Photon Energy;Energy (GeV); Entries per 0.1 GeV bin", 1000, 0.0, 100.0 );
    TH1Manager.at(id_phiHist2) = new MyTH1D("phiHist2","Photon Phi;Phi (rad); Entries per bin", 40, -PI, PI );
    TH1Manager.at(id_runHist) = new MyTH1D("runHist",";Run Number; Events per bin", 5000, 320500.5, 325500.5);
    TH1Manager.at(id_isdataHist) = new MyTH1D("isdataHist",";isData; Events per bin", 2, -0.5, 1.5);
    TH1Manager.at(id_dminHist) = new MyTH1D("dminHist",";dmin; Conversions per bin", 200, -1.0, 1.0);
    TH1Manager.at(id_dphiHist) = new MyTH1D("dphiHist",";dphi; Conversions per bin", 200, -0.8, 0.8);
    TH1Manager.at(id_dcotthetaHist) = new MyTH1D("dcotthetaHist",";dcottheta; Conversions per bin", 200, -0.5, 0.5);
    TH1Manager.at(id_mpairHist) = new MyTH1D("mpairHist",";mpair (GeV); Conversions per bin", 52, -0.005, 0.255);
// These two are only relevant if MC
    TH1Manager.at(id_nPUHist) = new MyTH1D("nPUHist",";nPU (MC truth); Events per bin", 100, -0.5, 99.5);
    TH1Manager.at(id_PUHist) = new MyTH1D("PUDistributionHist",";PU (MC); Events per bin", 400, -0.5, 99.5);
    TH1Manager.at(id_q0Hist) = new MyTH1D("q0Hist","; q0; Events per bin", 5, -2.5, 2.5);
    TH1Manager.at(id_q1Hist) = new MyTH1D("q1Hist","; q1; Events per bin", 5, -2.5, 2.5);
    TH1Manager.at(id_qtotHist) = new MyTH1D("qtotHist","; qtot; Events per bin", 5, -2.5, 2.5);

    TH1Manager.at(id_xplus0Hist) = new MyTH1D("xplus0", "Photon pT <= 1 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus1Hist) = new MyTH1D("xplus1", "Photon pT > 1 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus2Hist) = new MyTH1D("xplus2", "Photon pT > 2 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus4Hist) = new MyTH1D("xplus4", "Photon pT > 4 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus8Hist) = new MyTH1D("xplus8", "Photon pT > 8 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus16Hist) = new MyTH1D("xplus16", "Photon pT > 16 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);

    TH1Manager.at(id_alpha0Hist) = new MyTH1D("alpha0", "Photon pT <= 1 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha1Hist) = new MyTH1D("alpha1", "Photon pT > 1 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha2Hist) = new MyTH1D("alpha2", "Photon pT > 2 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha4Hist) = new MyTH1D("alpha4", "Photon pT > 4 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha8Hist) = new MyTH1D("alpha8", "Photon pT > 8 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha16Hist) = new MyTH1D("alpha16", "Photon pT > 16 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);

    TH1Manager.at(id_r0Hist) = new MyTH1D("r0", "Photon pT <= 1 GeV; Radius (cm); Entries per 0.1 cm bin", 250, 0.0, 25.0);
    TH1Manager.at(id_r1Hist) = new MyTH1D("r1", "Photon pT > 1 GeV; Radius (cm); Entries per 0.1 cm bin", 250,  0.0, 25.0);
    TH1Manager.at(id_r2Hist) = new MyTH1D("r2", "Photon pT > 2 GeV; Radius (cm); Entries per 0.1 cm bin", 250,  0.0, 25.0);
    TH1Manager.at(id_r4Hist) = new MyTH1D("r4", "Photon pT > 4 GeV; Radius (cm); Entries per 0.1 cm bin", 250,  0.0, 25.0);
    TH1Manager.at(id_r8Hist) = new MyTH1D("r8", "Photon pT > 8 GeV; Radius (cm); Entries per 0.1 cm bin", 250,  0.0, 25.0);
    TH1Manager.at(id_r16Hist) = new MyTH1D("r16", "Photon pT > 16 GeV; Radius (cm); Entries per 0.1 cm bin", 250, 0.0, 25.0);

    TH1Manager.at(id_r1Hist2) = new MyTH1D("r1p", "Photon pT [1,2] GeV; Radius (cm); Entries per 0.1 cm bin", 250,  0.0, 25.0);
    TH1Manager.at(id_r2Hist2) = new MyTH1D("r2p", "Photon pT [2,4] GeV; Radius (cm); Entries per 0.1 cm bin", 250,  0.0, 25.0);
    TH1Manager.at(id_r4Hist2) = new MyTH1D("r4p", "Photon pT [4,8] GeV; Radius (cm); Entries per 0.1 cm bin", 250,  0.0, 25.0);
    TH1Manager.at(id_r8Hist2) = new MyTH1D("r8p", "Photon pT [8,16] GeV; Radius (cm); Entries per 0.1 cm bin", 250,  0.0, 25.0);
    TH1Manager.at(id_rAsymmetricHist) = new MyTH1D("rAsymmetric", "Photon pT > 16 GeV (|#alpha| > 0.96); Radius (cm); Entries per 0.25 cm bin", 100, 0.0, 25.0);

    TH1Manager.at(id_conversionCandidateMassHist) = new MyTH1D("conversionCandidateMassHist","e+e- Pair Mass; Mass (GeV); Conversions per bin",200,0.0,1.0);
    TH1Manager.at(id_conversionCandidateMassHist2) = new MyTH1D("conversionCandidateMassHist2","e+e- Pair Mass; Mass (GeV); Conversions per bin",100,0.0,0.1);
    TH1Manager.at(id_lambdasCandidateMassHist) = new MyTH1D("lambdasCandidateMassHist","Lambda/Lambdabar Candidate Mass; Mass (GeV); Candidates per bin",225,1.05,1.50);
    TH1Manager.at(id_lambdasBkgdMassHist) = new MyTH1D("lambdasBkgdMassHist","Lambda/Lambdabar Candidate Mass; Mass (GeV); Candidates per bin",90,1.05,1.50);
    TH1Manager.at(id_lambdasBkgdMassHistR1) = new MyTH1D("lambdasBkgdMassHistR1","Lambda/Lambdabar Candidate Mass; Mass (GeV); Candidates per bin",90,1.05,1.50);
    TH1Manager.at(id_lambdasBkgdMassHistR2) = new MyTH1D("lambdasBkgdMassHistR2","Lambda/Lambdabar Candidate Mass; Mass (GeV); Candidates per bin",90,1.05,1.50);
    TH1Manager.at(id_lambdasBkgdMassHistR3) = new MyTH1D("lambdasBkgdMassHistR3","Lambda/Lambdabar Candidate Mass; Mass (GeV); Candidates per bin",90,1.05,1.50);
    TH1Manager.at(id_lambdasSignalMassHist) = new MyTH1D("lambdasSignalMassHist","Lambda/Lambdabar Candidate Mass; Mass (GeV); Candidates per bin",90,1.05,1.50);
    TH1Manager.at(id_lambdabarCandidateMassHist) = new MyTH1D("lambdabarCandidateMassHist","Lambdabar Candidate Mass; Mass (GeV); Candidates per bin",225,1.05,1.50);
    TH1Manager.at(id_lambdaCandidateMassHist) = new MyTH1D("lambdaCandidateMassHist","Lambda Candidate Mass; Mass (GeV); Candidates per bin",225,1.05,1.50);

    TH1Manager.at(id_KShortMassHist) = new MyTH1D("KShortMassHist","KShort Candidate Mass; Mass (GeV); Candidates per bin",250,0.25,0.75);
    TH1Manager.at(id_KShortBkgdMassHist) = new MyTH1D("KShortBkgdMassHist","KShort Candidate Mass; Mass (GeV); Candidates per bin",250,0.25,0.75);
    TH1Manager.at(id_KShortBkgdMassHistR1) = new MyTH1D("KShortBkgdMassHistR1","KShort Candidate Mass; Mass (GeV); Candidates per bin",250,0.25,0.75);
    TH1Manager.at(id_KShortBkgdMassHistR2) = new MyTH1D("KShortBkgdMassHistR2","KShort Candidate Mass; Mass (GeV); Candidates per bin",250,0.25,0.75);
    TH1Manager.at(id_KShortBkgdMassHistR3) = new MyTH1D("KShortBkgdMassHistR3","KShort Candidate Mass; Mass (GeV); Candidates per bin",250,0.25,0.75);
    
    TH1Manager.at(id_AP_pTminHist) = new MyTH1D("AP_pTminHist","Armenteros-Podolanski pT; Minimum pT (GeV); Conversions per bin",200,0.0,0.1);
    TH1Manager.at(id_AP_pTmaxHist) = new MyTH1D("AP_pTmaxHist","Armenteros-Podolanski pT; Maximum pT (GeV); Conversions per bin",200,0.0,0.1);
    TH1Manager.at(id_AP_pTaveHist) = new MyTH1D("AP_pTaveHist","Armenteros-Podolanski pT; Average pT (GeV); Conversions per bin",200,0.0,0.1);
    TH1Manager.at(id_AP_alphaHist) = new MyTH1D("AP_alphaHist","Armenteros-Podolanski alpha; (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-}); Conversions per bin",200,-1.0,1.0);

    TH1Manager.at(id_alphaBkgdHist) = new MyTH1D("alphaBkgdHist","Armenteros-Podolanski alpha; (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-}); Conversions per bin",200,-1.0,1.0);
    TH1Manager.at(id_alphaBkgdHistR1) = new MyTH1D("alphaBkgdHistR1","Armenteros-Podolanski alpha; (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-}); Conversions per bin",200,-1.0,1.0);
    TH1Manager.at(id_alphaBkgdHistR2) = new MyTH1D("alphaBkgdHistR2","Armenteros-Podolanski alpha; (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-}); Conversions per bin",200,-1.0,1.0);
    TH1Manager.at(id_alphaBkgdHistR3) = new MyTH1D("alphaBkgdHistR3","Armenteros-Podolanski alpha; (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-}); Conversions per bin",200,-1.0,1.0);
    TH1Manager.at(id_alphaSignalHist) = new MyTH1D("alphaSignalHist","Armenteros-Podolanski alpha; (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-}); Conversions per bin",200,-1.0,1.0);

    TH1Manager.at(id_nconvHist) = new MyTH1D("nconvHist", "Photon Conversions; Number of Conversions; Entries per bin", 40, -0.5, 39.5 );
    TH1Manager.at(id_nassignedHist) = new MyTH1D("nassignedHist", "Photon Conversions; Number of Conversions; Entries per bin", 40, -0.5, 39.5 );
    TH1Manager.at(id_nnonassignedHist) = new MyTH1D("nnonassignedHist", "Photon Conversions; Number of Conversions; Entries per bin", 40, -0.5, 39.5 );

// init TH2D
    TH2Manager.at(id_pxpyHist) = new MyTH2D("pxpyHist", "p_{X} vs p_{Y} Distribution;p_{X};p_{Y}", 200, -10., 10., 200, -10., 10.);
    TH2Manager.at(id_xyHist) = new MyTH2D("xyHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
    TH2Manager.at(id_xywideHist) = new MyTH2D("xywideHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_rphiHist) = new MyTH2D("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",250,0.0,25.0,40,-PI,PI);
    TH2Manager.at(id_rzHist) = new MyTH2D("rzHist", "Conversion Vertices in R-z per mm^{2} bin; PC_z (cm); R (cm)",200,-10.,10.,100,0.,10.);
    TH2Manager.at(id_rzHist2) = new MyTH2D("rzHist2", "Conversion Vertices in R-z per mm^{2} bin; |z| (cm); R (cm)",250,25.0,50.0,250,0.0,25.0);
    TH2Manager.at(id_rzHist3) = new MyTH2D("rzHist3", "Conversion Vertices in R-z per mm^{2} bin; |z| (cm); R (cm)",250,25.0,50.0,250,0.0,25.0);
    TH2Manager.at(id_rzHist4) = new MyTH2D("rzHist4", "Conversion Vertices in R-z per mm^{2} bin; |z| (cm); R (cm)",250,25.0,50.0,250,0.0,25.0);
    TH2Manager.at(id_xycutHist) = new MyTH2D("xycutHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
    TH2Manager.at(id_xywidecutHist) = new MyTH2D("xywidecutHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_xywidecutHist2) = new MyTH2D("xywidecutHist2", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
    TH2Manager.at(id_npv_rcutHist) = new MyTH2D("npv_rcutHist", "Conversion Vertices per mm; R (cm); nPV",250,0.0,25.0,100,-0.5,99.5);
    TH2Manager.at(id_npc_npvHist) = new MyTH2D("npc_npvHist", " ; nPV; nPC",100,-0.5,99.5,100,-0.5,99.5);
    TH2Manager.at(id_mgg2Hist) = new MyTH2D("mgg2Hist","Di-#gamma Mass;Mass (GeV); nPV", 400, 0.0, 1.0, 100, -0.5, 99.5 );
    TH2Manager.at(id_mggRCutHist) = new MyTH2D("mggRCutHist","Di-#gamma Mass;Mass (GeV); Radius (cm)", 400, 0.0, 1.0, 25, 0.0, 25.0 );
    TH2Manager.at(id_rhophiHist) = new MyTH2D("rhophiHist","Conversion Radius w.r.t Beam Pipe Center and Quality Cuts; R (cm); Phi (rad)",100,0.0,5.0,40,-PI,PI);
    TH2Manager.at(id_AP_pT_alphaHist) = new MyTH2D("AP_pT_alphaHist","Armenteros-Podolanski Plot; #alpha = (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-}); p_{T} (GeV)",200,-1.0,1.0,200,0.0,0.22);
    TH2Manager.at(id_pTalphaHist) = new MyTH2D("pTalphaHist","Photon Candidate pT vs #alpha; #alpha = (p_{L}^{+} - p_{L}^{-})/(p_{L}^{+} + p_{L}^{-}); p_{T} (photon) (GeV)",200,-1.0,1.0,100,0.0,25.0);
}

