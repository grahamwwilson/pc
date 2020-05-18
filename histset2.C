#ifndef HISTS
#define HISTS
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "ROOT/TThreadedObject.hxx"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "convsel.C"
#include "myweights.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
//#include <vector>
using MyTH1D = ROOT::TThreadedObject<TH1D>;
using MyTH2D = ROOT::TThreadedObject<TH2D>;

class histset2{
	
    public:
 	   double PI =4.0*atan(1.0);

       histset2();	
       void init(); 
       void setweightoption(); 

       void AnalyzeEntry(convsel& s); 
       //bookeeping enumeration: 
       //(if we do this we don't need to worry about hist pointer copies and merging)
       enum th1d_ids{id_ptHist, id_pzHist, id_numpcHist, id_numpvHist,
                       id_numpvWHist,
                       id_rerrHist, id_phierrHist, id_zerrHist,
                       id_r1dHist, id_r1dcutHist, 
                       id_r1dlowPUHist, id_r1dmedPUHist, id_r1dhiPUHist, 
                       id_r1dlowPUcutHist, id_r1dmedPUcutHist, id_r1dhiPUcutHist, 
                       id_r1dwideHist, id_r1dwidecutHist, id_r1dwidecutWHist,
                       id_r1dwidecutPSHist,
                       id_r1dwidelowPUHist, id_r1dwidemedPUHist, id_r1dwidehiPUHist, 
                       id_r1dwidelowPUcutHist, id_r1dwidemedPUcutHist, id_r1dwidehiPUcutHist, 
                       id_rhobpHist, id_rbpHist, id_mggHist, id_mggCutHist,
                       id_numnopcHist, id_numpvnopcHist, id_phiHist,
                       id_mggallHist, id_pfitHist, id_zHist, id_costhetaHist,
                       id_pTHist, id_EHist,
                       id_pTHist2, id_EHist2, id_phiHist2, id_runHist,
                       id_isdataHist, id_nPUHist, id_PUHist, id_wtHist,
                       id_wwtHist, id_numpvUWHist, 
                       id_dminHist, id_dphiHist, id_mpairHist, id_dcotthetaHist,
                       id_r1dwidecutDHist, id_r1dwidecutDDHist,
                       id_r1dwidecutNomHist, id_rnomHist,
                       id_xplus1Hist, id_xplus2Hist, id_xplus4Hist, id_xplus8Hist, id_xplus16Hist,
                       id_alpha1Hist, id_alpha2Hist, id_alpha4Hist, id_alpha8Hist, id_alpha16Hist,
                       id_q0Hist, id_q1Hist, id_qtotHist,
                       id_zcutHist, id_zcutHist2, id_rendcapHist,
                       id_conversionCandidateMassHist,
                       id_conversionCandidateMassHist2,
                       id_lambdaCandidateMassHist,id_lambdabarCandidateMassHist,id_lambdasCandidateMassHist,
                       id_lambdasBkgdMassHist,id_lambdasSignalMassHist,
                       id_lambdasBkgdMassHistR1,id_lambdasBkgdMassHistR2,id_lambdasBkgdMassHistR3,
                       id_KShortMassHist,id_KShortBkgdMassHistR1,id_KShortBkgdMassHistR2,
                       id_KShortBkgdMassHistR3,id_KShortBkgdMassHist,
                       id_AP_pTminHist, id_AP_pTmaxHist, id_AP_pTaveHist,
                       id_AP_alphaHist,
                       id_alphaBkgdHist, id_alphaSignalHist,
                       id_alphaBkgdHistR1, id_alphaBkgdHistR2, id_alphaBkgdHistR3,
                       numTH1Hist};
       enum th2d_ids{id_pxpyHist,
                     id_xyHist,
                     id_xywideHist,
                     id_rphiHist, 
                     id_rzHist, id_rzHist2, id_rzHist3, id_rzHist4,
                     id_xycutHist,
                     id_xywidecutHist, id_xywidecutHist2,
                     id_npv_rcutHist,
                     id_mgg2Hist,
                     id_mggRCutHist,
                     id_npc_npvHist, id_rhophiHist,
                     id_AP_pT_alphaHist,
                     numTH2Hist};

	   //make a big vector and load enumerated histograms onto the vector
       std::vector<MyTH1D*>  TH1Manager{};
       std::vector<MyTH2D*>  TH2Manager{};

	   //locate the histogram and perform pointer copying 
       void FillTH1(int index, double x, double w);
       void FillTH2(int index, double x, double y, double w);
	
       void WriteHist(); 
};

histset2::histset2(){

    std::vector<MyTH1D*>  Manager1(numTH1Hist);
    TH1Manager=Manager1;

    std::vector<MyTH2D*>  Manager2(numTH2Hist);
    TH2Manager=Manager2;

    init();
    setweightoption();

}

void histset2::setweightoption(){

	for(int i=0; i<numTH1Hist; i++){
        auto hptr = TH1Manager.at(i)->Get();
        hptr->Sumw2(kTRUE);
    }

	for(int i=0; i<numTH2Hist; i++){
        auto hptr = TH2Manager.at(i)->Get();
        hptr->Sumw2(kTRUE);
    }

}

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

/*    asym = new TH1D("asym", "Positron pT Fraction; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    asymcut0 = new TH1D("asymcut0", "Positron pT Fraction; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    asymcut1 = new TH1D("asymcut1", "Positron pT Fraction; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    asymcut2 = new TH1D("asymcut2", "Positron pT Fraction; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    asymcut3 = new TH1D("asymcut3", "Positron pT Fraction; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    asymcut4 = new TH1D("asymcut4", "Positron pT Fraction; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0); */

    TH1Manager.at(id_xplus1Hist) = new MyTH1D("xplus1", "Photon pT > 1 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus2Hist) = new MyTH1D("xplus2", "Photon pT > 2 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus4Hist) = new MyTH1D("xplus4", "Photon pT > 4 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus8Hist) = new MyTH1D("xplus8", "Photon pT > 8 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_xplus16Hist) = new MyTH1D("xplus16", "Photon pT > 16 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);

    TH1Manager.at(id_alpha1Hist) = new MyTH1D("alpha1", "Photon pT > 1 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha2Hist) = new MyTH1D("alpha2", "Photon pT > 2 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha4Hist) = new MyTH1D("alpha4", "Photon pT > 4 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha8Hist) = new MyTH1D("alpha8", "Photon pT > 8 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);
    TH1Manager.at(id_alpha16Hist) = new MyTH1D("alpha16", "Photon pT > 16 GeV; AP alpha; Entries per 0.02 bin", 100, -1.0, 1.0);

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
    

/*    asym2 = new TH1D("asym2", "Photon pT > 2 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    asym4 = new TH1D("asym4", "Photon pT > 4 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    asym8 = new TH1D("asym8", "Photon pT > 8 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    asym16 = new TH1D("asym16", "Photon pT > 16 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0); */

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
//                    id_AP_pT_alphaHist,
	TH2Manager.at(id_AP_pT_alphaHist) = new MyTH2D("AP_pT_alphaHist","Armenteros-Podolanski Plot; alpha; pT (GeV)",200,-1.0,1.0,200,0.0,0.22);
}

void histset2::FillTH1(int index, double x, double w=1.0){
	//we must make pointer copies for performance reasons when trying to fill a histogram
	auto myhist = TH1Manager.at(index)->Get();
	myhist->Fill(x,w);
}

void histset2::FillTH2(int index, double x, double y, double w=1.0){
	auto myhist = TH2Manager.at(index)->Get();
	myhist->Fill(x,y,w);
}

void histset2::WriteHist(){

	TFile* outfile = new TFile("Outfile.root", "RECREATE");

	for(int i=0; i<numTH1Hist; i++){
	//do a check for entries, merge isn't safe for empty histograms
        auto hptr = TH1Manager.at(i)->Get();
	    if(hptr->GetEntries() > 0){	
           auto histmerged = TH1Manager.at(i)->Merge();
           TH1D* h = (TH1D*) histmerged->Clone();
		   outfile->WriteObject(h, h->GetName() );
        }
        else{
           auto h = TH1Manager.at(i)->Get()->Clone();
           outfile->WriteObject(h, h->GetName() );
        }
	}

	for(int i=0; i<numTH2Hist; i++){
	//do a check for entries, merge isn't safe for empty histograms
        auto hptr = TH2Manager.at(i)->Get();
	    if(hptr->GetEntries() > 0){
           auto histmerged = TH2Manager.at(i)->Merge();
           TH2D* h = (TH2D*) histmerged->Clone();
		   outfile->WriteObject(h, h->GetName() );
        }
        else{
           auto h = TH2Manager.at(i)->Get()->Clone();
           outfile->WriteObject(h, h->GetName() );
        }
	}	
}

void histset2::AnalyzeEntry(convsel& s){

    #include "mylocaltree.h"     //All the variable incantations needed

	double px,py,pz;
    double px0p,py0p,pz0p;
    double px1p,py1p,pz1p;
    double pt0,pt1,tanl0,tanl1,phi0,phi1,qR0,qR1;
	double x,y,z;
    double x0p,y0p,z0p;
    double x1p,y1p,z1p;
    int nBefore0,nBefore1;
	double r,theta,phi;
    double rho,phip;
    double rps;
    double rnominal;
	
	const double RERRCUT = 0.25;
	const double COSTCUT = 0.85;
	const double ZCUT = 25.0;
	const double FITPROBCUT = 0.010;
    const double MASSCUT = 0.15;

    const double MASS_ELECTRON = 0.5109989461e-3;
    const double MASS_PION = 139.57061e-3;
//    const double MASS_KAON = 493.677e-3;
    const double MASS_PROTON = 938.272081e-3;

// Scale MC to data based on number 
// of events with at least 1 conversion
    const double NFACTOR = 3.68068361*1.00000128;

	double fitprob;

// We now have various "centers" to compare to for radial coordinates.

	//beam pipe displacement (in cm) from Anna's DPF2019 talk
    const double x0bpdata =  0.171;
    const double y0bpdata = -0.176;
    const double x0bpmc = 0.0;
    const double y0bpmc = 0.0;

    //BPIX center displacement (in cm) (Anna's 16-Dec-2019 talk page 3)
    const double x0data =  0.086;
    const double y0data = -0.102;
    const double x0mc = 0.0;
    const double y0mc = 0.0;

    //Pixel support displacement (in cm) (Anna's December talk page 4)
    const double x0psdata = -0.080;
    const double y0psdata = -0.318;
    const double x0psmc = 0.0;
    const double y0psmc = 0.0;

    //Also nominal - where all (x,y) coordinates are referred to (0,0).

    double x0,y0,x0bp,y0bp,x0ps,y0ps;

	//error calcs
	double vxx,vxy,vyy,vzz; //variances
	double varsum_r, varsum_phi; //intermediate calculation variables
	double rerr,phierr,zerr;//errors
	double sphi,cphi;
    
    double wtPU;

    FillTH1(id_runHist, runNumber);
    FillTH1(id_isdataHist, isRealData);
    if(isRealData){
       x0bp = x0bpdata;
       y0bp = y0bpdata;
       x0 = x0data;
       y0 = y0data;
       x0ps = x0psdata;
       y0ps = y0psdata;
       wtPU = 1.0;
    }
    else
    {
// MC  
       x0bp = x0bpmc;
       y0bp = y0bpmc;
       x0 = x0mc;
       y0 = y0mc;
       x0ps = x0psmc;
       y0ps = y0psmc;

       if(numberOfPV<100){
          wtPU = wt[numberOfPV]*NFACTOR;
       }
       else{
          wtPU = wt[100]*NFACTOR;
       }
       FillTH1(id_nPUHist, nMCPU);
    }

// Skip the MC events with no PCs for now.
    if(numberOfPC < 1){
       FillTH1(id_numnopcHist, numberOfPC);
       FillTH1(id_numpvnopcHist, numberOfPV);
       return;
    }
    else{     // not sure if this is needed - but may be need to make sure this histogram exists for data.
       FillTH1(id_numnopcHist, -1);
       FillTH1(id_numpvnopcHist, -1);
    }

// First check the weight distribution
    FillTH1(id_wtHist, wtPU);
    FillTH1(id_wwtHist, wtPU, wtPU);

	FillTH1(id_numpcHist, numberOfPC, wtPU);
	FillTH1(id_numpvHist, numberOfPV);
	FillTH1(id_numpvUWHist, numberOfPV);
	FillTH1(id_numpvWHist, numberOfPV, wtPU);
	FillTH2(id_npc_npvHist, numberOfPV, numberOfPC, wtPU);

// Pileup values (here assume first 12 are in-time, next 4 are out-of-time?)
/*
    if(!isRealData){
       int sum1 = 0;
       for(int i=0; i<16; i++){
           sum1 += mcpu[i];
       }
       FillTH1(id_PUHist, double(sum1)/16.0);
    }
*/

    std::vector<bool> vcuts;

// Probably a good idea to keep track of whether each conversion 
// candidate passes particular cuts, to ease later code 
// with gamma-gamma invariant mass
	for(int i=0; i<PC_x.GetSize(); i++){

        vcuts.push_back(false);
		x = PC_x[i];
		y = PC_y[i];
		z = PC_z[i];
        nBefore0 = PC_vTrack0_nBefore[i];
        nBefore1 = PC_vTrack1_nBefore[i];

// Swim each track from inner hit to the conversion vertex - see convsel.C code
        pt0 = sqrt(Tk0_px[i]*Tk0_px[i] + Tk0_py[i]*Tk0_py[i]);
        tanl0 = Tk0_pz[i]/pt0;
        phi0 = atan2(Tk0_py[i],Tk0_px[i]);
        qR0 = double(PC_vTrack0_charge[i])*100.0*pt0/(0.2998*3.80); //in cm
        double A0 =  2.0*qR0*( (Tk0_y[i]-y)*cos(phi0) - (Tk0_x[i]-x)*sin(phi0) - qR0 );
        double B0 = -2.0*qR0*( (Tk0_y[i]-y)*sin(phi0) + (Tk0_x[i]-x)*cos(phi0) );
        double alp0 = atan2(B0/A0,1.0);
        px0p = pt0*cos(phi0 + alp0);
        py0p = pt0*sin(phi0 + alp0);
        pz0p = Tk0_pz[i];
// Also swim the track position
        x0p = Tk0_x[i] + qR0*( (1.0-cos(alp0))*sin(phi0) - sin(alp0)*cos(phi0) );
        y0p = Tk0_y[i] - qR0*( (1.0-cos(alp0))*cos(phi0) + sin(alp0)*sin(phi0) );
        z0p = Tk0_z[i] - qR0*tanl0*alp0;
        ROOT::Math::XYZVector xyzv0 = ROOT::Math::XYZVector(px0p,py0p,pz0p);

// Swim each track from inner hit to the conversion vertex - see convsel.C code
        pt1 = sqrt(Tk1_px[i]*Tk1_px[i] + Tk1_py[i]*Tk1_py[i]);
        tanl1 = Tk1_pz[i]/pt1;
        phi1 = atan2(Tk1_py[i],Tk1_px[i]);
        qR1 = double(PC_vTrack1_charge[i])*100.0*pt1/(0.2998*3.80); //in cm
        double A1 =  2.0*qR1*( (Tk1_y[i]-y)*cos(phi1) - (Tk1_x[i]-x)*sin(phi1) - qR1 );
        double B1 = -2.0*qR1*( (Tk1_y[i]-y)*sin(phi1) + (Tk1_x[i]-x)*cos(phi1) );
        double alp1 = atan2(B1/A1,1.0);
        px1p = pt1*cos(phi1 + alp1);
        py1p = pt1*sin(phi1 + alp1);
        pz1p = Tk1_pz[i];
// Also swim the track position
        x1p = Tk1_x[i] + qR1*( (1.0-cos(alp1))*sin(phi1) - sin(alp1)*cos(phi1) );
        y1p = Tk1_y[i] - qR1*( (1.0-cos(alp1))*cos(phi1) + sin(alp1)*sin(phi1) );
        z1p = Tk1_z[i] - qR1*tanl1*alp1;
        ROOT::Math::XYZVector xyzv1 = ROOT::Math::XYZVector(px1p,py1p,pz1p);

        double px = PC_Px[i];
        double py = PC_Py[i];
        double pz = PC_Pz[i];
        double E  = PC_E[i];
        double pt = sqrt(px*px + py*py);
        ROOT::Math::XYZVector xyzv = ROOT::Math::XYZVector(px,py,pz);

// calculate asymmetry
        int q0 = PC_vTrack0_charge[i];
        int q1 = PC_vTrack1_charge[i];

        ROOT::Math::XYZVector vcross0 = xyzv.Cross(xyzv0);
        ROOT::Math::XYZVector vcross1 = xyzv.Cross(xyzv1);
        double APpT0 = sqrt(vcross0.Mag2())/sqrt(xyzv.Mag2());
        double APpT1 = sqrt(vcross1.Mag2())/sqrt(xyzv.Mag2());
        double pL0 = xyzv0.Dot(xyzv)/sqrt(xyzv.Mag2());
        double pL1 = xyzv1.Dot(xyzv)/sqrt(xyzv.Mag2());
        double APalpha = (pL0 - pL1)/(pL0 + pL1);
        if (q0<0) APalpha = -APalpha;

		r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
		phi = atan2(y-y0, x-x0);
        theta = atan2(pt,pz);
		rho  =  sqrt( (x-x0bp)*(x-x0bp) + (y-y0bp)*(y-y0bp)) ;
        phip =  atan2(y-y0bp, x-x0bp);
		rps = sqrt( (x-x0ps)*(x-x0ps) + (y-y0ps)*(y-y0ps) );
        rnominal = sqrt( x*x + y*y );
		vxx = PC_vtx_sigmaxx[i];
		vxy = PC_vtx_sigmaxy[i];
		vyy = PC_vtx_sigmayy[i];
		vzz = PC_vtx_sigmazz[i];
		cphi = cos(phi);
		sphi = sin(phi);
		// This is the correct one
		varsum_r   = cphi*cphi*vxx + 2.0*sphi*cphi*vxy + sphi*sphi*vyy;
		varsum_phi = sphi*sphi*vxx - 2.0*sphi*cphi*vxy + cphi*cphi*vyy;
		rerr = sqrt(varsum_r);
		phierr = sqrt(varsum_phi)/r;
		zerr = sqrt(vzz);
	 	fitprob = TMath::Prob(PC_vtx_chi2[i], 3);

        ROOT::Math::PxPyPzMVector v0,v0pi,v0p,v1,v1pi,v1p;
        v0 = ROOT::Math::PxPyPzMVector( px0p, py0p, pz0p, MASS_ELECTRON );
        v0pi = ROOT::Math::PxPyPzMVector( px0p, py0p, pz0p, MASS_PION );
        v0p  = ROOT::Math::PxPyPzMVector( px0p, py0p, pz0p, MASS_PROTON );
        v1   = ROOT::Math::PxPyPzMVector( px1p, py1p, pz1p, MASS_ELECTRON );
        v1pi = ROOT::Math::PxPyPzMVector( px1p, py1p, pz1p, MASS_PION );
        v1p  = ROOT::Math::PxPyPzMVector( px1p, py1p, pz1p, MASS_PROTON );
        ROOT::Math::PxPyPzMVector vpair, vpairpipi, vpairpip, vpairppi;
        vpair += v0;
        vpair += v1;
        vpairpipi += v0pi;
        vpairpipi += v1pi;
        vpairpip  += v0pi;
        vpairpip  += v1p;
        vpairppi  += v0p;
        vpairppi  += v1pi;

        bool region1 = (r>3.7) && (r<6.2);
        bool region2 = (r>7.6) && (r<10.3);
        bool region3 = (r>11.7) && (r<15.3);

// Apply fiducial cuts to all candidates
		if(abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
           FillTH1(id_r1dHist, r, wtPU);
           FillTH1(id_r1dwideHist, r, wtPU);
           FillTH1(id_rerrHist, rerr, wtPU);
           FillTH1(id_phierrHist, phierr, wtPU);
           FillTH1(id_zerrHist, zerr, wtPU);
           FillTH1(id_pfitHist, fitprob, wtPU);
           FillTH1(id_zHist, z, wtPU);
           FillTH1(id_costhetaHist, cos(theta), wtPU);      
    	   FillTH2(id_xyHist, x, y, wtPU);
           FillTH2(id_xywideHist, x, y, wtPU);
//           FillTH2(id_rphiHist, r, phi, wtPU);
           FillTH2(id_rzHist, z, r, wtPU);
           FillTH1(id_pTHist2, pt, wtPU);
           FillTH1(id_EHist2, E, wtPU);
           FillTH1(id_phiHist2, phi, wtPU);
/*    auto& PC_mpair = s.PC_pairInvariantMass;
    auto& PC_dcottheta = s.PC_pairCotThetaSeparation;
    auto& PC_dmin = s.PC_distOfMinimumApproach;
    auto& PC_dphi = s.PC_dPhiTracksAtVtx; */
           FillTH1(id_dminHist, PC_dmin[i], wtPU);
           FillTH1(id_dphiHist, PC_dphi[i], wtPU);
           if(PC_mpair[i]<0.25){
              FillTH1(id_mpairHist, PC_mpair[i], wtPU);
           }
           else{
//Fill overflow bin
              FillTH1(id_mpairHist, 0.2501, wtPU);
           }

           FillTH1(id_dcotthetaHist, PC_dcottheta[i], wtPU);
        }

		//make quality cuts
		if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT 
                && fitprob > FITPROBCUT && std::max(nBefore0,nBefore1)==0 ){

            vcuts[i] = true;

            if(region1||region2||region3){
               FillTH1(id_alphaBkgdHist, APalpha, wtPU);
               if(region1)FillTH1(id_alphaBkgdHistR1, APalpha, wtPU);
               if(region2)FillTH1(id_alphaBkgdHistR2, APalpha, wtPU);
               if(region3)FillTH1(id_alphaBkgdHistR3, APalpha, wtPU);
            }
            else{
               FillTH1(id_alphaSignalHist, APalpha, wtPU);
            }

            FillTH1(id_conversionCandidateMassHist, vpair.M(), wtPU);
            FillTH1(id_conversionCandidateMassHist2, vpair.M(), wtPU);

            FillTH1(id_KShortMassHist, vpairpipi.M(), wtPU);
            if(region1||region2||region3){
               FillTH1(id_KShortBkgdMassHist, vpairpipi.M(), wtPU);
               if(region1)FillTH1(id_KShortBkgdMassHistR1, vpairpipi.M(), wtPU);
               if(region2)FillTH1(id_KShortBkgdMassHistR2, vpairpipi.M(), wtPU);
               if(region3)FillTH1(id_KShortBkgdMassHistR3, vpairpipi.M(), wtPU);
            }            

// pip: Track 0 = pi, Track 1 = p
// ppi: Track 0 = p , Track 1 = pi
            if(APalpha >= 0.0){
// Lambda candidate with a proton ( Lambda -> pi- p ) only populates +ve alpha
               if(q1 == 1){
                  FillTH1(id_lambdaCandidateMassHist, vpairpip.M(), wtPU);
                  FillTH1(id_lambdasCandidateMassHist, vpairpip.M(), wtPU);
                  if(region1||region2||region3){
                     FillTH1(id_lambdasBkgdMassHist, vpairpip.M(), wtPU);
                     if(region1)FillTH1(id_lambdasBkgdMassHistR1, vpairpip.M(), wtPU);
                     if(region2)FillTH1(id_lambdasBkgdMassHistR2, vpairpip.M(), wtPU);
                     if(region3)FillTH1(id_lambdasBkgdMassHistR3, vpairpip.M(), wtPU);
                  }
                  else{
                     FillTH1(id_lambdasSignalMassHist, vpairpip.M(), wtPU);
                  }
               }
               else{
                  FillTH1(id_lambdaCandidateMassHist, vpairppi.M(), wtPU);
                  FillTH1(id_lambdasCandidateMassHist, vpairppi.M(), wtPU);
                  if(region1||region2||region3){
                     FillTH1(id_lambdasBkgdMassHist, vpairppi.M(), wtPU);
                     if(region1)FillTH1(id_lambdasBkgdMassHistR1, vpairppi.M(), wtPU);
                     if(region2)FillTH1(id_lambdasBkgdMassHistR2, vpairppi.M(), wtPU);
                     if(region3)FillTH1(id_lambdasBkgdMassHistR3, vpairppi.M(), wtPU);
                  }
                  else{
                     FillTH1(id_lambdasSignalMassHist, vpairppi.M(), wtPU);
                  }
               }
            }
            else{
// Lambdabar candidate with an antiproton ( Lambdabar -> pi+ pbar) only populates -ve alpha
               if(q0 == 1){
                  FillTH1(id_lambdabarCandidateMassHist, vpairpip.M(), wtPU);
                  FillTH1(id_lambdasCandidateMassHist, vpairpip.M(), wtPU);
                  if(region1||region2||region3){
                     FillTH1(id_lambdasBkgdMassHist, vpairpip.M(), wtPU);
                     if(region1)FillTH1(id_lambdasBkgdMassHistR1, vpairpip.M(), wtPU);
                     if(region2)FillTH1(id_lambdasBkgdMassHistR2, vpairpip.M(), wtPU);
                     if(region3)FillTH1(id_lambdasBkgdMassHistR3, vpairpip.M(), wtPU);
                  }
                  else{
                     FillTH1(id_lambdasSignalMassHist, vpairpip.M(), wtPU);
                  }
               }
               else{
                  FillTH1(id_lambdabarCandidateMassHist, vpairppi.M(), wtPU);
                  FillTH1(id_lambdasCandidateMassHist, vpairppi.M(), wtPU);
                  if(region1||region2||region3){
                     FillTH1(id_lambdasBkgdMassHist, vpairppi.M(), wtPU);
                     if(region1)FillTH1(id_lambdasBkgdMassHistR1, vpairppi.M(), wtPU);
                     if(region2)FillTH1(id_lambdasBkgdMassHistR2, vpairppi.M(), wtPU);
                     if(region3)FillTH1(id_lambdasBkgdMassHistR3, vpairppi.M(), wtPU);
                  }
                  else{
                     FillTH1(id_lambdasSignalMassHist, vpairppi.M(), wtPU);
                  }
               }
            }
            FillTH1(id_AP_pTminHist, std::min(APpT0,APpT1), wtPU);
            FillTH1(id_AP_pTmaxHist, std::max(APpT0,APpT1), wtPU);
            FillTH1(id_AP_pTaveHist, 0.5*(APpT0+APpT1), wtPU);
            FillTH1(id_AP_alphaHist, APalpha, wtPU);
            FillTH2(id_AP_pT_alphaHist, APalpha, std::min(APpT0,APpT1), wtPU);

			FillTH1(id_zcutHist, z, wtPU);
			FillTH1(id_r1dcutHist, r, wtPU);
			FillTH1(id_r1dwidecutHist, r);
			FillTH1(id_r1dwidecutWHist, r, wtPU);
			FillTH1(id_r1dwidecutPSHist, rps, wtPU);
			FillTH1(id_r1dwidecutNomHist, rnominal, wtPU);
            if(PC_dmin[i]>-998.0){
   			   FillTH1(id_r1dwidecutDHist, r, wtPU);
            }
            if(PC_dmin[i]>0.0){
   			   FillTH1(id_r1dwidecutDDHist, r, wtPU);
            }

            FillTH2(id_rphiHist, r, phi, wtPU);

			FillTH2(id_xycutHist, x, y, wtPU);
			FillTH2(id_xywidecutHist, x, y, wtPU);
            FillTH1(id_phiHist, phi, wtPU);	
			FillTH1(id_rhobpHist, rho, wtPU);
			FillTH1(id_rbpHist, r, wtPU);
            FillTH1(id_rnomHist, rnominal, wtPU);
			FillTH2(id_rhophiHist, rho, phip, wtPU);
            FillTH2(id_npv_rcutHist, r, numberOfPV, wtPU);
            FillTH1(id_pTHist, pt, wtPU);
            FillTH1(id_EHist,E, wtPU);

            int qtot = q0+q1;
// Check charges
            FillTH1(id_q0Hist, q0, wtPU);
            FillTH1(id_q1Hist, q1, wtPU);
            FillTH1(id_qtotHist, qtot, wtPU);
            double ptasym = (pt0-pt1)/(pt0+pt1);
            if (q0<0) ptasym = -ptasym;
            double xplus = (1.0 + ptasym)/2.0;
            if(pt>1.0)FillTH1(id_xplus1Hist, xplus, wtPU);
            if(pt>2.0)FillTH1(id_xplus2Hist, xplus, wtPU);
            if(pt>4.0)FillTH1(id_xplus4Hist, xplus, wtPU);
            if(pt>8.0)FillTH1(id_xplus8Hist, xplus, wtPU);
            if(pt>16.0)FillTH1(id_xplus16Hist, xplus, wtPU);

            if(pt>1.0)FillTH1(id_alpha1Hist, APalpha, wtPU);
            if(pt>2.0)FillTH1(id_alpha2Hist, APalpha, wtPU);
            if(pt>4.0)FillTH1(id_alpha4Hist, APalpha, wtPU);
            if(pt>8.0)FillTH1(id_alpha8Hist, APalpha, wtPU);
            if(pt>16.0)FillTH1(id_alpha16Hist, APalpha, wtPU);

		}
/// Endcap plots
		if( rerr < RERRCUT && abs(z) > ZCUT && abs(z) < 50.0 && abs(cos(theta)) < COSTCUT 
                && fitprob > FITPROBCUT ){
			FillTH2(id_xywidecutHist2, x, y, wtPU);
			FillTH2(id_rzHist2, abs(z), r, wtPU);
            if(z>0.0)FillTH2(id_rzHist3, abs(z), r, wtPU);
            if(z<0.0)FillTH2(id_rzHist4, abs(z), r, wtPU);
			FillTH1(id_rendcapHist, r, wtPU);
			FillTH1(id_zcutHist2, z, wtPU);
        }
	
		//pileup cuts
		if( numberOfPV <= 16){
			FillTH1(id_r1dlowPUHist, r, wtPU);
			FillTH1(id_r1dwidelowPUHist, r, wtPU) ;
			//low PU quality cuts
			if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
				FillTH1(id_r1dlowPUcutHist, r, wtPU);
				FillTH1(id_r1dwidelowPUcutHist, r, wtPU);
			}
		}
		if( numberOfPV > 16 && numberOfPV < 36){
			FillTH1(id_r1dmedPUHist, r, wtPU);
			FillTH1(id_r1dwidemedPUHist, r, wtPU);
			//low PU quality cuts
			if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
				FillTH1(id_r1dmedPUcutHist, r, wtPU);
				FillTH1(id_r1dwidemedPUcutHist, r, wtPU);
			}
		}
		if( numberOfPV >= 36){
			FillTH1(id_r1dhiPUHist, r, wtPU);
			FillTH1(id_r1dwidehiPUHist, r, wtPU);
			if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
				FillTH1(id_r1dhiPUcutHist, r, wtPU);
				FillTH1(id_r1dwidehiPUcutHist, r, wtPU);
			}
		}				
	}

	//gamma gamma stuff IS THIS LOOP CORRECT
    // Can also speed this up a bit ...
 	if (numberOfPC>=2){
        for(unsigned int i=0; i<numberOfPC-1; i++){
            double pxi = PC_Px[i];
            double pyi = PC_Py[i];
            double pzi = PC_Pz[i];
            double Ei = PC_E[i];
            double xi = PC_x[i];
            double yi = PC_y[i];
            double Ri = sqrt(xi*xi + yi*yi);
            for(unsigned int j=i+1; j<numberOfPC; j++){
                double pxj = PC_Px[j];
                double pyj = PC_Py[j];
                double pzj = PC_Pz[j];
                double Ej = PC_E[j];
                double pxgg = pxi + pxj;
                double pygg = pyi + pyj;
                double pzgg = pzi + pzj;
                double Egg  = Ei + Ej;
                double m12 = sqrt(Egg*Egg - pxgg*pxgg - pygg*pygg - pzgg*pzgg);
                double xj = PC_x[j];
                double yj = PC_y[j];
                double Rj = sqrt(xj*xj + yj*yj);
                if (abs(Ri-3.0)<1.0 && abs(Rj-3.0)<1.0&&std::min(Ei,Ej)>2.0){
                    FillTH1(id_mggHist, m12, wtPU);
                    if(vcuts[i]&&vcuts[j])FillTH2(id_mgg2Hist, m12, numberOfPV, wtPU);
                    if(vcuts[i]&&vcuts[j])FillTH1(id_mggCutHist, m12, wtPU);
                }
                if (vcuts[i]&&vcuts[j]&&std::min(Ei,Ej)>2.0){
                    FillTH1(id_mggallHist, m12, wtPU);
                    FillTH2(id_mggRCutHist, m12, Ri, wtPU);
                    FillTH2(id_mggRCutHist, m12, Rj, wtPU);
                }
            }
        }
    }//end numpc
		
/*
	for(int i=0; i<PC_vTrack_pt.GetSize(); i++){
        for(int j=0; j<PC_vTrack_pt[i].size(); j++){
			px = PC_vTrack_pt[i][j] * cos( PC_vTrack_phi[i][j] );
			py = PC_vTrack_pt[i][j] * sin( PC_vTrack_phi[i][j] );
			pz = PC_vTrack_pt[i][j] * sinh( PC_vTrack_eta[i][j] );
            FillTH1(id_ptHist, PC_vTrack_pt[i][j], wtPU/PC_vTrack_pt[i][j]);
			FillTH1(id_pzHist, pz, wtPU);
			FillTH2(id_pxpyHist, px, py, wtPU);
        }
    }
*/
//    cout << " vcuts length " << vcuts.size() << endl;
    vcuts.clear();
    
}
#endif
