#ifndef HISTS
#define HISTS
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "ROOT/TThreadedObject.hxx"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "myselector.C"
#include "myweights.h"
using MyTH1D = ROOT::TThreadedObject<TH1D>;
using MyTH2D = ROOT::TThreadedObject<TH2D>;

class histset{
	
    public:
 	   double PI =4.0*atan(1.0);

       histset();	
       void init(); 
       void setweightoption(); 

       void AnalyzeEntry(myselector& s); 
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
                       id_asym1Hist, id_asym2Hist, id_asym4Hist, id_asym8Hist, id_asym16Hist,
                       id_q0Hist, id_q1Hist, id_qtotHist,
                       id_zcutHist, id_zcutHist2, id_rendcapHist,
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
                     numTH2Hist};

	   //make a big vector and load enumerated histograms onto the vector
       std::vector<MyTH1D*>  TH1Manager{};
       std::vector<MyTH2D*>  TH2Manager{};

	   //locate the histogram and perform pointer copying 
       void FillTH1(int index, double x, double w);
       void FillTH2(int index, double x, double y, double w);
	
       void WriteHist(); 
};

histset::histset(){

    std::vector<MyTH1D*>  Manager1(numTH1Hist);
    TH1Manager=Manager1;

    std::vector<MyTH2D*>  Manager2(numTH2Hist);
    TH2Manager=Manager2;

    init();
    setweightoption();

}

void histset::setweightoption(){

	for(int i=0; i<numTH1Hist; i++){
        auto hptr = TH1Manager.at(i)->Get();
        hptr->Sumw2(kTRUE);
    }

	for(int i=0; i<numTH2Hist; i++){
        auto hptr = TH2Manager.at(i)->Get();
        hptr->Sumw2(kTRUE);
    }

}

void histset::init(){
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

    TH1Manager.at(id_asym1Hist) = new MyTH1D("asym1", "Photon pT > 1 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_asym2Hist) = new MyTH1D("asym2", "Photon pT > 2 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_asym4Hist) = new MyTH1D("asym4", "Photon pT > 4 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_asym8Hist) = new MyTH1D("asym8", "Photon pT > 8 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);
    TH1Manager.at(id_asym16Hist) = new MyTH1D("asym16", "Photon pT > 16 GeV; Positron pT Fraction; Entries per 0.01 bin", 100, 0.0, 1.0);

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
}

void histset::FillTH1(int index, double x, double w=1.0){
	//we must make pointer copies for performance reasons when trying to fill a histogram
	auto myhist = TH1Manager.at(index)->Get();
	myhist->Fill(x,w);
}

void histset::FillTH2(int index, double x, double y, double w=1.0){
	auto myhist = TH2Manager.at(index)->Get();
	myhist->Fill(x,y,w);
}

void histset::WriteHist(){

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

void histset::AnalyzeEntry(myselector& s){
   	
	//always make a local copy, if its a value dereference. 
    //if you dont do this scope/dereferencing will get really weird, clunky, and unmanageable
	//have to auto& or myreader will try to register copy of the readerarray pointer
	auto& PC_vTrack_pt = s.PC_vTrack_pt;
	auto& PC_vTrack_phi = s.PC_vTrack_phi;
	auto& PC_vTrack_eta = s.PC_vTrack_eta;

	auto& PC_vtx_sigmaxx = s.PC_vtx_sigmaxx;
	auto& PC_vtx_sigmaxy = s.PC_vtx_sigmaxy;
	auto& PC_vtx_sigmayy = s.PC_vtx_sigmayy;
	auto& PC_vtx_sigmazz = s.PC_vtx_sigmazz;

	auto& PC_x = s.PC_x;
	auto& PC_y = s.PC_y;
	auto& PC_z = s.PC_z;

// New variables
    auto& PC_mpair = s.PC_pairInvariantMass;
    auto& PC_dcottheta = s.PC_pairCotThetaSeparation;
    auto& PC_dmin = s.PC_distOfMinimumApproach;
    auto& PC_dphi = s.PC_dPhiTracksAtVtx;

    auto isRealData = *(s.isRealData);
    auto lumiSection = *(s.lumiSection);
    auto& mcpu = s.MC_PUInfo_numberOfInteractions; 
    auto runNumber = *(s.runNumber);
    auto nMCPU = *(s.numberOfMC_PUInfo);

	auto& PC_vtx_chi2 = s.PC_vtx_chi2;

	auto numberOfPC = *(s.numberOfPC);
	auto numberOfPV = *(s.numberOfPV);

	auto& PC_fitmomentumOut_pt = s.PC_fitmomentumOut_pt;
	auto& PC_fitmomentumOut_phi = s.PC_fitmomentumOut_phi;
	auto& PC_fitmomentumOut_theta = s.PC_fitmomentumOut_theta;
//	auto& PC_fitmomentumOut_mass = s.PC_fitmomentumOut_mass;

//Also needed for asymmetry
    auto& PC_vTrack_charge = s.PC_vTrack_charge;

	double px,py,pz;
	double x,y,z;
	double r,theta,phi;
    double rho,phip;
    double rps;
    double rnominal;
	
	const double RERRCUT = 0.25;
	const double COSTCUT = 0.85;
	const double ZCUT = 25.0;
	const double FITPROBCUT = 0.010;
    const double MASSCUT = 0.15;

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
    if(!isRealData){
       int sum1 = 0;
       for(int i=0; i<16; i++){
           sum1 += mcpu[i];
       }
       FillTH1(id_PUHist, double(sum1)/16.0);
    }

    std::vector<bool> vcuts;

// Probably a good idea to keep track of whether each conversion 
// candidate passes particular cuts, to ease later code 
// with gamma-gamma invariant mass
	for(int i=0; i<PC_x.GetSize(); i++){
        vcuts.push_back(false);
		x = PC_x[i];
		y = PC_y[i];
		z = PC_z[i];
        double px = PC_fitmomentumOut_pt[i]*cos(PC_fitmomentumOut_phi[i]);
        double py = PC_fitmomentumOut_pt[i]*sin(PC_fitmomentumOut_phi[i]);
        double pz = PC_fitmomentumOut_pt[i]/tan(PC_fitmomentumOut_theta[i]);
        double E = sqrt(px*px + py*py + pz*pz);
        double pt = PC_fitmomentumOut_pt[i]; 

		r = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) );
		phi = atan2(y-y0, x-x0);
		theta = PC_fitmomentumOut_theta[i];

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
		//make r plots with quality cuts and without

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
                && fitprob > FITPROBCUT ){
            vcuts[i] = true;
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
// calculate asymmetry
            double pt0 = PC_vTrack_pt[i][0];
            double pt1 = PC_vTrack_pt[i][1];
            int q0 = PC_vTrack_charge[i][0];
            int q1 = PC_vTrack_charge[i][1];
            int qtot = q0+q1;
// Check charges
            FillTH1(id_q0Hist, q0, wtPU);
            FillTH1(id_q1Hist, q1, wtPU);
            FillTH1(id_qtotHist, qtot, wtPU);
            double ptasym = pt0/(pt0+pt1);
            if (q0<0) ptasym = 1.0-ptasym;
            if(pt>1.0)FillTH1(id_asym1Hist, ptasym, wtPU);
            if(pt>2.0)FillTH1(id_asym2Hist, ptasym, wtPU);
            if(pt>4.0)FillTH1(id_asym4Hist, ptasym, wtPU);
            if(pt>8.0)FillTH1(id_asym8Hist, ptasym, wtPU);
            if(pt>16.0)FillTH1(id_asym16Hist, ptasym, wtPU);
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
            double pxi = PC_fitmomentumOut_pt[i]*cos(PC_fitmomentumOut_phi[i]);
            double pyi = PC_fitmomentumOut_pt[i]*sin(PC_fitmomentumOut_phi[i]);
            double pzi = PC_fitmomentumOut_pt[i]/tan(PC_fitmomentumOut_theta[i]);
            double Ei = sqrt(pxi*pxi + pyi*pyi + pzi*pzi);
            double xi = PC_x[i];
            double yi = PC_y[i];
            double Ri = sqrt(xi*xi + yi*yi);
            for(unsigned int j=i+1; j<numberOfPC; j++){
                double pxj = PC_fitmomentumOut_pt[j]*cos(PC_fitmomentumOut_phi[j]);
                double pyj = PC_fitmomentumOut_pt[j]*sin(PC_fitmomentumOut_phi[j]);
                double pzj = PC_fitmomentumOut_pt[j]/tan(PC_fitmomentumOut_theta[j]);
                double Ej = sqrt(pxj*pxj + pyj*pyj + pzj*pzj);
                double pxgg = pxi + pxj;
                double pygg = pyi + pyj;
                double pzgg = pzi + pzj;
                double Egg  = Ei + Ej;
                double m12 = sqrt(Egg*Egg - pxgg*pxgg - pygg*pygg - pzgg*pzgg);
                double xj = PC_x[j];
                double yj = PC_y[j];
                double Rj = sqrt(xj*xj + yj*yj);
                if (abs(Ri-3.0)<1.0 && abs(Rj-3.0)<1.0&&min(Ei,Ej)>2.0){
                    FillTH1(id_mggHist, m12, wtPU);
                    if(vcuts[i]&&vcuts[j])FillTH2(id_mgg2Hist, m12, numberOfPV, wtPU);
                    if(vcuts[i]&&vcuts[j])FillTH1(id_mggCutHist, m12, wtPU);
                }
                if (vcuts[i]&&vcuts[j]&&min(Ei,Ej)>2.0){
                    FillTH1(id_mggallHist, m12, wtPU);
                    FillTH2(id_mggRCutHist, m12, Ri, wtPU);
                    FillTH2(id_mggRCutHist, m12, Rj, wtPU);
                }
            }
        }
    }//end numpc
		
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
//    cout << " vcuts length " << vcuts.size() << endl;
    vcuts.clear();
    
}
#endif
