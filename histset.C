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
                       id_r1dwidelowPUHist, id_r1dwidemedPUHist, id_r1dwidehiPUHist, 
                       id_r1dwidelowPUcutHist, id_r1dwidemedPUcutHist, id_r1dwidehiPUcutHist, 
                       id_rhobpHist, id_mggHist, id_mggCutHist,
                       id_numnopcHist, id_numpvnopcHist, id_phiHist,
                       id_mggallHist, id_pfitHist, id_zHist, id_costhetaHist,
                       id_pTHist, id_EHist,
                       id_pTHist2, id_EHist2, id_phiHist2, id_runHist,
                       id_isdataHist, id_nPUHist, id_PUHist,
                       numTH1Hist};
       enum th2d_ids{id_pxpyHist,
                     id_xyHist,
                     id_xywideHist,
                     id_rphiHist, 
                     id_rzHist,
                     id_xycutHist,
                     id_xywidecutHist,
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

}
void histset::init(){
//init TH1D
    TH1Manager.at(id_ptHist) = new MyTH1D("ptHist", "p_{T} Distribution;p_{T};1/p_{T} dN/dp_{T}", 100, 0.0, 5.0);
    TH1Manager.at(id_pzHist) = new MyTH1D("pzHist", "p_{Z} Distribution;p_{Z};dN/dp_{Z}", 100, 0.0, 5.0);
	TH1Manager.at(id_numpcHist) = new MyTH1D("numpcHist", "Number of PC;;Entries per bin", 100,-0.5, 99.5);
	TH1Manager.at(id_numpvHist) = new MyTH1D("numpvHist", "Number of PV;;Entries per bin", 100,-0.5, 99.5);
	TH1Manager.at(id_numpvWHist) = new MyTH1D("numpvWHist", "Number of PV;;Entries per bin", 200,-0.5, 199.5);
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
	TH1Manager.at(id_r1dcutHist) = new MyTH1D("r1dcutHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 100, 0.0, 10.0);
	TH1Manager.at(id_r1dwidecutHist) = new MyTH1D("r1dwidecutHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
	TH1Manager.at(id_r1dwidecutWHist) = new MyTH1D("r1dwidecutWHist","Conversions Radius With Cuts; R (cm); Entries per 0.1 bin", 250, 0.0, 25.0);
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
// These two are only relevant if MC
	TH1Manager.at(id_nPUHist) = new MyTH1D("nPUHist",";nPU (MC truth); Events per bin", 100, -0.5, 99.5);
	TH1Manager.at(id_PUHist) = new MyTH1D("PUDistributionHist",";PU (MC); Events per bin", 400, -0.5, 99.5);

// init TH2D
	TH2Manager.at(id_pxpyHist) = new MyTH2D("pxpyHist", "p_{X} vs p_{Y} Distribution;p_{X};p_{Y}", 200, -10., 10., 200, -10., 10.);
	TH2Manager.at(id_xyHist) = new MyTH2D("xyHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
	TH2Manager.at(id_xywideHist) = new MyTH2D("xywideHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
	TH2Manager.at(id_rphiHist) = new MyTH2D("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",100,0.,10.,100,-3.2,3.2);
	TH2Manager.at(id_rzHist) = new MyTH2D("rzHist", "Conversion Vertices in R-z per mm^{2} bin; PC_z (cm); R (cm)",200,-10.,10.,100,0.,10.);
	TH2Manager.at(id_xycutHist) = new MyTH2D("xycutHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",200,-10.,10.,200,-10.,10.);
	TH2Manager.at(id_xywidecutHist) = new MyTH2D("xywidecutHist", "Conversion Vertices per mm^{2} bin; x (cm); y (cm)",500,-25.,25.,500,-25.,25.);
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
		auto histmerged = TH1Manager.at(i)->Merge();
		TH1D* h = (TH1D*) histmerged->Clone();
		outfile->WriteObject(h, h->GetName() );
	}

	for(int i=0; i<numTH2Hist; i++){
		auto histmerged = TH2Manager.at(i)->Merge();
		TH2D* h = (TH2D*) histmerged->Clone();
		outfile->WriteObject(h, h->GetName() );
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
    auto& PC_cottheta = s.PC_pairCotThetaSeparation;
    auto& PC_dist = s.PC_distOfMinimumApproach;
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

	double px,py,pz;
	double x,y,z;
	double r,theta,phi;
    double rho,phip;
	
	const double RERRCUT = 0.25;
	const double COSTCUT = 0.85;
	const double ZCUT = 25.0;
	const double FITPROBCUT = 0.010;

	double fitprob;

	//beam pipe displacement (in cm) from Anna's DPF2019 talk
    const double x0data =  0.171;
    const double y0data = -0.176;
    const double x0mc =  0.0;
    const double y0mc =  0.0;
    double x0,y0;

	//error calcs
	double vxx,vxy,vyy,vzz; //variances
	double varsum_r, varsum_phi; //intermediate calculation variables
	double rerr,phierr,zerr;//errors
	double sphi,cphi;
    
    double PUwt;

    FillTH1(id_runHist, runNumber);
    FillTH1(id_isdataHist, isRealData);
    if(isRealData){
       x0 = x0data;
       y0 = y0data;
       PUwt = 1.0;
//Fix
//       FillTH1(id_nPUHist, 1.0);
    }
    else
    {
// MC  
       x0 = x0mc;
       y0 = y0mc;
       if(numberOfPV<100){
          PUwt = wt[numberOfPV];
       }
       else{
          PUwt = wt[100];
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

	FillTH1(id_numpcHist, numberOfPC);
	FillTH1(id_numpvHist, numberOfPV);
	FillTH1(id_numpvWHist, numberOfPV, PUwt);
	FillTH2(id_npc_npvHist, numberOfPV, numberOfPC);

// Pileup values (here assume first 12 are in-time, next 4 are out-of-time?)
    if(!isRealData){
       int sum1 = 0;
       for(int i=0; i<16; i++){
           sum1 += mcpu[i];
       }
       FillTH1(id_PUHist, double(sum1)/16.0);
    }
    else{
// Fix for data
//       FillTH1(id_PUHist, 0.0);
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

		r = sqrt( x*x + y*y );
		phi = atan2(y, x);
		theta = PC_fitmomentumOut_theta[i];

		rho  =  sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0)) ;
        phip =  atan2(y-y0, x-x0);

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
           FillTH1(id_r1dHist, r);
           FillTH1(id_r1dwideHist, r);
           FillTH1(id_rerrHist, rerr);
           FillTH1(id_phierrHist, phierr);
           FillTH1(id_zerrHist, zerr);
           FillTH1(id_pfitHist, fitprob);
           FillTH1(id_zHist, z);
           FillTH1(id_costhetaHist, cos(theta));      
    	   FillTH2(id_xyHist, x, y);
           FillTH2(id_xywideHist, x, y);
           FillTH2(id_rphiHist, r, phi);
           FillTH2(id_rzHist, z, r);
           FillTH1(id_pTHist2, pt);
           FillTH1(id_EHist2,E);
           FillTH1(id_phiHist2, phi);
        }

		//make quality cuts
		if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
            vcuts[i] = true;
			FillTH1(id_r1dcutHist, r);
			FillTH1(id_r1dwidecutHist, r);
			FillTH1(id_r1dwidecutWHist, r, PUwt);
			FillTH2(id_xycutHist, x, y);
			FillTH2(id_xywidecutHist, x, y);
            FillTH1(id_phiHist, phi);	
			FillTH1(id_rhobpHist, rho);
			FillTH2(id_rhophiHist, rho, phip);
            FillTH2(id_npv_rcutHist, r, numberOfPV);
            FillTH1(id_pTHist, pt);
            FillTH1(id_EHist,E);
		}			
	
		//pileup cuts
		if( numberOfPV <= 16){
			FillTH1(id_r1dlowPUHist, r);
			FillTH1(id_r1dwidelowPUHist, r);
			//low PU quality cuts
			if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
				FillTH1(id_r1dlowPUcutHist, r);
				FillTH1(id_r1dwidelowPUcutHist, r);
			}
		}
		if( numberOfPV > 16 && numberOfPV < 36){
			FillTH1(id_r1dmedPUHist, r);
			FillTH1(id_r1dwidemedPUHist, r);
			//low PU quality cuts
			if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
				FillTH1(id_r1dmedPUcutHist, r);
				FillTH1(id_r1dwidemedPUcutHist, r);
			}
		}
		if( numberOfPV >= 36){
			FillTH1(id_r1dhiPUHist, r);
			FillTH1(id_r1dwidehiPUHist, r);
			if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
				FillTH1(id_r1dhiPUcutHist, r);
				FillTH1(id_r1dwidehiPUcutHist, r);
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
                    FillTH1(id_mggHist, m12);
                    if(vcuts[i]&&vcuts[j])FillTH2(id_mgg2Hist, m12, numberOfPV);
                    if(vcuts[i]&&vcuts[j])FillTH1(id_mggCutHist, m12);
                }
                if (vcuts[i]&&vcuts[j]&&min(Ei,Ej)>2.0){
                    FillTH1(id_mggallHist, m12);
                    FillTH2(id_mggRCutHist, m12, Ri);
                    FillTH2(id_mggRCutHist, m12, Rj);
                }

            }
        }
    }//end numpc
		
	for(int i=0; i<PC_vTrack_pt.GetSize(); i++){
        for(int j=0; j<PC_vTrack_pt[i].size(); j++){
			px = PC_vTrack_pt[i][j] * cos( PC_vTrack_phi[i][j] );
			py = PC_vTrack_pt[i][j] * sin( PC_vTrack_phi[i][j] );
			pz = PC_vTrack_pt[i][j] * sinh( PC_vTrack_eta[i][j] );
            FillTH1(id_ptHist, PC_vTrack_pt[i][j], 1./PC_vTrack_pt[i][j]);
			FillTH1(id_pzHist, pz);
			FillTH2(id_pxpyHist, px, py);
        }
    }
//    cout << " vcuts length " << vcuts.size() << endl;
    vcuts.clear();
    
}
#endif
