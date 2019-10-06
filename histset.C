#ifndef HISTS
#define HISTS


#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "ROOT/TThreadedObject.hxx"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "myselector.C"

class histset{
	
	public:
	   double pi =4.0*atan(1.0);

	   histset();	
	   void init(); 

	   void AnalyzeEntry(myselector& s); 
	 //  void AnalyzeEntry(mytreevalues& s);

	   //bookeeping enumeration: (if we do this we dont need to worry about hist ptr copies and merging)
	   enum th1d_index{ind_ptHist, ind_pzHist, ind_numpcHist,ind_rerrHist, ind_phierrHist, ind_zerrHist,ind_r1dHist, ind_r1dcutHist, ind_r1dlowPUHist, ind_r1dhiPUHist, ind_r1dlowPUcutHist, ind_r1dhiPUcutHist, ind_rhobpHist, ind_mgg1Hist, numTH1Hist};
	   enum th2d_index{ind_pxpyHist,ind_xyHist,ind_rphiHist, ind_rzHist, ind_xycutHist, numTH2Hist};


	   // make a big vector and load enumerated histograms onto the vector
	 std::vector<ROOT::TThreadedObject<TH1D>* >  TH1Manager{};
	 std::vector<ROOT::TThreadedObject<TH2D>* >  TH2Manager{};


	  //locate the histogram and perform ptr copying 
	  void FillTH1(int index, double x, double w);
	  void FillTH2(int index, double x, double y);
	
	  void WriteHist(); 
};

histset::histset(){

	 std::vector<ROOT::TThreadedObject<TH1D>*>  Manager1(numTH1Hist);
	TH1Manager=Manager1;

	 std::vector<ROOT::TThreadedObject<TH2D>*>  Manager2(numTH2Hist);
	TH2Manager=Manager2;

	init();

}
void histset::init(){
//init TH1D
	TH1Manager.at(ind_ptHist) = new ROOT::TThreadedObject<TH1D>("ptHist", "p_{T} Distribution;p_{T};dN/p_{T}dp_{T}", 100, 0, 5);
	TH1Manager.at(ind_pzHist) = new ROOT::TThreadedObject<TH1D>("pzHist", "p_{Z} Distribution;p_{Z};dN/dp_{Z}", 100, 0, 5);
	TH1Manager.at(ind_numpcHist) = new ROOT::TThreadedObject<TH1D>("numpcHist", "Number of PC", 51,-0.5, 50.5);

	TH1Manager.at(ind_rerrHist) = new ROOT::TThreadedObject<TH1D>("rerrHist", "Conversion Radial Error; #Delta R (cm); Entries Per 0.05 bin", 40, 0.0, 2.0);
	TH1Manager.at(ind_phierrHist) = new ROOT::TThreadedObject<TH1D>("phierrHist", "Conversion Azimuthal Error;#Delta #phi; Entries Per 0.002 bin",50, 0.0, 0.1);
	TH1Manager.at(ind_zerrHist) = new ROOT::TThreadedObject<TH1D>("zerrHist","Conversion Z Error;#Delta z (cm); Entries Per 0.1 bin", 50, 0.0, 5.0);

	TH1Manager.at(ind_r1dHist) = new ROOT::TThreadedObject<TH1D>("r1dHist","Conversion Radius No Cuts;R (cm);Entries Per 0.1 bin",100, 0.0, 10.0);
	TH1Manager.at(ind_r1dcutHist) = new ROOT::TThreadedObject<TH1D>("r1dcutHist","Conversions Radius With Cuts; R (cm); Entries Per 0.1 bin", 100, 0.0, 10.0);

	TH1Manager.at(ind_r1dlowPUHist) = new ROOT::TThreadedObject<TH1D>("r1dlowPUHist","Conversion Radius: No Quality Cuts, PV #leq 16;R (cm);Entries Per 0.1 bin",100,0.,10.);
	TH1Manager.at(ind_r1dhiPUHist) = new ROOT::TThreadedObject<TH1D>("r1dhiPUHist","Conversion Radius: No Quality Cuts, PV #geq 36;R (cm);Entries Per 0.1 bin",100,0.,10.);
	TH1Manager.at(ind_r1dlowPUcutHist) = new ROOT::TThreadedObject<TH1D>("r1dlowPUcutHist","Conversion Radius: Quality Cuts, PV #leq 16;R (cm);Entries Per 0.1 bin",100,0.,10.);
	TH1Manager.at(ind_r1dhiPUcutHist) = new ROOT::TThreadedObject<TH1D>("r1dhiPUcutHist","Conversion Radius: Quality Cuts, PV #geq 36;R (cm);Entries Per 0.1 bin",100,0.,10.);
	TH1Manager.at(ind_rhobpHist) = new ROOT::TThreadedObject<TH1D>("rhobpHist","Conversion Radius w.r.t Beam Pipe and Quality Cuts; R (cm); Entries Per 0.05 bin",100,0.,5.);

	TH1Manager.at(ind_mgg1Hist) = new ROOT::TThreadedObject<TH1D>("mgg1Hist","Di-#gamma Mass;Mass GeV; Entries Per 2.5 MeV bin", 100, 0.0, 0.25 );

// init TH2D
	TH2Manager.at(ind_pxpyHist) = new ROOT::TThreadedObject<TH2D>("pxpyHist", "p_{X} vs p_{Y} Distribution;p_{X};p_{Y}", 200, -10., 10., 200, -10., 10.);

	TH2Manager.at(ind_xyHist) = new ROOT::TThreadedObject<TH2D>("xyHist", "Conversion Vertices per mm^{2} bin; -PC_x (cm); PC_y (cm)",200,-10.,10.,200,-10.,10.);
	TH2Manager.at(ind_rphiHist) = new ROOT::TThreadedObject<TH2D>("rphiHist", "Conversion Vertices in R-#phi per mm*60mrad bin; R (cm); #phi",100,0.,10.,100,-3.2,3.2);
	TH2Manager.at(ind_rzHist) = new ROOT::TThreadedObject<TH2D>("rzHist", "Conversion Vertices in R-z per mm^{2} bin; PC_z (cm); R (cm)",200,-10.,10.,100,0.,10.);

	TH2Manager.at(ind_xycutHist) = new ROOT::TThreadedObject<TH2D>("xycutHist", "Conversion Vertices per mm^{2} bin; -PC_x (cm); PC_y (cm)",200,-10.,10.,200,-10.,10.);

}
void histset::FillTH1(int index, double x, double w=1){
	//we must make ptr copies for performance reasons when trying to fill a histogram
	auto myhist = TH1Manager.at(index)->Get();
	myhist->Fill(x,w);
}
void histset::FillTH2(int index, double x, double y){
	auto myhist = TH2Manager.at(index)->Get();
	myhist->Fill(x,y);
}
void histset::WriteHist(){

	//TFile* outfile = new TFile("MCPCBigSet.root","RECREATE");
	TFile* outfile = new TFile("PCBigSet.root", "RECREATE");

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
   	
	//always make a local copy, if its a value dereference.. if you dont do this scope/dereferencing will get really weird, clunky, and unmanageable
	//have to auto& or myreader will try to register copy of the readerarray ptr
	auto& PC_vTrack_pt = s.PC_vTrack_pt;
	auto& PC_vTrack_phi = s.PC_vTrack_phi;
	auto& PC_vTrack_eta = s.PC_vTrack_eta;

	auto& PC_vtx_sigmaxx = s.PC_vtx_sigmaxx;
	auto& PC_vtx_sigmaxy = s.PC_vtx_sigmaxy;
	auto& PC_vtx_sigmayy = s.PC_vtx_sigmayy;
	auto& PC_vtx_sigmazz = s.PC_vtx_sigmazz;

	auto& PC_fitmomentumOut_theta = s.PC_fitmomentumOut_theta;

	auto& PC_x = s.PC_x;
	auto& PC_y = s.PC_y;
	auto& PC_z = s.PC_z;

	auto& PC_vtx_chi2 = s.PC_vtx_chi2;

	auto numberOfPC = *(s.numberOfPC);
	auto numberOfPV = *(s.numberOfPV);

	double px,py,pz;
	double x,y,z;
	double r,theta,phi;
	
	const double RERRCUT = 0.25;
	const double COSTCUT = 0.85;
	const double ZCUT = 25.0;
	const double FITPROBCUT = 0.010;

	double fitprob;

	//beam pipe displacement
	double x0 =  0.171;
    	double y0 = -0.176;
   	double rho;


	//error calcs
	double vxx,vxy,vyy,vzz; //variances
	double varsum_r, varsum_phi; //intermediate calculation variables
	double rerr,phierr,zerr;//errors
	double sphi,cphi;

	FillTH1(ind_numpcHist, numberOfPC);


	for(int i=0; i<PC_x.GetSize(); i++){
		x = PC_x[i];
		y = PC_y[i];
		z = PC_z[i];

		r = sqrt( x*x + y*y );
		phi = atan2(y, x);
		theta = PC_fitmomentumOut_theta[i];
		
		FillTH2(ind_xyHist, -x, y);
		FillTH2(ind_rphiHist, r, phi);
		FillTH2(ind_rzHist, z, r);	


		vxx = PC_vtx_sigmaxx[i];
		vxy = PC_vtx_sigmaxy[i];
		vyy = PC_vtx_sigmayy[i];
		vzz = PC_vtx_sigmazz[i];
//BAD
//		varsum_r   = ( vxx/(x*x) ) + 2.0*( vxy/(x*y) ) + ( vyy/(y*y) );
//       		varsum_phi = ( vxx/(y*y) ) - 2.0*( vxy/(x*y) ) + ( vyy/(x*x) );

		cphi = cos(phi);
		sphi = sin(phi);

		// This is the correct one
		varsum_r   = cphi*cphi*vxx + 2.0*sphi*cphi*vxy + sphi*sphi*vyy;
		varsum_phi = sphi*sphi*vxx - 2.0*sphi*cphi*vxy + cphi*cphi*vyy;
		                 


//		phierr = sqrt(varsum_phi); 1/r factor?
//       	rerr = r*sqrt(varsum_r); no r factor?
		rerr = sqrt(varsum_r);
		phierr = sqrt(varsum_phi)/r;

		zerr = sqrt(vzz);

		FillTH1(ind_rerrHist, rerr);
		FillTH1(ind_phierrHist, phierr);
		FillTH1(ind_zerrHist, zerr);

			
	 	fitprob = TMath::Prob(PC_vtx_chi2[i], 3);
		//make r plots with quality cuts and without
		
		FillTH1(ind_r1dHist, r);

		//make quality cuts
		if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
			FillTH1(ind_r1dcutHist, r);
			FillTH2(ind_xycutHist, x,y);

			rho =  sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0)) ;	
			FillTH1(ind_rhobpHist, rho); 

		}			
	
		//pileup cuts
		if( numberOfPV <= 16){
			FillTH1(ind_r1dlowPUHist, r);
			//low PU quality cuts
			 if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
				FillTH1(ind_r1dlowPUcutHist, r);
			}
		}
		if( numberOfPV >= 36){
			FillTH1(ind_r1dhiPUHist, r);
			 if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT && fitprob > FITPROBCUT){
				FillTH1(ind_r1dhiPUcutHist,r);
			}
		}	
				
	}


	auto& PC_fitmomentumOut_pt = s.PC_fitmomentumOut_pt;
	auto& PC_fitmomentumOut_phi = s.PC_fitmomentumOut_phi;

	//gamma gamma stuff
 	if (numberOfPC>=2){
        for(unsigned int i=0; i<numberOfPC-1; i++){
            double pxi = PC_fitmomentumOut_pt[i]*cos(PC_fitmomentumOut_phi[i]);
            double pyi = PC_fitmomentumOut_pt[i]*sin(PC_fitmomentumOut_phi[i]);
            double pzi = PC_fitmomentumOut_pt[i]/tan(PC_fitmomentumOut_theta[i]);
            double Ei = sqrt(pxi*pxi + pyi*pyi + pzi*pzi);
            double xi = PC_x[i];
            double yi = PC_y[i];
            double Ri = sqrt(xi*xi + yi*yi);
            for(unsigned int j=i; j<numberOfPC; j++){
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
                    //mgg->Fill(m12);
                    FillTH1(ind_mgg1Hist, m12);
                    //hmgg->Fill(m12);
                }
                //hmggall->Fill(m12);
            }
        }
        }//end numpc
		


	for(int i=0; i<PC_vTrack_pt.GetSize(); i++){
                for(int j=0; j<PC_vTrack_pt[i].size(); j++){

			px = PC_vTrack_pt[i][j] * cos( PC_vTrack_phi[i][j] );
			py = PC_vTrack_pt[i][j] * sin( PC_vTrack_phi[i][j] );
			pz = PC_vTrack_pt[i][j] * sinh( PC_vTrack_eta[i][j] );


                        FillTH1(ind_ptHist, PC_vTrack_pt[i][j], 1./PC_vTrack_pt[i][j]);
			FillTH1(ind_pzHist, pz);
			
			FillTH2(ind_pxpyHist, px, py);
                }
        }
	
	

}
#endif
