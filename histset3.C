#ifndef HISTS
#define HISTS
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "ROOT/TThreadedObject.hxx"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "recosim.C"
#include "myweights.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <map>
#include <tuple>
#include <iomanip>
#include "Hungarian.h"
using MyTH1D = ROOT::TThreadedObject<TH1D>;
using MyTH2D = ROOT::TThreadedObject<TH2D>;

// struct for derived quantities of each conversion
#include  "mystruct.h"

class histset3{
    public:
       double PI =4.0*atan(1.0);
       histset3();
       void init();
       void setweightoption();
       void AnalyzeEntry(recosim& s);
       #include "histset2enums.h"
// make a big vector and load enumerated histograms onto the vector
       std::vector<MyTH1D*>  TH1Manager{};
       std::vector<MyTH2D*>  TH2Manager{};
// locate the histogram and perform pointer copying
       void FillTH1(int index, double x, double w);
       void FillTH2(int index, double x, double y, double w);
       void WriteHist();
};

histset3::histset3(){
    std::vector<MyTH1D*>  Manager1(numTH1Hist);
    TH1Manager=Manager1;
    std::vector<MyTH2D*>  Manager2(numTH2Hist);
    TH2Manager=Manager2;
    init();
    setweightoption();
}

void histset3::setweightoption(){
    for(int i=0; i<numTH1Hist; i++){
        auto hptr = TH1Manager.at(i)->Get();
        hptr->Sumw2(kTRUE);
    }
    for(int i=0; i<numTH2Hist; i++){
        auto hptr = TH2Manager.at(i)->Get();
        hptr->Sumw2(kTRUE);
    }
}

#include "histset3init.h"     //Put the histogram declarations in one separate file

void histset3::FillTH1(int index, double x, double w=1.0){
	//we must make pointer copies for performance reasons when trying to fill a histogram
	auto myhist = TH1Manager.at(index)->Get();
	myhist->Fill(x,w);
}

void histset3::FillTH2(int index, double x, double y, double w=1.0){
	auto myhist = TH2Manager.at(index)->Get();
	myhist->Fill(x,y,w);
}

void histset3::WriteHist(){

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

void histset3::AnalyzeEntry(recosim& s){

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

    const bool lpr = false;       // print flag
    const bool lreduce = true;   // do problem reduction
    const bool lassign = true;   // Do assignment problem
    const bool lpr2 = false;      // another print flag
	
    const double RERRCUT = 0.25;
    const double COSTCUT = 0.85;
    const double ZCUT = 25.0;
    const double FITPROBCUT = 0.010;
    const double MASSCUT = 0.15;

    const double MASS_ELECTRON = 0.5109989461e-3;
    const double MASS_PION = 139.57061e-3;
    const double MASS_KAON = 493.677e-3;
    const double MASS_PROTON = 938.272081e-3;

// Scale MC to data based on number 
// of events with at least 1 conversion
    const double NFACTOR = 3.68068361*1.00000128;
    double fitprob;
// We now have various "centers" to compare to for radial coordinates.
// beam pipe displacement (in cm) from Anna's DPF2019 talk
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
    double eta;
    
    double wtPU;

    std::vector<bool> vcuts;
    std::vector<double> PosTkInfo;
    std::vector<double> NegTkInfo;
    std::vector<double> PosPt;
    std::vector<double> NegPt;     
    std::vector<int> vcandidate;
    std::map<double, int> mNeg;
    std::multimap<double, int> mmNeg;
    std::map<double, int> mPos;
    std::multimap<double, int> mmPos;
    std::multimap<int, double> mmCandidate;
    std::set <std::pair<double,double> > trkPair;

//    if(eventNumber==772775){
//    if(eventNumber==773069){
//    if(eventNumber==764002){
//    if(eventNumber==764127){
    if(eventNumber!=-1){
//       if(lpr)std::cout << "Skipping this event" << std::endl;
//       goto SKIP;
//    }

/*    if(lpr)std::cout << " " << std::endl;
    if(lpr)std::cout << " ------------BoE-------------- " << eventNumber << std::endl;
    if(lpr)std::cout << " " << std::endl; */

    std::cout << " " << std::endl;
    std::cout << " ------------BoE-------------- " << eventNumber << std::endl;
    std::cout << " " << std::endl;


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
//       return;
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

// Probably a good idea to keep track of whether each conversion 
// candidate passes particular cuts, to ease later code 
// with gamma-gamma invariant mass

    int nassigned = 0;
  
// Vector of derived quantities managed using a struct
    std::vector<mytuple> tup;

	for(int i=0; i<PC_x.GetSize(); i++){
// push-back new mytuple struct with default constructor
        tup.push_back(mytuple());

        vcuts.push_back(false);
        PosTkInfo.push_back(-1.0);
        NegTkInfo.push_back(-1.0);
        PosPt.push_back(-1.0);
        NegPt.push_back(-1.0);
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
        eta = -log(tan(theta/2.0));
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

// This logic is a bit confused. It would likely be less confusing if it was reliably -ve/+ve.
// Let's try using a deque. I could have a deque of vectors that is ordered by -ve then +ve.
        std::deque < std::vector <double> > Mom;
        std::vector <double> momentum0 = {px0p, py0p, pz0p};
        std::vector <double> momentum1 = {px1p, py1p, pz1p};
// Fill Mom with the -ve then +ve track 3-momenta
        Mom.push_back(momentum0);
        if(q1>0){
           Mom.push_back(momentum1);
        }
        else{
           Mom.push_front(momentum1);
        }

        ROOT::Math::PxPyPzMVector vNege,vNegpi,vNegk,vNegp;
        ROOT::Math::PxPyPzMVector vPose,vPospi,vPosk,vPosp;
        vNege  = ROOT::Math::PxPyPzMVector( Mom[0][0], Mom[0][1], Mom[0][2], MASS_ELECTRON );
        vNegpi = ROOT::Math::PxPyPzMVector( Mom[0][0], Mom[0][1], Mom[0][2], MASS_PION );
        vNegk  = ROOT::Math::PxPyPzMVector( Mom[0][0], Mom[0][1], Mom[0][2], MASS_KAON );
        vNegp  = ROOT::Math::PxPyPzMVector( Mom[0][0], Mom[0][1], Mom[0][2], MASS_PROTON );
        vPose  = ROOT::Math::PxPyPzMVector( Mom[1][0], Mom[1][1], Mom[1][2], MASS_ELECTRON );
        vPospi = ROOT::Math::PxPyPzMVector( Mom[1][0], Mom[1][1], Mom[1][2], MASS_PION );
        vPosk  = ROOT::Math::PxPyPzMVector( Mom[1][0], Mom[1][1], Mom[1][2], MASS_KAON );
        vPosp  = ROOT::Math::PxPyPzMVector( Mom[1][0], Mom[1][1], Mom[1][2], MASS_PROTON );

        ROOT::Math::PxPyPzMVector vpair, vpairpipi, vpairkk, vpairpip, vpairppi, vpairLambda, vpairALambda;
        vpair += vNege;
        vpair += vPose;
        vpairpipi += vNegpi;
        vpairpipi += vPospi;
        vpairkk += vNegk;
        vpairkk += vPosk;
        vpairLambda += vPosp;
        vpairLambda += vNegpi;
        vpairALambda += vNegp;
        vpairALambda += vPospi;

// Make primary quality cuts
		if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT 
                && fitprob > FITPROBCUT && std::max(nBefore0,nBefore1)==0 ){
            vcuts[i] = true;
// Keep track of charge-signed chisq/dof of constituent tracks to later identify conversion 
// candidates that use the same track (sign by charge of track)
            if(q0 == 1 && q1 == -1){
               PosTkInfo[i] =  Tk0_chi2[i]/double(Tk0_ndof[i]);
               NegTkInfo[i] = -Tk1_chi2[i]/double(Tk1_ndof[i]);
               PosPt[i] = sqrt(Tk0_px[i]*Tk0_px[i] + Tk0_py[i]*Tk0_py[i]);
               NegPt[i] = sqrt(Tk1_px[i]*Tk1_px[i] + Tk1_py[i]*Tk1_py[i]);
//               PosTkInfo[i] =  abs(Tk0_sd0[i])*Tk0_chi2[i]/double(Tk0_ndof[i]);
//               NegTkInfo[i] = -abs(Tk1_sd0[i])*Tk1_chi2[i]/double(Tk1_ndof[i]);
            }
            else{
               PosTkInfo[i] =  Tk1_chi2[i]/double(Tk1_ndof[i]);
               NegTkInfo[i] = -Tk0_chi2[i]/double(Tk0_ndof[i]);
               PosPt[i] = sqrt(Tk1_px[i]*Tk1_px[i] + Tk1_py[i]*Tk1_py[i]);
               NegPt[i] = sqrt(Tk0_px[i]*Tk0_px[i] + Tk0_py[i]*Tk0_py[i]);
//               PosTkInfo[i] =  abs(Tk1_sd0[i])*Tk1_chi2[i]/double(Tk1_ndof[i]);
//               NegTkInfo[i] = -abs(Tk0_sd0[i])*Tk0_chi2[i]/double(Tk0_ndof[i]);
            }
// Need to make sure that each track pair is distinguishable from pairs already selected
            if(!trkPair.insert(std::make_pair(NegTkInfo[i], PosTkInfo[i])).second){
               if(lpr)std::cout << "INDISTINGUISHABLE EDGE " << i << " ignored in event " << eventNumber << std::endl;
// good idea to keep a tally in some histogram bin
               vcuts[i] = false;
            }
            else{
// Fill STL containers that help define the "matching problem"
               vcandidate.push_back(i);
// For the key/value pair use charge-signed chi2/dof as key, and conversion index as edge id for value
               mNeg.insert(std::make_pair(NegTkInfo[i],i));
               mmNeg.insert(std::make_pair(NegTkInfo[i],i));
               mPos.insert(std::make_pair(PosTkInfo[i],i));
               mmPos.insert(std::make_pair(PosTkInfo[i],i));
// Fill tuple with relevant derived quantities
               tup[i].radius = r;
               tup[i].rerr = rerr;
               tup[i].pfit = fitprob;
               tup[i].alpha = APalpha;
               tup[i].qpt0 = APpT0;
               tup[i].qpt1 = APpT1;
               tup[i].Pos = PosTkInfo[i];
               tup[i].Neg = NegTkInfo[i];
               tup[i].phierr = phierr;
               tup[i].zerr = zerr;
               tup[i].mass = {vpair.M(), vpairpipi.M(), vpairkk.M(), vpairLambda.M(), vpairALambda.M()};
// Check
               if(lpr)std::cout << "Masses: " << tup[i].mass[0] << " " << tup[i].mass[1] << " " 
                                << tup[i].mass[2] << " " << tup[i].mass[3] << " " << tup[i].mass[4] << std::endl;
               tup[i].rho = rho;
               tup[i].rps = rps;
               tup[i].rnominal = rnominal;
               tup[i].phi = phi;
               tup[i].pt = pt;
               tup[i].E = E;
               tup[i].phip = phip;
               double ptasym = (pt0-pt1)/(pt0+pt1);
               if (q0<0) ptasym = -ptasym;
               double xplus = (1.0 + ptasym)/2.0;
               tup[i].xplus = xplus;
               tup[i].theta = theta;
               tup[i].ptPos = PosPt[i];
               tup[i].ptNeg = NegPt[i];
               tup[i].etaphys = eta;
               tup[i].id = i;
               tup[i].x = x;
               tup[i].y = y;
               tup[i].z = z;
// for technical checks
               FillTH1( id_r1dHist2, tup[i].radius );
            }
        }
// Also include correlation plot between APalpha and ptasym.
	}

// Characterize our matching problem for this event
// Note this is unnecessary if there are no duplicates, ie. n- = n+ = nedges.
    std::vector<int> vsel;   // Vector for indices of selected conversions after removing duplicates

    if(std::min(mNeg.size(), mPos.size()) < vcandidate.size() && lassign ){
// We actually have an assignment problem to worry about and we want to worry about it
       if(lpr){
          std::cout << " " << std::endl;
          std::cout << "Event " << eventNumber << " numberOfPC " << numberOfPC << std::endl;
          std::cout << "nedges = " << vcandidate.size() << std::endl;
          std::cout << "N- = " << mmNeg.size() << std::endl;
          std::cout << "N+ = " << mmPos.size() << std::endl;
          std::cout << "n- = " << mNeg.size() << std::endl;
          std::cout << "n+ = " << mPos.size() << std::endl;
          std::cout << "Target maximum possible cardinality of solution = " << std::min(mNeg.size(),mPos.size()) << std::endl;
// Prior solution with no arbitration
          std::cout << "Prior edge solution (no arbitration at all) : ";
          for(int i=0; i<numberOfPC; ++i){
              if(vcuts[i])std::cout << std::setw(3) << i; 
          }
          std::cout << std::endl;
       }
// Let's do some more characterization of the ambiguity complexity by investigating the multimaps
// with a view to removing the non-ambiguities.
// Note that in case of count==1 from the multimap, the map value is the unique edge ID for this polarity of track 
       if(lreduce){
          for (auto i = mNeg.begin(); i!= mNeg.end(); ++i) {
              auto tkInfo = i->first;
              auto tkEdge = i->second;
              auto negCount = mmNeg.count(tkInfo);
              if(lpr){
                  std::cout << negCount << " edge(s) for e- with first key,value=(" 
                            << tkInfo << "," << tkEdge << ") [edges: ";
// Example from http://www.cplusplus.com/reference/map/multimap/count/
// equal_range returns a pair with lower_bound and upper_bound positions.
                  for (auto it=mmNeg.equal_range(tkInfo).first; it!=mmNeg.equal_range(tkInfo).second; ++it){
                       std::cout << ' ' << (*it).second;
                  }
                  std::cout << " ] " << std::endl;
              }
              if(negCount==1){
// Make candidate multimap when only one possibility for this electron. The key is the edge id. 
                 mmCandidate.insert(std::make_pair(tkEdge, tkInfo));
              }
          }
          for (auto i = mPos.begin(); i!= mPos.end(); ++i) {
              auto tkInfo = i->first;
              auto tkEdge = i->second;
              auto posCount = mmPos.count(tkInfo);
              if(lpr){
                  std::cout << posCount << " edge(s) for e+ with first key,value=(" 
                            << tkInfo << "," << tkEdge << ") [edges: ";
                  for (auto it=mmPos.equal_range(tkInfo).first; it!=mmPos.equal_range(tkInfo).second; ++it){
                       std::cout << ' ' << (*it).second;
                  }
                  std::cout << " ] " << std::endl;
              }
              if(posCount==1){
                 mmCandidate.insert(std::make_pair(tkEdge, tkInfo));
              }
          }
          if(lpr)std::cout << "mmCandidate multimap size " << mmCandidate.size() << std::endl;
          for (auto i = mmCandidate.begin(); i!= mmCandidate.end(); ++i) {
               auto tkEdge = i->first;
               auto tkInfo = i->second;
               if(lpr)std::cout << "mmCandidate key= " << tkEdge << " multiplicity " 
                                << mmCandidate.count(tkEdge) << " (value " << tkInfo << " ) " << std::endl;
               if(mmCandidate.count(tkEdge) == 2 ){
// There is an electron-positron pairing where the degree of each vertex is 1, and the same edge is incident on both.
// So if we add this edge to the selected candidates we can erase this one and its constituents from the assignment problem.
                  if(tkInfo<0.0){
                     vsel.push_back(tkEdge);
                     vcandidate.erase(remove(vcandidate.begin(),vcandidate.end(),tkEdge),vcandidate.end());
                     mNeg.erase(tkInfo);
                  }
                  else{
                     mPos.erase(tkInfo);
                  }
               }
          }
          if(lpr){
              std::cout << "Already selected " << vsel.size() << " non-ambiguous pairings " << std::endl;
              std::cout << "After reduction, nedges = " << vcandidate.size() << std::endl;
              std::cout << "n- = " << mNeg.size() << std::endl;
              std::cout << "n+ = " << mPos.size() << std::endl;

          }
      }  // end of lreduce clause

// Set up the (n- * n+) cost matrix as input to the Hungarian Algorithm after reduction.
// (this works for both nRows <= nCols and nRows > nCols)

    std::vector< std::vector <double> > costMatrix;
    std::vector< std::vector <int> > edgeMatrix;
    int irow = -1;
    for ( auto i = mNeg.begin(); i != mNeg.end(); ++i ){        // n- rows for each -ve track
         irow++;
         std::vector<double> v;
         std::vector<int> e;
         int nfound = 0;
         for ( auto j = mPos.begin(); j != mPos.end(); ++j ){   // n+ cols for each +ve track
// Go through all the edge candidates and see if there is an edge corresponding to this pairing.
             double negInfo = i->first;
             double posInfo = j->first;
             bool found = false;
             for (auto iter = vcandidate.begin(); iter != vcandidate.end(); ++iter) {
                 unsigned int k = *iter;
                 if(NegTkInfo[k] == negInfo && PosTkInfo[k] == posInfo ) {
// Found matching edge
                    if(lpr)std::cout << "Found matching edge " << k << std::endl;
// Add cost value of corresponding chi-squared value for 1 dof.
                    v.push_back(TMath::ChisquareQuantile(1.0-tup[k].pfit,1.0));
                    e.push_back(k);
                    found=true;
                    nfound++;
                 }
             }
             if(!found){
                v.push_back(10000.0);
                e.push_back(-1);
             }
         }
         if(lpr)std::cout << "irow " << irow << " nedges = " << nfound << std::endl;
         costMatrix.push_back(v);
         edgeMatrix.push_back(e);
    }
// Debug printing
    if(lpr){
       for (auto irow = costMatrix.begin(); irow != costMatrix.end(); ++irow) {
           std::cout << "Row weights:   ";
           for (auto pos = irow->begin(); pos != irow->end(); ++pos) {
               std::cout << std::setw(10) << *pos << " ";
           }
           std::cout << std::endl;
       }
       for (auto irow = edgeMatrix.begin(); irow != edgeMatrix.end(); ++irow) {
           std::cout << "Edges      :   ";
           for (auto pos = irow->begin(); pos != irow->end(); ++pos) {
               std::cout << std::setw(10) << *pos << " ";
           }
           std::cout << std::endl;
       }
    }

// Use Hungarian Algorithm to find the minimum cost assignment of e- to e+. 
// The total cost will be the total chi-squared of all assignments (each with 1 dof).
// Need to take care also of the case where there is no one-sided perfect matching, 
// and algorithmically, the extra assignments are assigned to the fictional 
// high-cost edges with nominal weight of 10000.
	HungarianAlgorithm HungAlgo;
	std::vector<int> assignment;
	double cost = HungAlgo.Solve(costMatrix, assignment);
    if(lpr)std::cout << "costMatrix.size() " << costMatrix.size() << std::endl;
	for (unsigned int x = 0; x < costMatrix.size(); x++){
// Need to check if this row is assigned (may not be if nRows > nCols)
        if(assignment[x] >= 0){
           if(edgeMatrix[x][assignment[x]] == -1){
              cost-=10000.0;
           }
           else{
              vsel.push_back(edgeMatrix[x][assignment[x]]);
              nassigned++;
           }
           if(lpr)std::cout << x << "," << assignment[x] 
                     << " " << costMatrix[x][assignment[x]] << " " 
                     << edgeMatrix[x][assignment[x]] << std::endl;
        }
    }
    if(nassigned > 0){
       if(lpr)std::cout << "Assigned " << nassigned << " initially ambiguous pairings" << std::endl;
       double psel = TMath::Prob(cost, nassigned);
	   if(lpr)std::cout << "Minimized total chisq: " << cost << " ( " << nassigned << " ) " << " p-value " << psel <<std::endl;
    }
    if(vsel.size() > 0){
       std::sort(vsel.begin(), vsel.end());
       if(lpr){
          std::cout << "Selected conversions (" << vsel.size() << ")" ;
          for(int i=0; i<vsel.size(); ++i){
              std::cout << std::setw(3) << " " << vsel[i]; 
          }
          std::cout << std::endl;
       }
    }
 }
  else{
// No ambiguities - so no assignment problem to solve 
     if(vcandidate.size() > 0){
        vsel = vcandidate;
        std::sort(vsel.begin(), vsel.end());
        if(lpr){
           std::cout << " " << std::endl;
           std::cout << "Selected conversions S: " << eventNumber << "  (" << vsel.size() << ")";
           for(int i=0; i<vsel.size(); ++i){
              std::cout << std::setw(3) << " " << vsel[i]; 
           }
           std::cout << std::endl;
        }
     }
  }

// Now with duplicates eliminated histograms for selected candidates
    FillTH1( id_nconvHist, vsel.size() );
    FillTH1( id_nassignedHist, nassigned );
    FillTH1( id_nnonassignedHist, vsel.size() - nassigned );

    for(unsigned int j=0; j<vsel.size(); ++j){

        int i = vsel[j];
        r = tup[i].radius;
        bool region1 = (r>3.7) && (r<6.2);
        bool region2 = (r>7.6) && (r<10.3);
        bool region3 = (r>11.7) && (r<15.3);

        FillTH1( id_r1dHist4, tup[i].radius );
        if( numberOfPV <= 16)FillTH1(id_r1dwidelowPUcutHist, tup[i].radius, wtPU);
        if( numberOfPV > 16 && numberOfPV < 36)FillTH1(id_r1dwidemedPUcutHist, tup[i].radius, wtPU);
        if( numberOfPV >= 36)FillTH1(id_r1dwidehiPUcutHist, tup[i].radius, wtPU);

        FillTH1(id_AP_pTminHist, std::min(tup[i].qpt0,tup[i].qpt1), wtPU);
        FillTH1(id_AP_pTmaxHist, std::max(tup[i].qpt0,tup[i].qpt1), wtPU);
        FillTH1(id_AP_pTaveHist, 0.5*(tup[i].qpt0+tup[i].qpt1), wtPU);
        FillTH1(id_AP_alphaHist, tup[i].alpha, wtPU);
        FillTH2(id_AP_pT_alphaHist, tup[i].alpha, std::min(tup[i].qpt0,tup[i].qpt1), wtPU);

        FillTH1(id_zcutHist, PC_z[i], wtPU);
        FillTH1(id_zPVHist, PC_zPV[i], wtPU);
        FillTH1(id_r1dcutHist, r, wtPU);
        FillTH1(id_r1dwidecutHist, r);
        FillTH1(id_r1dwidecutWHist, r, wtPU);

        FillTH1(id_r1dwidecutPSHist, tup[i].rps, wtPU);
        FillTH1(id_r1dwidecutNomHist, tup[i].rnominal, wtPU);
        if(PC_dmin[i]>-998.0){
   		   FillTH1(id_r1dwidecutDHist, r, wtPU);
        }
        if(PC_dmin[i]>0.0){
   		   FillTH1(id_r1dwidecutDDHist, r, wtPU);
        }
        FillTH2(id_rphiHist, r, tup[i].phi, wtPU);
		FillTH2(id_xycutHist, PC_x[i], PC_y[i], wtPU);
		FillTH2(id_xywidecutHist, PC_x[i], PC_y[i], wtPU);
        FillTH1(id_phiHist, tup[i].phi, wtPU);	
		FillTH1(id_rhobpHist, tup[i].rho, wtPU);
		FillTH1(id_rbpHist, r, wtPU);
        FillTH1(id_rnomHist, tup[i].rnominal, wtPU);
		FillTH2(id_rhophiHist, tup[i].rho, tup[i].phip, wtPU);
        FillTH2(id_npv_rcutHist, r, numberOfPV, wtPU);
        FillTH1(id_pTHist, tup[i].pt, wtPU);
        FillTH1(id_EHist, tup[i].E, wtPU);

// Other histograms
        FillTH1(id_rerrHist, tup[i].rerr, wtPU);
        FillTH1(id_phierrHist, tup[i].phierr, wtPU);
        FillTH1(id_zerrHist, tup[i].zerr, wtPU);
        FillTH1(id_pfitHist, tup[i].pfit, wtPU);
        FillTH1(id_zHist, PC_z[i], wtPU);
        FillTH1(id_costhetaHist, cos(tup[i].theta), wtPU);
        FillTH1(id_etaHist, tup[i].etaphys, wtPU);    
    	FillTH2(id_xyHist, PC_x[i], PC_y[i], wtPU);
        FillTH2(id_xywideHist, PC_x[i], PC_y[i], wtPU);
        FillTH2(id_rzHist, PC_z[i], r, wtPU);
        FillTH1(id_pTHist2, tup[i].pt, wtPU);
        FillTH1(id_EHist2, tup[i].E, wtPU);
        FillTH1(id_phiHist2, tup[i].phi, wtPU);
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

        if(tup[i].pt<=1.0)FillTH1(id_xplus0Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>1.0)FillTH1(id_xplus1Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>2.0)FillTH1(id_xplus2Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>4.0)FillTH1(id_xplus4Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>8.0)FillTH1(id_xplus8Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>16.0)FillTH1(id_xplus16Hist, tup[i].xplus, wtPU);

        if(tup[i].pt<=1.0)FillTH1(id_alpha0Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>1.0)FillTH1(id_alpha1Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>2.0)FillTH1(id_alpha2Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>4.0)FillTH1(id_alpha4Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>8.0)FillTH1(id_alpha8Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>16.0)FillTH1(id_alpha16Hist, tup[i].alpha, wtPU);
        FillTH2(id_pTalphaHist, tup[i].alpha, tup[i].pt, wtPU);

        if(tup[i].pt>1.0)FillTH1(id_r1Hist, r, wtPU);
        if(tup[i].pt>2.0)FillTH1(id_r2Hist, r, wtPU);
        if(tup[i].pt>4.0)FillTH1(id_r4Hist, r, wtPU);
        if(tup[i].pt>8.0)FillTH1(id_r8Hist, r, wtPU);

        if(tup[i].pt<=1.0)FillTH1(id_r0Hist, r, wtPU);
        if(tup[i].pt>1.0&&tup[i].pt<=2.0)FillTH1(id_r1Hist2, r, wtPU);
        if(tup[i].pt>2.0&&tup[i].pt<=4.0)FillTH1(id_r2Hist2, r, wtPU);
        if(tup[i].pt>4.0&&tup[i].pt<=8.0)FillTH1(id_r4Hist2, r, wtPU);
        if(tup[i].pt>8.0&&tup[i].pt<=16.0)FillTH1(id_r8Hist2, r, wtPU);
        if(tup[i].pt>16.0)FillTH1(id_r16Hist, r, wtPU);

        if(tup[i].pt<=1.0)FillTH1(id_m0Hist, tup[i].mass[0], wtPU);
        if(tup[i].pt>1.0&&tup[i].pt<=2.0)FillTH1(id_m1Hist2, tup[i].mass[0], wtPU);
        if(tup[i].pt>2.0&&tup[i].pt<=4.0)FillTH1(id_m2Hist2, tup[i].mass[0], wtPU);
        if(tup[i].pt>4.0&&tup[i].pt<=8.0)FillTH1(id_m4Hist2, tup[i].mass[0], wtPU);
        if(tup[i].pt>8.0&&tup[i].pt<=16.0)FillTH1(id_m8Hist2, tup[i].mass[0], wtPU);
        if(tup[i].pt>16.0)FillTH1(id_m16Hist, tup[i].mass[0], wtPU);

        if(tup[i].pt>16.0 && abs(tup[i].alpha)>0.96){
           FillTH1(id_rAsymmetricHist, r, wtPU);
           FillTH1(id_convMassAnomaly, tup[i].mass[0], wtPU);
        }

        FillTH1(id_minptHist, std::min(tup[i].ptPos,tup[i].ptNeg), wtPU);
        FillTH1(id_maxptHist, std::max(tup[i].ptPos,tup[i].ptNeg), wtPU);

        if(region1||region2||region3){
           FillTH1(id_alphaBkgdHist, tup[i].alpha, wtPU);
           if(region1)FillTH1(id_alphaBkgdHistR1, tup[i].alpha, wtPU);
           if(region2)FillTH1(id_alphaBkgdHistR2, tup[i].alpha, wtPU);
           if(region3)FillTH1(id_alphaBkgdHistR3, tup[i].alpha, wtPU);
        }
        else{
           FillTH1(id_alphaSignalHist, tup[i].alpha, wtPU);
        }
        FillTH1(id_conversionCandidateMassHist, tup[i].mass[0], wtPU);
        FillTH1(id_conversionCandidateMassHist2, tup[i].mass[0], wtPU);
        FillTH1(id_conversionCandidateMassHist3, tup[i].mass[0], wtPU);
        FillTH1(id_KShortMassHist, tup[i].mass[1], wtPU);
        if(region1||region2||region3){
           FillTH1(id_KShortBkgdMassHist, tup[i].mass[1], wtPU);
           if(region1)FillTH1(id_KShortBkgdMassHistR1, tup[i].mass[1], wtPU);
           if(region2)FillTH1(id_KShortBkgdMassHistR2, tup[i].mass[1], wtPU);
           if(region3)FillTH1(id_KShortBkgdMassHistR3, tup[i].mass[1], wtPU);
           if(region1)FillTH1(id_convMassR1, tup[i].mass[0], wtPU);
           if(region2)FillTH1(id_convMassR2, tup[i].mass[0], wtPU);
           if(region3)FillTH1(id_convMassR3, tup[i].mass[0], wtPU);
        }
        else{
           FillTH1(id_convMassSignal, tup[i].mass[0], wtPU);
        }
        if(tup[i].alpha >= 0.0){
// Lambda candidate with a proton ( Lambda -> pi- p ) only populates +ve alpha
           FillTH1(id_lambdaCandidateMassHist, tup[i].mass[3], wtPU);
           FillTH1(id_lambdasCandidateMassHist, tup[i].mass[3], wtPU);
           if(region1||region2||region3){
              FillTH1(id_lambdasBkgdMassHist, tup[i].mass[3], wtPU);
              if(region1)FillTH1(id_lambdasBkgdMassHistR1, tup[i].mass[3], wtPU);
              if(region2)FillTH1(id_lambdasBkgdMassHistR2, tup[i].mass[3], wtPU);
              if(region3)FillTH1(id_lambdasBkgdMassHistR3, tup[i].mass[3], wtPU);
           }
           else{
              FillTH1(id_lambdasSignalMassHist, tup[i].mass[3], wtPU);
           }
        }
        else{
// Lambdabar candidate with an antiproton ( Lambdabar -> pi+ pbar) only populates -ve alpha
           FillTH1(id_lambdabarCandidateMassHist, tup[i].mass[4], wtPU);
           FillTH1(id_lambdasCandidateMassHist, tup[i].mass[4], wtPU);
           if(region1||region2||region3){
              FillTH1(id_lambdasBkgdMassHist, tup[i].mass[4], wtPU);
              if(region1)FillTH1(id_lambdasBkgdMassHistR1, tup[i].mass[4], wtPU);
              if(region2)FillTH1(id_lambdasBkgdMassHistR2, tup[i].mass[4], wtPU);
              if(region3)FillTH1(id_lambdasBkgdMassHistR3, tup[i].mass[4], wtPU);
           }
           else{
              FillTH1(id_lambdasSignalMassHist, tup[i].mass[4], wtPU);
           }
        } 

// DEBUG stuff
        if(lpr){
//           if(abs(tup[i].alpha) > 0.96){
//              std::cout << " BIGALPHA: Event " << eventNumber << std::endl;
              std::cout << " alpha: " << tup[i].alpha << std::endl;
              std::cout << " pt: " << tup[i].pt << std::endl;
              std::cout << " qpt0: " << tup[i].qpt0 << std::endl;
              std::cout << " qpt1: " << tup[i].qpt1 << std::endl;
              std::cout << " pfit: " << tup[i].pfit << std::endl;
              std::cout << " mee : " << tup[i].mass[0] << std::endl;
              std::cout << " radius : " << tup[i].radius << std::endl;
              std::cout << " pT- : " << tup[i].ptNeg << std::endl;
              std::cout << " pT+ : " << tup[i].ptPos << std::endl;
              std::cout << " (x,y,z): " << tup[i].x << " " << tup[i].y << " " << tup[i].z << std::endl;
//           }
        }

    }  // End of selected candidate loop

//   for(int i=0; i<Conv_vtxdl.GetSize(); i++){
//	vtxdl->Fill(Conv_vtxdl[i]);	
//   }

    auto& Conv_vtxdl = s.Conv_vtxdl;
    auto& Conv_Tk0_Idx = s.Conv_Tk0_Idx;
    auto& Conv_Tk1_Idx = s.Conv_Tk1_Idx;
    auto& Conv_convVtxIdx = s.Conv_convVtxIdx;
	auto nSimTrk = *(s.nSimTrk);
	auto nSimVtx = *(s.nSimVtx);

    auto& SimVtx_x = s.SimVtx_x;
    auto& SimVtx_y = s.SimVtx_y;
    auto& SimVtx_z = s.SimVtx_z;
    auto& SimVtx_tof = s.SimVtx_tof;
    auto& SimVtx_processType = s.SimVtx_processType;
    auto& SimVtx_simtrk_parent_tid = s.SimVtx_simtrk_parent_tid;

    auto& SimTrk_charge = s.SimTrk_charge;
    auto& SimTrk_eta = s.SimTrk_eta;
    auto& SimTrk_phi = s.SimTrk_phi;
    auto& SimTrk_pt = s.SimTrk_pt;
    auto& SimTrk_pdgId = s.SimTrk_pdgId;
    auto& SimTrk_simvtx_Idx = s.SimTrk_simvtx_Idx;
    auto& SimTrk_trackId = s.SimTrk_trackId;

// Maybe make some data structure for our own goals?

    if (vsel.size()>=1) {
       std::cout << "nSimVtx: " << nSimVtx << " nSimTrk: " << nSimTrk << std::endl; 
// Check SimVtx and SimTrk info
       for (unsigned int j=0; j<vsel.size(); j++){
           int i = vsel[j];
           std::cout << "Distance: " << Conv_vtxdl[i] << " " << Conv_convVtxIdx[i] << " " << Conv_Tk0_Idx[i] << " " << Conv_Tk1_Idx[i] << std::endl; 
       }
//       for (unsigned int i=0; i<nSimVtx; i++){
           unsigned int i = Conv_convVtxIdx[vsel[0]];
           std::cout << "SimVtx " << i << " " << SimVtx_x[i] << " " << SimVtx_y[i] << " " << SimVtx_z[i] << " " 
                     << SimVtx_tof[i] << " " << SimVtx_processType[i] << " " << SimVtx_simtrk_parent_tid[i] << std::endl;
           if(SimVtx_processType[i]==14){
              std::cout << "CONVERSION PROCESS " << SimVtx_processType[i] << std::endl;
           }
           else{
              std::cout << "SOMETHING-ELSE " << SimVtx_processType[i] << std::endl;
           }

//       }
       for (unsigned int j=0; j<nSimTrk; j++){

// Search first for the parent
           if(SimTrk_trackId[j]==SimVtx_simtrk_parent_tid[i])
           std::cout << "Vertex SimTrk PARENT " << j << " " << SimTrk_charge[j] << " " << SimTrk_eta[j] << " " << SimTrk_phi[j] 
                                       << " " << SimTrk_pt[j] << " " << SimTrk_pt[j]*cosh(SimTrk_eta[j]) << " " << SimTrk_pdgId[j] << " " << SimTrk_simvtx_Idx[j] << " " << SimTrk_trackId[j] << std::endl;
//           if(abs(SimTrk_pdgId[j])==11&&SimTrk_simvtx_Idx[j]==i)
           if(SimTrk_simvtx_Idx[j]==i)
           std::cout << "Vertex SimTrk Daughter " << j << " " << SimTrk_charge[j] << " " << SimTrk_eta[j] << " " << SimTrk_phi[j] 
                                       << " " << SimTrk_pt[j] << " " << SimTrk_pt[j]*cosh(SimTrk_eta[j]) << " " << SimTrk_pdgId[j] << " " << SimTrk_simvtx_Idx[j] << " " << SimTrk_trackId[j] << std::endl;
       }
    }

// Loop over all SimVtx
//    if(vsel.size()>=1){
    if(lpr2){
       std::cout << " " << std::endl;
       std::cout << "All SimVtx elements" << std::endl;
    }
    for (unsigned int i=0; i<nSimVtx; i++){
         if(lpr2){
         std::cout << "SimVtx " << i << " " << SimVtx_x[i] << " " << SimVtx_y[i] << " " << SimVtx_z[i] << " " 
                   << SimVtx_tof[i] << " " << SimVtx_processType[i] << " " << SimVtx_simtrk_parent_tid[i] << std::endl;
         if(SimVtx_processType[i]==14){
              std::cout << "Process: CONVERSION " << SimVtx_processType[i] << std::endl;
         }
         else if(SimVtx_processType[i]==0){
              std::cout << "Process: PRIMARY VERTEX " << SimVtx_processType[i] << std::endl;
         }
         else if(SimVtx_processType[i]==201){
              std::cout << "Process: DECAY " << SimVtx_processType[i] << std::endl;
         }
         else if(SimVtx_processType[i]==3){
              std::cout << "Process: BREMSSTRAHLUNG " << SimVtx_processType[i] << std::endl;
         }
         else if(SimVtx_processType[i]==121){
              std::cout << "Process: HADRONIC INTERACTION " << SimVtx_processType[i] << std::endl;
         }
         else{
              std::cout << "Process: SOMETHING-ELSE " << SimVtx_processType[i] << std::endl;
         }
         }
         FillTH1(id_SimVtxProcess, SimVtx_processType[i]);
         if(lpr2){
         for (unsigned int j=0; j<nSimTrk; j++){
// Search first for the parent SimTrk
             if(SimTrk_trackId[j]==SimVtx_simtrk_parent_tid[i])
             std::cout << "Vertex SimTrk PARENT " << j << " " << SimTrk_charge[j] << " " << SimTrk_eta[j] << " " << SimTrk_phi[j] 
                                        << " " << SimTrk_pt[j] << " " << SimTrk_pt[j]*cosh(SimTrk_eta[j]) << " " << SimTrk_pdgId[j] << " " << SimTrk_simvtx_Idx[j] << " " << SimTrk_trackId[j] << std::endl;
             if(SimTrk_simvtx_Idx[j]==i)
             std::cout << "Vertex SimTrk Daughter " << j << " " << SimTrk_charge[j] << " " << SimTrk_eta[j] << " " << SimTrk_phi[j] 
                                        << " " << SimTrk_pt[j] << " " << SimTrk_pt[j]*cosh(SimTrk_eta[j]) << " " << SimTrk_pdgId[j] << " " << SimTrk_simvtx_Idx[j] << " " << SimTrk_trackId[j] << std::endl;
         }
         }
    }

    int nSimPhotons = 0;
    for (unsigned int j=0; j<nSimTrk; j++){
       if(SimTrk_pdgId[j]==22){
          if(lpr2){
             std::cout << "SimTrk Photon " << j << " " << SimTrk_charge[j] << " " << SimTrk_eta[j] << " " << SimTrk_phi[j] 
                                           << " " << SimTrk_pt[j] << " " << SimTrk_pt[j]*cosh(SimTrk_eta[j]) 
                                           << " " << SimTrk_pdgId[j] << " " << SimTrk_simvtx_Idx[j] << " " << SimTrk_trackId[j] << std::endl;
             std::cout << "SimTrkPhotonProcess " << SimVtx_processType[SimTrk_simvtx_Idx[j]] << std::endl;
          }
          FillTH1(id_SimTrkPhotonProcess, SimVtx_processType[SimTrk_simvtx_Idx[j]]);
          nSimPhotons++;
       }
    }
    std::cout << "nSimPhotons " << nSimPhotons << std::endl;

//    }
    
	//gamma gamma stuff
 	if (vsel.size()>=2){
        for(unsigned int ii=0; ii<vsel.size()-1; ++ii){
            int i = vsel[ii];
            double pxi = PC_Px[i];
            double pyi = PC_Py[i];
            double pzi = PC_Pz[i];
            double Ei = PC_E[i];
            double xi = PC_x[i];
            double yi = PC_y[i];
            double Ri = sqrt(xi*xi + yi*yi);
            for(unsigned int jj=ii+1; jj<vsel.size(); ++jj){
                int j = vsel[jj];
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
                    FillTH2(id_mgg2Hist, m12, numberOfPV, wtPU);
                    FillTH1(id_mggCutHist, m12, wtPU);
                }
                if (std::min(Ei,Ej)>2.0){
                    FillTH1(id_mggallHist, m12, wtPU);
                    FillTH2(id_mggRCutHist, m12, Ri, wtPU);
                    FillTH2(id_mggRCutHist, m12, Rj, wtPU);
                }
            }
        }
    }	
/*  for(int i=0; i<PC_vTrack_pt.GetSize(); i++){
        for(int j=0; j<PC_vTrack_pt[i].size(); j++){
			px = PC_vTrack_pt[i][j] * cos( PC_vTrack_phi[i][j] );
			py = PC_vTrack_pt[i][j] * sin( PC_vTrack_phi[i][j] );
			pz = PC_vTrack_pt[i][j] * sinh( PC_vTrack_eta[i][j] );
            FillTH1(id_ptHist, PC_vTrack_pt[i][j], wtPU/PC_vTrack_pt[i][j]);
			FillTH1(id_pzHist, pz, wtPU);
			FillTH2(id_pxpyHist, px, py, wtPU);
        }
    } */
    vcuts.clear(); 

// Loop over all SimTrks and fill photon information 
/*    TH1Manager.at(id_SimPhotonPt) = new MyTH1D("SimPhotonPt", "Sim Photons; Photon pT; Entries per bin", 1000, 0.0, 10.0 );
    TH1Manager.at(id_SimPhotonEta) = new MyTH1D("SimPhotonEta", "Sim Photons; Photon eta; Entries per bin", 200, -10.0, 10.0 );
    TH1Manager.at(id_SimPhotonPhi) = new MyTH1D("SimPhotonPhi", "Sim Photons; Photon phi; Entries per bin", 200, -PI, PI );
    TH1Manager.at(id_SimPhotons) = new MyTH1D("SimPhotons", "Sim Photons; Number of photons; Entries per bin", 1000, -0.5, 999.5 );*/

    vector<int> vSimPhotons;
    int nSimPhoton = 0;
    std::cout << "vSimPhotons ";
    for (unsigned int i=0; i<nSimTrk; i++){
       vSimPhotons.push_back(0);
       if(SimTrk_pdgId[i]==22){
          vSimPhotons[i] += 1;    // Flag SimTrk as a SimPhoton
          nSimPhoton++;
          unsigned int ivtx = SimTrk_simvtx_Idx[i];
// There sometimes seem to be non-statistical fluctuations in the 
// pT histogram. This look like a numerical precision issue. 
// Also maybe some cuts on the photon origin in space ?
          double xv = SimVtx_x[ivtx]; double yv = SimVtx_y[ivtx]; double zv = SimVtx_z[ivtx];
          double distv = sqrt(xv*xv +yv*yv + zv*zv);
          double distxy = sqrt(xv*xv + yv*yv);
          if(abs(SimTrk_eta[i]) < 1.2 && SimTrk_pt[i] > 1.000 && abs(zv) < 12.5 && distxy < 1.0){
// For each SimPhoton - we need to keep track of whether it passes these acceptance cuts
             vSimPhotons[i] += 2;  // Flag SimPhoton as accepted
             FillTH1(id_SimPhotonPt, SimTrk_pt[i]);
             FillTH1(id_SimPhotonEta, SimTrk_eta[i]);
             FillTH1(id_SimPhotonPhi, SimTrk_phi[i]);
             FillTH1(id_SimPhotonDistance, distv);
             FillTH1(id_SimPhotonDistance2, distv);
             FillTH1(id_SimPhotonRadialDistance, distxy);
             FillTH1(id_SimPhotonRadialDistance2, distxy);
             FillTH1(id_SimPhotonzVertex, zv);
          }
       }
       if(vSimPhotons[i] > 0)std::cout << " " << vSimPhotons[i];
    }
    std::cout << endl;
    FillTH1(id_SimPhotons, nSimPhoton);
    
    std::map<unsigned int, std::tuple<unsigned int, unsigned int, unsigned int, double, double, double, double>> mapIndices;
    std::map<unsigned int, std::tuple<unsigned int, unsigned int, unsigned int, double, double, double, double>>::iterator mitr;
    std::cout << "vCnvPhotons ";
    vector<int> vConvertedPhotons(nSimTrk);
    for (unsigned int i=0; i<nSimTrk; i++){
       vConvertedPhotons.push_back(0);
       if(SimTrk_pdgId[i]==22){
// Loop over all SimVtx's trying to find one that has this SimPhoton 
// as the parent
          for (unsigned int j=0; j<nSimVtx; j++){
             if(SimVtx_simtrk_parent_tid[j]==SimTrk_trackId[i]){
                vConvertedPhotons[i] += 1;
                if(SimVtx_processType[j] == 14){
                   vConvertedPhotons[i] +=2;
// Now that we have identified it as a conversion, find the 
// resulting daughter tracks so that we can make kinematic cuts on the 
// energy partitioning
                   vector<unsigned int> vDaughters;
                   int chargesum = 0;
                   int nelectrons = 0;
                   for (unsigned int k=0; k<nSimTrk; k++){
                      if(SimTrk_simvtx_Idx[k]==j){
                         chargesum += int(SimTrk_charge[k]);
                         if(abs(SimTrk_pdgId[k]) == 11)nelectrons +=1;
                         vDaughters.push_back(k);
                      }
                   }
// Now some cuts
                   double pxPos, pyPos, pzPos, pPos, pxNeg, pyNeg, pzNeg, pNeg;
                   if(vDaughters.size()==2 && nelectrons==2 && chargesum ==0){
                      vConvertedPhotons[i] +=4;
// Also check energy matching and direction
                      double pxParent = double(SimTrk_pt[i])*cos(double(SimTrk_phi[i]));
                      double pyParent = double(SimTrk_pt[i])*sin(double(SimTrk_phi[i]));
                      double pzParent = double(SimTrk_pt[i])*sinh(double(SimTrk_eta[i]));
                      double pParent = sqrt(pxParent*pxParent + pyParent*pyParent + pzParent*pzParent);
                      unsigned int k0 = vDaughters[0];
                      unsigned int k1 = vDaughters[1];
                      unsigned int kNeg, kPos;
                      if(SimTrk_pdgId[k0] == 11){  // k0 is the electron
                         pxNeg = double(SimTrk_pt[k0])*cos(double(SimTrk_phi[k0]));
                         pyNeg = double(SimTrk_pt[k0])*sin(double(SimTrk_phi[k0]));
                         pzNeg = double(SimTrk_pt[k0])*sinh(double(SimTrk_eta[k0]));
                         pxPos = double(SimTrk_pt[k1])*cos(double(SimTrk_phi[k1]));
                         pyPos = double(SimTrk_pt[k1])*sin(double(SimTrk_phi[k1]));
                         pzPos = double(SimTrk_pt[k1])*sinh(double(SimTrk_eta[k1]));
                         kNeg = k0;
                         kPos = k1;
                      }
                      else{                        // k1 is the electron
                         pxNeg = double(SimTrk_pt[k1])*cos(double(SimTrk_phi[k1]));
                         pyNeg = double(SimTrk_pt[k1])*sin(double(SimTrk_phi[k1]));
                         pzNeg = double(SimTrk_pt[k1])*sinh(double(SimTrk_eta[k1]));
                         pxPos = double(SimTrk_pt[k0])*cos(double(SimTrk_phi[k0]));
                         pyPos = double(SimTrk_pt[k0])*sin(double(SimTrk_phi[k0]));
                         pzPos = double(SimTrk_pt[k0])*sinh(double(SimTrk_eta[k0]));
                         kNeg = k1;
                         kPos = k0;
                      }
                      pNeg = sqrt(pxNeg*pxNeg + pyNeg*pyNeg + pzNeg*pzNeg);
                      pPos = sqrt(pxPos*pxPos + pyPos*pyPos + pzPos*pzPos);
                      double pxtot = pxNeg + pxPos;
                      double pytot = pyNeg + pyPos;
                      double pztot = pzNeg + pzPos;
                      double ptot = sqrt(pxtot*pxtot + pytot*pytot + pztot*pztot);
//                      std::cout << "Parent    " << pxParent << " " << pyParent << " " << pzParent << " " << pParent << std::endl;
//                      std::cout << "Daughters " << pNeg << " " << pPos << " " << pxtot << " " << pytot << " " << pztot << " " << ptot << std::endl;
// should histogram ptot/pParent and cosangle
//                      std::cout << "Momentum Ratio " << ptot/pParent << std::endl;
                      double cosangle = (pxtot*pxParent + pytot*pyParent + pztot*pzParent)/(ptot*pParent);   // pointing between parent photon and daughter system
//                      std::cout << "cosangle " << cosangle << std::endl;
// Also calculate invariant mass of di-electron pair
                      double melsq = MASS_ELECTRON*MASS_ELECTRON;
// Let's assume very relativistic to avoid precision issues?
                      double eNeg = sqrt(pNeg*pNeg + melsq);
                      double ePos = sqrt(pPos*pPos + melsq);
                      double cosNegPos = (pxNeg*pxPos + pyNeg*pyPos + pzNeg*pzPos)/(pNeg*pPos);
                      double msq = 2.0*melsq + 2.0*eNeg*ePos - 2.0*pNeg*pPos*cosNegPos;
//                      double msq = 2.0*melsq + std::max(0.0,2.0*pNeg*pPos*(1.0-cosNegPos));    // Add max function to clean up numerical issues
//                      std::cout << "cosNegPos " << cosNegPos << std::endl;
//                      std::cout << "pair mass  " << sqrt(msq) << std::endl;
// Add some histograms
                      double Eratio = (eNeg+ePos)/pParent;
//                      std::cout << "Eratio " << Eratio << std::endl;
//                      std::cout << "xPlus " << ePos/(ePos+eNeg) << std::endl;
                      double xPlusLocal = ePos/(ePos+eNeg);
                      FillTH1(id_SimConvPhotonERatio, Eratio);
                      if(Eratio>0.98)FillTH1(id_SimConvPhotonERatioCut, Eratio);
                      FillTH1(id_SimConvPhotonCosAngle, cosangle);
                      FillTH1(id_SimConvPhotonMass, sqrt(msq));
                      FillTH1(id_SimConvPhotonxPlus, ePos/(ePos+eNeg));
                      if(Eratio>0.98)vConvertedPhotons[i] +=8;
                      double xv = SimVtx_x[j]; double yv = SimVtx_y[j]; double zv = SimVtx_z[j];
                      double distv = sqrt(xv*xv +yv*yv + zv*zv);
                      double distxy = sqrt(xv*xv + yv*yv);
                      if(vConvertedPhotons[i] == 15)FillTH1(id_SimConvPhotonR, distxy);
                      if(vConvertedPhotons[i] == 15 && vSimPhotons[i]==3)FillTH1(id_SimConvPhotonAccR, distxy);
// Add this instance to the map
                      if(vConvertedPhotons[i] == 15){
// map with [key = SimPhotonTrk i, value = tuple(SimVtx j, SimTrk kNeg, SimTrk kPos, x, y, z, xPlus)]
                         mapIndices.emplace(i, std::make_tuple(j, kNeg, kPos, double(SimVtx_x[j]), double(SimVtx_y[j]), double(SimVtx_z[j]), xPlusLocal));  // should add some check on this being unique
                      }
                   }
                }
             }
          }
       }
       if(vSimPhotons[i] > 0)std::cout << " " << vConvertedPhotons[i];
    }
    std::cout << endl;
    std::cout <<"Accepted SimPhotons : " << std::endl;
    for (unsigned int i=0; i<nSimTrk; i++){
        FillTH1(id_SimPhotonId, vSimPhotons[i]);
        FillTH1(id_SimConvPhotonId, vConvertedPhotons[i]);
        if(vSimPhotons[i]==3&&vConvertedPhotons[i]==15){
// Extract map info with SimConvPhoton info
           mitr = mapIndices.find(i);
           if(mitr != mapIndices.end()){
              auto mytuple = mitr->second;
              auto x = std::get<3>(mytuple);
              auto y = std::get<4>(mytuple);
              double Rlocal = sqrt(x*x + y*y);
              FillTH1(id_EffSimConvPhotonsDen, Rlocal);
           }
        }
    }
    std::cout << endl;
   // We want to find the best match if any amongst flagged SimConvPhotons and selected reconstructed conversions.
   // Maybe only bother with Accepted SimPhotons to start with.
    for(unsigned int icandidate=0; icandidate<vsel.size(); ++icandidate){
// First loop over reco conversion candidates
        int iconv = vsel[icandidate];
        double x = PC_x[iconv]; double y = PC_y[iconv]; double z = PC_z[iconv];
        std::cout << "Reco candidate " << icandidate << " " << iconv << " " << x << " " << y << " " << z << " " << PC_E[iconv] << std::endl;
        double dmax = 2.0;  
        double Rmatch;
        int imatch = -1;
// Loop over SimTrks trying to match to SimConvPhotons passing the acceptance cuts above
        for(unsigned int i=0; i<nSimTrk; ++i){
            if(vSimPhotons[i]==3 && vConvertedPhotons[i]==15){
// Extract map info with SimConvPhoton info
               mitr = mapIndices.find(i);
               if(mitr != mapIndices.end()){
                  auto mytuple = mitr->second;
                  auto xi = std::get<3>(mytuple);
                  auto yi = std::get<4>(mytuple);
                  auto zi = std::get<5>(mytuple);
                  double dist = sqrt((x-xi)*(x-xi) + (y-yi)*(y-yi) + (z-zi)*(z-zi));
                  std::cout << "SimConvPhoton " << i << " " << xi << " " << yi << " " << zi << " " << " dist (cm) = " << dist << std::endl;
                  if(dist < dmax){
                     dmax = dist;
                     imatch = i;
                     Rmatch = sqrt(xi*xi + yi*yi);
                  }
               }
            }
        }
        if(imatch >= 0){
           std::cout << "imatch set to " << imatch << std::endl;
           FillTH1(id_EffSimConvPhotonsNum, Rmatch);
// TODO histogram match distance and other things ...
        }
        else{
           std::cout << "no match found " << std::endl;
        }
    }
}  // End of event number selection if block
}  // End of sub-program
#endif
