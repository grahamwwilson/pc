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
#include <map>
#include <iomanip>
#include "Hungarian.h"
using MyTH1D = ROOT::TThreadedObject<TH1D>;
using MyTH2D = ROOT::TThreadedObject<TH2D>;

// struct for derived quantities of each conversion
#include  "mystruct.h"

class histset2{
    public:
       double PI =4.0*atan(1.0);
       histset2();	
       void init(); 
       void setweightoption(); 
       void AnalyzeEntry(convsel& s);
       #include "histset2enums.h" 
// make a big vector and load enumerated histograms onto the vector
       std::vector<MyTH1D*>  TH1Manager{};
       std::vector<MyTH2D*>  TH2Manager{};
// locate the histogram and perform pointer copying 
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

#include "histset2init.h"     //Put the histogram declarations in one separate file

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

    const bool lpr = false;      // print flag
    const bool lreduce = true;   // do problem reduction
	
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
    
    double wtPU;

    if(lpr)std::cout << " " << std::endl;
    if(lpr)std::cout << " ------------BoE-------------- " << std::endl;
    if(lpr)std::cout << " " << std::endl;

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

    std::vector<bool> vcuts;
    std::vector<double> PosTkInfo;
    std::vector<double> NegTkInfo;    
    std::vector<int> vcandidate;
    std::map<double, int> mNeg;
    std::multimap<double, int> mmNeg;
    std::map<double, int> mPos;
    std::multimap<double, int> mmPos;
    std::multimap<int, double> mmCandidate;

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
            }
            else{
               PosTkInfo[i] =  Tk1_chi2[i]/double(Tk1_ndof[i]);
               NegTkInfo[i] = -Tk0_chi2[i]/double(Tk0_ndof[i]);
            }
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
//            tup[i].mpair = vpair.M();
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
// for technical checks
            FillTH1( id_r1dHist2, tup[i].radius );
        }
// Also include correlation plot between APalpha and ptasym.
	}

// Characterize our matching problem for this event
// Note this is unnecessary if there are no duplicates, ie. n- = n+ = nedges.
    std::vector<int> vsel;   // Vector for indices of selected conversions after removing duplicates

    if(std::min(mNeg.size(), mPos.size()) < vcandidate.size()){
// We actually have an assignment problem to worry about
       if(lpr){
          std::cout << " " << std::endl;
          std::cout << "Event " << eventNumber << std::endl;
          std::cout << "nedges = " << vcandidate.size() << std::endl;
          std::cout << "N- = " << mmNeg.size() << std::endl;
          std::cout << "N+ = " << mmPos.size() << std::endl;
          std::cout << "n- = " << mNeg.size() << std::endl;
          std::cout << "n+ = " << mPos.size() << std::endl;
          std::cout << "Target maximum possible cardinality of solution = " << std::min(mNeg.size(),mPos.size()) << std::endl;
// Prior solution with no arbitration
          std::cout << "Prior edge solution: ";
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

        if(tup[i].pt>1.0)FillTH1(id_xplus1Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>2.0)FillTH1(id_xplus2Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>4.0)FillTH1(id_xplus4Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>8.0)FillTH1(id_xplus8Hist, tup[i].xplus, wtPU);
        if(tup[i].pt>16.0)FillTH1(id_xplus16Hist, tup[i].xplus, wtPU);

        if(tup[i].pt>1.0)FillTH1(id_alpha1Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>2.0)FillTH1(id_alpha2Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>4.0)FillTH1(id_alpha4Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>8.0)FillTH1(id_alpha8Hist, tup[i].alpha, wtPU);
        if(tup[i].pt>16.0)FillTH1(id_alpha16Hist, tup[i].alpha, wtPU);

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
        FillTH1(id_KShortMassHist, tup[i].mass[1], wtPU);
        if(region1||region2||region3){
           FillTH1(id_KShortBkgdMassHist, tup[i].mass[1], wtPU);
           if(region1)FillTH1(id_KShortBkgdMassHistR1, tup[i].mass[1], wtPU);
           if(region2)FillTH1(id_KShortBkgdMassHistR2, tup[i].mass[1], wtPU);
           if(region3)FillTH1(id_KShortBkgdMassHistR3, tup[i].mass[1], wtPU);
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

    }  // End of selected candidate loop
    
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
}
#endif
