#define convsel_cxx
// The class definition in convsel.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("convsel.C")
// root> T->Process("convsel.C","some options")
// root> T->Process("convsel.C+")
//


#include "convsel.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include "TMath.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
using namespace ROOT::Math;

void convsel::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void convsel::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t convsel::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
//typedef ROOT::math::PxPyPzMVector LorentzVector;

   fReader.SetLocalEntry(entry);

  for(int i=0; i< *nConv; i++){
    xyhist_25.Fill(Conv_vtx_X[i],Conv_vtx_Y[i]);
    xyhist_10.Fill(Conv_vtx_X[i],Conv_vtx_Y[i]);
  }

double pt1,eta1,phi1,pt2,eta2,phi2;
double px1,py1,pz1;
double px2,py2,pz2;
int q1,q2;
double fitprob;
double ip1,ip2;
double chisq1,chisq2;
int ndof1,ndof2;

const double RERRCUT = 0.25;
const double COSTCUT = 0.85;
const double ZCUT = 25.0;
const double FITPROBCUT = 0.01;
const double MASS_ELECTRON = 0.5109989461e-3;
const double MASS_PION = 139.57061e-3;
const double MASS_KAON = 493.677e-3;
const double MASS_PROTON = 938.272081e-3;

double theta;
double x,y,z,r;

double x0=0.0;
double y0=0.0;
double vxx,vxy,vyy;
double phi,cphi,sphi,varsum_r,rerr;
double px,py,pz,pt;
double ptasym,xplus;
double alp1,alp2;
double tanl1,tanl2;
double f1,f2;
double A1,A2,B1,B2;
double qR1,qR2;
double x1,y1,z1;
double x2,y2,z2;
double px1p,py1p,pz1p;
double px2p,py2p,pz2p;
double x1p,y1p,z1p;
double x2p,y2p,z2p;
double dsqxy,dsqxyz;
double xc,yc,zc,tphi,dc;
double lamx,lamy;
double pgamma;

int goodConvs = 0;

std::cout << "Event " << *event << "   nPV " << *nPV << std::endl;

std::vector<double> tkvector;

double PVX,PVY,PVZ;
// let's just use the first one
PVX = PV_X[0];
PVY = PV_Y[0];
PVZ = PV_Z[0];

for(int i=0; i<*nPV; i++){
/*   std::cout << " PV " << i << " " << PV_X[i] << " " << PV_Y[i] << " " << PV_Z[i] << std::endl;
   std::cout << " PV chi-squared/ndof " << PV_chi2[i] << " " << PV_ndof[i] << " " << PV_normalizedChi2[i] << std::endl;
   std::cout << " PV sigma " << sqrt(PV_cov_00[i]) << " " << sqrt(PV_cov_11[i]) << " " << sqrt(PV_cov_22[i]) << std::endl;
   std::cout << " PV cov   " << PV_cov_01[i] << " " << PV_cov_02[i] << " " << PV_cov_12[i] << std::endl;
   std::cout << " " << std::endl;
*/
}

//Begin scan with DR matching for conv tracks
for(int i=0; i<*nConv; i++){
    x = Conv_vtx_X[i];
    y = Conv_vtx_Y[i];
    z = Conv_vtx_Z[i];
    r = sqrt(x*x+y*y);

    ip1 = Conv_tracksSigned_d0_Tk0[i];
    chisq1 = Conv_Tk0_chi2[i];
    ndof1 = Conv_Tk0_ndof[i];
    pt1 = Conv_Tk0_pt[i];
    eta1 = Conv_Tk0_eta[i];
    phi1 = Conv_Tk0_phi[i];  // THIS LOOKS BUGGY ?
    q1 = Conv_Tk0_charge[i];
    px1 = Conv_tracksPin_Px_Tk0[i];
    py1 = Conv_tracksPin_Py_Tk0[i];
    pz1 = Conv_tracksPin_Pz_Tk0[i];
    x1 = Conv_tracksInnerPosition_X_Tk0[i];
    y1 = Conv_tracksInnerPosition_Y_Tk0[i];
    z1 = Conv_tracksInnerPosition_Z_Tk0[i];
    tanl1 = pz1/sqrt(px1*px1 + py1*py1);
    f1 = atan2(py1,px1);
    qR1 = 100.0*double(q1)*pt1/(0.2998*3.80);   // in cm
    A1 = 2.0*qR1*( (y1-y)*cos(f1) - (x1-x)*sin(f1) - qR1 );
    B1 =-2.0*qR1*( (y1-y)*sin(f1) + (x1-x)*cos(f1) );
    alp1 = atan2(B1/A1,1.0);
    px1p = pt1*cos(f1 + alp1);
    py1p = pt1*sin(f1 + alp1);
    pz1p = pz1;
// Also swim the track position
    x1p = x1 + qR1*( (1.0-cos(alp1))*sin(f1) - sin(alp1)*cos(f1) );
    y1p = y1 - qR1*( (1.0-cos(alp1))*cos(f1) + sin(alp1)*sin(f1) );
    z1p = z1 - qR1*tanl1*alp1;
    XYZVector xyzv1 = XYZVector(px1p,py1p,pz1p);
    XYZVector xyzvt1 = XYZVector(px1p,py1p,0.0);

    ip2 = Conv_tracksSigned_d0_Tk1[i];
    chisq2 = Conv_Tk1_chi2[i];
    ndof2 = Conv_Tk1_ndof[i];
    pt2= Conv_Tk1_pt[i];
    eta2 = Conv_Tk1_eta[i];
    phi2 = Conv_Tk1_phi[i];  // THIS LOOKS BUGGY ?  Maybe it is phi0 ?
    q2 = Conv_Tk1_charge[i];
    px2 = Conv_tracksPin_Px_Tk1[i];
    py2 = Conv_tracksPin_Py_Tk1[i];
    pz2 = Conv_tracksPin_Pz_Tk1[i];
    x2 = Conv_tracksInnerPosition_X_Tk1[i];
    y2 = Conv_tracksInnerPosition_Y_Tk1[i];
    z2 = Conv_tracksInnerPosition_Z_Tk1[i];
    tanl2 = pz2/sqrt(px2*px2 + py2*py2);
    f2 = atan2(py2,px2);

    qR2 = 100.0*double(q2)*pt2/(0.2998*3.80);   // in cm
    A2 = 2.0*qR2*( (y2-y)*cos(f2) - (x2-x)*sin(f2) - qR2 );
    B2 =-2.0*qR2*( (y2-y)*sin(f2) + (x2-x)*cos(f2) );
    alp2 = atan2(B2/A2,1.0);
    px2p = pt2*cos(f2 + alp2);
    py2p = pt2*sin(f2 + alp2);
    pz2p = pz2;
// Also swim the track position
    x2p = x2 + qR2*( (1.0-cos(alp2))*sin(f2) - sin(alp2)*cos(f2) );
    y2p = y2 - qR2*( (1.0-cos(alp2))*cos(f2) + sin(alp2)*sin(f2) );
    z2p = z2 - qR2*tanl2*alp2;
    XYZVector xyzv2 = XYZVector(px2p,py2p,pz2p);
    XYZVector xyzvt2 = XYZVector(px2p,py2p,0.0);

    dsqxy = (x1p-x2p)*(x1p-x2p) + (y1p-y2p)*(y1p-y2p);
    dsqxyz = dsqxy + (z1p-z2p)*(z1p-z2p);

//Normal analysis cuts
    fitprob = TMath::Prob(Conv_vtx_chi2[i], 3);

	vxx = Conv_vtx_cov_00[i];
	vxy = Conv_vtx_cov_01[i];
	vyy = Conv_vtx_cov_11[i];
    phi = atan2(y-y0, x-x0);
	cphi = cos(phi);
	sphi = sin(phi);
	varsum_r   = cphi*cphi*vxx + 2.0*sphi*cphi*vxy + sphi*sphi*vyy;
	rerr = sqrt(varsum_r);
    px = Conv_refittedPair4Momentum_Px[i];
    py = Conv_refittedPair4Momentum_Py[i];
    pz = Conv_refittedPair4Momentum_Pz[i];
    pgamma = sqrt(px*px+py*py+pz*pz);

    XYZVector xyzv = XYZVector(px,py,pz);

    pt = sqrt(px*px + py*py);
    theta = atan2(pt,pz);
    tphi = py/px;
// Try to find the dca of the photon to the PV in the r-phi plane.
 
    xc = ((PVX + x*tphi*tphi) - (y - PVY)*tphi)/(1.0 + tphi*tphi);
    yc = y + (xc-x)*tphi;
    dc = sqrt( (PVX-xc)*(PVX-xc) + (PVY-yc)*(PVY-yc) );
// Equation of a line
    lamx = (xc - x)/(px/pgamma);
    lamy = (yc - y)/(py/pgamma);
    zc = z + 0.5*(lamx+lamy)*pz/pgamma;

    ptasym = (pt1-pt2)/(pt1+pt2);
    if (q1<0) ptasym = -ptasym;
    xplus = 0.5*(1.0 + ptasym);

	if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT 
        && fitprob > FITPROBCUT && std::max(Conv_nHitsBeforeVtx_Tk0[i],Conv_nHitsBeforeVtx_Tk1[i])==0 ){

    asymHist.Fill(ptasym);
    xasymHist.Fill(xplus);

//	if( rerr < RERRCUT && abs(z) < ZCUT && abs(cos(theta)) < COSTCUT 
//        && fitprob > FITPROBCUT ){

        goodConvs++;
    std::cout << " " << std::endl;
	std::cout << "Conversion " << i << " pT " << pt << " algo " << Conv_algo[i] << " A(pT): " << ptasym << " x+: " << xplus << std::endl;
    std::cout << " xyzv : " << xyzv.X() << " " << xyzv.Y() << " " << xyzv.Z() << std::endl;
    std::cout << " xyzv1: " << xyzv1.X() << " " << xyzv1.Y() << " " << xyzv1.Z() << std::endl;
    std::cout << " xyzv2: " << xyzv2.X() << " " << xyzv2.Y() << " " << xyzv2.Z() << std::endl;
    std::cout << " Sum  : " << xyzv1.X()+xyzv2.X() << " " << xyzv1.Y()+xyzv2.Y() << " " << xyzv1.Z()+xyzv2.Z() << std::endl;

    XYZVector vcross1 = xyzv.Cross(xyzv1);
    XYZVector vcross2 = xyzv.Cross(xyzv2);
    double appt1 = sqrt(vcross1.Mag2())/sqrt(xyzv.Mag2());
    double appt2 = sqrt(vcross2.Mag2())/sqrt(xyzv.Mag2());
    double pL1 = xyzv1.Dot(xyzv)/sqrt(xyzv.Mag2());
    double pL2 = xyzv2.Dot(xyzv)/sqrt(xyzv.Mag2());
    double asympL = (pL1 - pL2)/(pL1 + pL2);
    if (q1<0) asympL = -asympL;

// It may be better to choose the track with the smaller AP pT for estimating AP pT??
// Or the one with the higher pT? 
// - or at least the track whose innermost hit is closest to the vertex
// (or in principle the input tracks should have been propagated to the vertex ...)
// For the pL analysis one could apply a similar approach in terms of the best measured one.
    std::cout << "AP analysis " << appt1 << " " << appt2 << " " << pL1 << " " << pL2 << " " << asympL << std::endl;
    minpTHist.Fill(std::min(appt1,appt2));
    minpTHist2.Fill(std::min(appt1,appt2));
    maxpTHist.Fill(std::max(appt1,appt2));

    std::cout << "Cut quantities: rerr, z, cos(theta), fitprob: " << rerr << " " << z << " " << cos(theta) << " " << fitprob << std::endl;
	std::cout << "c vtx (x,y,z) : " << Conv_vtx_X[i] << " " << Conv_vtx_Y[i] << " " << Conv_vtx_Z[i] 
                                      << " R: " << sqrt(Conv_vtx_X[i]*Conv_vtx_X[i] + Conv_vtx_Y[i]*Conv_vtx_Y[i]) 
                                      << " phi: " << atan2(Conv_vtx_Y[i], Conv_vtx_X[i]) << std::endl;
    std::cout << "c vtx z error : " << sqrt(Conv_vtx_cov_22[i]) << std::endl;
    std::cout << "First PV        " << PVX << " " << PVY << " " << PVZ << std::endl;
    std::cout << "2d photon pca : " << xc << " " << yc << " " << zc << " " << tphi << " " << dc << std::endl;
    std::cout << "lamx,lamy,z,PVZ " << lamx << " " << lamy << " " << zc << " " << PVZ << std::endl;
    dcaHist.Fill(dc);

    std::cout << "Conversion quality : " << Conv_nHitsBeforeVtx_Tk0[i] << " " << Conv_nHitsBeforeVtx_Tk1[i] << " " << Conv_nSharedHits[i] << " " << Conv_vtx_chi2[i] << " " << fitprob << std::endl;
    std::cout << "Pair 3-vector      : " << Conv_pairMomentum_Px[i] << " " << Conv_pairMomentum_Py[i] << " " << Conv_pairMomentum_Pz[i] 
              << " phi: " << atan2(Conv_pairMomentum_Py[i], Conv_pairMomentum_Px[i]) << std::endl;
    std::cout << "Refit 5-vector     : " << Conv_refittedPair4Momentum_E[i] << " " << Conv_refittedPair4Momentum_M[i] << " " 
              << Conv_refittedPair4Momentum_Px[i] << " " << Conv_refittedPair4Momentum_Py[i] << " " << Conv_refittedPair4Momentum_Pz[i] << 
              " phi: " <<  atan2(Conv_refittedPair4Momentum_Py[i], Conv_refittedPair4Momentum_Px[i]) <<  std::endl;
    fittedPairMassHist.Fill(Conv_refittedPair4Momentum_M[i]);
    if(Conv_refittedPair4Momentum_M[i]<0.0)std::cout << "SCREAM - negative mass " << std::endl;
    double hash1= double(q1)*chisq1/double(ndof1);
    double hash2= double(q2)*chisq2/double(ndof2);
    minDofHist.Fill(std::min(ndof1,ndof2));
    maxDofHist.Fill(std::max(ndof1,ndof2));
    maxNormalizedChi2Hist.Fill(std::max(abs(hash1),abs(hash2)));
    std::cout << "Tk0 (pt,eta,phi0,ip,chisq,ndof,charge) " << pt1 << " " << eta1 << " " << phi1 << " " << ip1 << " " << chisq1 << " " << ndof1 << " " << q1 << " " << hash1 << std::endl;
	std::cout << "Tk1 (pt,eta,phi0,ip,chisq,ndof,charge) " << pt2 << " " << eta2 << " " << phi2 << " " << ip2 << " " << chisq2 << " " << ndof2 << " " << q2 << " " << hash2 << std::endl;
    std::cout << "Pair quantities (dMin,dPhi,dCotTh,m,zOfPV,dlsig0,dlsig1) " 
              << Conv_distOfMinimumApproach[i] << " "
              << Conv_dPhiTracksAtVtx[i] << " " 
              << Conv_pairCotThetaSeparation[i] << " " << Conv_pairInvariantMass[i] << " " 
              << Conv_zOfPrimaryVertexFromTracks[i] << " "
              << Conv_dlClosestHitToVtx_sig_Tk0[i] << " " << Conv_dlClosestHitToVtx_sig_Tk1[i] << std::endl;
    std::cout << "Tk0 pinner (px,py,pz|pt,phi) " << Conv_tracksPin_Px_Tk0[i] << " " << Conv_tracksPin_Py_Tk0[i] << " " << Conv_tracksPin_Pz_Tk0[i] 
              << " | " << sqrt(pow(Conv_tracksPin_Px_Tk0[i],2) + pow(Conv_tracksPin_Py_Tk0[i],2)) << " " << atan2(Conv_tracksPin_Py_Tk0[i], Conv_tracksPin_Px_Tk0[i] ) << std::endl;
    std::cout << "Tk0 pouter (px,py,pz|pt,phi) " << Conv_tracksPout_Px_Tk0[i] << " " << Conv_tracksPout_Py_Tk0[i] << " " << Conv_tracksPout_Pz_Tk0[i] 
              << " | " << sqrt(pow(Conv_tracksPout_Px_Tk0[i],2) + pow(Conv_tracksPout_Py_Tk0[i],2)) << " " << atan2(Conv_tracksPout_Py_Tk0[i], Conv_tracksPout_Px_Tk0[i] ) << std::endl;
    std::cout << "Tk1 pinner (px,py,pz|pt,phi) " << Conv_tracksPin_Px_Tk1[i] << " " << Conv_tracksPin_Py_Tk1[i] << " " << Conv_tracksPin_Pz_Tk1[i] 
              << " | " << sqrt(pow(Conv_tracksPin_Px_Tk1[i],2) + pow(Conv_tracksPin_Py_Tk1[i],2)) << " " << atan2(Conv_tracksPin_Py_Tk1[i], Conv_tracksPin_Px_Tk1[i] ) << std::endl;
    std::cout << "Tk1 pouter (px,py,pz|pt,phi) " << Conv_tracksPout_Px_Tk1[i] << " " << Conv_tracksPout_Py_Tk1[i] << " " << Conv_tracksPout_Pz_Tk1[i] 
              << " | " << sqrt(pow(Conv_tracksPout_Px_Tk1[i],2) + pow(Conv_tracksPout_Py_Tk1[i],2)) << " " << atan2(Conv_tracksPout_Py_Tk1[i], Conv_tracksPout_Px_Tk1[i] ) << std::endl;
    std::cout << "Pair pinner (px,py,pz) " << Conv_tracksPin_Px_Tk0[i] + Conv_tracksPin_Px_Tk1[i] << " "  << Conv_tracksPin_Py_Tk0[i] + Conv_tracksPin_Py_Tk1[i] << " " 
                                           << Conv_tracksPin_Pz_Tk0[i] + Conv_tracksPin_Pz_Tk1[i] << std::endl;
    std::cout << "Tk0 alpha, pvertex, rpca " << alp1 << " " << px1p << " " << py1p << " " << pz1p << " " << x1p << " " << y1p << " " << z1p << std::endl;
    std::cout << "Tk1 alpha, pvertex, rpca " << alp2 << " " << px2p << " " << py2p << " " << pz2p << " " << x2p << " " << y2p << " " << z2p << std::endl;
    std::cout << "Miss distances at rpca " << sqrt(dsqxy) << " " << z1p-z2p << " " << sqrt(dsqxyz) << " " << std::endl;
    std::cout << "Angles at rpca (3-d): " << acos(xyzv1.Dot(xyzv2)/sqrt(xyzv1.Mag2()*xyzv2.Mag2())) << std::endl;
    std::cout << "Angles at rpca (2-d): " << acos(xyzvt1.Dot(xyzvt2)/sqrt(xyzvt1.Mag2()*xyzvt2.Mag2())) << std::endl;
    std::cout << "tanl difference:      " << tanl1 - tanl2 << std::endl;

    PxPyPzMVector v0,v0pi,v0p; 
//    v0   = PxPyPzMVector( Conv_tracksPin_Px_Tk0[i], Conv_tracksPin_Py_Tk0[i], Conv_tracksPin_Pz_Tk0[i], MASS_ELECTRON);
//    v0pi = PxPyPzMVector( Conv_tracksPin_Px_Tk0[i], Conv_tracksPin_Py_Tk0[i], Conv_tracksPin_Pz_Tk0[i], MASS_PION);
//    v0p  = PxPyPzMVector( Conv_tracksPin_Px_Tk0[i], Conv_tracksPin_Py_Tk0[i], Conv_tracksPin_Pz_Tk0[i], MASS_PROTON);
    v0   = PxPyPzMVector( px1p, py1p, pz1p, MASS_ELECTRON );
    v0pi = PxPyPzMVector( px1p, py1p, pz1p, MASS_PION );
    v0p  = PxPyPzMVector( px1p, py1p, pz1p, MASS_PROTON );
    PxPyPzMVector v1,v1pi,v1p; 
//    v1   = PxPyPzMVector( Conv_tracksPin_Px_Tk1[i], Conv_tracksPin_Py_Tk1[i], Conv_tracksPin_Pz_Tk1[i], MASS_ELECTRON);
//    v1pi = PxPyPzMVector( Conv_tracksPin_Px_Tk1[i], Conv_tracksPin_Py_Tk1[i], Conv_tracksPin_Pz_Tk1[i], MASS_PION);
//    v1p  = PxPyPzMVector( Conv_tracksPin_Px_Tk1[i], Conv_tracksPin_Py_Tk1[i], Conv_tracksPin_Pz_Tk1[i], MASS_PROTON);
    v1   = PxPyPzMVector( px2p, py2p, pz2p, MASS_ELECTRON );
    v1pi = PxPyPzMVector( px2p, py2p, pz2p, MASS_PION );
    v1p  = PxPyPzMVector( px2p, py2p, pz2p, MASS_PROTON );

    PxPyPzMVector vpair, vpairpipi, vpairpip, vpairppi;
    vpair += v0;
    vpair += v1;
    vpairpipi += v0pi;
    vpairpipi += v1pi;
    vpairpip  += v0pi;
    vpairpip  += v1p;
    vpairppi  += v0p;
    vpairppi  += v1pi;
   
//    v0.Print();
//    v1.Print();
    std::cout << "v0 eta,phi,theta,cotq: " << v0.Eta() << " " << v0.Phi() << " " << v0.Theta() << " " << 1.0/tan(v0.Theta()) << std::endl;
    std::cout << "v1 eta,phi,theta,cotq: " << v1.Eta() << " " << v1.Phi() << " " << v1.Theta() << " " << 1.0/tan(v1.Theta()) << std::endl;
    std::cout << "vpair (px,py,pz,m,E) : " << vpair.Px() << " " << vpair.Py() << " " << vpair.Pz() << " " << vpair.M() << " " << vpair.E() << std::endl;
    std::cout << "Masses (ee, pipi, pip, ppi): " << vpair.M() << " " << vpairpipi.M() << " " << vpairpip.M() << " " << vpairppi.M() << std::endl;

    XYZVector s0 = v0.Vect();
    XYZVector s1 = v1.Vect();
    double cosangle = s0.Dot(s1)/sqrt(s0.Mag2()*s1.Mag2());
    std::cout << "Dot product, cos(0,1), angle, dcotq, angleE " << s0.Dot(s1) << " " << cosangle << " " << acos(cosangle) << " " 
              << (1.0/tan(v0.Theta())) - (1.0/tan(v1.Theta())) << " " << acos(cosangle)*vpair.E() << std::endl;
    angleEHist.Fill(acos(cosangle)*vpair.E());
//    if(vpair.M()<0.05){
//    if(acos(cosangle)*vpair.E()>0.05){
// Restrict to conversions that either do or do not reconstruct well as photon conversions
       lambdaCandidateMassHist.Fill(vpairpip.M());
       lambdaCandidateMassHist.Fill(vpairppi.M());
//    }
    photonCandidateMassHist.Fill(vpair.M());
    photonCandidateMassHist2.Fill(vpair.M());
    std::cout << "Tk0 inner position " << Conv_tracksInnerPosition_X_Tk0[i] << " " << Conv_tracksInnerPosition_Y_Tk0[i] << " " << Conv_tracksInnerPosition_Z_Tk0[i] << std::endl;
    std::cout << "Tk0 aflq-criteria  " << Conv_Tk0_algo[i] << " " << Conv_Tk0_found[i] << " " << Conv_Tk0_lost[i] << " " << Conv_Tk0_quality[i] << std::endl; 
    std::cout << "Tk1 inner position " << Conv_tracksInnerPosition_X_Tk1[i] << " " << Conv_tracksInnerPosition_Y_Tk1[i] << " " << Conv_tracksInnerPosition_Z_Tk1[i] << std::endl;
    std::cout << "Tk1 aflq-criteria  " << Conv_Tk1_algo[i] << " " << Conv_Tk1_found[i] << " " << Conv_Tk1_lost[i] << " " << Conv_Tk1_quality[i] << std::endl; 
    double deltaX = Conv_tracksInnerPosition_X_Tk1[i] - Conv_tracksInnerPosition_X_Tk0[i];
    double deltaY = Conv_tracksInnerPosition_Y_Tk1[i] - Conv_tracksInnerPosition_Y_Tk0[i];
    double deltaZ = Conv_tracksInnerPosition_Z_Tk1[i] - Conv_tracksInnerPosition_Z_Tk0[i];
    double dL = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    double dR = sqrt(deltaX*deltaX + deltaY*deltaY);
    std::cout << "Inner hits: dL, dR " << dL << " " << dR << std::endl;
    dLHist.Fill(dL);
    dRHist.Fill(dR);

// Add to track vector if hash values not already in the vector
    bool lfound = false;
    for (int j = 0; j < tkvector.size(); j++){
       if(tkvector[j] == hash1) lfound = true;
    }
    if(!lfound)tkvector.push_back(hash1);
    lfound = false;
    for (int j = 0; j < tkvector.size(); j++){
       if(tkvector[j] == hash2) lfound = true;
    }
    if(!lfound)tkvector.push_back(hash2);

    }
}
    std::cout << " " << std::endl;
    std::cout << "Numbers : " << goodConvs << " " << tkvector.size() << std::endl;
    int nduplicates = 2*goodConvs - tkvector.size();
    if(nduplicates > 0)std::cout << "Duplicate conversions: " << nduplicates << std::endl;

    convHist.Fill(goodConvs);
    gconvHist.Fill(goodConvs - nduplicates);
    duplicatesHist.Fill(nduplicates);

   return kTRUE;
}

void convsel::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void convsel::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  TCanvas* c1 = new TCanvas();
  xyhist_25.Draw("COLZ");
  TCanvas* c2 = new TCanvas();
  xyhist_10.Draw("COLZ");
  TCanvas* c4 = new TCanvas();
  convHist.Draw("hist");
  TCanvas* c5 = new TCanvas();
  gconvHist.Draw("hist");
  TCanvas* c6 = new TCanvas();
  duplicatesHist.Draw("hist");
  TCanvas* c7 = new TCanvas();
  dLHist.Draw("hist");
  TCanvas* c8 = new TCanvas();
  dRHist.Draw("hist");
  TCanvas* c9 = new TCanvas();
  fittedPairMassHist.Draw("hist");
  TCanvas* c10 = new TCanvas();
  angleEHist.Draw("hist");
  TCanvas* c11 = new TCanvas();
  lambdaCandidateMassHist.Draw("hist");
  TCanvas* c12 = new TCanvas();
  minDofHist.Draw("hist");
  TCanvas* c13 = new TCanvas();
  maxDofHist.Draw("hist");
  TCanvas* c14 = new TCanvas();
  maxNormalizedChi2Hist.Draw("hist");
  TCanvas* c15 = new TCanvas();
  xasymHist.Draw("hist");
  TCanvas* c16 = new TCanvas();
  asymHist.Draw("hist");
  TCanvas* c17 = new TCanvas();
  minpTHist.Draw("hist");
  TCanvas* c18 = new TCanvas();
  minpTHist2.Draw("hist");
  TCanvas* c19 = new TCanvas();
  maxpTHist.Draw("hist");
  TCanvas* c20 = new TCanvas();
  photonCandidateMassHist.Draw("hist");
  TCanvas* c21 = new TCanvas();
  photonCandidateMassHist2.Draw("hist");
  TCanvas* c22 = new TCanvas();
  dcaHist.Draw("hist");

  
}
