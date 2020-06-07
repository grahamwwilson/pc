void plot3ratioMC(string histtype="r1dHist2", float ymax=3000000.0, 
          float ymin=0.1, float xlmin = 0.33, float ylmin = 0.66, 
          bool logscale = false, string xtitle="Radius (cm)", string ytitle = "Cands (new) / Cands (old) per 0.1cm bin")
{

TCanvas *c1 = new TCanvas("c1","multipads",800,600);
if(logscale)c1->SetLogy(1);

//TFile *fd = new TFile("../PC_DataHPC.root");
//TFile *fm = new TFile("../PC_MCHPC.root");

//TFile *fd = new TFile("../PC_DataHPC_R0.root");
//TFile *fm = new TFile("../PC_MCHPC_R0.root");

TFile *fd = new TFile("../PC_mc.root");
//TFile *fm = new TFile("../PC_mc.root");

TH1D * hd = (TH1D*)fd->Get(histtype.c_str());
TH1D * hm = (TH1D*)fd->Get("r1dwidecutHist");

//TH1D * hdpv = (TH1D*)fd->Get("numpvWHist");
//TH1D * hmpv = (TH1D*)fm->Get("numpvWHist");

hd->Print();
cout << hd->Integral() << endl;
hm->Print();
cout << hm->Integral() << endl;
double myratio=hd->Integral()/hm->Integral();
cout << "myratio: " << myratio << " " << 1.0/myratio << endl;

//double scale = double(hd->GetEntries())/double(hm->GetEntries());
double scale = myratio;
//cout << "Scale = " << scale << endl;

/*
vector<double> weights;
// Find total number of interactions from PVs 
long int npvdtotal = 0;
long int npvmtotal = 0;
for (int kmult=1; kmult<=99;kmult++){
    int i = kmult + 1;
    int nbindcontent = hdpv->GetBinContent(i);
    int nbinmcontent = hmpv->GetBinContent(i);
    cout << "PV multiplicity " << kmult << " Bin number " << i << " " 
         << nbindcontent << " " << nbinmcontent 
         << " Ratio " << double(nbindcontent)/(scale*double(nbinmcontent)) <<endl;
    npvdtotal += kmult*nbindcontent;
    npvmtotal += kmult*nbinmcontent;
    cout << "Running suma " << npvdtotal << " " << npvmtotal << endl;
    weights.push_back( double(nbindcontent)/(scale*double(nbinmcontent)) );
}
double ratio = double(npvdtotal)/double(npvmtotal);
cout << "Ratio = " << ratio << endl;
*/

hd->GetYaxis()->SetTitleOffset(1.4);
hd->SetMaximum(ymax);
hd->SetMinimum(ymin);
hd->SetStats(kFALSE);
hm->SetStats(kFALSE);

hd->Divide(hm,hd);
hd->SetMaximum(1.05);
hd->SetMinimum(0.75);

c1->SetTicks(1,1);

//scale=ratio;  // Overwrite with PV based one.

//hm->Scale(scale);

//hd->GetXaxis()->SetTitle("Rho (cm)");
//hd->GetYaxis()->SetTitle("Conversions per bin");

hd->GetXaxis()->SetTitle(xtitle.c_str());
hd->GetYaxis()->SetTitle(ytitle.c_str());

hd->SetTitle(" ");

hd->SetLineColor(kBlue);
hd->SetLineWidth(2);
hd->Draw("ehist");

hm->SetLineWidth(2);
hm->SetLineColor(kRed);
//hm->Draw("histsame");

hd->Print();
hm->Print();

// Also need a legend ...
const float dx=0.28;
const float dy=0.20;
//float xlmin = 0.58;
//float ylmin = 0.66;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
leg->SetLineWidth(2);
leg->SetTextFont(62);
leg->SetTextSize(0.04);
//leg->SetHeader("16 < n_{PV} < 36","C");
leg->AddEntry(hd,"2018C Data (All candidates)","l");
leg->AddEntry(hm,"2018C Data (Duplicates removed)","l");
leg->SetBorderSize(0);                          // Include border
//leg->Draw();

   TLatex * tex = new TLatex(0.18,0.91,"CMS Work in progress");
   tex->SetNDC();
   tex->SetTextFont(61);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex * tex2 = new TLatex(0.20,0.19,"2017 J/psi MC - effect of duplicate removal");
   tex2->SetNDC();
   tex2->SetTextFont(61);
   tex2->SetTextSize(0.03);
   tex2->SetLineWidth(2);
   tex2->Draw();

// Create include file with the histogram data
/*
   ofstream outfile;
   string outfilename = "myweights.h2";
   outfile.open(outfilename);
   outfile << "//" << endl;
   outfile << "// " << outfilename << endl;
   outfile << "// include file for weights" << endl
;
   outfile << "//" << endl;
   outfile << "vector<double> wt{" << endl;
   for (int i=0;i<20;i++){
      for (int j=0;j<5;j++){
         int ielement = 5*i + j;
         outfile << setw(11) << weights[ielement];
         if(ielement!=98)outfile << " , ";
      }
      outfile << endl;
   }
   outfile << "};" << endl;

   outfile.close();
*/

if(logscale){
c1->Print((histtype+"_Log.pdf").c_str());
}
else{
c1->Print((histtype+"ratioMC.pdf").c_str());
}

}
