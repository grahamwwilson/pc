void myplotr(string histtype="npv_rcutHist", float ymax=150000.0, float ymin=40.0)
{

TCanvas *c1 = new TCanvas("c1","multipads",800,600);

TFile *fd = new TFile("../PC_DataAllHPC.root");
TFile *fm = new TFile("../PC_MCHPC.root");

TH2D * hd2d = (TH2D*)fd->Get(histtype.c_str());
TH2D * hm2d = (TH2D*)fm->Get(histtype.c_str());

// Create two TH1D's for the projection.

TH1D* hd = new TH1D("hd", ";R (cm)", 250, 0.0, 25.0);
TH1D* hm = new TH1D("hm", ";R (cm)", 250, 0.0, 25.0);

TH1D * hdpv = (TH1D*)fd->Get("numpvHist");
TH1D * hmpv = (TH1D*)fm->Get("numpvHist");

// Now loop over the 2d histogram filling the 1d histogram

int ipu = 27;

for (int i=1; i<=250; i++){
    hd->Fill( (double(i)-0.5)/10.0, double(hd2d->GetBinContent(i,ipu)) );
    hm->Fill( (double(i)-0.5)/10.0, double(hm2d->GetBinContent(i,ipu)) );
} 

// Find total number of interactions from PVs 
long int npvdtotal = 0;
long int npvmtotal = 0;
for (int kmult=26; kmult<=26;kmult++){
    int i = kmult + 1;
    int nbindcontent = hdpv->GetBinContent(i);
    int nbinmcontent = hmpv->GetBinContent(i);
    cout << "PV multiplicity " << kmult << " Bin number " << i << " " << nbindcontent << " " << nbinmcontent << endl;
    npvdtotal += kmult*nbindcontent;
    npvmtotal += kmult*nbinmcontent;
    cout << "Running suma " << npvdtotal << " " << npvmtotal << endl;
}
double ratio = double(npvdtotal)/double(npvmtotal);
cout << "Ratio = " << ratio << endl;

hd->GetYaxis()->SetTitleOffset(1.4);
hd->SetMaximum(ymax);
hd->SetMinimum(ymin);
hd->SetStats(kFALSE);
hm->SetStats(kFALSE);

c1->SetTicks(1,1);

/*
double scale = 3.19688/2.48492;
cout << "Scale = " << scale << endl;
*/
double scale=ratio;  // Overwrite with PV based one.
hm->Scale(scale);

hd->GetYaxis()->SetTitle("Conversion vertices per mm");
hd->SetTitle(" ");

hd->SetLineColor(kBlue);
hd->SetLineWidth(2);
hd->Draw("hist");

hm->SetLineWidth(2);
hm->SetLineColor(kRed);
hm->Draw("histsame");

hd->Print();
hm->Print();

// Also need a legend ...
const float dx=0.28;
const float dy=0.20;
float xlmin = 0.58;
float ylmin = 0.66;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
leg->SetLineWidth(2);
leg->SetTextFont(62);
leg->SetTextSize(0.05);
leg->SetHeader("n_{PV} = 26","C");
leg->AddEntry(hd,"2018 Data","l");
leg->AddEntry(hm,"2018 MC","l");
leg->SetBorderSize(0);                          // Include border
leg->Draw();

   TLatex * tex = new TLatex(0.11,0.91,"CMS Work in progress");
   tex->SetNDC();
   tex->SetTextFont(61);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();

}
