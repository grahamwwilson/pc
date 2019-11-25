void mymedian(string histtype="pTHist", float ymax=3000000.0, float ymin=900.0)
{

TCanvas *c1 = new TCanvas("c1","multipads",800,600);

TFile *fd = new TFile("../PC_dataHPC.root");
TFile *fm = new TFile("../PC_mcHPC.root");

TH1D * hd = (TH1D*)fd->Get(histtype.c_str());
TH1D * hm = (TH1D*)fm->Get(histtype.c_str());

int n = hd->GetXaxis()->GetNbins();  
std::vector<double>  x(n);
hd->GetXaxis()->GetCenter( &x[0] );
const double * y = hd->GetArray(); 
// exclude underflow/overflows from bin content array y
cout << "Median = " << TMath::Median(n, &x[0], &y[1]) << endl;

//cout << Median(&hd) << endl;
//cout << Median(&hm) << endl;


TH1D * hdpv = (TH1D*)fd->Get("numpvHist");
TH1D * hmpv = (TH1D*)fm->Get("numpvHist");

// Find total number of interactions from PVs 
long int npvdtotal = 0;
long int npvmtotal = 0;
for (int kmult=17; kmult<=35;kmult++){
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

double scale = 3.19688/2.48492;
cout << "Scale = " << scale << endl;
scale=ratio;  // Overwrite with PV based one.
hm->Scale(scale);

hd->GetYaxis()->SetTitle("Conversion vertices per mm");
hd->SetTitle(" ");

hd->SetLineColor(kBlue);
hd->SetLineWidth(2);
hd->Draw("ehist");

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
leg->SetHeader("16 < n_{PV} < 36","C");
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
