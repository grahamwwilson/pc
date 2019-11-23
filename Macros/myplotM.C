void myplotM(string histtype="mgg1Hist", float ymax=0.5e5, float ymin=1.0e4)
{

TCanvas *c1 = new TCanvas("c1","multipads",800,600);

TFile *fd = new TFile("PC_DataAllHPC.root");
TFile *fm = new TFile("PC_MCHPC.root");

TH1D * hd = (TH1D*)fd->Get(histtype.c_str());
TH1D * hm = (TH1D*)fm->Get(histtype.c_str());

hd->GetYaxis()->SetTitleOffset(1.4);
hd->SetMaximum(ymax);
hd->SetMinimum(ymin);
hd->SetStats(kFALSE);
hm->SetStats(kFALSE);

c1->SetTicks(1,1);

double scale = 2.94407/2.05164;
//double scale = 1.0;
hm->Scale(scale);

hd->GetYaxis()->SetTitle("Di-photons per bin");
hd->SetTitle(" ");

hd->SetLineColor(kBlue);
hd->SetLineWidth(2);
hd->Draw("ehist");

hm->SetLineWidth(2);
hm->SetLineColor(kRed);
hm->Draw("histsame");

hd->Print();
hm->Print();

cout << hd->Integral(17,100) << endl;
cout << hm->Integral(17,100) << endl;

cout << hd->Integral(81,100) << endl;
cout << hm->Integral(81,100) << endl;


// Also need a legend ...
const float dx=0.28;
const float dy=0.20;
float xlmin = 0.58;
float ylmin = 0.66;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
leg->SetLineWidth(2);
leg->SetTextFont(62);
leg->SetTextSize(0.05);
leg->SetHeader("All n_{PV}","C");
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
