void myplotlow(string histtype="r1dwidelowPUcutHist", float ymax=300000.0, float ymin=90.0)
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

double scale = 3.68575/1.741423;
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
leg->SetHeader("n_{PV} < 17","C");
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
