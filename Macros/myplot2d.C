void myplot2d(string histtype="xywideHist", float ymax=3000000.0, float ymin=900.0)
{

TCanvas *c1 = new TCanvas("c1","multipads",900,900);

TFile *fd = new TFile("../PC_DataHPC.root");
TFile *fm = new TFile("../PC_MCHPC.root");

TH1D * hd = (TH1D*)fd->Get(histtype.c_str());
//TH1D * hd = (TH1D*)fm->Get(histtype.c_str());

//hd->GetYaxis()->SetTitleOffset(1.4);

//hd->SetMaximum(ymax);
//hd->SetMinimum(ymin);
hd->SetStats(kFALSE);
//hm->SetStats(kFALSE);

c1->SetTicks(1,1);

//double scale = 3.19688/2.48492;
//hm->Scale(scale);

hd->SetTitle(" ");
//hd->SetTitle(" ");

//hd->SetLineColor(kBlue);
hd->Draw("colz");

/*
hm->SetLineColor(kRed);
hm->Draw("histsame");

hd->Print();
hm->Print();
*/

// Also need a legend ...
/*
const float dx=0.20;
const float dy=0.08;
float xlmin = 0.58;
float ylmin = 0.66;
TLegend* leg = new TLegend(xlmin,ylmin,xlmin+dx,ylmin+dy);
leg->SetHeader("CMS Work in progress","C");
//leg->AddEntry(hd,"2018 Data","l");
//leg->AddEntry(hm,"2018 MC","l");
leg->SetBorderSize(0);                          // Include border
leg->Draw();
*/

   TLatex * tex = new TLatex(0.11,0.91,"CMS Work in progress");
   tex->SetNDC();
   tex->SetTextFont(61);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();

/*
   tex = new TLatex(0.14,0.72,"Work in progress");
   tex->SetNDC();
   tex->SetTextFont(61);
   tex->SetTextSize(0.040);
   tex->SetLineWidth(2);
   tex->Draw();
*/

   c1->SetRightMargin(0.12);





}
