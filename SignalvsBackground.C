#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include <TLatex.h>
#include "tdrstyle_mod22.C"
#include "TBox.h"

void SignalvsBackground() {
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  TFile *file1 = new TFile("Winter24_TTtoLNu2Q.root", "READ");

  TH2D* h2_control    = (TH2D*)file1->Get("h2_W_ptpair_cr");
  TH2D* h2_background = (TH2D*)file1->Get("h2_W_ptpair_bkg");
  TH2D* h2_signal     = (TH2D*)file1->Get("h2_W_ptpair_sig");
  
  TH1D *h1 = tdrHist("h1","Mass (GeV)",15,250.,"ptpair",15,230);
  //lumi_13TeV = "Run2";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);

  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetLabelSize(0.045);
  h1->GetYaxis()->SetLabelSize(0.045);


// ——— Custom axis ranges ———
// User can modify these to zoom in on a region of interest
double xMin = 15.0, xMax = 230.0;
double yMin = 15.0, yMax = 249.0;
// apply to the 2D “signal” axis, which is used to draw the frame
h2_control->GetXaxis()->SetRangeUser(xMin, xMax);
h2_control->GetYaxis()->SetRangeUser(yMin, yMax);
// ————————————————

  // Find maxima to normalize alpha
  double maxSig = h2_signal->GetMaximum();
  double maxBkg = h2_background->GetMaximum();
  double maxCR  = h2_control->GetMaximum();
  int nx = h2_control->GetNbinsX(), ny = h2_control->GetNbinsY();

  auto paint = [&](TH2D* h, Int_t color, double max, double globalAlpha){
    for(int ix=1; ix<=nx; ++ix){
      double x1 = h->GetXaxis()->GetBinLowEdge(ix);
      double x2 = h->GetXaxis()->GetBinUpEdge(ix);
      if (x2 < xMin || x1 > xMax) continue;
      for(int iy=1; iy<=ny; ++iy){
        double y1 = h->GetYaxis()->GetBinLowEdge(iy);
        double y2 = h->GetYaxis()->GetBinUpEdge(iy);
        if (y2 < yMin || y1 > yMax) continue;
        double val = h->GetBinContent(ix,iy);
        if(val<=0) continue;
        double alpha = std::min(val/max*globalAlpha, 1.0);
        TBox *box = new TBox(x1,y1,x2,y2);
        box->SetFillColorAlpha(color, alpha);
        box->SetLineColor(color);
        box->Draw("same");
      }
    }
  };

  // Paint control (green), background (blue), and signal (red)
  paint(h2_control,  kGreen+2, maxCR,  1);
  paint(h2_background, kBlue-4,  maxBkg, 1);
  paint(h2_signal,    kRed+1,   maxSig, 0.7);
  
  TLegend *legend = tdrLeg(0.72, 0.75-0.015*3, 0.9, 0.9);
  TH1D* dummySignal = new TH1D("dummySignal", "", 1, 0, 1);
  dummySignal->SetFillColorAlpha(kRed+1, 0.7);
  dummySignal->SetLineColor(kRed+1);

  TH1D* dummyBackground = new TH1D("dummyBackground", "", 1, 0, 1);
  dummyBackground->SetFillColorAlpha(kBlue-4, 0.7);
  dummyBackground->SetLineColor(kBlue-4);

  TH1D* dummyControl = new TH1D("dummyControl", "", 1, 0, 1);
  dummyControl->SetFillColorAlpha(kGreen+2, 1);
  dummyControl->SetLineColor(kGreen+2);

  legend->AddEntry(dummySignal,    "Signal",    "f");
  legend->AddEntry(dummyBackground,"Background","f");
  legend->AddEntry(dummyControl,   "Control",   "f");
  legend->SetTextSize(0.03);

  c1->RedrawAxis();

  TFile* fout = new TFile("Combined.root", "RECREATE");
  fout->cd();
  // write the canvas so that the overlaid colors are preserved
  c1->Write("c_SignalvsBackground");
  fout->Close();
  c1->Update();
  c1->SaveAs("pdf/SignalvsBackground.pdf");   // Save as PDF
}