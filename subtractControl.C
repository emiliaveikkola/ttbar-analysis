// subtractControl.C
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>
#include <algorithm>
#include <TROOT.h>

void subtractControl(){
  // 1) Open input and create output
  TFile *fin  = new TFile("Winter24_TTtoLNu2Q.root", "READ");
  TFile *fout = new TFile("sigback.root", "RECREATE");
  
  // 2) Retrieve the four profiles & TH2Ds
  TH2D*     h2_sig    = (TH2D*)fin->Get("h2_Wmass_ptpair_qq");
  TH2D*     h2_bkg1   = (TH2D*)fin->Get("h2_Wmass_ptpair_gg");
  TH2D*     h2_bkg2   = (TH2D*)fin->Get("h2_Wmass_ptpair_qg");
  TH2D*     h2_bkg3   = (TH2D*)fin->Get("h2_Wmass_ptpair_qq_others");
  
  
  // Threshold: only keep bins with more than 100 entries
  const double threshold = 100;
  // Apply threshold to all three histograms
  for (auto hist : {h2_sig, h2_bkg1, h2_bkg2, h2_bkg3}) {
      int nx = hist->GetNbinsX();
      int ny = hist->GetNbinsY();
      for (int ix = 1; ix <= nx; ++ix) {
          for (int iy = 1; iy <= ny; ++iy) {
              if (hist->GetBinContent(ix, iy) <= threshold) {
                  hist->SetBinContent(ix, iy, 0);
              }
          }
      }
  }

 const double threshold3 = 500;
  // Apply threshold to all three histograms
  for (auto hist : {h2_sig, h2_bkg2, h2_bkg3}) {
      int nx = hist->GetNbinsX();
      int ny = hist->GetNbinsY();
      for (int ix = 1; ix <= nx; ++ix) {
          for (int iy = 1; iy <= ny; ++iy) {
              if (hist->GetBinContent(ix, iy) <= threshold3) {
                  hist->SetBinContent(ix, iy, 0);
              }
          }
      }
  }
  const double threshold4 = 1000;
    for (auto hist : {h2_sig, h2_bkg2}) {
      int nx = hist->GetNbinsX();
      int ny = hist->GetNbinsY();
      for (int ix = 1; ix <= nx; ++ix) {
          for (int iy = 1; iy <= ny; ++iy) {
              if (hist->GetBinContent(ix, iy) <= threshold4) {
                  hist->SetBinContent(ix, iy, 0);
              }
          }
      }
  }
/*
    const double threshold5 = 3000;
    for (auto hist : {h2_sig}) {
      int nx = hist->GetNbinsX();
      int ny = hist->GetNbinsY();
      for (int ix = 1; ix <= nx; ++ix) {
          for (int iy = 1; iy <= ny; ++iy) {
              if (hist->GetBinContent(ix, iy) <= threshold5) {
                  hist->SetBinContent(ix, iy, 0);
              }
          }
      }
  }
  */
    // Print all non-zero bins
    for (auto hist : {h2_sig, h2_bkg1, h2_bkg2, h2_bkg3}) {
        std::cout << "Non-zero bins for " << hist->GetName() << ":\n";
        int nx = hist->GetNbinsX();
        int ny = hist->GetNbinsY();
        for (int ix = 1; ix <= nx; ++ix) {
            for (int iy = 1; iy <= ny; ++iy) {
                double content = hist->GetBinContent(ix, iy);
                if (content != 0) {
                    double xlow = hist->GetXaxis()->GetBinLowEdge(ix);
                    double xup = hist->GetXaxis()->GetBinUpEdge(ix);
                    double ylow = hist->GetYaxis()->GetBinLowEdge(iy);
                    double yup = hist->GetYaxis()->GetBinUpEdge(iy);
                    std::cout << "  " << hist->GetName()
                              << " bin(" << ix << "," << iy << "), x=[" << xlow << "," << xup << "]"
                              << ", y=[" << ylow << "," << yup << "] = "
                              << content << "\n";
                }
            }
        }
    }

      // Combine histograms into unified signal and background
  TH2D* h2_signal = (TH2D*)h2_sig->Clone("h2_Wmass_ptpair_signal");
  TH2D* h2_background = (TH2D*)h2_bkg1->Clone("h2_Wmass_ptpair_background");
  h2_background->Reset();
  h2_background->Add(h2_bkg1);
  h2_background->Add(h2_bkg2);
  h2_background->Add(h2_bkg3);

      for (auto hist : {h2_signal, h2_background}) {
        std::cout << "Non-zero bins for " << hist->GetName() << ":\n";
        int nx = hist->GetNbinsX();
        int ny = hist->GetNbinsY();
        for (int ix = 1; ix <= nx; ++ix) {
            for (int iy = 1; iy <= ny; ++iy) {
                double content = hist->GetBinContent(ix, iy);
                if (content != 0) {
                    double xlow = hist->GetXaxis()->GetBinLowEdge(ix);
                    double xup = hist->GetXaxis()->GetBinUpEdge(ix);
                    double ylow = hist->GetYaxis()->GetBinLowEdge(iy);
                    double yup = hist->GetYaxis()->GetBinUpEdge(iy);
                    std::cout << "  " << hist->GetName()
                              << " bin(" << ix << "," << iy << "), x=[" << xlow << "," << xup << "]"
                              << ", y=[" << ylow << "," << yup << "] = "
                              << content << "\n";
                }
            }
        }
    }

  // Write unified histograms
  h2_signal->Write();
  h2_background->Write();

  h2_sig->Write();
  h2_bkg1->Write();
  h2_bkg2->Write();
  h2_bkg3->Write();
  fout->Close();
  
}