#include "../tdrstyle_mod22.C"

#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TGraphErrors.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

void plot_Wmean_vs_era()
{
    setTDRStyle();

    std::vector<std::string> files = {
        "closure/V2/Muon_Run2024C_ReReco_V11M_Golden_e7.root",
        "closure/V2/Muon_Run2024D_ReReco_V11M_Golden_e7.root",
        "closure/V2/Muon_Run2024E_ReReco_V11M_Golden_e7.root",

        "closure/V2/Muon_Run2024F_nib1_Prompt_V11M_Golden_e7.root",
        "closure/V2/Muon_Run2024F_nib2_Prompt_V11M_Golden_e7.root",
        "closure/V2/Muon_Run2024F_nib3_Prompt_V11M_Golden_e7.root",
        "closure/V2/Muon_Run2024F_Prompt_V11M_Golden_e7.root",

        "closure/V2/Muon_Run2024G_nib1_Prompt_V11M_Golden_e7.root",
        "closure/V2/Muon_Run2024G_nib2_Prompt_V11M_Golden_e7.root",
        "closure/V2/Muon_Run2024G_Prompt_V11M_Golden_e7.root",

        "closure/V2/Muon_Run2024H_Prompt_V11M_Golden_e7.root",
        "closure/V2/Muon_Run2024I_Prompt_V11M_Golden_e7.root",

        "closure/V2/Muon_Run2025C_Prompt_V5M_Golden_e7.root",
        "closure/V2/Muon_Run2025D_Prompt_V5M_Golden_e7.root",
        "closure/V2/Muon_Run2025E_Prompt_V5M_Golden_e7.root",
        "closure/V2/Muon_Run2025F_Prompt_V5M_Golden_e7.root",
        "closure/V2/Muon_Run2025G_Prompt_V5M_Golden_e7.root",

        "closure/V2/Muon_Run2026B_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root",
        "closure/V2/Muon_Run2026C_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root",
        "closure/V2/Muon_Run2026D_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root"
    };

    std::vector<std::string> labels = {
        "2024C",
        "2024D",
        "2024E",

        "2024F nib1",
        "2024F nib2",
        "2024F nib3",
        "2024F",

        "2024G nib1",
        "2024G nib2",
        "2024G",

        "2024H",
        "2024I",

        "2025C",
        "2025D",
        "2025E",
        "2025F",
        "2025G",

        "2026B",
        "2026C",
        "2026D"
    };

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ex;
    std::vector<double> ey;

    for (size_t i = 0; i < files.size(); ++i) {
        TFile* f = TFile::Open(files[i].c_str(), "READ");

        if (!f || f->IsZombie()) {
            std::cerr << "[WARNING] Could not open file: " << files[i] << std::endl;
            continue;
        }

        TProfile* p = dynamic_cast<TProfile*>(f->Get("prof_W_inWindow"));

        if (!p) {
            std::cerr << "[WARNING] Could not find prof_W_inWindow in: "
                      << files[i] << std::endl;
            f->Close();
            continue;
        }

        // Mean of the profile y-values over bins with entries.
        double sumw = 0.0;
        double sumwy = 0.0;

        for (int ibin = 1; ibin <= p->GetNbinsX(); ++ibin) {
            double n = p->GetBinEntries(ibin);
            if (n <= 0) continue;

            double yi = p->GetBinContent(ibin);

            sumw  += n;
            sumwy += n * yi;
        }

        if (sumw <= 0) {
            std::cerr << "[WARNING] Profile has no filled bins in: "
                      << files[i] << std::endl;
            f->Close();
            continue;
        }

        double mean = sumwy / sumw;

        // Approximate uncertainty from profile RMS / sqrt(N effective).
        // This is useful for quick comparison plots.
        double variance = 0.0;

        for (int ibin = 1; ibin <= p->GetNbinsX(); ++ibin) {
            double n = p->GetBinEntries(ibin);
            if (n <= 0) continue;

            double yi = p->GetBinContent(ibin);
            variance += n * (yi - mean) * (yi - mean);
        }

        double err = 0.0;
        if (sumw > 1) {
            err = std::sqrt(variance / (sumw * (sumw - 1.0)));
        }

        x.push_back(x.size() + 1.0);
        y.push_back(mean);
        ex.push_back(0.0);
        ey.push_back(err);

        std::cout << labels[i] << "  mean = " << mean
                  << " +/- " << err
                  << "  entries = " << sumw << std::endl;

        f->Close();
    }

    if (x.empty()) {
        std::cerr << "[ERROR] No valid points were found." << std::endl;
        return;
    }

    const int n = x.size();


    double xmin = 0.5;
    double xmax = n + 0.5;

    // Adjust these depending on your result.
    double ymin = 70.0;
    double ymax = 90.0;

    TH1D* h = tdrHist("h_Wmean_vs_era",
                      "Mean of prof_W_inWindow",
                      ymin, ymax,
                      "Year / era",
                      xmin, xmax);

    h->GetXaxis()->SetNdivisions(n, false);

    for (int i = 0; i < n; ++i) {
        h->GetXaxis()->SetBinLabel(h->FindBin(x[i]), labels[i].c_str());
    }

    h->GetXaxis()->LabelsOption("v");
    h->GetXaxis()->SetLabelSize(0.035);

    lumi_136TeV = "";
    extraText = "Private work";

    TCanvas* c = tdrCanvas("c_Wmean_vs_era", h, 8, 11, kRectangular);
std::vector<int> colors = {
    kBlack,
    kRed + 1,
    kBlue + 1,
    kGreen + 2,
    kMagenta + 1,
    kCyan + 2,
    kOrange + 7,
    kViolet + 1,
    kAzure + 1,
    kSpring + 5,
    kPink + 7,
    kTeal + 3,
    kGray + 2,
    kRed - 4,
    kBlue - 4,
    kGreen - 6,
    kMagenta - 4,
    kOrange - 3,
    kCyan - 6,
    kViolet - 6
};

std::vector<TGraphErrors*> graphs;

for (int i = 0; i < n; ++i) {
    TGraphErrors* gr = new TGraphErrors(1);

    gr->SetPoint(0, x[i], y[i]);
    gr->SetPointError(0, ex[i], ey[i]);

    int col = colors[i % colors.size()];
    int marker = 20 + (i % 10);

    gr->SetMarkerSize(1.2);
    gr->SetLineWidth(2);

    tdrDraw(gr, "Pz", marker, col, kSolid, col);

    graphs.push_back(gr);
}

    c->RedrawAxis();

    c->SaveAs("Wmean_vs_era.pdf");

    std::cout << "[OK] Saved Wmean_vs_era.pdf and Wmean_vs_era.png" << std::endl;
}