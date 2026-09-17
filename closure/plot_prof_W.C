#include "../tdrstyle_mod22.C"

#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TAxis.h>

#include <iostream>
#include <vector>
#include <string>

void plot_prof_W()
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

    std::vector<int> colors = {
        kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
        kCyan + 2, kOrange + 7, kViolet + 1, kAzure + 1,
        kSpring + 5, kPink + 7, kTeal + 3, kGray + 2,
        kRed - 4, kBlue - 4, kGreen - 6, kMagenta - 4,
        kOrange - 3, kCyan - 6, kViolet - 6
    };

    std::vector<TFile*> openedFiles;
    std::vector<TProfile*> profiles;

     // Change these ranges if needed.
    // If prof_W_inWindow is a W-mass-like profile, 70--90 is probably better.
    double xmin = 0.0;
    double xmax = 250.0;
    double ymin = 40;
    double ymax = 100;

    // If prof_W_inWindow is W mass-like, use e.g. 70--90 instead of 0--1.2.
    TH1D *h = tdrHist("h_prof_W",
                      "prof_W_inWindow",
                      ymin, ymax,
                      "ptpair",
                      xmin, xmax);
    lumi_136TeV = "";
    extraText = "Private work";
    TCanvas *c = tdrCanvas("c_prof_W", h, 8, 11, kRectangular);

    TLegend *leg = tdrLeg(0.40, 0.15, 0.88, 0.5);
    leg->SetNColumns(2);
    leg->SetTextSize(0.022);

    for (size_t i = 0; i < files.size(); ++i) {
        const std::string& path = files[i];

        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "[WARNING] Could not open file: " << path << std::endl;
            continue;
        }

        TProfile* p = dynamic_cast<TProfile*>(f->Get("prof_W_inWindow_ptpair"));
        if (!p) {
            std::cerr << "[WARNING] Could not find prof_W_inWindow_ptpair in: "
                      << path << std::endl;
            f->Close();
            continue;
        }

        TProfile* pc = dynamic_cast<TProfile*>(p->Clone(Form("prof_W_inWindow_ptpair_%zu", i)));
        pc->SetDirectory(nullptr);

        openedFiles.push_back(f);
        profiles.push_back(pc);

        int col = colors[i % colors.size()];
        int marker = 20 + (i % 10);

        pc->SetLineColor(col);
        pc->SetMarkerColor(col);
        pc->SetMarkerStyle(marker);
        pc->SetMarkerSize(0.7);
        pc->SetLineWidth(2);

        // TDR drawing
        tdrDraw(pc, "E1 SAME", marker, col, kSolid, col, 1001, col);

        TString label = path;
        label.ReplaceAll("closure/V2/", "");
        label.ReplaceAll("Muon_Run", "");
        label.ReplaceAll("_Golden_e7.root", "");
        label.ReplaceAll("_MLEnhancedGolden_Latest1.6._e7.root", "");
        label.ReplaceAll("_Prompt_", " Prompt ");
        label.ReplaceAll("_ReReco_", " ReReco ");

        leg->AddEntry(pc, label, "lep");
    }

    if (profiles.empty()) {
        std::cerr << "[ERROR] No profiles were found. Check file paths and object name." << std::endl;
        return;
    }

    leg->Draw();

    lumi_136TeV = "";
    extraText = "Private work";
    CMS_lumi(c, 0, 0);

    c->SaveAs("prof_W_inWindow.pdf");

    std::cout << "[OK] Saved prof_W_inWindow.pdf and prof_W_inWindow.png" << std::endl;
}