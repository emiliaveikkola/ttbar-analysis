#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TLatex.h>
#include <vector>
#include <string>
#include <iostream>
#include "tdrstyle_mod22.C"

void resComparison() {
    setTDRStyle();

    const TString histName = "pjes";

    std::vector<TString> files = {
        "Summer24_TTtoLNu2Q_V9M_pt15_eta2.4_fixedres.root",
        "Summer24_TTtoLNu2Q_V9M_pt30_eta2.5_fixedres.root",
        "Summer24_TTtoLNu2Q_V9M_pt30_eta2.5_oldres.root",
        "Summer24_TTtoLNu2Q_V9M_pt30_eta2.5_fixedres_muonTight.root",
        "Summer24_TTtoLNu2Q_V9M_pt30_eta2.5_fixedres_muonTight_arithmeticptpair.root"
        // ,"file5.root"
    };

    std::vector<TString> labels = {
        "pt15_eta2.4_newJER",
        "pt30_eta2.5_newJER",
        "pt30_eta2.5_oldJER",
        "pt30_eta2.5_newJER_muonTight",
        "arit_pt30_eta2.5_newJER_muonTight"
        // ,"Sample 5"
    };

    if (files.size() != labels.size()) {
        std::cerr << "files and labels must have the same size" << std::endl;
        return;
    }

    std::vector<int> colors = {
        kBlack,
        kRed+1,
        kBlue+1,
        kGreen+2,
        kMagenta+1
    };

    std::vector<int> markers = {
        kFullCircle,
        kFullSquare,
        kFullTriangleUp,
        kFullDiamond,
        kOpenCircle
    };

    std::vector<TFile*> openedFiles;
    std::vector<TProfile*> hists;

    for (size_t i = 0; i < files.size(); ++i) {
        TFile* f = TFile::Open(files[i], "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Could not open file: " << files[i] << std::endl;
            return;
        }
        openedFiles.push_back(f);

        TProfile* h = dynamic_cast<TProfile*>(f->Get(histName));
        if (!h) {
            std::cerr << "Could not find " << histName << " in file: " << files[i] << std::endl;
            return;
        }

        TProfile* hc = static_cast<TProfile*>(h->Clone(Form("%s_clone_%zu", histName.Data(), i)));
        hc->SetDirectory(nullptr);
        hc->SetLineColor(colors[i % colors.size()]);
        hc->SetMarkerColor(colors[i % colors.size()]);
        hc->SetMarkerStyle(markers[i % markers.size()]);
        hc->SetMarkerSize(1.0);
        hc->SetLineWidth(2);

        hists.push_back(hc);
    }
    TH1D *h = tdrHist("h", "pjes", 0.9, 0.92, "p_{T,pair}", 0, 250);
    TCanvas* c1 = tdrCanvas("c1", h, 2024, 0, kSquare);

    tdrDraw(hists[0], "P", markers[0], colors[0], kSolid, colors[0], 0, 0);

    TLegend* leg = tdrLeg(0.40, 0.70, 0.88, 0.88);
    leg->SetTextSize(0.03);
    leg->AddEntry(hists[0], labels[0], "pl");

    for (size_t i = 1; i < hists.size(); ++i) {
        tdrDraw(hists[i], "P SAME", markers[i % markers.size()], colors[i % colors.size()],
                kSolid, colors[i % colors.size()], 0, 0);
        leg->AddEntry(hists[i], labels[i], "pl");
    }

    c1->SaveAs("resComparison.pdf");

    for (auto* f : openedFiles) {
        if (f) f->Close();
    }
}