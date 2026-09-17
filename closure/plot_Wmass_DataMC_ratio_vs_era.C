#include "../tdrstyle_mod22.C"

#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TBox.h>
#include <TLatex.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

struct Sample {
    std::string dataFile;
    std::string label;
    std::string mcFile;
    std::string mcGroup;
};

struct MeanResult {
    double mean = 0.0;
    double err = 0.0;
    double entries = 0.0;
};

struct Point {
    std::string label;
    std::string mcGroup;
    double ratio = 0.0;
    double ratioErr = 0.0;
};

bool GetProfileMean(const std::string& fileName,
                    const std::string& profName,
                    MeanResult& result)
{
    TFile* f = TFile::Open(fileName.c_str(), "READ");

    if (!f || f->IsZombie()) {
        std::cerr << "[WARNING] Could not open file: " << fileName << std::endl;
        return false;
    }

    TProfile* p = dynamic_cast<TProfile*>(f->Get(profName.c_str()));

    if (!p) {
        std::cerr << "[WARNING] Could not find " << profName
                  << " in file: " << fileName << std::endl;
        f->Close();
        return false;
    }

    double sumw = 0.0;
    double sumwy = 0.0;

    for (int ibin = 1; ibin <= p->GetNbinsX(); ++ibin) {
        double n = p->GetBinEntries(ibin);
        if (n <= 0) continue;

        double y = p->GetBinContent(ibin);

        sumw  += n;
        sumwy += n * y;
    }

    if (sumw <= 0.0) {
        std::cerr << "[WARNING] Empty profile in file: " << fileName << std::endl;
        f->Close();
        return false;
    }

    result.mean = sumwy / sumw;
    result.entries = sumw;

    double variance = 0.0;

    for (int ibin = 1; ibin <= p->GetNbinsX(); ++ibin) {
        double n = p->GetBinEntries(ibin);
        if (n <= 0) continue;

        double y = p->GetBinContent(ibin);
        variance += n * (y - result.mean) * (y - result.mean);
    }

    result.err = 0.0;
    if (sumw > 1.0) {
        result.err = std::sqrt(variance / (sumw * (sumw - 1.0)));
    }

    f->Close();
    return true;
}

MeanResult CombineMeans(const std::vector<MeanResult>& results)
{
    MeanResult combined;

    double sumw = 0.0;
    double sumwy = 0.0;

    for (const auto& r : results) {
        if (r.entries <= 0.0) continue;

        sumw  += r.entries;
        sumwy += r.entries * r.mean;
    }

    if (sumw <= 0.0) return combined;

    combined.mean = sumwy / sumw;
    combined.entries = sumw;

    double errNumerator2 = 0.0;

    for (const auto& r : results) {
        if (r.entries <= 0.0) continue;

        errNumerator2 += std::pow(r.entries * r.err, 2);
    }

    combined.err = std::sqrt(errNumerator2) / sumw;

    return combined;
}

Point MakeRatioPoint(const std::string& label,
                     const std::string& mcGroup,
                     const MeanResult& data,
                     const MeanResult& mc)
{
    Point p;
    p.label = label;
    p.mcGroup = mcGroup;

    if (mc.mean == 0.0 || data.mean == 0.0) {
        p.ratio = 0.0;
        p.ratioErr = 0.0;
        return p;
    }

    p.ratio = data.mean / mc.mean;

    const double relDataErr = data.err / data.mean;
    const double relMCErr   = mc.err / mc.mean;

    p.ratioErr = p.ratio * std::sqrt(relDataErr * relDataErr +
                                     relMCErr   * relMCErr);

    return p;
}

int ColorForMCGroup(const std::string& mcGroup)
{
    if (mcGroup == "2024 MC")   return kBlue - 10;
    if (mcGroup == "2025 MC")   return kGreen - 10;
    if (mcGroup == "2026BD MC") return kOrange - 9;
    if (mcGroup == "2026C MC")  return kViolet - 9;

    return kGray;
}

void DrawMCGroupShading(const std::vector<Point>& points,
                        TH1D* h,
                        double ymin,
                        double ymax)
{
    if (points.empty()) return;

    std::vector<TBox*> boxes;

    int start = 0;
    std::string currentGroup = points[0].mcGroup;

    for (int i = 1; i <= (int)points.size(); ++i) {
        bool endGroup = false;

        if (i == (int)points.size()) {
            endGroup = true;
        }
        else if (points[i].mcGroup != currentGroup) {
            endGroup = true;
        }

        if (endGroup) {
            double x1 = start + 0.5;
            double x2 = i + 0.5;

            TBox* box = new TBox(x1, ymin, x2, ymax);
            box->SetFillColorAlpha(ColorForMCGroup(currentGroup), 0.28);
            box->SetLineColor(0);
            box->Draw("SAME");
            boxes.push_back(box);

            start = i;
            if (i < (int)points.size()) {
                currentGroup = points[i].mcGroup;
            }
        }
    }

    h->Draw("AXIS SAME");

    start = 0;
    currentGroup = points[0].mcGroup;

    for (int i = 1; i <= (int)points.size(); ++i) {
        bool endGroup = false;

        if (i == (int)points.size()) {
            endGroup = true;
        }
        else if (points[i].mcGroup != currentGroup) {
            endGroup = true;
        }

        if (endGroup) {
            double xCenter = 0.5 * ((start + 1.0) + i);

            start = i;
            if (i < (int)points.size()) {
                currentGroup = points[i].mcGroup;
            }
        }
    }
}

void plot_Wmass_DataMC_ratio_vs_era()
{
    setTDRStyle();

    const std::string profName = "prof_W_inWindow_ptpair";

    const std::string mc2024 =
        "closure/V2/Summer24_TTtoLNu2Q_JMENano_V11M_JER2024nib_e7.root";

    const std::string mc2025 =
        "closure/V2/Summer24_TTtoLNu2Q_JMENano_V5M_JER2025CDEFG_e7.root";

    const std::string mc2026BD =
        "closure/V2/Summer24_TTtoLNu2Q_JMENano_V2M_JER2026BD_e7.root";

    const std::string mc2026C =
        "closure/V2/Summer24_TTtoLNu2Q_JMENano_V2M_JER2026C_e7.root";

    std::vector<Sample> samples = {
        {"closure/V2/Muon_Run2024C_ReReco_V11M_Golden_e7.root",
         "2024C", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024D_ReReco_V11M_Golden_e7.root",
         "2024D", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024E_ReReco_V11M_Golden_e7.root",
         "2024E", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024F_nib1_Prompt_V11M_Golden_e7.root",
         "2024F nib1", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024F_nib2_Prompt_V11M_Golden_e7.root",
         "2024F nib2", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024F_nib3_Prompt_V11M_Golden_e7.root",
         "2024F nib3", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024F_Prompt_V11M_Golden_e7.root",
         "2024F", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024G_nib1_Prompt_V11M_Golden_e7.root",
         "2024G nib1", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024G_nib2_Prompt_V11M_Golden_e7.root",
         "2024G nib2", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024G_Prompt_V11M_Golden_e7.root",
         "2024G", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024H_Prompt_V11M_Golden_e7.root",
         "2024H", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2024I_Prompt_V11M_Golden_e7.root",
         "2024I", mc2024, "2024 MC"},

        {"closure/V2/Muon_Run2025C_Prompt_V5M_Golden_e7.root",
         "2025C", mc2025, "2025 MC"},

        {"closure/V2/Muon_Run2025D_Prompt_V5M_Golden_e7.root",
         "2025D", mc2025, "2025 MC"},

        {"closure/V2/Muon_Run2025E_Prompt_V5M_Golden_e7.root",
         "2025E", mc2025, "2025 MC"},

        {"closure/V2/Muon_Run2025F_Prompt_V5M_Golden_e7.root",
         "2025F", mc2025, "2025 MC"},

        {"closure/V2/Muon_Run2025G_Prompt_V5M_Golden_e7.root",
         "2025G", mc2025, "2025 MC"},

        {"closure/V2/Muon_Run2026B_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root",
         "2026B", mc2026BD, "2026BD MC"},

        {"closure/V2/Muon_Run2026D_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root",
         "2026D", mc2026BD, "2026BD MC"},

        {"closure/V2/Muon_Run2026C_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root",
         "2026C", mc2026C, "2026C MC"}
    };

    std::vector<int> pointColors = {
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
        kViolet - 6,
        kBlack,
        kRed + 2,
        kBlue + 2,
        kGreen + 3
    };

    std::vector<Point> points;

    std::vector<MeanResult> data2024Results;
    std::vector<MeanResult> data2025Results;
    std::vector<MeanResult> data2026BDResults;
    std::vector<MeanResult> data2026CResults;
    std::vector<MeanResult> data2026AllResults;

    MeanResult mc2024Result;
    MeanResult mc2025Result;
    MeanResult mc2026BDResult;
    MeanResult mc2026CResult;

    const bool okMC2024 =
        GetProfileMean(mc2024, profName, mc2024Result);

    const bool okMC2025 =
        GetProfileMean(mc2025, profName, mc2025Result);

    const bool okMC2026BD =
        GetProfileMean(mc2026BD, profName, mc2026BDResult);

    const bool okMC2026C =
        GetProfileMean(mc2026C, profName, mc2026CResult);

    MeanResult mc2026CombinedResult =
        CombineMeans({mc2026BDResult, mc2026CResult});

    for (size_t i = 0; i < samples.size(); ++i) {
        MeanResult dataResult;
        MeanResult mcResult;

        const bool okData =
            GetProfileMean(samples[i].dataFile, profName, dataResult);

        const bool okMC =
            GetProfileMean(samples[i].mcFile, profName, mcResult);

        if (!okData || !okMC || mcResult.mean == 0.0) {
            std::cerr << "[WARNING] Skipping " << samples[i].label << std::endl;
            continue;
        }

        Point p = MakeRatioPoint(samples[i].label,
                                 samples[i].mcGroup,
                                 dataResult,
                                 mcResult);

        points.push_back(p);

        if (samples[i].label.rfind("2024", 0) == 0) {
            data2024Results.push_back(dataResult);
        }
        else if (samples[i].label.rfind("2025", 0) == 0) {
            data2025Results.push_back(dataResult);
        }
        else if (samples[i].label == "2026B" || samples[i].label == "2026D") {
            data2026BDResults.push_back(dataResult);
            data2026AllResults.push_back(dataResult);
        }
        else if (samples[i].label == "2026C") {
            data2026CResults.push_back(dataResult);
            data2026AllResults.push_back(dataResult);
        }

        std::cout << samples[i].label
                  << "  data mean = " << dataResult.mean
                  << " +/- " << dataResult.err
                  << "  MC mean = " << mcResult.mean
                  << " +/- " << mcResult.err
                  << "  Data/MC = " << p.ratio
                  << " +/- " << p.ratioErr
                  << "  divided by " << samples[i].mcGroup
                  << std::endl;
    }

    MeanResult data2024Full = CombineMeans(data2024Results);
    MeanResult data2025Full = CombineMeans(data2025Results);
    MeanResult data2026BDFull = CombineMeans(data2026BDResults);
    MeanResult data2026CFull = CombineMeans(data2026CResults);
    MeanResult data2026Full = CombineMeans(data2026AllResults);

    if (okMC2024 && data2024Full.entries > 0.0) {
        Point p = MakeRatioPoint("2024",
                                 "2024 MC",
                                 data2024Full,
                                 mc2024Result);
        points.push_back(p);
    }

    if (okMC2025 && data2025Full.entries > 0.0) {
        Point p = MakeRatioPoint("2025",
                                 "2025 MC",
                                 data2025Full,
                                 mc2025Result);
        points.push_back(p);
    }

    if (okMC2026BD && data2026BDFull.entries > 0.0) {
        Point p = MakeRatioPoint("2026BD",
                                 "2026BD MC",
                                 data2026BDFull,
                                 mc2026BDResult);
        points.push_back(p);
    }

    if (okMC2026C && data2026CFull.entries > 0.0) {
        Point p = MakeRatioPoint("2026C",
                                 "2026C MC",
                                 data2026CFull,
                                 mc2026CResult);
        points.push_back(p);
    }

    if (okMC2026BD && okMC2026C && data2026Full.entries > 0.0) {
        Point p = MakeRatioPoint("2026",
                                 "2026 MC",
                                 data2026Full,
                                 mc2026CombinedResult);
        points.push_back(p);
    }

    if (points.empty()) {
        std::cerr << "[ERROR] No valid Data/MC ratio points found." << std::endl;
        return;
    }

    const int n = points.size();

    double xmin = 0.5;
    double xmax = n + 0.5;

    double ymin = 0.95;
    double ymax = 1.03;

    TH1D* h = tdrHist("h_Wmass_DataMC_ratio_vs_era",
                      "Data / MC  #LTm_{W}#GT",
                      ymin, ymax,
                      "Year / era",
                      xmin, xmax);

    h->GetXaxis()->SetNdivisions(n, false);

    for (int i = 0; i < n; ++i) {
        h->GetXaxis()->SetBinLabel(h->FindBin(i + 1.0),
                                   points[i].label.c_str());
    }

    h->GetXaxis()->LabelsOption("v");
    h->GetXaxis()->SetLabelSize(0.027);

    lumi_136TeV = "";
    extraText = "Private work";

    TCanvas* c = tdrCanvas("c_Wmass_DataMC_ratio_vs_era",
                           h,
                           8,
                           11,
                           kRectangular);

    DrawMCGroupShading(points, h, ymin, ymax);

    TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray + 2);
    line->SetLineWidth(2);
    line->Draw("SAME");

    std::vector<TGraphErrors*> graphs;

    for (int i = 0; i < n; ++i) {
        TGraphErrors* gr = new TGraphErrors(1);

        gr->SetPoint(0, i + 1.0, points[i].ratio);
        gr->SetPointError(0, 0.0, points[i].ratioErr);

        int col = pointColors[i % pointColors.size()];
        int marker = 20 + (i % 10);

        if (points[i].label.find("/") != std::string::npos) {
            marker = kOpenSquare;
            col = kBlack;
            gr->SetMarkerSize(1.6);
        }
        else {
            gr->SetMarkerSize(1.2);
        }

        gr->SetLineWidth(2);

        tdrDraw(gr, "Pz", marker, col, kSolid, col);

        graphs.push_back(gr);
    }

    TLegend* leg = tdrLeg(0.18, 0.65, 0.88, 0.75);
    leg->SetNColumns(2);
    leg->SetTextSize(0.025);

    TBox* leg2024 = new TBox();
    leg2024->SetFillColorAlpha(ColorForMCGroup("2024 MC"), 0.28);
    leg2024->SetLineColor(0);

    TBox* leg2025 = new TBox();
    leg2025->SetFillColorAlpha(ColorForMCGroup("2025 MC"), 0.28);
    leg2025->SetLineColor(0);

    TBox* leg2026BD = new TBox();
    leg2026BD->SetFillColorAlpha(ColorForMCGroup("2026BD MC"), 0.28);
    leg2026BD->SetLineColor(0);

    TBox* leg2026C = new TBox();
    leg2026C->SetFillColorAlpha(ColorForMCGroup("2026C MC"), 0.28);
    leg2026C->SetLineColor(0);

    leg->AddEntry(leg2024, "2024 MC", "f");
    leg->AddEntry(leg2025, "2025 MC", "f");
    leg->AddEntry(leg2026BD, " 2026BD MC", "f");
    leg->AddEntry(leg2026C, "2026C MC", "f");

    leg->Draw();

    c->RedrawAxis();

    c->SaveAs("Wmass_DataMC_ratio_vs_era.pdf");

    std::cout << "[OK] Saved Wmass_DataMC_ratio_vs_era.pdf and .png" << std::endl;
}