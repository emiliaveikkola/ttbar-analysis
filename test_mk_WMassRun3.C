#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorWrapper.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"

#include "WMassRun3.h"

#include "TSystem.h"
#include "TChain.h"
#include "TString.h"

#include <cstdio>
#include <iostream>
#include <memory>
#include <string>

R__LOAD_LIBRARY(WMassRun3_C.so);
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrectorWrapper_cc)
R__LOAD_LIBRARY(WMassRun3_C)

void AddFilesFromDAS(const std::string& dataset, TChain* chain, int maxFiles = -1) {
    std::string cmd = "dasgoclient -query='file dataset=" + dataset + "'";
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);

    if (!pipe) {
        std::cerr << "ERROR: failed to run DAS query\n";
        return;
    }

    char buffer[4096];
    int nAdded = 0;
    while (fgets(buffer, sizeof(buffer), pipe.get())) {
        std::string lfn(buffer);

        while (!lfn.empty() && (lfn.back() == '\n' || lfn.back() == '\r' || lfn.back() == ' ' || lfn.back() == '\t')) {
            lfn.pop_back();
        }
        if (lfn.empty()) continue;

        std::string pfn = "root://cms-xrd-global.cern.ch/" + lfn;
        std::cout << "Adding file: " << pfn << std::endl;
        chain->Add(pfn.c_str());

        ++nAdded;
        if (maxFiles > 0 && nAdded >= maxFiles) break;
    }

    std::cout << "Total files added: " << nAdded << std::endl;
}

void mk_WMassRun3() {
    gErrorIgnoreLevel = 9999;
    TChain *c = new TChain("Events");

    AddFilesFromDAS("/Muon0/Run2026B-PromptReco-v1/AOD", c);
    AddFilesFromDAS("/Muon1/Run2026B-PromptReco-v1/AOD", c);
    AddFilesFromDAS("/Muon2/Run2026B-PromptReco-v1/AOD", c);
    AddFilesFromDAS("/Muon3/Run2026B-PromptReco-v1/AOD", c);

    std::cout << "Chain entries = " << c->GetEntries() << std::endl;

    WMassRun3 s(c);
    s.Loop();
}