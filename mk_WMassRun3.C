#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "WMassRun3.h"

#include "TSystem.h"

#include <fstream>
#include <string>

R__LOAD_LIBRARY(WMassRun3_C.so);
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector_cc)

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty_cc)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty_cc)

R__LOAD_LIBRARY(WMassRun3_C)

void mk_WMassRun3(){
  TChain *c = new TChain("Events");
  string filename;
  ifstream fin("input_files/test_files.txt");
  while (fin >> filename) { c->AddFile(filename.c_str()); }
  WMassRun3 s(c);
  s.Loop();
}