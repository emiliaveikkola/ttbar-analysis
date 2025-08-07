// Run this script to compile CondFormats libraries. After this can easily run 
// root -l -b -q mk_WMassRun3.C
// using R__LOAD_LIBRARY to load *.so
{
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrectorWrapper.cc+");

  gROOT->ProcessLine(".L WMassRun3.C+g");
}