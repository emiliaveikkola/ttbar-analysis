#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorWrapper.h"

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

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrectorWrapper_cc)

R__LOAD_LIBRARY(WMassRun3_C)
/*
void mk_WMassRun3(){
  gErrorIgnoreLevel = 9999;
  TChain *c = new TChain("Events");
  string filename;
  //ifstream fin("input_files/test_files.txt");
  ifstream fin("input_files/Skimmed_Wqqm_Muon.txt");
  while (fin >> filename) { c->AddFile(filename.c_str()); }
  WMassRun3 s(c);
  s.Loop();
}
*/

#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TChain.h>
#include <TString.h>
#include <iostream>
/*
// Recursive function to add .root files from a directory and its subdirectories.
void AddFilesRecursively(const TString &dirPath, TChain *chain) {
   TSystemDirectory dir(dirPath, dirPath);
   TList *files = dir.GetListOfFiles();
   if (!files) return;
   TIter next(files);
   TSystemFile *file;
   while ((file = (TSystemFile*)next())) {
      TString fname = file->GetName();
      if (fname == "." || fname == "..") continue;
      TString fullPath = dirPath + "/" + fname;
      if (file->IsDirectory()) {
         // Recurse into subdirectory
         AddFilesRecursively(fullPath, chain);
      } else if (fname.EndsWith(".root")) {
         std::cout << "Adding file: " << fullPath.Data() << std::endl;
         chain->AddFile(fullPath);
      }
   }
}
*/

void AddFilesRecursively(const TString &dirPath, TChain *chain) {
   TSystemDirectory dir(dirPath, dirPath);
   TList *files = dir.GetListOfFiles();
   if (!files) return;
   TIter next(files);
   TSystemFile *file;
   while ((file = (TSystemFile*)next())) {
      TString fname = file->GetName();
      if (fname == "." || fname == "..") continue;
      TString fullPath = dirPath + "/" + fname;
      if (file->IsDirectory()) {
         // Always recurse into directories
         AddFilesRecursively(fullPath, chain);
      } else if (fname.EndsWith(".root") && fullPath.Contains("date-24Sep2025_time-114102_commit-53ad610")) { 
         //MC: date-14May2025_time-124123_commit-23dc684
         //CDE DataReprocessing: date-08May2025_time-160559_commit-23dc684
         //FGHI Prompt24: date-20May2025_time-131924_commit-23dc684
         //MCSummer24: date-16Jun2025_time-094328_commit-23dc684
         //2025C: date-04Aug2025_time-123721_commit-23dc684
         //2025D: date-04Aug2025_time-123721_commit-23dc684
         //2025E: date-11Sep2025_time-082745_commit-e998a03
         //2025F: date-24Sep2025_time-114102_commit-53ad610
         //MCWinter25: date-17Sep2025_time-102845_commit-6ad4a3d
         std::cout << "Adding file: " << fullPath.Data() << std::endl;
         chain->AddFile(fullPath);
      }
   }
}

void mk_WMassRun3() {
   gErrorIgnoreLevel = 9999;
   TChain *c = new TChain("Events");

   // Main directory containing subdirectories like 2024A, 2024B, etc.
   //TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2024/DataReprocessing/";
   //TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2024/MC/";
   //TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2024/Data/";
   //TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2024/MCSummer24/";
   //TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2025/Data/2025C/";
   //TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2025/Data/2025D/";
   //TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2025/Data/2025E/";
   TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2025/Data/2025F/";
   //TString mainDirectory = "/eos/user/e/eveikkol/Skim/Wqqm/2025/MCWinter25/";

   // Recursively add all ROOT files found under the main directory.
   AddFilesRecursively(mainDirectory, c);

   // Create instance of your analysis class and run the loop.
   WMassRun3 s(c);
   s.Loop();
}