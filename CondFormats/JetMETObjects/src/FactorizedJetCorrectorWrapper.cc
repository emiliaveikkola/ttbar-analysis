#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include "TSystemDirectory.h"
#include "TList.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
#include <glob.h>    // for runtime directory listing when no JEC files found

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorWrapper.h"

// Helper trim function
std::string trim(const std::string &s) {
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start)) ++start;
  auto end = s.end();
  do { --end; } while (std::distance(start, end) > 0 && std::isspace(*end));
  return std::string(start, end+1);
}

void FactorizedJetCorrectorWrapper::addJECset(std::string id, int verbosity) {

  // Base directory with JEC text files and fibs.txt
  TString baseDir("CondFormats/JetMETObjects/data/");

  // Build a list of directories to scan.
  std::vector<TString> dirsToScan;
  dirsToScan.push_back(baseDir);

  // Check if a subfolder with the name of the id exists.
  TString subDir = baseDir + TString(id.c_str()) + "/";
  if (!gSystem->AccessPathName(subDir)) {  // exists if AccessPathName returns false
    dirsToScan.push_back(subDir);
    std::cout << "Found subfolder: " << subDir << std::endl;
  }
  /*
  TSystemDirectory dir("dataDir", dataDir);
  TList* files = dir.GetListOfFiles();
  if (!files) {
    std::cerr << "Directory " << dataDir << " not found!" << std::endl;
    return;
  }
  */
  // Map: era (e.g. "2024G_nib1") -> (correction type -> file path)
  std::map<std::string, std::map<std::string, TString>> eraFiles;

  // Loop over all directories.
  for (auto dirPath : dirsToScan) {
    TSystemDirectory dir("dir", dirPath);
    TList* files = dir.GetListOfFiles();
    if (!files) {
      std::cerr << "Directory " << dirPath << " not found!" << std::endl;
      continue;
    }
  
    TIter next(files);
    TObject* obj = nullptr;
    while ((obj = next())) {
      TString fileName = obj->GetName();
      if (fileName == "." || fileName == "..") continue;
      // Only process the desired jet type.
      if (!fileName.EndsWith("_AK4PFPuppi.txt")) continue;
      
      // Tokenize by '_' assuming:
      // prefixTag: e.g. "Prompt24"
      // runTag:   e.g. "Run2024G"
      // nibTag:   e.g. "nib1"
      // versionTag: e.g. "V6M"
      // dataTag:  "DATA" or "MC"
      // corrType: correction type ("L1FastJet", "L2Relative", "L2L3Residual", or "L2Residual")
      TObjArray* tokens = fileName.Tokenize("_");
      TString prefixTag, runTag, nibTag, versionTag, dataTag, corrType;
      if (tokens->GetEntries() < 7) { tokens->Delete(); continue; } // e.g.,
      prefixTag = ((TObjString*)tokens->At(0))->GetString(); // "Prompt24"
      runTag = ((TObjString*)tokens->At(1))->GetString(); // "Run2024G"
      nibTag = ((TObjString*)tokens->At(2))->GetString(); // "nib1"
      versionTag = ((TObjString*)tokens->At(3))->GetString(); // "V9M"
      dataTag = ((TObjString*)tokens->At(4))->GetString(); // "Data"
      corrType = ((TObjString*)tokens->At(5))->GetString(); // "L2Relative"
 
      TString expectedId = Form("%s_%s_%s", prefixTag.Data(), versionTag.Data(), dataTag.Data());
      if (expectedId != id) { tokens->Delete(); continue; }
      //TString tok2 = ((TObjString*)tokens->At(2))->GetString(); // "nib1"
      //TString tok3 = ((TObjString*)tokens->At(3))->GetString(); // "V6M"
      //TString tok4 = ((TObjString*)tokens->At(4))->GetString(); // "DATA","MC"
      //TString tok5 = ((TObjString*)tokens->At(5))->GetString(); // correction type
    
      // Reconstruct the file id (e.g. "Prompt24_V6M_DATA") and check.
      //TString fileId = tok0 + "_" + tok3 + "_" + tok4;
      //if (fileId != id) { tokens->Delete(); continue; }
    
      // Define era key: remove "Run" from tok1 and append "_" + tok2, e.g. "2024G_nib1"
      /*TString era = tok1; 
      era.Remove(0, 3); 
      era += "_" + tok2; 
      //changed this for 2025 without fibs file (below)
      */

      // Define era key: remove "Run" from tok1; ignore token2 (nib) if the year is 2025
      TString runStr = runTag;
      runStr.Remove(0, 3);  // e.g. "2024G" or "2025C"
      TString era;
      era = runStr + "_" + nibTag;
      
      // Save the full file path for this correction.
      //eraFiles[std::string(era.Data())][std::string(tok5.Data())] = baseDir + fileName;
      eraFiles[std::string(era.Data())][std::string(corrType.Data())] = dirPath + fileName;
      tokens->Delete();
    } // end file loop
  } // end dir loop

  if (eraFiles.empty()) {
      // No JEC files matched—report which id and folder, and list all .txt files
      std::cerr << ".ERROR: no JEC files found for id '" << id
                << "' in directory '" << subDir << "'" << std::endl;
      std::cerr << "Available .txt files in " << subDir << ":" << std::endl;
      glob_t globbuf;
      std::string pattern = std::string(subDir.Data()) + "*.txt";
      if (glob(pattern.c_str(), 0, nullptr, &globbuf) == 0) {
          for (size_t i = 0; i < globbuf.gl_pathc; ++i) {
              std::cerr << "  " << globbuf.gl_pathv[i] << std::endl;
          }
          globfree(&globbuf);
      } else {
          std::cerr << "  (failed to list files in directory)" << std::endl;
      }
      exit(1);
  }
  // Determine which fibs.txt file to use.
  TString fibsFileSub = baseDir + TString(id.c_str()) + "/fibs.txt";
  TString fibsFileBase = baseDir + "fibs.txt";
  TString fibsFile;
  if (!gSystem->AccessPathName(fibsFileSub)) {
    fibsFile = fibsFileSub;
    if (verbosity>1)
      std::cout << "Using fibs file from subdirectory: " << fibsFile.Data() << std::endl;
  } else {
    fibsFile = fibsFileBase;
    if (verbosity>1)
      std::cout << "Using fibs file from base directory: " << fibsFile.Data() << std::endl;
  }
  
  // Read fibs.txt to build nib -> run-range map.
  std::map<std::string, std::pair<int,int>> nibRunRanges;
  //TString fibsFile = baseDir + "fibs.txt";
  std::ifstream ifs(fibsFile.Data());
  if (!ifs) {
    std::cerr << "Could not open fibs.txt: " << fibsFile.Data() << std::endl;
    return;
  }
  std::string line;
  // Skip header line
  std::getline(ifs, line);
  while (std::getline(ifs, line)) {
    if (line.empty()) continue;
    std::istringstream iss(line);
    std::string runRangeStr, nameStr;
    if (!std::getline(iss, runRangeStr, '|')) continue;
    if (!std::getline(iss, nameStr, '|')) continue;
    runRangeStr = trim(runRangeStr);
    nameStr = trim(nameStr);
    // Remove surrounding brackets from run range, e.g. "[355374, 355769]"
    if (runRangeStr.front() == '[' && runRangeStr.back() == ']')
      runRangeStr = runRangeStr.substr(1, runRangeStr.size()-2);
    std::istringstream riss(runRangeStr);
    std::string runMinStr, runMaxStr;
    if (!std::getline(riss, runMinStr, ',')) continue;
    if (!std::getline(riss, runMaxStr)) continue;
    int runMin = std::stoi(trim(runMinStr));
    int runMax = std::stoi(trim(runMaxStr));
    
    // Parse nib from name field. Example: "2022B-nib1-fib1" -> nib key "2022B_nib1"
    std::istringstream niss(nameStr);
    std::string part;
    std::vector<std::string> parts;
    while (std::getline(niss, part, '-'))
      parts.push_back(trim(part));
    if (parts.size() < 2) continue;
    std::string nibKey = parts[0] + "_" + parts[1];
    
    // Update nib run-range: take overall min and max among fibs.
    if (nibRunRanges.find(nibKey) == nibRunRanges.end()) {
      nibRunRanges[nibKey] = std::make_pair(runMin, runMax);
    } else {
      auto &range = nibRunRanges[nibKey];
      range.first = std::min(range.first, runMin);
      range.second = std::max(range.second, runMax);
    }
  }
  ifs.close();
  
  // Loop over eras and build JEC chains.
  for (auto &eraPair : eraFiles) {
    std::string eraKey = eraPair.first;  // e.g. "2024G_nib1"
    auto &corrMap = eraPair.second;
    std::vector<JetCorrectorParameters> params;
    
    // Add L1FastJet if available.
    if (corrMap.find("L1FastJet") != corrMap.end()) {
      params.push_back(JetCorrectorParameters(corrMap["L1FastJet"].Data()));
      if (verbosity>1)
	std::cout << " +" << corrMap["L1FastJet"] << std::endl;
    }
    
    // L2Relative is always required.
    if (corrMap.find("L2Relative") != corrMap.end()) {
      params.push_back(JetCorrectorParameters(corrMap["L2Relative"].Data()));
      if (verbosity>1)
	std::cout <<  " +" << corrMap["L2Relative"] << std::endl;
    }
    else {
      std::cerr << "Missing L2Relative correction for era " << eraKey << std::endl;
      //for (auto p : params) delete p;
      continue;
    }
    
    // Determine if this is DATA (id contains "DATA").
    bool isData = TString(id.c_str()).Contains("DATA");
    if (isData) {
      // For DATA, add additional residual: prefer L2L3Residual over L2Residual.
      if (corrMap.find("L2L3Residual") != corrMap.end()) {
        params.push_back(JetCorrectorParameters(corrMap["L2L3Residual"].Data()));
	if (verbosity>1)
	  std::cout <<  " +" << corrMap["L2L3Residual"] << std::endl;
      }
      else if (corrMap.find("L2Residual") != corrMap.end()) {
        params.push_back(JetCorrectorParameters(corrMap["L2Residual"].Data()));
	if (verbosity>1)
	  std::cout <<  " +" << corrMap["L2Residual"] << std::endl;
      }
      else {
        std::cerr << "Missing residual correction (L2L3Residual or L2Residual) for DATA era " << eraKey << std::endl;
        //for (auto p : params) delete p;
        continue;
      }
    }
    
    // Create the FactorizedJetCorrector from the parameters.
    FactorizedJetCorrector* fJEC = new FactorizedJetCorrector(params);
    
    // Find the run-range for this era using the nib key.
    if (nibRunRanges.find(eraKey) == nibRunRanges.end()) {
      std::cerr << "No run range found for era " << eraKey << " in fibs.txt" << std::endl;
      // Optionally delete fJEC and its parameters
      continue;
    }
    int fRunMin = nibRunRanges[eraKey].first;
    int fRunMax = nibRunRanges[eraKey].second;
    mJECmap[std::make_pair(fRunMin, fRunMax)] = fJEC;
    if (verbosity>0)
      std::cout << "Added JEC chain for era " << eraKey 
		<< " with run range (" << fRunMin << ", " << fRunMax << ")" << std::endl;
  }
}
