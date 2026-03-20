#define WMassRun3_cxx
#include "WMassRun3.h"
#include <TH2.h>
#include <TH2D.h>
#include <TH3.h>
#include <TH3D.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TLorentzVector.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <vector>
#include "TMath.h"
#include <random>
#include <algorithm>
#include <TLegend.h>
#include <TColor.h>
#include <TStopwatch.h>
#include "tdrstyle_mod22.C"
#include <set>
#include <unordered_map>
#include <cstdint>
#include "THnSparse.h"

// --- JetMETObjects (standalone JEC/JER) ---
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"


bool _gh_debug = false;


bool WMassRun3::LoadLumiAvgPU(const std::string& lumifile)
{
  PrintInfo(string("Processing LoadLumiAvgPU() with ") + lumifile + "...",true);


  // Check lumi against the list of good runs
  const int a_goodruns[] = {};
  const int ngoodruns = sizeof(a_goodruns)/sizeof(a_goodruns[0]);
  set<int> goodruns;
  if (ngoodruns>0) { // This is an old remnant
    for (int runidx = 0; runidx != ngoodruns; ++runidx)
      goodruns.insert(a_goodruns[runidx]);


    for (auto runit = goodruns.begin(); runit != goodruns.end(); ++runit) cout << *runit << ", ";
    cout << endl;
  }
  set<pair<int, int> > nolums;


  ifstream f(lumifile, ios::in);
  if (!f.is_open()) return false;
  float secLS = 2.3310e+01;
  string s;
  int rn, fill, ls, ifoo;
  float del, rec, avgpu;
  char sfoo[512];
  char shlt[512];
  bool getsuccess1 = static_cast<bool>(getline(f, s, '\n'));
  if (!getsuccess1) return false;
  PrintInfo(string("\nstring: ") + s + " !",true);


  // HOX: the lumi file format has been changing. Change the conditions when needed.
  if (s!="#Data tag : online , Norm tag: None") return false;


  bool getsuccess2 = static_cast<bool>(getline(f, s, '\n'));
  if (!getsuccess2) return false;
  PrintInfo(string("\nstring: ") + s + " !",true);
  if (s!="#run:fill,ls,time,hltpath,delivered(/fb),recorded(/fb),avgpu,source") return false;


  int nls(0);
  double lumsum(0);
  double lumsum_good(0);
  double lumsum_json(0);
  bool skip(false);
  while (getline(f, s, '\n')) {

    // Skip if wrong number of arguments
    if (sscanf(s.c_str(),"%d:%d,%d:%d,%d/%d/%d %d:%d:%d,%[^,],%f,%f,%f,%s",
        &rn,&fill,&ls,&ifoo,&ifoo,&ifoo,&ifoo,&ifoo,&ifoo,&ifoo,shlt,&del,&rec,&avgpu,sfoo)!=15)
      skip=true;


    if (_gh_debug) PrintInfo(Form("Run %d ls %d lumi %f/pb",rn,ls,rec*1000.),true);


    if (skip) { // The user should know if this happens
      if (skip) PrintInfo(string("Skipping line (effects the recorded lumi):\n")+s,true);
      skip = false;
      continue;
    }


    if (_avgpu[rn][ls]!=0) return false;

    // brilcalc returns lumi in units of /fb here; convert to /pb
    double lum = rec*1000.;

    if (lum==0 and goodruns.find(rn)!=goodruns.end() and _json[rn][ls]==1)
      nolums.insert(pair<int, int>(rn,ls));


    _avgpu[rn][ls] = avgpu;

    lumsum += lum;
    if (goodruns.find(rn)!=goodruns.end())
      lumsum_good += lum;
    if (_json[rn][ls])
      lumsum_json += lum;
    ++nls;
    if (nls>100000000) return false;
    if (s.size()>0 && s[0]=='#') continue;
  }

    // If no explicit good-run list is provided, treat all loaded lumi as "good runs"
  if (ngoodruns==0) lumsum_good = lumsum;


  PrintInfo(Form("Called LoadLumiAvgPU() with %s:\nLoaded %lu runs with %d lumi sections containing %f"
                 " pb-1 of data,\n of which %f pb-1 is in good runs (%f%%)\nThis corresponds to %f"
                 " hours of data-taking\nThe JSON file contains %f pb-1 (%f%%)",
                 lumifile.c_str(),_avgpu.size(),nls,lumsum,lumsum_good,
                 (lumsum>0 ? 100.*lumsum_good/lumsum : 0.),nls*secLS/3600,lumsum_json,(lumsum>0 ? 100.*lumsum_json/lumsum : 0.)),true);


  // Report any empty lumi section
  if (nolums.size()!=0) {
    PrintInfo(Form("Warning, found %lu non-normalizable LS:",nolums.size()),true);
    for (auto lumit = nolums.begin(); lumit != nolums.end(); ++lumit) {
      cout << " ["<<lumit->first<<","<<lumit->second;
      auto lumit2 = lumit;
      ++lumit2;
      if (lumit2->first!=lumit->first or lumit2->second!=lumit->second+1) cout << "]";
      else {
        for (int lumadd = 0; lumit2!=nolums.end() and lumit2->first==lumit->first and
                             lumit2->second==lumit->second+lumadd+1; ++lumadd, ++lumit2) {};
        lumit = --lumit2;
        cout << "-" << lumit->second << "]";
      }
    } // for lumit
    cout << endl;
  } // nolums
  return true;
} // LoadLumiAvgPU




// Helper function to retrieve FactorizedJetCorrector 
FactorizedJetCorrector *getFJC(string l1="", string l2="", string res="",
			       string path="") {

  // Set default jet algo                                                       
  if (l1!="" && !(TString(l1.c_str()).Contains("_AK")))
    l1 += "_AK4PFPuppi";
  if (l2!="" && !(TString(l2.c_str()).Contains("_AK")))
    l2 += "_AK4PFPuppi";
  if (res!="" && !(TString(res.c_str()).Contains("_AK")))
    res += "_AK4PFPuppi";

  // Set default path
  if (path=="") path = "CondFormats/JetMETObjects/data";
  const char *cd = path.c_str();
  const char *cl1 = l1.c_str();
  const char *cl2 = l2.c_str();
  const char *cres = res.c_str();
  string s("");

  vector<JetCorrectorParameters> v;
  if (l1!=""){
    s = Form("%s/%s.txt",cd,cl1);
    cout << s << endl << flush;
    JetCorrectorParameters *pl1 = new JetCorrectorParameters(s);
    v.push_back(*pl1);
  }
  if (l2!="") {
    s = Form("%s/%s.txt",cd,cl2);
    cout << s << endl << flush;
    JetCorrectorParameters *pl2 = new JetCorrectorParameters(s);
    v.push_back(*pl2);
  }
  if (res!="") {
    s = Form("%s/%s.txt",cd,cres);
    cout << s << endl << flush;
    JetCorrectorParameters *pres = new JetCorrectorParameters(s);
    v.push_back(*pres);
  }
  FactorizedJetCorrector *jec2024 = new FactorizedJetCorrector(v);

  return jec2024;
} // getJFC 


static bool isUpDownPair(int f1, int f2) {
   // up-type: u(2), c(4);  down-type: d(1), s(3), b(5)
   bool upDown = ( (abs(f1)==2 || abs(f1)==4) && (abs(f2)==1 || abs(f2)==3 || abs(f2)==5) );
   bool downUp = ( (abs(f2)==2 || abs(f2)==4) && (abs(f1)==1 || abs(f1)==3 || abs(f1)==5) );
   return upDown || downUp;
}


// Function to check if a flavor pair is one of the six allowed quark–antiquark combinations
bool qqPairs(int flav1, int flav2) {
    // Work with absolute values
    int f1 = std::abs(flav1);
    int f2 = std::abs(flav2);
    // cs
    if ((f1 == 4 && f2 == 3) || (f1 == 3 && f2 == 4)) return true;
    // ud
    if ((f1 == 2 && f2 == 1) || (f1 == 1 && f2 == 2)) return true;
    // cd
    if ((f1 == 4 && f2 == 1) || (f1 == 1 && f2 == 4)) return true;
    // us
    if ((f1 == 2 && f2 == 3) || (f1 == 3 && f2 == 2)) return true;
    // cb
    if ((f1 == 4 && f2 == 5) || (f1 == 5 && f2 == 4)) return true;
    // ub
    if ((f1 == 2 && f2 == 5) || (f1 == 5 && f2 == 2)) return true;
    return false;
}

// Simple, dependency-free DeltaR for (eta,phi)
static inline double deltaPhi(double p1, double p2) {
   double dphi = std::fmod(p1 - p2, 2.0*TMath::Pi());
   if (dphi >  TMath::Pi()) dphi -= 2.0*TMath::Pi();
   if (dphi < -TMath::Pi()) dphi += 2.0*TMath::Pi();
   return dphi;
}
static inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
   const double dphi = deltaPhi(phi1, phi2);
   const double deta = eta1 - eta2;
   return std::sqrt(deta*deta + dphi*dphi);
}

std::map<UInt_t, Long64_t> runCountAfterSplit;

void WMassRun3::Loop()
{
//   In a ROOT session, you can do:
//      root> .L WMassRun3.C
//      root> WMassRun3 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   // /Muo0/Run2024G-PromptReco-v1/NANOAOD

   // Disable all branches
   fChain->SetBranchStatus("*",0);

   // Enable only the branches you need
   fChain->SetBranchStatus("HLT_IsoMu24",1);
   //fChain->SetBranchStatus("",1);
   fChain->SetBranchStatus("nJet",1);
   fChain->SetBranchStatus("Jet_pt",1);
   fChain->SetBranchStatus("Jet_eta",1);
   fChain->SetBranchStatus("Jet_phi",1);
   fChain->SetBranchStatus("Jet_mass",1);
   // Enable the muon branches you need
   fChain->SetBranchStatus("nMuon", 1);
   fChain->SetBranchStatus("Muon_charge", 1);
   fChain->SetBranchStatus("Muon_pt", 1);
   fChain->SetBranchStatus("Muon_eta", 1);
   fChain->SetBranchStatus("Muon_phi", 1);
   fChain->SetBranchStatus("Muon_mass", 1);
   fChain->SetBranchStatus("Muon_jetIdx", 1);

   fChain->SetBranchStatus("Muon_looseId", 1);
   fChain->SetBranchStatus("Muon_mediumId", 1);
   fChain->SetBranchStatus("Muon_tightId", 1);

   fChain->SetBranchStatus("Muon_pfIsoId", 1);
   fChain->SetBranchStatus("Muon_puppiIsoId", 1);
   fChain->SetBranchStatus("Muon_softId", 1);
   
   fChain->SetBranchStatus("run",1);
   fChain->SetBranchStatus("luminosityBlock",1);
   fChain->SetBranchStatus("event",1);

   fChain->SetBranchStatus("Flag_goodVertices", 1);
   fChain->SetBranchStatus("Flag_globalSuperTightHalo2016Filter", 1);
   fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);
   fChain->SetBranchStatus("Flag_BadPFMuonFilter", 1);
   fChain->SetBranchStatus("Flag_BadPFMuonDzFilter", 1);
   fChain->SetBranchStatus("Flag_hfNoisyHitsFilter", 1);
   fChain->SetBranchStatus("Flag_eeBadScFilter", 1);
   fChain->SetBranchStatus("Flag_ecalBadCalibFilter", 1);

   fChain->SetBranchStatus("RawPuppiMET_pt", 1);
   fChain->SetBranchStatus("RawPuppiMET_phi", 1);

   fChain->SetBranchStatus("Jet_btagUParTAK4B", 1);
   fChain->SetBranchStatus("Jet_muonIdx1", 1);
   fChain->SetBranchStatus("Jet_muonIdx2", 1);
   
   fChain->SetBranchStatus("Jet_chEmEF",1);
   fChain->SetBranchStatus("Jet_chHEF" ,1);
   fChain->SetBranchStatus("Jet_chMultiplicity" ,1);
   fChain->SetBranchStatus("Jet_neMultiplicity" ,1);
   fChain->SetBranchStatus("Jet_muEF"  ,1);
   fChain->SetBranchStatus("Jet_neEmEF",1);
   fChain->SetBranchStatus("Jet_neHEF" ,1);
   fChain->SetBranchStatus("Jet_rawFactor", 1);
   // --- PU & UE branches ---
   fChain->SetBranchStatus("PV_npvs", 1);
   fChain->SetBranchStatus("PV_npvsGood", 1);
   fChain->SetBranchStatus("Rho_fixedGridRhoFastjetAll", 1);

   fChain->SetBranchStatus("bunchCrossing", 1);

   TDirectory *curdir = gDirectory;
   // Create the output file based on a condition
   TFile* fout;

   //fout = new TFile("Winter24_TTtoLNu2Q.root", "RECREATE");
   //fout = new TFile("Muon_Run2024CDE_Reprocessing_Wfinal.root", "RECREATE");
   //fout = new TFile("Muon_Run2024FGHI_Prompt.root", "RECREATE");
   //fout = new TFile("Muon_Run2024CDE_ReReco.root", "RECREATE");
   //fout = new TFile("Muon_Run2024FGHI_ReReco.root", "RECREATE");
   //fout = new TFile("Summer24_TTtoLNu2Q.root", "RECREATE");
   // Year-based configuration
   static const int runYear = 2026; // set to 2024 or 2025
   static const std::string runEra = "B"; // for 2025: "A", "B", etc.
   std::string jsonFile, jecMCSet, jecDataSet, outputFile, jetVetoMap, JERSFvsPt, JERres;

   bool TopPtDependentMass = false;
   bool isMC = false;
   // JER smearing (JER SF)
   bool smearJets = false;
   
   bool useJERSFvsPt = false; // new file format
   int smearNMax = 3;
   std::uint32_t _seed;
   std::mt19937 _mersennetwister;
   if (isMC && smearJets) {
      _seed = 4;
		_mersennetwister = std::mt19937(_seed);
   }

   string jerpath(""), jerpathsf("");

   if (isMC && (runYear == 2024)) {
      if (TopPtDependentMass) {
         jecMCSet   = "RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI"; //Winter24Run3_V1_MC_L2Relative_AK4PUPPI //RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI
         jetVetoMap = "jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root"; //jet_veto_maps/Winter24Prompt24/Winter24Prompt24_2024BCDEFGHI.root //jetvetoReReco2024_V9M.root
         // JER Scale Factors (pt-dependent via FactorizedJetCorrector)
         JERSFvsPt = "ReReco24_V10M_MC/ReReco24_2024_nib_JRV10M_MC_SF_AK4PFPuppi"; //"ReReco24_2024_nib_JRV10M_MC_SF_AK4PFPuppi";
         // JER resolution (Pt Resolution txt)
         JERres   = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt";              
         outputFile = "Summer24_TTtoLNu2Q_TopPtDependentMass.root";
      } else {
         jecMCSet   = "RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI"; //Winter24Run3_V1_MC_L2Relative_AK4PUPPI //RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI
         jetVetoMap = "jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root"; //jet_veto_maps/Winter24Prompt24/Winter24Prompt24_2024BCDEFGHI.root //jetvetoReReco2024_V9M.root
         JERSFvsPt = "ReReco24_V10M_MC/ReReco24_2024_nib_JRV10M_MC_SF_AK4PFPuppi"; //"ReReco24_2024_nib_JRV10M_MC_SF_AK4PFPuppi";
         // JER resolution (Pt Resolution txt)
         JERres   = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt";              
         outputFile = "Summer24_TTtoLNu2Q_V9M_FSR.root";
      }
   }
   else if (isMC && (runYear == 2025)) {
      if (TopPtDependentMass) {
         jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI"; //Winter24Run3_V1_MC_L2Relative_AK4PUPPI //RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI
         jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //jet_veto_maps/Winter24Prompt24/Winter24Prompt24_2024BCDEFGHI.root //jetvetoReReco2024_V9M.root
         // JER Scale Factors (pt-dependent via FactorizedJetCorrector)
         JERSFvsPt = "Prompt25_V2M_MC/Prompt25_2025CDE_JRV2M_MC_SF_AK4PFPuppi"; //"ReReco24_2024_nib_JRV10M_MC_SF_AK4PFPuppi";
         // JER resolution (Pt Resolution txt)
         JERres   = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt";              
         outputFile = "Winter25_TTtoLNu2Q_TopPtDependentMass.root";
      } else {
         jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI"; //Winter24Run3_V1_MC_L2Relative_AK4PUPPI //RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI
         jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //jet_veto_maps/Winter24Prompt24/Winter24Prompt24_2024BCDEFGHI.root //jetvetoReReco2024_V9M.root
         // JER Scale Factors (pt-dependent via FactorizedJetCorrector)
         JERSFvsPt = "Prompt25_V3M_MC/Prompt25_2025CDEFG_JRV3M_MC_SF_AK4PFPuppi.txt"; //"ReReco24_2024_nib_JRV10M_MC_SF_AK4PFPuppi";
         // JER resolution (Pt Resolution txt)
         JERres   = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt";
         outputFile = "Winter25_TTtoLNu2Q_V3M.root"; //"Winter25_TTtoLNu2Q.root";
      }
   }
   else if (runYear == 2024) {
       jsonFile   = "Cert_Collisions2024_378981_386951_Golden.json";
       jecMCSet   = "RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI"; //Winter24Run3_V1_MC_L2Relative_AK4PUPPI //RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI
       jecDataSet = "ReReco24_V9M_DATA"; //Reprocessing24_V8M_DATA //Prompt24_V8M_DATA //ReReco24_V9M_DATA
       jetVetoMap = "jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root"; //jet_veto_maps/Winter24Prompt24/Winter24Prompt24_2024BCDEFGHI.root //jetvetoReReco2024_V9M.root
       if (runEra == "CDE") {
         outputFile = "Muon_Run2024CDE_ReReco_V9M.root";
       } else if (runEra == "FGHI"){
         outputFile = "Muon_Run2024FGHI_Prompt_V9M.root";
       }
   } else if (runYear == 2025) {
       if (runEra == "C") {
           jsonFile   = "Cert_Collisions2025_391658_398860_Golden.json";
           jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI";
           jecDataSet = "Prompt25_V3M_DATA/Prompt25_Run2025C_V3M_DATA_L2L3Residual_AK4PFPuppi"; //"Prompt25_V1M_DATA";
           jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //"jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root";               
           outputFile = "Muon_Run2025C_Prompt_V3M_test.root";
       } else if (runEra == "C_TrkRadDamage") {
           jsonFile   = "Cert_Collisions2025_391658_398860_Golden.json";
           jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI";
           jecDataSet = "Prompt25_V3M_DATA/Prompt25_Run2025C_V3M_DATA_L2L3Residual_AK4PFPuppi"; //"Prompt25_V1M_DATA";
           jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //"jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root";               
           outputFile = "Muon_Run2025C_TrkRadDamage_Prompt_V3M.root";
       } else if (runEra == "D") {
           jsonFile   = "Cert_Collisions2025_391658_398860_Golden.json"; //"Cert_Collisions2025D_daily_dials_12-08-2025.json"; //"Collisions25_13p6TeV_391658_395372_DCSOnly_TkPx.json";
           jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI";
           jecDataSet = "Prompt25_V3M_DATA/Prompt25_Run2025D_V3M_DATA_L2L3Residual_AK4PFPuppi"; //"Prompt25_V1M_DATA";
           jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //"jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root";               
           outputFile = "Muon_Run2025D_Prompt_V3M.root";
       } else if (runEra == "E") {
           jsonFile   = "Cert_Collisions2025_391658_398860_Golden.json"; //"Cert_Collisions2025D_daily_dials_12-08-2025.json"; //"Collisions25_13p6TeV_391658_395372_DCSOnly_TkPx.json";
           jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI";
           jecDataSet = "Prompt25_V3M_DATA/Prompt25_Run2025E_V3M_DATA_L2L3Residual_AK4PFPuppi"; //"Prompt25_V1M_DATA";
           jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //"jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root";               
           outputFile = "Muon_Run2025E_Prompt_V3M.root";
       } else if (runEra == "F") {
           jsonFile   = "Cert_Collisions2025_391658_398860_Golden.json"; //"Cert_Collisions2025D_daily_dials_12-08-2025.json"; //"Collisions25_13p6TeV_391658_395372_DCSOnly_TkPx.json";
           jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI";
           jecDataSet = "Prompt25_V3M_DATA/Prompt25_Run2025F_V3M_DATA_L2L3Residual_AK4PFPuppi"; //"Prompt25_V1M_DATA";
           jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //"jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root";               
           outputFile = "Muon_Run2025F_Prompt_V3M_test.root";
      } else if (runEra == "G") {
           jsonFile   = "Cert_Collisions2025_391658_398860_Golden.json"; //"Cert_Collisions2025D_daily_dials_12-08-2025.json"; //"Collisions25_13p6TeV_391658_395372_DCSOnly_TkPx.json";
           jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI";
           jecDataSet = "Prompt25_V3M_DATA/Prompt25_Run2025G_V3M_DATA_L2L3Residual_AK4PFPuppi"; //"Prompt25_V1M_DATA";
           jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //"jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root";               
           outputFile = "Muon_Run2025G_Prompt_V3M.root";
      } else if (runEra == "CDEFG") {
           jsonFile   = "Cert_Collisions2025_391658_398860_Golden.json"; //"Cert_Collisions2025D_daily_dials_12-08-2025.json"; //"Collisions25_13p6TeV_391658_395372_DCSOnly_TkPx.json";
           jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI";
           jecDataSet = "Prompt25_V3M_DATA/Prompt25_Run2025CDEFG_V3M_DATA_L2L3Residual_AK4PFPuppi"; //"Prompt25_V1M_DATA";
           jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //"jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root";               
           outputFile = "Muon_Run2025CDEFG_Prompt_V3M.root";     
       } else {
           std::cerr << "Unsupported runEra: " << runEra << std::endl;
           return;
       }
   } else if (runYear == 2026) {
         jsonFile   = "CombinedJSONS_DCSCMLE_Runs_401624to402040.json";//"Collisions26_13p6TeV_401623_401961_DCSOnly_TkPx.json"; //"Cert_Collisions2025D_daily_dials_12-08-2025.json"; //"Collisions25_13p6TeV_391658_395372_DCSOnly_TkPx.json";
         //jecDataSet = "Prompt25_V3M_DATA/Prompt25_Run2025G_V3M_DATA_L2L3Residual_AK4PFPuppi"; //"Prompt25_V1M_DATA";
         jetVetoMap = "jet_veto_maps/jetveto2025CDEFG_V3M.root"; //"jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root";
         if (runEra == "A") {
         jecMCSet   = "Winter25Run3_V1_MC_L2Relative_AK4PUPPI";
         outputFile = "Muon_Run2026A_Prompt_L2L32025G.root";
      } else if (runEra == "Bnib1") {
         jecMCSet   = "Run3Winter26_PhiDependent_L2Relative_AK4PUPPI";
         outputFile = "Muon_Run2026Bnib1_Prompt_CombinedJSON.root";
      } else if (runEra == "Bnib2") {
         jecMCSet   = "Run3Winter26_PhiDependent_L2Relative_AK4PUPPI";
        outputFile = "Muon_Run2026Bnib2_Prompt_CombinedJSON.root";               
      } else if (runEra == "B") {
         jecMCSet   = "Run3Winter26_PhiDependent_L2Relative_AK4PUPPI";
        outputFile = "Muon_Run2026B_Prompt_CombinedJSON.root";      
      }
   } else {
       std::cerr << "Unsupported runYear: " << runYear << std::endl;
       return;
   }
   fout = new TFile(outputFile.c_str(), "RECREATE");
   fout->cd();
   // Separate output for run vs BX map (filled once, no cuts)
   TFile* frunbx = nullptr;

   if (!isMC) {
      // Separate output for run vs BX map (filled once, no cuts)
      frunbx = new TFile("runbx_map.root", "RECREATE");
      fout->cd();
   }

   TFile* fmap = nullptr;
   if (isMC) {
      fmap = new TFile("GenJetPartonPtMap_updated.root", "RECREATE");
      fout->cd();
   }

   if (isMC) {
      fChain->SetBranchStatus("Jet_hadronFlavour", 1);
      fChain->SetBranchStatus("Jet_partonFlavour", 1);
      fChain->SetBranchStatus("genWeight", 1);
      fChain->SetBranchStatus("GenVtx_z", 1);
		fChain->SetBranchStatus("PV_z", 1);

      fChain->SetBranchStatus("Pileup_nTrueInt", 1); // μ (true interactions) in MC
      fChain->SetBranchStatus("GenPart_pt", 1);
      fChain->SetBranchStatus("GenPart_eta", 1);
      fChain->SetBranchStatus("GenPart_phi", 1);
      fChain->SetBranchStatus("GenPart_mass", 1);
      fChain->SetBranchStatus("nGenPart", 1);
      fChain->SetBranchStatus("GenPart_pdgId", 1);
      fChain->SetBranchStatus("GenPart_statusFlags", 1);
      fChain->SetBranchStatus("GenPart_status", 1);
      fChain->SetBranchStatus("GenPart_genPartIdxMother", 1);

      fChain->SetBranchStatus("Jet_genJetIdx", 1);
      fChain->SetBranchStatus("nGenJet", 1);
      fChain->SetBranchStatus("GenJet_pt", 1);
      fChain->SetBranchStatus("GenJet_eta", 1);
      fChain->SetBranchStatus("GenJet_phi", 1);
      fChain->SetBranchStatus("GenJet_mass", 1);
      fChain->SetBranchStatus("GenJet_partonFlavour", 1);
      fChain->SetBranchStatus("Rho_fixedGridRhoFastjetAll", 1);

      // PS weights (w_var / w_nominal): [0] isr.murfac=2.0; [1] fsr.murfac=2.0; [2] isr.murfac=0.5; [3] fsr.murfac=0.5
      fChain->SetBranchStatus("nPSWeight", 1);
      fChain->SetBranchStatus("PSWeight", 1);
      
   }

   TH1::SetDefaultSumw2();

   double xmax = 450000.5; //389000.5
   double xmin = 355000.5;
   double histnx = xmax-xmin;

   double vx[] = {15, 20, 25, 30, 35, 40, 50, 60, 75, 90, 110, 130, 175, 230,
         300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750,
         2100, 2500, 3000};
   const int nx = sizeof(vx)/sizeof(vx[0])-1;

   

      // --- Control plots (explicit only W and Top masses) ---
      // --- JER control: pre-JER vs post-JER (smeared) ---
      TH1D* h_ctrl_Wmass_preJER    = new TH1D("h_ctrl_Wmass_preJER",    "; W Mass pre-JER (GeV); Events", 180, 0, 900);
      TH1D* h_ctrl_Topmass_preJER  = new TH1D("h_ctrl_Topmass_preJER",  "; Top Mass pre-JER (GeV); Events", 180, 0, 900);

      TH1D* h_Wmass = new TH1D("h_Wmass", "; W Mass (GeV); Events", 180, 0, 900);
      TH1D* h_Wpt = new TH1D("h_Wpt", "; Pt (GeV); Events", 80, 0, 400);
      TH1D* h_Weta = new TH1D("h_Weta", "; W #eta; Events", 50, -5, 5);

      TH2D* h_Wmass_pt = new TH2D("h_Wmass_pt", "; W Pt (GeV) ;W Mass (GeV); Events", 80, 0, 400, 180, 0, 900);
      TProfile *p_Wmass_vs_Wpt = new TProfile("p_Wmass_vs_Wpt", ";W p_{T} (GeV);Mean W mass (GeV)", 80, 0, 400);

      TH1D* h_singleJetWmass = new TH1D("h_singleJetWmass", "; W Mass (GeV); Events", 180, 0, 900);
      TH1D* h_singleJetWpt = new TH1D("h_singleJetWpt", "; Pt (GeV); Events", 80, 0, 400);
      TProfile *p_singleJetWmass_vs_Wpt = new TProfile("p_singleJetWmass_vs_Wpt", ";W p_{T} (GeV);Mean W mass (GeV)", 80, 0, 400);

      TH1D* h_Topmass = new TH1D("h_Topmass", "; Top Mass (GeV); Events", 180, 0, 900);
      TH1D* h_Toppt = new TH1D("h_Toppt", "; Pt (GeV); Events", 80, 0, 400);
      TH2D* h_Topmass_pt = new TH2D("h_Topmass_pt", "; Top Pt (GeV) ;Top Mass (GeV); Events", 80, 0, 400, 180, 0, 900);
      TProfile *p_Topmass_vs_Toppt = new TProfile("p_Topmass_vs_Toppt", ";Top p_{T} (GeV);Mean Top mass (GeV)", 80, 0, 400);

      TProfile *p_Topmass_vs_TopptTrue = new TProfile("p_Topmass_vs_TopptTrue", ";p_{T}^{true}(t_{had}) (GeV);#LT m_{qqb} #GT (GeV)", 200, 0, 1000);
      TProfile *p_Topmass_vs_Toppt_hadronic = new TProfile("p_Topmass_vs_Toppt_hadronic", ";p_{T}^{had}(t_{had}) (GeV);#LT m_{qqb} #GT (GeV)", 200, 0, 1000);

      TProfile *prof_top = new TProfile("prof_top",";Run;Top mass;",histnx,xmin,xmax);
      TProfile *prof_top_inWindow = new TProfile("prof_top_inWindow",";Run;Top mass;",histnx,xmin,xmax);
      TProfile *prof_top_improved = new TProfile("prof_top_improved",";Run;Top mass (improved);",histnx,xmin,xmax);
      TProfile *prof_W = new TProfile("prof_W",";Run;W mass;",histnx,xmin,xmax);
      TProfile *prof_W_inWindow = new TProfile("prof_W_inWindow",";Run;W mass;",histnx,xmin,xmax);

      TProfile *prof_L2L3Res = new TProfile("prof_L2L3Res", ";Run;L2L3Res factor;", histnx, xmin, xmax);
      TProfile *prof_L2L3 = new TProfile("prof_L2L3", ";Run;L2L3 factor;", histnx, xmin, xmax);
      TProfile *prof_L1 = new TProfile("prof_L1", ";Run;L1 factor;", histnx, xmin, xmax);
      TProfile *prof_corr = new TProfile("prof_corr", ";Run;corr;", histnx, xmin, xmax);

      TProfile* pres = new TProfile("pres", ";pt_pair;1/(L2L3Res);", nx, vx);
      TProfile* pjes = new TProfile("pjes", ";pt_pair;1/(L2);",       nx, vx);

      TProfile *prof_L2L3Res_ptpair = new TProfile("prof_L2L3Res_ptpair", ";pt_pair;L2L3Res factor;", nx, vx);
      TProfile *prof_L2L3_ptpair = new TProfile("prof_L2L3_ptpair", ";pt_pair;L2L3 factor;", nx, vx);
      TProfile *prof_L1_ptpair = new TProfile("prof_L1_ptpair", ";pt_pair;L1 factor;", nx, vx);
      TProfile *prof_corr_ptpair = new TProfile("prof_corr_ptpair", ";pt_pair;corr;", nx, vx);

      //TH1D *h1_test = new TH1D("h1_test", ";pt;N", 100,0,400);
      TH1D* h_Wmass_inWindow = new TH1D("h_Wmass_inWindow", "; W Mass in Window (GeV); Events", 180, 0, 900);
      TH2D* h2_Wjets = new TH2D("h2_Wjets", "; flav1;flav2;", 12,-5.5,6.5,12,-5.5,6.5);

      TH1D* h_Jet1Pt = new TH1D("h_Jet1Pt", "; Pt (GeV); Events", 80, 0, 400);
      TH1D* h_Jet2Pt = new TH1D("h_Jet2Pt", "; Pt (GeV); Events", 80, 0, 400);

      TH1D* h_Top1mass = new TH1D("h_Top1mass", "; Top Mass b1 (GeV); Events", 180, 0, 900);
      TH1D* h_Top2mass = new TH1D("h_Top2mass", "; Top Mass b2 (GeV); Events", 180, 0, 900);
      TH1D* h_Topmass_improved = new TH1D("h_Topmass_improved", "; Top Mass in Window (GeV); Events", 180, 0, 900);
      TH1D* h_Topmass_inWindow = new TH1D("h_Topmass_inWindow", "; Top Mass in Window (GeV); Events", 180, 0, 900);
      TH1D* h_Topmass_inWindow_improved = new TH1D("h_Topmass_inWindow_improved", "; Top Mass in Window (GeV); Events", 180, 0, 900);
      
      TH2D* h2_Wmass_ptpair = new TH2D("h2_Wmass_ptpair", "; ptpair; W Mass (GeV)", nx, vx, 180, 0, 900);
      TH2D* h2_Wmass_inWindow_ptpair = new TH2D("h2_Wmass_inWindow_ptpair", "; ptpair; W Mass (GeV)", nx, vx, 180, 0, 900);
      TH2D* h2_Wmass_ptpair_improved = new TH2D("h2_Wmass_ptpair_improved", "; ptpair; W Mass (GeV)", nx, vx, 180, 0, 900);
      TH2D* h2_Wmass_inWindow_ptpair_improved = new TH2D("h2_Wmass_inWindow_ptpair_improved", "; ptpair; W Mass (GeV)", nx, vx, 180, 0, 900);

      // Diagonal‐window histograms for Mass ≃ 2*ptpair + 4
      TH1D*  h_Wmass_diagWindow               = new TH1D("h_Wmass_diagWindow",               "; W Mass (GeV) diag window; Events", 180, 0, 900);
      TH1D*  h_ptpair_diagWindow              = new TH1D("h_ptpair_diagWindow",              "; ptpair diag window; Events",          nx, vx);
      TH2D*  h2_Wmass_ptpair_diagWindow       = new TH2D("h2_Wmass_ptpair_diagWindow",       "; ptpair; W Mass (GeV); Events",        nx, vx, 180, 0, 900);
      TProfile* prof_W_diagWindow_ptpair      = new TProfile("prof_W_diagWindow_ptpair",      ";pt_pair;W mass;",                       nx, vx);
      TProfile* prof_ptpair_diagWindow_mass   = new TProfile("prof_ptpair_diagWindow_mass",   ";W mass; Mean ptpair (GeV)",             180, 0, 900);

      TProfile *prof_W_ptpair = new TProfile("prof_W_ptpair",";pt_pair;W mass;", nx, vx);
      TProfile *prof_W_ptpair_improved = new TProfile("prof_W_ptpair_improved",";pt_pair;W mass;", nx, vx);
      TProfile *prof_W_inWindow_ptpair = new TProfile("prof_W_inWindow_ptpair",";pt_pair;W mass;", nx, vx);
      TProfile *prof_W_inWindow_ptpair_improved = new TProfile("prof_W_inWindow_ptpair_improved",";pt_pair;W mass;", nx, vx);

      // --- PU plots for UE density response (HDM unclustered-energy correction) ---
      // x-axis: μ (MC: Pileup_nTrueInt; Data: PV_npvsGood as proxy)
      TProfile* pnpvvsmu        = new TProfile("pnpvvsmu",        ";#mu; NPV (all vertices)",    200, 0, 200);
      TProfile* pnpvGoodvsmu    = new TProfile("pnpvGoodvsmu",    ";#mu; NPV^{good}",            200, 0, 200);
      TProfile* prhovsmu        = new TProfile("prhovsmu",        ";#mu; #rho (GeV)",            200, 0, 200);
      // Self-consistency / bookkeeping: <mu> vs mu (weighted)
      TProfile* pmuvsmu         = new TProfile("pmuvsmu",         ";#mu; #LT#mu#GT",              200, 0, 200);

      TH1D* h_muPU = new TH1D("h_muPU", ";#mu; Events", 200, 0, 200);

      // --- PU plots for UE density response (HDM unclustered-energy correction) ---
      // x-axis: μ (MC: Pileup_nTrueInt; Data: PV_npvsGood as proxy)
      TProfile* pnpvvsmu_after        = new TProfile("pnpvvsmu_after",        ";#mu; NPV (all vertices)",    200, 0, 200);
      TProfile* pnpvGoodvsmu_after    = new TProfile("pnpvGoodvsmu_after",    ";#mu; NPV^{good}",            200, 0, 200);
      TProfile* prhovsmu_after        = new TProfile("prhovsmu_after",        ";#mu; #rho (GeV)",            200, 0, 200);
      // Self-consistency / bookkeeping: <mu> vs mu (weighted)
      TProfile* pmuvsmu_after         = new TProfile("pmuvsmu_after",         ";#mu; #LT#mu#GT",              200, 0, 200);

      TH1D* h_muPU_after = new TH1D("h_muPU_after", ";#mu; Events", 200, 0, 200);

      // --- region classified W mass profiles ---
      TH2D *h2_W_ptpair_sig  = new TH2D("h2_W_ptpair_sig",  ";pt_pair;W mass (signal)",       nx, vx, 180, 0, 900);
      TH2D *h2_W_ptpair_bkg = new TH2D("h2_W_ptpair_bkg", ";pt_pair;W mass (background1)",  nx, vx, 180, 0, 900);
      TH2D *h2_W_ptpair_cr   = new TH2D("h2_W_ptpair_cr",   ";pt_pair;W mass (control)",      nx, vx, 180, 0, 900);

      TH1D* h_ptpair = new TH1D("h_ptpair",";#sqrt{p_{T,1}p_{T,2}};Events", nx, vx);
      TH1D* h_ptpair_improved = new TH1D("h_ptpair_improved",";#sqrt{p_{T,1}p_{T,2}};Events", nx, vx);

      TH1D* h_ptpair_inWindow = new TH1D("h_ptpair_inWindow",";#sqrt{p_{T,1}p_{T,2}};Events", nx, vx);
      TH1D* h_ptpair_inWindow_improved = new TH1D("h_ptpair_inWindow_improved",";#sqrt{p_{T,1}p_{T,2}};Events", nx, vx);
      
      // --- begin gg/qg/qq ptpair histograms ---
      // gg category
      TH2D*    h2_Wmass_ptpair_gg                   = new TH2D("h2_Wmass_ptpair_gg",                   "; ptpair (gg); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_ptpair_gg                    = new TProfile("prof_W_ptpair_gg",                    ";pt_pair (gg);W mass;", nx, vx);
      TH1D*    h_ptpair_gg                          = new TH1D("h_ptpair_gg",                          ";#sqrt{p_{T,1}p_{T,2}} (gg);Events", nx, vx);
      TH2D*    h2_Wmass_ptpair_improved_gg          = new TH2D("h2_Wmass_ptpair_improved_gg",          "; ptpair (gg); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_ptpair_improved_gg           = new TProfile("prof_W_ptpair_improved_gg",           ";pt_pair (gg);W mass;", nx, vx);
      TH1D*    h_ptpair_improved_gg                 = new TH1D("h_ptpair_improved_gg",                 ";#sqrt{p_{T,1}p_{T,2}} (gg);Events", nx, vx);
      TH2D*    h2_Wmass_inWindow_ptpair_gg          = new TH2D("h2_Wmass_inWindow_ptpair_gg",          "; ptpair (gg); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_inWindow_ptpair_gg           = new TProfile("prof_W_inWindow_ptpair_gg",           ";pt_pair (gg);W mass;", nx, vx);
      TH1D*    h_ptpair_inWindow_gg                 = new TH1D("h_ptpair_inWindow_gg",                 ";#sqrt{p_{T,1}p_{T,2}} (gg);Events", nx, vx);
      TH2D*    h2_Wmass_inWindow_ptpair_improved_gg = new TH2D("h2_Wmass_inWindow_ptpair_improved_gg", "; ptpair (gg); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_inWindow_ptpair_improved_gg  = new TProfile("prof_W_inWindow_ptpair_improved_gg",  ";pt_pair (gg);W mass;", nx, vx);
      TH1D*    h_ptpair_inWindow_improved_gg        = new TH1D("h_ptpair_inWindow_improved_gg",        ";#sqrt{p_{T,1}p_{T,2}} (gg);Events", nx, vx);

      // qg category
      TH2D*    h2_Wmass_ptpair_qg                   = new TH2D("h2_Wmass_ptpair_qg",                   "; ptpair (qg); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_ptpair_qg                    = new TProfile("prof_W_ptpair_qg",                    ";pt_pair (qg);W mass;", nx, vx);
      TH1D*    h_ptpair_qg                          = new TH1D("h_ptpair_qg",                          ";#sqrt{p_{T,1}p_{T,2}} (qg);Events", nx, vx);
      TH2D*    h2_Wmass_ptpair_improved_qg          = new TH2D("h2_Wmass_ptpair_improved_qg",          "; ptpair (qg); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_ptpair_improved_qg           = new TProfile("prof_W_ptpair_improved_qg",           ";pt_pair (qg);W mass;", nx, vx);
      TH1D*    h_ptpair_improved_qg                 = new TH1D("h_ptpair_improved_qg",                 ";#sqrt{p_{T,1}p_{T,2}} (qg);Events", nx, vx);
      TH2D*    h2_Wmass_inWindow_ptpair_qg          = new TH2D("h2_Wmass_inWindow_ptpair_qg",          "; ptpair (qg); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_inWindow_ptpair_qg           = new TProfile("prof_W_inWindow_ptpair_qg",           ";pt_pair (qg);W mass;", nx, vx);
      TH1D*    h_ptpair_inWindow_qg                 = new TH1D("h_ptpair_inWindow_qg",                 ";#sqrt{p_{T,1}p_{T,2}} (qg);Events", nx, vx);
      TH2D*    h2_Wmass_inWindow_ptpair_improved_qg = new TH2D("h2_Wmass_inWindow_ptpair_improved_qg", "; ptpair (qg); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_inWindow_ptpair_improved_qg  = new TProfile("prof_W_inWindow_ptpair_improved_qg",  ";pt_pair (qg);W mass;", nx, vx);
      TH1D*    h_ptpair_inWindow_improved_qg        = new TH1D("h_ptpair_inWindow_improved_qg",        ";#sqrt{p_{T,1}p_{T,2}} (qg);Events", nx, vx);

      // qq category
      TH2D*    h2_Wmass_ptpair_qq                   = new TH2D("h2_Wmass_ptpair_qq",                   "; ptpair (qq); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_ptpair_qq                    = new TProfile("prof_W_ptpair_qq",                    ";pt_pair (qq);W mass;", nx, vx);
      TH1D*    h_ptpair_qq                          = new TH1D("h_ptpair_qq",                          ";#sqrt{p_{T,1}p_{T,2}} (qq);Events", nx, vx);
      TH2D*    h2_Wmass_ptpair_improved_qq          = new TH2D("h2_Wmass_ptpair_improved_qq",          "; ptpair (qq); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_ptpair_improved_qq           = new TProfile("prof_W_ptpair_improved_qq",           ";pt_pair (qq);W mass;", nx, vx);
      TH1D*    h_ptpair_improved_qq                 = new TH1D("h_ptpair_improved_qq",                 ";#sqrt{p_{T,1}p_{T,2}} (qq);Events", nx, vx);
      TH2D*    h2_Wmass_inWindow_ptpair_qq          = new TH2D("h2_Wmass_inWindow_ptpair_qq",          "; ptpair (qq); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_inWindow_ptpair_qq           = new TProfile("prof_W_inWindow_ptpair_qq",           ";pt_pair (qq);W mass;", nx, vx);
      TH1D*    h_ptpair_inWindow_qq                 = new TH1D("h_ptpair_inWindow_qq",                 ";#sqrt{p_{T,1}p_{T,2}} (qq);Events", nx, vx);
      TH2D*    h2_Wmass_inWindow_ptpair_improved_qq = new TH2D("h2_Wmass_inWindow_ptpair_improved_qq", "; ptpair (qq); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_inWindow_ptpair_improved_qq  = new TProfile("prof_W_inWindow_ptpair_improved_qq",  ";pt_pair (qq);W mass;", nx, vx);
      TH1D*    h_ptpair_inWindow_improved_qq        = new TH1D("h_ptpair_inWindow_improved_qq",        ";#sqrt{p_{T,1}p_{T,2}} (qq);Events", nx, vx);

      // qq_others category
      TH2D*    h2_Wmass_ptpair_qq_others                   = new TH2D("h2_Wmass_ptpair_qq_others",                   "; ptpair (qq_others); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_ptpair_qq_others                    = new TProfile("prof_W_ptpair_qq_others",                    ";pt_pair (qq_others);W mass;", nx, vx);
      TH1D*    h_ptpair_qq_others                          = new TH1D("h_ptpair_qq_others",                          ";#sqrt{p_{T,1}p_{T,2}} (qq_others);Events", nx, vx);
      TH2D*    h2_Wmass_ptpair_improved_qq_others          = new TH2D("h2_Wmass_ptpair_improved_qq_others",          "; ptpair (qq_others); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_ptpair_improved_qq_others           = new TProfile("prof_W_ptpair_improved_qq_others",           ";pt_pair (qq_others);W mass;", nx, vx);
      TH1D*    h_ptpair_improved_qq_others                 = new TH1D("h_ptpair_improved_qq_others",                 ";#sqrt{p_{T,1}p_{T,2}} (qq_others);Events", nx, vx);
      TH2D*    h2_Wmass_inWindow_ptpair_qq_others          = new TH2D("h2_Wmass_inWindow_ptpair_qq_others",          "; ptpair (qq_others); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_inWindow_ptpair_qq_others           = new TProfile("prof_W_inWindow_ptpair_qq_others",           ";pt_pair (qq_others);W mass;", nx, vx);
      TH1D*    h_ptpair_inWindow_qq_others                 = new TH1D("h_ptpair_inWindow_qq_others",                 ";#sqrt{p_{T,1}p_{T,2}} (qq_others);Events", nx, vx);
      TH2D*    h2_Wmass_inWindow_ptpair_improved_qq_others = new TH2D("h2_Wmass_inWindow_ptpair_improved_qq_others", "; ptpair (qq_others); W Mass (GeV)", nx, vx, 180, 0, 900);
      TProfile* prof_W_inWindow_ptpair_improved_qq_others  = new TProfile("prof_W_inWindow_ptpair_improved_qq_others",  ";pt_pair (qq_others);W mass;", nx, vx);
      TH1D*    h_ptpair_inWindow_improved_qq_others        = new TH1D("h_ptpair_inWindow_improved_qq_others",        ";#sqrt{p_{T,1}p_{T,2}} (qq_others);Events", nx, vx);
      // --- end gg/qg/qq ptpair histograms ---


      TH1D* h_mbl1 = new TH1D("h_mbl1", ";; Events", 80, 0, 400);
      TH1D* h_mbl2 = new TH1D("h_mbl2", ";; Events", 80, 0, 400);
      TH1D* h_mbl = new TH1D("h_mbl", ";; Events", 80, 0, 400);
      TH1D* h_mbl_red = new TH1D("h_mbl_red", ";; Events", 100, 0, 5);
      TH1D* h_mbl_red_inWindow = new TH1D("h_mbl_red_inWindow", ";; Events", 100, 0, 5);
      TH1D* h_mbl_inWindow = new TH1D("h_mbl_inWindow", ";; Events", 80, 0, 400);
      TH1D* h_rbq = new TH1D("h_rbq", ";; Events", 100, 0, 5);
      TH1D* h_rbq_inWindow = new TH1D("h_rbq_inWindow", ";; Events", 100, 0, 5);
      TH1D* h_rbl = new TH1D("h_rbl", ";; Events", 100, 0, 5);
      TH1D* h_rbl_inWindow = new TH1D("h_rbl_inWindow", ";; Events", 100, 0, 5);

      // Sum of generator weights (for later normalization)
      TH1D* h_sumw = new TH1D("h_sumw", ";bin;sumw", 1, 0, 1);


      // Sparse (run, BX) event-count map
      const Int_t runbx_ndim = 2;
      Int_t    runbx_nbins[runbx_ndim] = {55000, 3600};
      Double_t runbx_xmin[runbx_ndim]  = {355000., 0.};
      Double_t runbx_xmax[runbx_ndim]  = {410000., 3600.};

      THnSparseF *h2runbx = new THnSparseF(
         "h2runbx", ";Run number;Bunch crossing",
         runbx_ndim, runbx_nbins, runbx_xmin, runbx_xmax
      );

     // --- JER control plots (same binning and names as DijetHistosFill.C) ---
   // pT,avp (HDM) binning (GeV)
   double vptd[] = {15, 20, 25, 30, 35, 40, 50, 60, 75, 90, 110, 130, 175, 230,
                  300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750,
                  2100, 2500, 3000};
   const int nptd = int(sizeof(vptd)/sizeof(vptd[0])) - 1;

      // Parton/Gen matching control profiles
      // x-axis is the matched parton pT (from GenPart) for a given GenJet
      TProfile* p_partonPt_vs_partonPt        = new TProfile("p_partonPt_vs_partonPt",        ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV)",        nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt        = new TProfile("p_genJetPt_vs_partonPt",        ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV)",        nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt = new TProfile("p_genJetPtOverPartonPt_vs_partonPt", ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT", nptd, vptd);

      // FSR up/down (PSWeight[1]/PSWeight[3]) versions for Gen/parton studies
      TProfile* p_partonPt_vs_partonPt_fsrUp        = new TProfile("p_partonPt_vs_partonPt_fsrUp",        ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV) [FSR Up]",        nptd, vptd);
      TProfile* p_partonPt_vs_partonPt_fsrDown      = new TProfile("p_partonPt_vs_partonPt_fsrDown",      ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV) [FSR Down]",      nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt_fsrUp        = new TProfile("p_genJetPt_vs_partonPt_fsrUp",        ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV) [FSR Up]",        nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt_fsrDown      = new TProfile("p_genJetPt_vs_partonPt_fsrDown",      ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV) [FSR Down]",      nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_fsrUp   = new TProfile("p_genJetPtOverPartonPt_vs_partonPt_fsrUp",   ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [FSR Up]",   nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_fsrDown = new TProfile("p_genJetPtOverPartonPt_vs_partonPt_fsrDown", ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [FSR Down]", nptd, vptd);

      // Category splits for Gen/parton studies (keep inclusive profiles above)
      // b: |flav|==5, light: |flav| in {1,2,3,4}
      TProfile* p_partonPt_vs_partonPt_b        = new TProfile("p_partonPt_vs_partonPt_b",        ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV) [b]",        nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt_b        = new TProfile("p_genJetPt_vs_partonPt_b",        ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV) [b]",        nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_b = new TProfile("p_genJetPtOverPartonPt_vs_partonPt_b", ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [b]", nptd, vptd);

      TProfile* p_partonPt_vs_partonPt_light        = new TProfile("p_partonPt_vs_partonPt_light",        ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV) [light]",        nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt_light        = new TProfile("p_genJetPt_vs_partonPt_light",        ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV) [light]",        nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_light = new TProfile("p_genJetPtOverPartonPt_vs_partonPt_light", ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [light]", nptd, vptd);

      // FSR up/down category splits
      TProfile* p_partonPt_vs_partonPt_b_fsrUp        = new TProfile("p_partonPt_vs_partonPt_b_fsrUp",        ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV) [b, FSR Up]",        nptd, vptd);
      TProfile* p_partonPt_vs_partonPt_b_fsrDown      = new TProfile("p_partonPt_vs_partonPt_b_fsrDown",      ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV) [b, FSR Down]",      nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt_b_fsrUp        = new TProfile("p_genJetPt_vs_partonPt_b_fsrUp",        ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV) [b, FSR Up]",        nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt_b_fsrDown      = new TProfile("p_genJetPt_vs_partonPt_b_fsrDown",      ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV) [b, FSR Down]",      nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_b_fsrUp   = new TProfile("p_genJetPtOverPartonPt_vs_partonPt_b_fsrUp",   ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [b, FSR Up]",   nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_b_fsrDown = new TProfile("p_genJetPtOverPartonPt_vs_partonPt_b_fsrDown", ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [b, FSR Down]", nptd, vptd);

      TProfile* p_partonPt_vs_partonPt_light_fsrUp        = new TProfile("p_partonPt_vs_partonPt_light_fsrUp",        ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV) [light, FSR Up]",        nptd, vptd);
      TProfile* p_partonPt_vs_partonPt_light_fsrDown      = new TProfile("p_partonPt_vs_partonPt_light_fsrDown",      ";parton p_{T} (GeV);#LT parton p_{T} #GT (GeV) [light, FSR Down]",      nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt_light_fsrUp        = new TProfile("p_genJetPt_vs_partonPt_light_fsrUp",        ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV) [light, FSR Up]",        nptd, vptd);
      TProfile* p_genJetPt_vs_partonPt_light_fsrDown      = new TProfile("p_genJetPt_vs_partonPt_light_fsrDown",      ";parton p_{T} (GeV);#LT GenJet p_{T} #GT (GeV) [light, FSR Down]",      nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_light_fsrUp   = new TProfile("p_genJetPtOverPartonPt_vs_partonPt_light_fsrUp",   ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [light, FSR Up]",   nptd, vptd);
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_light_fsrDown = new TProfile("p_genJetPtOverPartonPt_vs_partonPt_light_fsrDown", ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [light, FSR Down]", nptd, vptd);

            // Selection-category splits (based on chosen jets in analysis selection)
      // selB: the two b-tag jets selected in the event
      // selLight: the W(light) jets selected in the event
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_selB = new TProfile(
         "p_genJetPtOverPartonPt_vs_partonPt_selB",
         ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [sel b]",
         nptd, vptd
      );
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_selLight = new TProfile(
         "p_genJetPtOverPartonPt_vs_partonPt_selLight",
         ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [sel light]",
         nptd, vptd
      );

      // FSR up/down for selection-category splits
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_selB_fsrUp = new TProfile(
         "p_genJetPtOverPartonPt_vs_partonPt_selB_fsrUp",
         ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [sel b, FSR Up]",
         nptd, vptd
      );
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_selB_fsrDown = new TProfile(
         "p_genJetPtOverPartonPt_vs_partonPt_selB_fsrDown",
         ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [sel b, FSR Down]",
         nptd, vptd
      );
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrUp = new TProfile(
         "p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrUp",
         ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [sel light, FSR Up]",
         nptd, vptd
      );
      TProfile* p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrDown = new TProfile(
         "p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrDown",
         ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [sel light, FSR Down]",
         nptd, vptd
      );

      // region thresholds use Rjet=0.4 (AK4)
TProfile* p_genJetPtOverPartonPt_vs_partonPt_dlt0p5R =
  new TProfile("p_genJetPtOverPartonPt_vs_partonPt_dlt0p5R",
               ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [d<0.5R]",
               nptd, vptd);

TProfile* p_genJetPtOverPartonPt_vs_partonPt_d0p5to2R =
  new TProfile("p_genJetPtOverPartonPt_vs_partonPt_d0p5to2R",
               ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [0.5R<d<2R]",
               nptd, vptd);

TProfile* p_genJetPtOverPartonPt_vs_partonPt_dgt2R =
  new TProfile("p_genJetPtOverPartonPt_vs_partonPt_dgt2R",
               ";parton p_{T} (GeV);#LT GenJet p_{T} / parton p_{T} #GT [d>2R]",
               nptd, vptd);


   // |eta| binning
   double vxd[]  = {0.00, 0.13, 0.26, 0.39, 0.52, 0.65, 0.78, 0.91, 1.04, 1.17, 1.30,
                  1.44, 1.57, 1.74, 1.93, 2.10, 2.30, 2.50};
   const int nxd = int(sizeof(vxd)/sizeof(vxd[0])) - 1;

   // p_reco/p_gen binning
   int reco_nedges = 101;
   double minValue = 0.0;
   double maxValue = 2.0;
   double stepSize = (maxValue - minValue) / (reco_nedges - 1);
      // Create the vector
   std::vector<double> vres(reco_nedges);

   // Populate the vector with values
   for (int i = 0; i < reco_nedges; ++i) {
         vres[i] = minValue + i * stepSize;
   }
   const int nres = sizeof(vres) / sizeof(vres[0]) - 1;


   fout->cd();
   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
   auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
   std::cout << std::ctime(&now) << std::endl<< flush;  
   cout << "Processing " << nentries << " events" << endl << flush;
   TStopwatch t;
   t.Start();
   const int nlap = 1000;
   const int nlap2 = 80000;
   //nentries = 100000;

   auto passRunEraSelection = [&](UInt_t runNumber) {
   if (isMC) return true;
   if (runEra == "Bnib1") return (runNumber == 401844u || runNumber == 401848u);
   if (runEra == "Bnib2") return (runNumber != 401844u && runNumber != 401848u);
   return true;
};

   TLorentzVector p4jet, lhe, jet, jet2, jetn;
   TLorentzVector met, met1, metn, metu, metnu, rawmet, corrmet, rawgam;
   TLorentzVector jeti, corrjets, rawjet, rawjets, rcjet, rcjets, rcoffsets;
   TLorentzVector geni, genjet, genjet2;
   TLorentzVector p4, p4g;
   TLorentzVector lj1_p4, lj2_p4, W_p4, singleJetW_p4, top1_p4, top2_p4, b1_p4, b2_p4, muon_p4, b1l_p4, b2l_p4, toph_p4, mbl_p4, bl_p4, bh_p4, toph_improved_p4, W_improved_p4;

   int _ntot(0), _nevents(0), _nbadevents_json(0), _nbadevents_trigger(0), _nbadevents_lumi(0);
   int _nbadevents_veto(0);
   int _nPass_3jets = 0;
   int _nPass_1btag = 0;
   int _nPass_3jets_1btag = 0;


   
   // JSON valinta (Photon + jet code) (Bettina)

   // -------------------------------------------------
   // JSON + lumi handling (DATA ONLY)
   // -------------------------------------------------
   if (!isMC) {

      std::cout << "[MODE] Data mode → loading JSON and lumi CSV" << std::endl;

      // JSON
      std::cout << "[JSON] Loading JSON: " << jsonFile << std::endl;
      LoadJSON(jsonFile.c_str());

      // Lumi CSV
      //source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
      //brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json --datatag online -i Collisions26_13p6TeV_Latest.json --hltpath HLT_IsoMu24_v*  -u /fb --minBiasXsec 75300 -o lumi_26_IsoMu24.csv
      std::string lumifile;
      if (runYear == 2024)
         lumifile = "lumi_378981_386951_IsoMu24.csv";
      else if (runYear == 2025)
         lumifile = "lumi_391658_398860_Golden_IsoMu24.csv";
      //else if (runYear == 2026)
         //lumifile = "lumi_26_IsoMu24.csv";

      if (lumifile.empty()) {
         std::cout << "[LUMI][ERROR] lumifile is empty in data mode!" << std::endl;
      } else {
         std::cout << "[LUMI] Loading lumi CSV: " << lumifile << std::endl;
         bool ok = LoadLumiAvgPU(lumifile);
         if (!ok)
            std::cout << "[LUMI][ERROR] LoadLumiAvgPU FAILED for " << lumifile << std::endl;
         else
            std::cout << "[LUMI] Lumi CSV loaded successfully" << std::endl;
      }

   } else {
      std::cout << "[MODE] MC mode → JSON and lumi CSV not used" << std::endl;
   }

    // Jetti korjaukset + tyyppi 1 Met korjaukset (JEC) //Winter24Prompt24 (Bettina)

   //FactorizedJetCorrector *jec(0);
   //jec = getFJC("", "Winter24Prompt24_RunG_V2_DATA_L2Relative_AK4PFPuppi", "Winter24Prompt24_RunG_V2_DATA_L2L3Residual_AK4PFPuppi");
   //jec = getFJC("", "Winter24Run3_V1_MC_L2Relative_AK4PUPPI", "Prompt24_Run2024G_V7M_DATA_L2L3Residual_AK4PFPuppi");

   // Pointers for JEC
   FactorizedJetCorrectorWrapper *jec2024(0);
   FactorizedJetCorrector *jec(0), *jersfvspt(0);;

   // Initialize JEC according to mode
   //if (isMC) {
   //   jec = getFJC("", "RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI", ""); //Winter24Run3_V1_MC_L2Relative_AK4PUPPI //RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI
   //   assert(jec);
   //} else {
   //   // Data: instantiate and assign directly to outer jec pointer
   //   jec = new FactorizedJetCorrectorWrapper();
   //   jec->addJECset("ReReco24_V9M_DATA"); //Reprocessing24_V8M_DATA //Prompt24_V8M_DATA //ReReco24_V9M
   //   assert(jec);
   //}
   // Dynamic JEC initialization
   if (isMC) {
       jec = getFJC("", jecMCSet, "");
       assert(jec);
       if (smearJets) {
         jerpathsf = Form("CondFormats/JetMETObjects/data/%s.txt", JERSFvsPt.c_str());
         jersfvspt = getFJC("", JERSFvsPt, "");
         // JER resolution (Pt Resolution txt)
         jerpath   = Form("CondFormats/JetMETObjects/data/%s", JERres.c_str());
         // Tell the smearing code to use the pt-dependent SF provider
         useJERSFvsPt = true;
       }
   } else if ((runYear == 2025) || (runYear == 2026)){
      jec = getFJC("", jecMCSet, jecDataSet);
      assert(jec);
   } else { //2024 with fibs.txt
       jec2024 = new FactorizedJetCorrectorWrapper();
       jec2024->addJECset(jecDataSet);
       assert(jec2024);
   }
   //FactorizedJetCorrectorWrapper *jec = new FactorizedJetCorrectorWrapper();
   //jec->addJECset("Prompt24_V8M_DATA");
   //jec->addJECset("Reprocessing24_V8M_DATA");
   //assert(jec);



    TFile *fjv(0);
    //fjv = new TFile("jet_veto_maps/Summer24ReReco/jetvetoReReco2024_V9M.root","READ");
    fjv = new TFile(jetVetoMap.c_str(),"READ");
    //jet_veto_maps/Winter24Prompt24/Winter24Prompt24_2024BCDEFGHI.root //jetvetoReReco2024_V9M.root
    if (!fjv) cout << "Jetvetomap file not found"  << endl << flush;
    assert(fjv);

    TH2D *h2jv = 0;
    TH2D *bpixjv = 0;
    TH2D *fpixjv = 0;
    h2jv = (TH2D*)fjv->Get("jetvetomap");
    bpixjv = (TH2D*)fjv->Get("jetvetomap_bpix"); //loading the bpix vetomap for all '24 stuff
    fpixjv = (TH2D*)fjv->Get("jetvetomap_fpix"); //loading the fpix vetomap for all '24 stuff
    if (!h2jv) cout << "Jetvetomap histo not found" << endl << flush;
    assert(h2jv);

   std::vector<int> goodMuonIndices;
   goodMuonIndices.reserve(nMuon);

   std::vector<int> bjetIndices;
   bjetIndices.reserve(nJet);

   std::vector<int> lightJetIndices;
   lightJetIndices.reserve(nJet);

   std::vector<int> gluonJetIndices;
   gluonJetIndices.reserve(nJet);

   // Smear JER
	JME::JetResolution *jer(0);
	JME::JetResolutionScaleFactor *jersf(0);

   if (isMC && smearJets)
   {
   std::cout << jerpath << std::endl << std::flush;
   if (!useJERSFvsPt)
      std::cout << jerpathsf << std::endl << std::flush;

   if (jerpath == "" || (jerpathsf == "" && !useJERSFvsPt))
      std::cout << "Missing JER file paths" << std::endl << std::flush;

   assert(jerpath != "");
   assert(jerpathsf != "" || useJERSFvsPt);
   assert(jersfvspt || !useJERSFvsPt);

   jer = new JME::JetResolution(jerpath.c_str());
   if (!useJERSFvsPt)
      jersf = new JME::JetResolutionScaleFactor(jerpathsf.c_str());

   if (!jer || (!jersf && !useJERSFvsPt) || (!jersfvspt && useJERSFvsPt))
      std::cout << "Missing JER files" << std::endl << std::flush;
   }

   //nentries = 1000000;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (!passRunEraSelection(run)) continue;
      runCountAfterSplit[run]++;

      if (jentry%nlap==0) {
         cout << "." << flush;
      }
      if (jentry%nlap2==0 && jentry!=0) {
         double time = t.RealTime();
      if (time>0) cout << Form("\n\%1.0f ev/s\n",nlap2/time) << flush;
         t.Reset();
         t.Start();
      }

      // ------------------------------------------------------------
      // Build GenJet_partonPt by ΔR-matching GenPart -> GenJet
      // Strategy:
      //   For each GenJet: find all GenParts with ΔR<0.2,
      //   pick the one with the largest pT, compare its flavour to
      //   GenJet_partonFlavour, and if equal store its pT.
      // Notes:
      //   - We consider quarks (|pdgId|=1..5) and gluons (pdgId=21).
      //   - If no match or flavour mismatch, the value stays -1.
      // ------------------------------------------------------------

      // Jet radius used for "in-cone" collection (AK4 => 0.4)
      const double Rjet = 0.4;

      std::vector<double> GenJet_partonPt;
      std::vector<double> GenJet_minPartonDR;
      std::vector<char>   GenJet_passSel; // 0/1: at least one candidate passed W-mother+lastCopy+quark
      if (isMC) {
         GenJet_partonPt.assign(nGenJet, -1.0f);
         GenJet_minPartonDR.assign(nGenJet, 999.0);
         GenJet_passSel.assign(nGenJet, 0);

         // Nominal and FSR PS weights for Gen/parton studies
         double w_fsrUp = genWeight * PSWeight[3]; // fsr.murfac=0.5
         double w_fsrDown  = genWeight * PSWeight[1]; // fsr.murfac=2.0

         // ===== Loop 1: for each GenJet, find best core parton (dr<0.2) and compute min DR between in-cone partons =====
         for (int igj = 0; igj < nGenJet; ++igj) {

            double gj_eta = GenJet_eta[igj];
            double gj_phi = GenJet_phi[igj];
            int    gj_flav = std::abs(GenJet_partonFlavour[igj]);

            int bestIdx = -1;
            double bestPt = -1.0;

            std::vector<int> conePartons;
            conePartons.clear();

            for (int igp = 0; igp < nGenPart; ++igp) {

               // keep only quarks
               if (abs(GenPart_pdgId[igp]) < 1 || abs(GenPart_pdgId[igp]) > 5) continue;

               // require isLastCopy (bit 13)
               if ( (GenPart_statusFlags[igp] & (1u << 13)) == 0u ) continue;

               // require valid mother
               if (GenPart_genPartIdxMother[igp] < 0 || GenPart_genPartIdxMother[igp] >= nGenPart) continue;

               // mother must be a W
               if (abs(GenPart_pdgId[GenPart_genPartIdxMother[igp]]) != 24) continue;

               // At least one GenPart passed your selection for this GenJet
               GenJet_passSel[igj] = 1;

               double dr = deltaR(gj_eta, gj_phi, GenPart_eta[igp], GenPart_phi[igp]);
               
               // collect partons in the full jet cone for multi-parton structure
               if (dr < Rjet) conePartons.push_back(igp);

               // keep your "best match" for partonPt mapping (core match)
               if (dr >= 0.2) continue;

               double pt = GenPart_pt[igp];
               if (pt > bestPt) {
                  bestPt = pt;
                  bestIdx = igp;
               }

               if (conePartons.size() >= 2) {
                  double dmin = 999.0;
                  for (size_t ia = 0; ia < conePartons.size(); ++ia) {
                     const int iA = conePartons[ia];
                     for (size_t ib = ia + 1; ib < conePartons.size(); ++ib) {
                        const int iB = conePartons[ib];
                        const double d = deltaR(GenPart_eta[iA], GenPart_phi[iA],
                                                GenPart_eta[iB], GenPart_phi[iB]);
                        if (d < dmin) dmin = d;
                     }
                  }
                  GenJet_minPartonDR[igj] = dmin;
                  } else {
                  // interpret as "0/1 selected parton in cone"
                  GenJet_minPartonDR[igj] = 999.0;
               }

            if (bestIdx >= 0) {
               int bestFlav = std::abs(GenPart_pdgId[bestIdx]);
               // Compare to GenJet_partonFlavour (NanoAOD convention is usually 1..5 or 21)
               if (bestFlav == gj_flav) {
                  GenJet_partonPt[igj] = bestPt;
               }
            }
         }
      } // end igj loop

         // ===== Loop 2: fill profiles ONLY for jets that passed your selection and have a core match =====
         for (int igj = 0; igj < nGenJet; ++igj) {

            // require: at least one candidate passed quark+lastCopy+mother=W checks
            if (!GenJet_passSel[igj]) continue;

            const double partonPt = GenJet_partonPt[igj];
                
            if (partonPt <= 0.0) continue; // require also the core match dr<0.2 + flavour match

            const double genPt = GenJet_pt[igj];
            p_partonPt_vs_partonPt->Fill(partonPt, partonPt, genWeight);
            p_genJetPt_vs_partonPt->Fill(partonPt, genPt, genWeight);
            p_genJetPtOverPartonPt_vs_partonPt->Fill(partonPt, genPt / partonPt, genWeight);

            // FSR up/down
            p_partonPt_vs_partonPt_fsrUp->Fill(partonPt, partonPt, w_fsrUp);
            p_partonPt_vs_partonPt_fsrDown->Fill(partonPt, partonPt, w_fsrDown);
            
            p_genJetPt_vs_partonPt_fsrUp->Fill(partonPt, genPt, w_fsrUp);
            p_genJetPt_vs_partonPt_fsrDown->Fill(partonPt, genPt, w_fsrDown);

            p_genJetPtOverPartonPt_vs_partonPt_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
            p_genJetPtOverPartonPt_vs_partonPt_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);

            const double dmin = GenJet_minPartonDR[igj];
            const double thr1 = 0.5 * Rjet;
            const double thr2 = 2.0 * Rjet;

            const bool reg1 = (dmin < thr1);
            const bool reg2 = (dmin >= thr1 && dmin < thr2);
            const bool reg3 = (dmin >= thr2); // includes dmin=999 (0/1 parton in cone)

            if (reg1) p_genJetPtOverPartonPt_vs_partonPt_dlt0p5R->Fill(partonPt, genPt/partonPt, genWeight);
            if (reg2) p_genJetPtOverPartonPt_vs_partonPt_d0p5to2R->Fill(partonPt, genPt/partonPt, genWeight);
            if (reg3) p_genJetPtOverPartonPt_vs_partonPt_dgt2R->Fill(partonPt, genPt/partonPt, genWeight);
            

            // Category splits by GenJet_partonFlavour
            bool isBJet = (std::abs(GenJet_partonFlavour[igj]) == 5);
            bool isLightJet = (std::abs(GenJet_partonFlavour[igj]) >= 1 && std::abs(GenJet_partonFlavour[igj]) <= 4);

            if (isBJet) {
               p_partonPt_vs_partonPt_b->Fill(partonPt, partonPt, genWeight);
               p_genJetPt_vs_partonPt_b->Fill(partonPt, genPt, genWeight);
               p_genJetPtOverPartonPt_vs_partonPt_b->Fill(partonPt, genPt / partonPt, genWeight);

               p_partonPt_vs_partonPt_b_fsrUp->Fill(partonPt, partonPt, w_fsrUp);
               p_partonPt_vs_partonPt_b_fsrDown->Fill(partonPt, partonPt, w_fsrDown);
               p_genJetPt_vs_partonPt_b_fsrUp->Fill(partonPt, genPt, w_fsrUp);
               p_genJetPt_vs_partonPt_b_fsrDown->Fill(partonPt, genPt, w_fsrDown);
               p_genJetPtOverPartonPt_vs_partonPt_b_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
               p_genJetPtOverPartonPt_vs_partonPt_b_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);
            }

            if (isLightJet) {
               p_partonPt_vs_partonPt_light->Fill(partonPt, partonPt, genWeight);
               p_genJetPt_vs_partonPt_light->Fill(partonPt, genPt, genWeight);
               p_genJetPtOverPartonPt_vs_partonPt_light->Fill(partonPt, genPt / partonPt, genWeight);

               p_partonPt_vs_partonPt_light_fsrUp->Fill(partonPt, partonPt, w_fsrUp);
               p_partonPt_vs_partonPt_light_fsrDown->Fill(partonPt, partonPt, w_fsrDown);
               p_genJetPt_vs_partonPt_light_fsrUp->Fill(partonPt, genPt, w_fsrUp);
               p_genJetPt_vs_partonPt_light_fsrDown->Fill(partonPt, genPt, w_fsrDown);
               p_genJetPtOverPartonPt_vs_partonPt_light_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
               p_genJetPtOverPartonPt_vs_partonPt_light_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);
            }
            
         }
      }


      // Fill run vs BX with NO cuts
      {
      Double_t x[2];
      x[0] = run;
      x[1] = bunchCrossing;
      h2runbx->Fill(x, 1.0);
      }

      // Trigger IsoMu24
      bool pass_trig = (HLT_IsoMu24);
      if (isMC) { 
         ++_nevents;
      } else {
         // JSON filter
         if (_json[run][luminosityBlock] == 0) {
            ++_nbadevents_json;
            continue;
         }

         // avgPU filter
         if (runYear != 2026) {
               if (_avgpu[run][luminosityBlock] == 0) {
               ++_nbadevents_lumi;
               continue;
            }
         }
         ++_nevents;
      }
         double w = (isMC ? (PSWeight[3]*genWeight) : 1);    //in case of MC set w to genWeight, otherwise (data) leave it 1
         if (isMC) h_sumw->Fill(0.5, genWeight);

         double Pileup_nTrue = Pileup_nTrueInt;

         // --- PU profiles (filled after basic filters) ---
         double muPU = 0.0;
         if (isMC)
            muPU = Pileup_nTrue;
         else
            muPU = _avgpu[run][luminosityBlock];

         double NPV  = PV_npvs;
         double NPV_Good = PV_npvsGood;
         double rho_Fastjet  = Rho_fixedGridRhoFastjetAll;

         h_muPU->Fill(muPU, w);
         pnpvvsmu->Fill(muPU, NPV, w);
         pnpvGoodvsmu->Fill(muPU, NPV_Good, w);
         prhovsmu->Fill(muPU, rho_Fastjet, w);
         pmuvsmu->Fill(muPU, muPU, w);

         // Select leading jets. Just exclude muon, don't apply JetID yet
         static const int nJetMax = 200;
         Float_t         Jet_resFactor[nJetMax]; // Custom addition
         Float_t         Jet_deltaJES[nJetMax]; // Custom addition
         Float_t         Jet_CF[nJetMax]; // Custom addition
         Float_t         Jet_pt_corr_preJER[nJetMax];   // NEW: corrected pT before any JER
         Float_t         Jet_mass_corr_preJER[nJetMax]; // NEW: corrected mass before any JER
         int iJet(-1), iJet2(-1), nJets(0);
         double djes(1), jes(1), res(1);
         jet.SetPtEtaPhiM(0,0,0,0);
         jet2.SetPtEtaPhiM(0,0,0,0);
         jetn.SetPtEtaPhiM(0,0,0,0);
         // Also calculate corrected type-I chsMET and HDM inputs
         corrjets.SetPtEtaPhiM(0,0,0,0);
         rawjets.SetPtEtaPhiM(0,0,0,0);
         rcjets.SetPtEtaPhiM(0,0,0,0);
         rcoffsets.SetPtEtaPhiM(0,0,0,0);
         // Jet loop for JEC
         for (int i = 0; i != nJet; ++i) {

            // Redo JEC on the fly (should be no previous use of corrected jets)
            if (jec2024!=0 || jec!=0) {

               double rawJetPt = Jet_pt[i] * (1.0 - Jet_rawFactor[i]);
               double rawJetMass = Jet_mass[i] * (1.0 - Jet_rawFactor[i]);
               if (isMC || (runYear != 2024)) {
                  jec->setJetPt(rawJetPt);
                  jec->setJetEta(Jet_eta[i]);
                  jec->setJetPhi(Jet_phi[i]);
               } else {
                  jec2024->setRun(run);
                  jec2024->setJetPt(rawJetPt);
                  jec2024->setJetEta(Jet_eta[i]);
                  jec2024->setJetPhi(Jet_phi[i]);
               }
               //double corr = jec->getCorrection();
               vector<float> v;
               if (isMC || (runYear != 2024)) {v = jec->getSubCorrections();
               } else {v = jec2024->getSubCorrections();}
               double corr = v.back();
               double res = (v.size()>1 ? v[v.size()-1]/v[v.size()-2] : 1.);

               double L2 = 1.0, L2L3Res = 1.0, L1 = 1.0;
               if (v.size() >= 1)          // always true here
                  L2 = v[0];
               if (v.size() >= 2)          // true if residual file present
                  L2L3Res = v.back() / v[v.size()-2];   // = L2L3Residual
               if (v.size() >= 3)          // only if you later add an L1 file
                  L1 = v[v.size()-3];

               prof_L2L3Res->Fill(run, L2L3Res);
               prof_L2L3   ->Fill(run, L2);        // this is just L2Relative here
               prof_L1     ->Fill(run, L1);        // stays 1.0 with two-level chain
               prof_corr   ->Fill(run, v.back());  // full factor
               //Jet_RES[i] = 1./res;
               Jet_deltaJES[i] = (1./corr) / (1.0 - Jet_rawFactor[i]);
               Jet_pt[i] = corr * rawJetPt;
               Jet_mass[i] = corr * rawJetMass;
               Jet_rawFactor[i] = (1.0 - 1.0/corr);
               Jet_resFactor[i] = (1.0 - 1.0/res);
               // Cache post-JEC, pre-JER values for JER control comparisons
               Jet_pt_corr_preJER[i]   = Jet_pt[i];
               Jet_mass_corr_preJER[i] = Jet_mass[i];
            }

            jeti.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

            if (iJet==-1) { // Leading jet for balance
               iJet = i;
               jet = jeti;
            }
            else { // Subleading jets 
               jetn += jeti;
               if (iJet2==-1) { // First subleading jet for alpha
                  iJet2 = i;
                  jet2 = jeti;
               }
            }
         } // for i in nJet 

         int njet = nJet;
         if (isMC && smearJets) {
            for (int i = 0; i != njet; ++i) {
               Jet_CF[i] = 1.;
               if (i < smearNMax) {
                  // Retrieve genJet and calculate dR
                  double dR(999);
                  p4.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);


                  if (Jet_genJetIdx[i] >= 0 ){ // && Jet_genJetIdx[i] < nGenJet)
                     int j = Jet_genJetIdx[i];
                     p4g.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j],
                                    GenJet_mass[j]);
                     dR = p4g.DeltaR(p4);
                  }
                  else
                     p4g.SetPtEtaPhiM(0, 0, 0, 0);


                  // Rename variables to keep naming as in jetphys/IOV.h.
                  double jPt = Jet_pt[i];
                  double jEta = Jet_eta[i];
                  double rho = Rho_fixedGridRhoFastjetAll;
                  double jE = p4.E();
                  double jPtGen = p4g.Pt();
                  // Set constants
                  double MIN_JET_ENERGY = 0.01; // TBD


                  // Some problems with the code below:
                  // 1) JER should us primarily genPt, secondary recoPt
                  // 2) relDpt  should evaluate vs genPt to avoid <1/x> != 1/<x> bias
                  // 3) For (JME)NANO, should also check DR of genJet
                  // Probably small impact except for the last, which I add


                  // The method presented here can be found in https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
                  // and the corresponding code in https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
                  assert(jer);
                  double Reso = jer->getResolution({{JME::Binning::JetPt, jPt}, {JME::Binning::JetEta, jEta}, {JME::Binning::Rho, rho}});
                  double SF(1);
                  if (useJERSFvsPt && jersfvspt) {
                     jersfvspt->setJetEta(jEta);
                     jersfvspt->setJetPt(jPt);
                     jersfvspt->setRho(rho);
                     SF = jersfvspt->getCorrection();
                  }
                  else if (!useJERSFvsPt && jersf){
                     SF = jersf->getScaleFactor({{JME::Binning::JetEta, jEta}, {JME::Binning::Rho, rho}}, Variation::NOMINAL);
                  //cout << "JER SF from jersf path" << endl;
                  }
                  else {
                     cout << "No JER SF available" << endl << flush;
                     assert(false);
                  }


                  // Case 0: by default the JER correction factor is equal to 1
                  double CF = 1.;
                  // We see if the gen jet meets our requirements
                  bool condPt = (jPtGen > MIN_JET_ENERGY && dR < 0.2);
                  double relDPt = condPt ? (jPt - jPtGen) / jPt : 0.0;
                  bool condPtReso = fabs(relDPt) < 3 * Reso;
                  if (condPt and condPtReso) {
                     // Case 1: we have a "good" gen jet matched to the reco jet (indicated by positive gen jet pt)
                     CF += (SF - 1.) * relDPt;
                  }
                  else if (SF > 1) {
                     // Case 2: we don't have a gen jet. Smear jet pt using a random gaussian variation
                     double sigma = Reso * std::sqrt(SF * SF - 1);
                     std::normal_distribution<> d(0, sigma);
                     CF += d(_mersennetwister);
                  }


                  // Negative or too small smearFactor. Safety precautions.
                  double CFLimit = MIN_JET_ENERGY / jE;
                  if (CF < CFLimit)
                     CF = CFLimit;


                  double origPt = Jet_pt[i];
                  double origJetMass = Jet_mass[i];
                  Jet_pt[i] = CF * origPt;
                  Jet_mass[i] = CF * origJetMass;
                  Jet_CF[i] = CF;
                  // Jet_smearFactor[i] = (1.0 - 1.0/CF);
               } // i<smearNMax
            }   // for njet
         }     // JER smearing

         /*for (int i = 0; i != nJet; ++i) {
               cout <<"CF: " << Jet_CF[i] << ", Jet_pt_preJER: " << Jet_pt_corr_preJER[i] << ", Jet_mass_preJER: " <<
               Jet_mass_corr_preJER[i] << ", Jet_pt: " << Jet_pt[i] << ", Jet_mass: " <<
               Jet_mass[i] << endl;
         }*/
         // Met filtterit (Noise filter 2024):

         // Run3: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Run_3_recommendations
         bool pass_filt = (//(isRun3 && Flag_METFilters>0) ||
            (//isRun3 &&
            Flag_goodVertices &&
            Flag_globalSuperTightHalo2016Filter &&
            Flag_EcalDeadCellTriggerPrimitiveFilter &&
            Flag_BadPFMuonFilter &&
            Flag_BadPFMuonDzFilter &&
            Flag_hfNoisyHitsFilter &&
            Flag_eeBadScFilter &&
            Flag_ecalBadCalibFilter));
         //) || isRun3; // pass_filt
   

         //JetVetoMap
         bool pass_jetveto = true;
         if (true) { // jet veto
            int i1 = h2jv->GetXaxis()->FindBin(jet.Eta());
            int j1 = h2jv->GetYaxis()->FindBin(jet.Phi());
            if (h2jv->GetBinContent(i1,j1)>0 or bpixjv->GetBinContent(i1,j1)>0 or fpixjv->GetBinContent(i1,j1)>0) {
               ++_nbadevents_veto;
               pass_jetveto = false;
            }
         } // jet veto
         // Met leikkaus? met+ sum pt_raw - sum pt_corr (>15 GeV for corrected) https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#Type_I_Correction_Propagation_of
         // ei tehdä
         rawmet.SetPtEtaPhiM(RawPuppiMET_pt, 0, RawPuppiMET_phi, 0);

         bool pass_basic = (pass_trig && pass_filt && pass_jetveto);

         // ========== Independent check: At least 3 jets, 1 b-tag ==========
         int nSelectedJets = 0;
         int nBTaggedJets = 0;
         double btagThreshold_skim = 0.0849; // same threshold as in skimming
         for (int j = 0; j < nJet; j++) {
               if (Jet_pt[j] > 15. && std::abs(Jet_eta[j]) < 2.4) {
                  ++nSelectedJets;
                  if (Jet_btagUParTAK4B[j] > btagThreshold_skim) {
                     ++nBTaggedJets;
                  }
               }
         }
         if (nSelectedJets >= 3) ++_nPass_3jets;
         if (nBTaggedJets >= 1) ++_nPass_1btag;
         if (nSelectedJets >= 3 && nBTaggedJets >= 1) ++_nPass_3jets_1btag;

        
         if (pass_basic){
               // valitse lepton = muon ja jet_index 

               // -------------------------------------------------
               // 1) Find and select muons; store their indices (iMu)
               // -------------------------------------------------

               goodMuonIndices.clear();

               for (int iMu = 0; iMu < nMuon; iMu++) {
                  double mu_pt  = Muon_pt[iMu];
                  double mu_eta = Muon_eta[iMu];
                  bool passMuonId = Muon_looseId[iMu];//Muon_tightId[iMu];
                  bool passMuonIsoId = (Muon_puppiIsoId[iMu] >= 1); //==3
                  if (passMuonIsoId && passMuonId && (mu_pt > 15.) && (fabs(mu_eta) < 2.4)){
                     // This muon passes our selection, so store the muon index (iMu)
                     goodMuonIndices.push_back(iMu);
                  }
               } // end muon loop

               if (goodMuonIndices.size() != 1) continue;

               // 2 johtavaa b jettiä (exclude muon)

               // -------------------------------------------------
               // 2) Select b-jets (exclude muon)
               // -------------------------------------------------

               bjetIndices.clear();
               lightJetIndices.clear();
               gluonJetIndices.clear();

               // Example threshold for "medium" btagThreshold
               double btagThreshold = 0.4319; //0.0849; //0.8482; //https://indico.cern.ch/event/1364748/contributions/5742611/subcontributions/458214/attachments/2778231/4854062/OH_btagSF_tnp_2022preEE_v0.pdf (slide 15/17)

               for (int j = 0; j < nJet; j++) {
                  // Exclude if matched to a selected muon
                  //if (Jet_muonIdx1[j] == goodMuonIndices[0]) continue;
                  //if (Jet_muonIdx2[j] == goodMuonIndices[0]) continue;
                  if (j == Muon_jetIdx[goodMuonIndices[0]]) continue;
                 
                  // If this jet passes the b-tag threshold => b-jet
                  if ((Jet_btagUParTAK4B[j] > btagThreshold) && (Jet_pt[j]  > 15.) && (fabs(Jet_eta[j]) < 2.4)) {
                     bjetIndices.push_back(j);
                  } 
                  else if (lightJetIndices.size() < 2 && (Jet_pt[j]  > 15.) && (fabs(Jet_eta[j]) < 2.4)){
                     //bool passesPt = (Jet_pt[j]  > 15.);
                     //bool passesEta  = (fabs(Jet_eta[j]) < 2.4); //myös blle^ > 30GeV jos pt alle 30GeV mutta b -> ei light!
                     // otherwise, treat it as a "light jet" 
                     lightJetIndices.push_back(j); //
                  }
                  else if ((Jet_pt[j]  > 15.)){
                     gluonJetIndices.push_back(j);
                  }
               }
               // 3) Now pick the leading two from each set (if available)
               if (bjetIndices.size() != 2) {
                  // Not enough b-jets => skip event
                  continue;
               }
               // Must have at least one light jet
               if (lightJetIndices.empty()) {
                  continue;
               }
               // If exactly one light jet and it's not "heavy", skip
               if (lightJetIndices.size() == 1 && Jet_mass[lightJetIndices[0]] < 60.) {
                  continue;
               }

               // The leading two b-jets
               int leadBJet   = bjetIndices[0];
               int secondBJet = bjetIndices[1];

               // The leading light jet (and duplicate if only one)
               int leadLightJet = lightJetIndices[0];
               int secondLightJet = -1;
               if (lightJetIndices.size() > 1) {
                   secondLightJet = lightJetIndices[1];
               }

               // -------------------------------------------------
               // 4) Lastly, check Jet ID for the chosen jets
               // -------------------------------------------------
               // In NanoAOD: bit 2 = tight, bit 3 = tightLeptonVeto
               bool passJetID = true;

               //if ((Jet_jetId[j] < 4)) continue;
/*
               // Macro-like function to test a single jet’s ID
               auto checkJetID = [&](int j) {
                  // If bit 2 is not set => fails tight
                  return ((Jet_jetId[j] & 4) != 0); //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV#Recommendations_for_the_13_6_AN1 JetID = (abs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0.01 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.99 );
               };
*/
               // Define a lambda that checks the new jet ID criteria for jet index j.
               auto checkJetID = [&](int j) -> bool {
                  // Apply the tight jet ID based on the jet's eta.
                  bool Jet_passJetIdTight = false;
                  if (abs(Jet_eta[j]) <= 2.6)
                     Jet_passJetIdTight = (Jet_neHEF[j] < 0.99) && (Jet_neEmEF[j] < 0.9) && (Jet_chMultiplicity[j]+Jet_neMultiplicity[j] > 1) && (Jet_chHEF[j] > 0.01) &&
                                          (Jet_chMultiplicity[j] > 0);
                  else if (abs(Jet_eta[j]) > 2.6 && abs(Jet_eta[j]) <= 2.7)
                     Jet_passJetIdTight = (Jet_neHEF[j] < 0.90) && (Jet_neEmEF[j] < 0.99);
                  else if (abs(Jet_eta[j]) > 2.7 && abs(Jet_eta[j]) <= 3.0)
                     Jet_passJetIdTight = (Jet_neHEF[j] < 0.99);
                  else if (abs(Jet_eta[j]) > 3.0)
                     Jet_passJetIdTight = (Jet_neMultiplicity[j] >= 2) && (Jet_neEmEF[j] < 0.4);

                  // Apply the additional lepton veto criteria.
                  bool Jet_passJetIdTightLepVeto = false;
                  if (abs(Jet_eta[j]) <= 2.7)
                     Jet_passJetIdTightLepVeto = Jet_passJetIdTight && (Jet_muEF[j] < 0.8) &&
                                                (Jet_chEmEF[j] < 0.8);
                  else
                     Jet_passJetIdTightLepVeto = Jet_passJetIdTight;

                  return Jet_passJetIdTightLepVeto;
               };


               // Check the two b-jets
               if (!checkJetID(leadBJet))   passJetID = false;
               if (!checkJetID(secondBJet)) passJetID = false;

               // Check the two light jets
               if (!checkJetID(leadLightJet))   passJetID = false;
               if (secondLightJet >= 0) {if (!checkJetID(secondLightJet)) passJetID = false;}

               if (!passJetID) {
                  // If any chosen jet fails ID => discard the event
                  continue;
               }

               // If you reach here, you have:
               //   - 1 selected muon
               //   - 2 highest-pt b-jets
               //   - 2 highest-pt light jets
               //   - All pass Jet ID
               // => Fill histograms or do further analysis

               // Build the 4-vectors for the two selected light jets
               lj1_p4.SetPtEtaPhiM(
                  Jet_pt[leadLightJet],
                  Jet_eta[leadLightJet],
                  Jet_phi[leadLightJet],
                  Jet_mass[leadLightJet]
               );
               if (secondLightJet >= 0) {
                  lj2_p4.SetPtEtaPhiM(
                     Jet_pt[secondLightJet],
                     Jet_eta[secondLightJet],
                     Jet_phi[secondLightJet],
                     Jet_mass[secondLightJet]
                  );
               }
               b1_p4.SetPtEtaPhiM(
                  Jet_pt[leadBJet],
                  Jet_eta[leadBJet],
                  Jet_phi[leadBJet],
                  Jet_mass[leadBJet]
               );

               b2_p4.SetPtEtaPhiM(
                  Jet_pt[secondBJet],
                  Jet_eta[secondBJet],
                  Jet_phi[secondBJet],
                  Jet_mass[secondBJet]
               );

               muon_p4.SetPtEtaPhiM(
                  Muon_pt[goodMuonIndices[0]],
                  Muon_eta[goodMuonIndices[0]],
                  Muon_phi[goodMuonIndices[0]],
                  Muon_mass[goodMuonIndices[0]]
               );


               bool heavyFallback = false;
               // Determine if one jet is “heavy” (mass > 60 GeV)
               if (lj1_p4.M() > 60.0 || lj2_p4.M() > 60.0) {
                  heavyFallback = true;
                  // Use the heavier jet alone for W
                  if (lj1_p4.M() > lj2_p4.M()) {
                     singleJetW_p4 = lj1_p4;
                  } else {
                     singleJetW_p4 = lj2_p4;
                  }
                  W_p4 = singleJetW_p4;
               } else {
                  // No heavy jet: form W from both light jets
                  W_p4 = lj1_p4 + lj2_p4;
               }


               // ---------- Selection-category fills (no auto, no function) ----------
               if (isMC) {
               // Nominal and FSR PS weights for Gen/parton studies
               double w_fsrDown = genWeight * PSWeight[3]; // fsr.murfac=0.5
               double w_fsrUp   = genWeight * PSWeight[1]; // fsr.murfac=2.0

                  // ---- selB: leadBJet ----
                  {
                     int genIdx = Jet_genJetIdx[leadBJet];
                     if (genIdx >= 0 && genIdx < nGenJet) {
                        double partonPt = GenJet_partonPt[genIdx];
                        if (partonPt > 0) {
                           double genPt = GenJet_pt[genIdx];

                           p_genJetPtOverPartonPt_vs_partonPt_selB->Fill(partonPt, genPt / partonPt, genWeight);
                           p_genJetPtOverPartonPt_vs_partonPt_selB_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
                           p_genJetPtOverPartonPt_vs_partonPt_selB_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);
                        }
                     }
                  }

                  // ---- selB: secondBJet ----
                  {
                     int genIdx = Jet_genJetIdx[secondBJet];
                     if (genIdx >= 0 && genIdx < nGenJet) {
                        double partonPt = GenJet_partonPt[genIdx];
                        if (partonPt > 0) {
                           double genPt = GenJet_pt[genIdx];

                           p_genJetPtOverPartonPt_vs_partonPt_selB->Fill(partonPt, genPt / partonPt, genWeight);
                           p_genJetPtOverPartonPt_vs_partonPt_selB_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
                           p_genJetPtOverPartonPt_vs_partonPt_selB_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);
                        }
                     }
                  }

                  // ---- selLight: choose indices according to your W definition ----
                  if (secondLightJet < 0) {
                     // only one light jet
                     int lj = leadLightJet;
                     int genIdx = Jet_genJetIdx[lj];
                     if (genIdx >= 0 && genIdx < nGenJet) {
                        double partonPt = GenJet_partonPt[genIdx];
                        if (partonPt > 0) {
                           double genPt = GenJet_pt[genIdx];

                           p_genJetPtOverPartonPt_vs_partonPt_selLight->Fill(partonPt, genPt / partonPt, genWeight);
                           p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
                           p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);
                        }
                     }
                  }
                  else if (heavyFallback) {
                     // only the heavier-mass jet contributes to W
                     int heavyIdx = (lj1_p4.M() > lj2_p4.M()) ? leadLightJet : secondLightJet;
                     int genIdx = Jet_genJetIdx[heavyIdx];
                     if (genIdx >= 0 && genIdx < nGenJet) {
                        double partonPt = GenJet_partonPt[genIdx];
                        if (partonPt > 0) {
                           double genPt = GenJet_pt[genIdx];

                           p_genJetPtOverPartonPt_vs_partonPt_selLight->Fill(partonPt, genPt / partonPt, genWeight);
                           p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
                           p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);
                        }
                     }
                  }
                  else {
                     // both light jets contribute to W
                     // leadLightJet
                     {
                        int genIdx = Jet_genJetIdx[leadLightJet];
                        if (genIdx >= 0 && genIdx < nGenJet) {
                           double partonPt = GenJet_partonPt[genIdx];
                           if (partonPt > 0) {
                              double genPt = GenJet_pt[genIdx];

                              p_genJetPtOverPartonPt_vs_partonPt_selLight->Fill(partonPt, genPt / partonPt, genWeight);
                              p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
                              p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);
                           }
                        }
                     }
                     // secondLightJet
                     {
                        int genIdx = Jet_genJetIdx[secondLightJet];
                        if (genIdx >= 0 && genIdx < nGenJet) {
                           double partonPt = GenJet_partonPt[genIdx];
                           if (partonPt > 0) {
                              double genPt = GenJet_pt[genIdx];

                              p_genJetPtOverPartonPt_vs_partonPt_selLight->Fill(partonPt, genPt / partonPt, genWeight);
                              p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrUp->Fill(partonPt, genPt / partonPt, w_fsrUp);
                              p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrDown->Fill(partonPt, genPt / partonPt, w_fsrDown);
                           }
                        }
                     }
                  }
               }
               // ---------- end selection-category fills ----------

               double topHad_pt_true = -1.0;

               if (isMC) {
                  int idxTopHad = -1;

                  // 1) find top quarks
                  for (int ig = 0; ig < nGenPart; ++ig) {
                     if (std::abs(GenPart_pdgId[ig]) != 6) continue;

                     // 2) decide if this top is hadronic: look at its W decay
                     bool hasLeptonDaughter = false;

                     for (int j = 0; j < nGenPart; ++j) {

                        int mother_full = GenPart_genPartIdxMother[j];
                        // decode: keep only lower 16 bits; negative stays -1
                        int mother_idx  = (mother_full >= 0) ? (mother_full & 0xFFFF) : -1;

                        if (mother_idx != ig) continue;
                        int pdg = std::abs(GenPart_pdgId[j]);
                        if (pdg == 11 || pdg == 13 || pdg == 15) {
                           hasLeptonDaughter = true;
                           break;
                        }
                     }

                     // if no lepton among daughters → treat as hadronic top
                     if (!hasLeptonDaughter) {
                           idxTopHad = ig;
                           break;
                     }
                  }

                  if (idxTopHad >= 0) {
                     TLorentzVector topHad_true_p4;
                     topHad_true_p4.SetPtEtaPhiM(
                           GenPart_pt[idxTopHad],
                           GenPart_eta[idxTopHad],
                           GenPart_phi[idxTopHad],
                           GenPart_mass[idxTopHad]
                     );
                     topHad_pt_true = topHad_true_p4.Pt();
                  }          
               }

               // ===== Build pre-JER versions for resolution-SF control =====
               TLorentzVector lj1_p4_preJER, lj2_p4_preJER, b1_p4_preJER, b2_p4_preJER, b1l_p4_preJER, b2l_p4_preJER, bh_p4_preJER, bl_p4_preJER;
               TLorentzVector W_p4_preJER, singleJetW_p4_preJER, toph_p4_preJER, top1_p4_preJER, top2_p4_preJER; 

               lj1_p4_preJER.SetPtEtaPhiM(
                  Jet_pt_corr_preJER[leadLightJet],
                  Jet_eta[leadLightJet],
                  Jet_phi[leadLightJet],
                  Jet_mass_corr_preJER[leadLightJet]
               );
               if (secondLightJet >= 0) {
                  lj2_p4_preJER.SetPtEtaPhiM(
                     Jet_pt_corr_preJER[secondLightJet],
                     Jet_eta[secondLightJet],
                     Jet_phi[secondLightJet],
                     Jet_mass_corr_preJER[secondLightJet]
                  );
               }
               b1_p4_preJER.SetPtEtaPhiM(
                  Jet_pt_corr_preJER[leadBJet],
                  Jet_eta[leadBJet],
                  Jet_phi[leadBJet],
                  Jet_mass_corr_preJER[leadBJet]
               );
               b2_p4_preJER.SetPtEtaPhiM(
                  Jet_pt_corr_preJER[secondBJet],
                  Jet_eta[secondBJet],
                  Jet_phi[secondBJet],
                  Jet_mass_corr_preJER[secondBJet]
               );

               // Mirror the same W-building logic using the smeared-jet heavyFallback decision
               if (secondLightJet < 0) {
                  W_p4_preJER = lj1_p4_preJER;
               } else if (heavyFallback) {
                  // choose heavier *smeared* jet, but use corresponding pre-JER vector
                  if (lj1_p4.M() > lj2_p4.M()) singleJetW_p4_preJER = lj1_p4_preJER;
                  else                         singleJetW_p4_preJER = lj2_p4_preJER;
                  W_p4_preJER = singleJetW_p4_preJER;
               } else {
                  W_p4_preJER = lj1_p4_preJER + lj2_p4_preJER;
               }
               double W_mass_preJER   = W_p4_preJER.M();

               top1_p4_preJER = b1_p4_preJER + W_p4_preJER;
               top2_p4_preJER = b2_p4_preJER + W_p4_preJER;
               
               double Top1_mass_preJER = top1_p4_preJER.M();
               double Top2_mass_preJER = top2_p4_preJER.M();

               bh_p4_preJER = (fabs(Top1_mass_preJER-172.5)< fabs(Top2_mass_preJER-172.5) ? b1_p4_preJER : b2_p4_preJER);
               bl_p4_preJER = (fabs(Top1_mass_preJER-172.5)< fabs(Top2_mass_preJER-172.5) ? b2_p4_preJER : b1_p4_preJER);

               toph_p4_preJER = bh_p4_preJER + W_p4_preJER;
               double Top_mass_preJER = toph_p4.M();

               // Control plots (pre-JER baselines)
               if (!TopPtDependentMass){
                  h_ctrl_Wmass_preJER->Fill(W_mass_preJER, w);
                  h_ctrl_Topmass_preJER->Fill(Top_mass_preJER, w);
               }
               // ===== Build pre-JER versions for resolution-SF control =====

               top1_p4 = b1_p4 + W_p4;
               double Top1_mass = top1_p4.M();  // invariant mass in GeV

               top2_p4 = b2_p4 + W_p4;
               double Top2_mass = top2_p4.M();  // invariant mass in GeV

               double W_mass = W_p4.M();  // invariant mass in GeV
               double W_pt = W_p4.Pt();   // invariant pt in GeV
               // Fill W eta distribution
               
               double singleJetW_mass = singleJetW_p4.M();

               W_improved_p4 = W_p4 * (80.4/W_p4.M());
               double W_mass_improved = W_improved_p4.M();

               b1l_p4 = b1_p4 + muon_p4;
               double mbl1 =  b1l_p4.M();

               b2l_p4 = b2_p4 + muon_p4;
               double mbl2 =  b2l_p4.M();

               bh_p4 = (fabs(Top1_mass-172.5)< fabs(Top2_mass-172.5) ? b1_p4 : b2_p4);
               bl_p4 = (fabs(Top1_mass-172.5)< fabs(Top2_mass-172.5) ? b2_p4 : b1_p4);

               toph_p4 = bh_p4 + W_p4;
               double Top_mass = toph_p4.M();  // invariant mass in GeV

               if (isMC && topHad_pt_true > 0.0) {
                  p_Topmass_vs_TopptTrue->Fill(topHad_pt_true, Top_mass, w);
               }
               
               p_Topmass_vs_Toppt_hadronic->Fill(toph_p4.Pt(), Top_mass, w);

               toph_improved_p4 = bh_p4 + W_improved_p4;
               double Top_mass_improved = toph_improved_p4.M();

               TLorentzVector Top_final;
               // 2) compute scale to PDG top (e.g. 172.5 GeV)
               Top_final = toph_improved_p4 * (172.5 / toph_improved_p4.M());
               TLorentzVector W_final;
               W_final = Top_final - bh_p4;

               // 3) apply that scale to W’s 3-vector only
               
               W_final.SetPxPyPzE(
                  W_final.Px(),
                  W_final.Py(),
                  W_final.Pz(),
                  W_improved_p4.E()   // keep energy fixed
               );

               // 4) rebuild the “final” top
               TLorentzVector top_final = bh_p4 + W_final;
               double Top_mass_final = top_final.M();   // should now sit at ≃ pdgTop
               double W_mass_final = W_final.M();
               mbl_p4 = bl_p4 + muon_p4;
               double mbl = mbl_p4.M();

               double mbl_red = mbl/Top_mass;

               prof_top->Fill(run, Top_mass, w); 
               prof_top_improved->Fill(run, Top_mass_improved, w); 
               prof_W->Fill(run, W_mass, w);

               if (isMC)
                  muPU = Pileup_nTrue;
               else
                  muPU = _avgpu[run][luminosityBlock];
               
               NPV  = PV_npvs;
               NPV_Good = PV_npvsGood;
               rho_Fastjet  = Rho_fixedGridRhoFastjetAll;

               h_muPU_after->Fill(muPU, w);
               pnpvvsmu_after->Fill(muPU, NPV, w);
               pnpvGoodvsmu_after->Fill(muPU, NPV_Good, w);
               prhovsmu_after->Fill(muPU, rho_Fastjet, w);
               pmuvsmu_after->Fill(muPU, muPU, w);

               int flav1 = Jet_partonFlavour[leadLightJet];
               int flav2 = Jet_partonFlavour[secondLightJet];

               // Fill the histogram
               if (!TopPtDependentMass){
                  h_Wmass->Fill(W_mass, w);
                  h_Wmass_pt->Fill(W_pt, W_mass, w);
                  h_singleJetWmass->Fill(singleJetW_mass, w);
                  h_Wpt->Fill(W_pt, w);
                  h_Weta->Fill(W_p4.Eta(), w);
               }

               // Fill W mass / pt and profile
               if (TopPtDependentMass) {
                  if (W_mass_preJER > 30.0 && W_mass_preJER < 130.0) h_ctrl_Wmass_preJER->Fill(W_mass_preJER, w);
                  if (singleJetW_mass > 30.0 && singleJetW_mass < 130.0){
                     // --- Single-jet W (heavy-jet fallback) ---
                     if (heavyFallback) {
                        double singleJetW_pt = singleJetW_p4.Pt();
                        h_singleJetWpt->Fill(singleJetW_pt, w);
                        p_singleJetWmass_vs_Wpt->Fill(singleJetW_pt, singleJetW_mass, w);
                     }
                  }
                  if (W_mass > 30.0 && W_mass < 130.0){   
                     h_Wmass->Fill(W_mass, w);
                     h_Weta->Fill(W_p4.Eta(), w);
                     h_Wmass_pt->Fill(W_pt, W_mass, w);
                     h_singleJetWmass->Fill(singleJetW_mass, w);
                     h_Wpt->Fill(W_pt, w);
                     p_Wmass_vs_Wpt->Fill(W_pt, W_mass, w);
                  }

                  h_ctrl_Topmass_preJER->Fill(Top_mass_preJER, w);
                  // --- Top mass vs Top pT ---
                  double Top_pt = toph_p4.Pt();  // b_h + W
                  h_Topmass->Fill(Top_mass, w);
                  h_Toppt->Fill(Top_pt, w);
                  h_Topmass_pt->Fill(Top_pt, Top_mass, w);
                  p_Topmass_vs_Toppt->Fill(Top_pt, Top_mass, w);

               }

               if ((isMC || (runYear != 2024)) && secondLightJet >= 0) {
                  int raw1 = Jet_partonFlavour[leadLightJet];
                  int raw2 = Jet_partonFlavour[secondLightJet];
                  bool isQuark1 = abs(raw1) >= 1 && abs(raw1) <= 5;
                  bool isQuark2 = abs(raw2) >= 1 && abs(raw2) <= 5;
                  bool isGluon1 = (raw1 == 21);
                  bool isGluon2 = (raw2 == 21);

                  if (
                     // antiquark + quark (only up-type + down-type)
                     (isQuark1 && isQuark2 && raw1*raw2 < 0 && isUpDownPair(flav1, flav2)) ||
                     // quark + gluon
                     (isQuark1 && isGluon2) ||
                     (isGluon1 && isQuark2)  ||
                     // gluon + gluon
                     (isGluon1 && isGluon2)
                  ) {
                     int f1 = min(raw1, 6);
                     int f2 = min(raw2, 6);
                     h2_Wjets->Fill(f1, f2, w);
                  }
               }
               h_Jet1Pt->Fill(lj1_p4.Pt(), w);
               if (!heavyFallback && secondLightJet >= 0) {
                   h_Jet2Pt->Fill(lj2_p4.Pt(), w);
               }

               // histograms: top1 = b1 + lj1 + lj2, top2 = b2 + lj1+ lj2
               h_Top2mass->Fill(Top2_mass, w);
               h_Top1mass->Fill(Top1_mass, w);
               h_Topmass->Fill(Top_mass, w);
               h_Topmass_improved->Fill(Top_mass_improved, w);
                // massb&lep = mbl
               h_mbl1->Fill(mbl1, w);
               h_mbl2->Fill(mbl2, w);
               h_mbl->Fill(mbl, w);
               h_mbl_red->Fill(mbl_red, w);

               double rbq = (b1_p4.Pt()+b2_p4.Pt())/(lj1_p4.Pt()+lj2_p4.Pt());
               double rbl = (b1_p4.Pt()+b2_p4.Pt())/(2*muon_p4.Pt());
               // pt: b1+b2/(lj1+lj2) = rbq
               if (secondLightJet >= 0) {h_rbq->Fill(rbq, w);}
               // rbl: b1+b2/(2*l1)
               h_rbl->Fill(rbl, w);

               double pt_pair = TMath::Sqrt( lj1_p4.Pt() * lj2_p4.Pt());
               double pt_pair_improved = (80.4/W_p4.M()) * pt_pair;
               
               if (secondLightJet >= 0) {
                  if (!heavyFallback) {
                     double rawJetPt1 = Jet_pt[leadLightJet] * (1.0 - Jet_rawFactor[leadLightJet]);
                     if (isMC || (runYear != 2024)){
                        jec->setJetPt(rawJetPt1);
                        jec->setJetEta(Jet_eta[leadLightJet]);
                        jec->setJetPhi(Jet_phi[leadLightJet]);
                     } else {
                        jec2024->setRun(run);
                        jec2024->setJetPt(rawJetPt1);
                        jec2024->setJetEta(Jet_eta[leadLightJet]);
                        jec2024->setJetPhi(Jet_phi[leadLightJet]);
                     }

                     vector<float> v1;
                     if (isMC || (runYear != 2024)) {v1 = jec->getSubCorrections();
                     } else {v1 = jec2024->getSubCorrections();}

                     if (fabs(lj1_p4.Eta()) < 1.3 && fabs(lj2_p4.Eta())  < 1.3) { //&& pt_pair > 53. && pt_pair < 107.){
                        h2_Wmass_ptpair->Fill(pt_pair, W_mass, w);
                        h2_Wmass_ptpair_improved->Fill(pt_pair_improved, W_mass, w);

                        prof_W_ptpair->Fill(pt_pair, W_mass, w);
                        prof_W_ptpair_improved->Fill(pt_pair_improved, W_mass, w);

                        struct Range { double xlo,xhi,ylo,yhi; };
                        static const Range sigBins[] = {
                           // tighter condition
                           {30,35,60,65}, {30,35,65,70}, {35,40,65,70}, {35,40,70,75}, {35,40,75,80},
                           {40,50,60,65}, {40,50,65,70}, {40,50,70,75}, {40,50,75,80}, {40,50,80,85},
                           {40,50,85,90}, {40,50,90,95}, {40,50,95,100},
                           {50,60,65,70}, {50,60,70,75}, {50,60,75,80}, {50,60,80,85}, {50,60,85,90}, {50,60,90,95}, {50,60,95,100},
                           {60,75,65,70}, {60,75,70,75}, {60,75,75,80}, {60,75,80,85}, {60,75,85,90}, {60,75,90,95},
                           {75,90,75,80}, {75,90,80,85}, {75,90,85,90},
                           {90,110,80,85}, {90,110,85,90}
                           // looser condition
                           /*
                           {25,30,50,55}, {25,30,55,60}, {25,30,60,65},
                           {30,35,50,55}, {30,35,55,60}, {30,35,60,65}, {30,35,65,70}, {30,35,70,75}, {30,35,75,80},
                           {35,40,50,55}, {35,40,55,60}, {35,40,60,65}, {35,40,65,70}, {35,40,70,75}, {35,40,75,80}, {35,40,80,85}, {35,40,85,90},
                           {40,50,45,50}, {40,50,50,55}, {40,50,55,60}, {40,50,60,65}, {40,50,65,70}, {40,50,70,75}, {40,50,75,80}, {40,50,80,85},
                           {40,50,85,90}, {40,50,90,95}, {40,50,95,100}, {40,50,100,105},
                           {50,60,50,55}, {50,60,55,60}, {50,60,60,65}, {50,60,65,70}, {50,60,70,75}, {50,60,75,80}, {50,60,80,85},
                           {50,60,85,90}, {50,60,90,95}, {50,60,95,100}, {50,60,100,105}, {50,60,105,110}, {50,60,110,115},
                           {60,75,55,60}, {60,75,60,65}, {60,75,65,70}, {60,75,70,75}, {60,75,75,80}, {60,75,80,85}, {60,75,85,90},
                           {60,75,90,95}, {60,75,95,100}, {60,75,100,105}, {60,75,105,110},
                           {75,90,65,70}, {75,90,70,75}, {75,90,75,80}, {75,90,80,85}, {75,90,85,90}, {75,90,90,95}, {75,90,95,100},
                           {90,110,70,75}, {90,110,75,80}, {90,110,80,85}, {90,110,85,90}, {90,110,90,95}, {90,110,95,100},
                           {110,130,75,80}, {110,130,80,85}, {110,130,85,90}, {110,130,90,95},
                           {130,175,80,85}, {130,175,85,90}, {130,175,90,95}*/
                        };
                        static const Range bkgBins[] = {
                        {20,25,35,40}, {20,25,40,45}, {20,25,45,50},
                        {25,30,45,50}, {25,30,50,55}, {25,30,55,60},
                        {30,35,25,30}, {30,35,40,45}, {30,35,45,50}, {30,35,50,55}, {30,35,55,60}, {30,35,60,65}, {30,35,65,70}, {30,35,70,75},
                        {35,40,25,30}, {35,40,30,35}, {35,40,35,40}, {35,40,40,45}, {35,40,45,50}, {35,40,50,55}, {35,40,55,60}, {35,40,60,65}, {35,40,65,70}, {35,40,70,75}, {35,40,75,80}, {35,40,80,85},
                        {40,50,25,30}, {40,50,30,35}, {40,50,35,40}, {40,50,40,45}, {40,50,45,50}, {40,50,50,55}, {40,50,55,60}, {40,50,60,65}, {40,50,65,70}, {40,50,70,75}, {40,50,75,80}, {40,50,80,85},
                        {40,50,85,90}, {40,50,90,95}, {40,50,95,100}, {40,50,100,105}, {40,50,105,110},
                        {50,60,30,35}, {50,60,35,40}, {50,60,40,45}, {50,60,45,50}, {50,60,50,55}, {50,60,55,60}, {50,60,60,65}, {50,60,65,70}, {50,60,70,75}, {50,60,75,80}, {50,60,80,85},
                        {50,60,85,90}, {50,60,90,95}, {50,60,95,100}, {50,60,100,105}, {50,60,105,110}, {50,60,110,115}, {50,60,115,120}, {50,60,120,125}, {50,60,125,130}, {50,60,130,135},
                        {60,75,35,40}, {60,75,40,45}, {60,75,45,50}, {60,75,50,55}, {60,75,55,60}, {60,75,60,65}, {60,75,65,70}, {60,75,70,75}, {60,75,75,80}, {60,75,80,85}, {60,75,85,90},
                        {60,75,90,95}, {60,75,95,100}, {60,75,100,105}, {60,75,105,110}, {60,75,110,115}, {60,75,115,120}, {60,75,120,125}, {60,75,125,130}, {60,75,130,135}, {60,75,135,140},
                        {60,75,140,145}, {60,75,145,150}, {60,75,150,155}, {60,75,155,160},
                        {75,90,55,60}, {75,90,60,65}, {75,90,145,150}, {75,90,150,155}, {75,90,155,160}, {75,90,160,165}, {75,90,165,170}, {75,90,170,175}, {75,90,175,180}
                        };

                        bool isSignal = false;
                        for (auto &r : sigBins) {
                        if (pt_pair >= r.xlo && pt_pair < r.xhi
                           && W_mass  >= r.ylo && W_mass  < r.yhi) {
                           isSignal = true;
                           break;
                        }
                        }

                        bool isBackground = false;
                        for (auto &r : bkgBins) {
                        if (pt_pair >= r.xlo && pt_pair < r.xhi
                           && W_mass  >= r.ylo && W_mass  < r.yhi) {
                           isBackground = true;
                           break;
                        }
                        }

                        if (isSignal || isBackground){
                           if (isSignal) {
                           h2_W_ptpair_sig->Fill(pt_pair, W_mass, w);
                           } 
                           if (isBackground) {
                              h2_W_ptpair_bkg->Fill(pt_pair, W_mass, w);
                           }
                        } else {
                           h2_W_ptpair_cr->Fill(pt_pair, W_mass, w);
                        }
                        
                        h_ptpair->Fill(pt_pair, w);
                        h_ptpair_improved->Fill(pt_pair_improved, w);
                        if (isMC || (runYear != 2024)){
                           bool isG1 = (flav1 == 21), isG2 = (flav2 == 21);
                           if      (isG1 && isG2) {
                              h2_Wmass_ptpair_gg->Fill(pt_pair, W_mass, w);
                              prof_W_ptpair_gg->Fill(pt_pair, W_mass, w);
                              h_ptpair_gg->Fill(pt_pair, w);
                              h2_Wmass_ptpair_improved_gg->Fill(pt_pair_improved, W_mass, w);
                              prof_W_ptpair_improved_gg->Fill(pt_pair_improved, W_mass, w);
                              h_ptpair_improved_gg->Fill(pt_pair_improved, w);
                           }
                           else if (isG1 != isG2) {
                              h2_Wmass_ptpair_qg->Fill(pt_pair, W_mass, w);
                              prof_W_ptpair_qg->Fill(pt_pair, W_mass, w);
                              h_ptpair_qg->Fill(pt_pair, w);
                              h2_Wmass_ptpair_improved_qg->Fill(pt_pair_improved, W_mass, w);
                              prof_W_ptpair_improved_qg->Fill(pt_pair_improved, W_mass, w);
                              h_ptpair_improved_qg->Fill(pt_pair_improved, w);
                           }
                           else if (qqPairs(flav1, flav2)) {
                              h2_Wmass_ptpair_qq->Fill(pt_pair, W_mass, w);
                              prof_W_ptpair_qq->Fill(pt_pair, W_mass, w);
                              h_ptpair_qq->Fill(pt_pair, w);
                              h2_Wmass_ptpair_improved_qq->Fill(pt_pair_improved, W_mass, w);
                              prof_W_ptpair_improved_qq->Fill(pt_pair_improved, W_mass, w);
                              h_ptpair_improved_qq->Fill(pt_pair_improved, w);
                           }
                           else {
                              h2_Wmass_ptpair_qq_others->Fill(pt_pair, W_mass, w);
                              prof_W_ptpair_qq_others->Fill(pt_pair, W_mass, w);
                              h_ptpair_qq_others->Fill(pt_pair, w);
                              h2_Wmass_ptpair_improved_qq_others->Fill(pt_pair_improved, W_mass, w);
                              prof_W_ptpair_improved_qq_others->Fill(pt_pair_improved, W_mass, w);
                              h_ptpair_improved_qq_others->Fill(pt_pair_improved, w);
                           }
                        }
                        // Geometric average of the two light‑jet factors (to match pt_pair definition)
                        double rawJetPt2 = Jet_pt[secondLightJet] * (1.0 - Jet_rawFactor[secondLightJet]);
                        if (isMC || (runYear != 2024)){
                           jec->setJetPt(rawJetPt2);
                           jec->setJetEta(Jet_eta[secondLightJet]);
                           jec->setJetPhi(Jet_phi[secondLightJet]);
                        } else {
                           jec2024->setRun(run);
                           jec2024->setJetPt(rawJetPt2);
                           jec2024->setJetEta(Jet_eta[secondLightJet]);
                           jec2024->setJetPhi(Jet_phi[secondLightJet]);
                        }
                        vector<float> v2;
                        if (isMC || (runYear != 2024)) {v2 = jec->getSubCorrections();
                        } else {v2 = jec2024->getSubCorrections();}

                        /* ---- derive geometric-mean sub-factors ------------------------- */
                        double L1_1 = 1.0, L2_1 = 1.0, Res_1 = 1.0;
                        if (v1.size() >= 1) L2_1  = v1[0];
                        if (v1.size() >= 2) Res_1 = v1.back() / v1[v1.size()-2];
                        if (v1.size() >= 3) L1_1  = v1[v1.size()-3];

                        double L1_2 = 1.0, L2_2 = 1.0, Res_2 = 1.0;
                        if (v2.size() >= 1) L2_2  = v2[0];
                        if (v2.size() >= 2) Res_2 = v2.back() / v2[v2.size()-2];
                        if (v2.size() >= 3) L1_2  = v2[v2.size()-3];

                        double L1_g      = TMath::Sqrt(L1_1 * L1_2);
                        double L2_g      = TMath::Sqrt(L2_1 * L2_2);          // η-shape
                        double Res_g     = TMath::Sqrt(Res_1 * Res_2);        // residual
                        double corrAvg_g = TMath::Sqrt(v1.back() * v2.back()); // full factor
                        /* ---------------------------------------------------------------- */
                        prof_L2L3Res_ptpair->Fill(pt_pair, Res_g, w);
                        prof_L2L3_ptpair->Fill(pt_pair, L2_g, w);
                        prof_L1_ptpair->Fill(pt_pair, L1_g, w);
                        prof_corr_ptpair->Fill(pt_pair, corrAvg_g, w);
                        if (Res_g > 0.0) pres->Fill(pt_pair, 1.0/Res_g, w);
                        if (L2_g  > 0.0) pjes->Fill(pt_pair, 1.0/L2_g, w);
                     }
                     //double top_mass = (fabs(Top1_mass-172.5)< fabs(Top2_mass-172.5) ? Top1_mass : Top2_mass);
                  }
               }
               //if (fabs(W_mass-80.4) < 20. && fabs(top_mass-172.5) < 40.){
               if ((W_mass > 40. && W_mass < 140.) && (Top_mass > 120. && Top_mass < 230.)){
                  h_Wmass_inWindow->Fill(W_mass, w);
                  prof_top_inWindow->Fill(run, Top_mass, w); 
                  prof_W_inWindow->Fill(run, W_mass, w);
                  h_Topmass_inWindow->Fill(Top_mass, w);
                  h_mbl_inWindow->Fill(mbl, w);
                  h_mbl_red_inWindow->Fill(mbl_red, w);
                  h_rbl_inWindow->Fill(rbl, w);
                  if (secondLightJet >= 0) {h_rbq_inWindow->Fill(rbq, w);}
                  
                  if (!heavyFallback && secondLightJet >= 0) {
                     if (fabs(lj1_p4.Eta()) < 1.3 && fabs(lj2_p4.Eta())  < 1.3) { //&& pt_pair > 53. && pt_pair < 107.){
                        h2_Wmass_inWindow_ptpair->Fill(pt_pair, W_mass, w);
                        prof_W_inWindow_ptpair->Fill(pt_pair, W_mass, w);
                        h2_Wmass_inWindow_ptpair_improved->Fill(pt_pair_improved, W_mass, w);
                        prof_W_inWindow_ptpair_improved->Fill(pt_pair_improved, W_mass, w);
                        h_ptpair_inWindow->Fill(pt_pair, w);
                        h_ptpair_inWindow_improved->Fill(pt_pair_improved, w);
                        if (isMC || (runYear != 2024)){
                           bool isG1 = (flav1 == 21), isG2 = (flav2 == 21);
                           if      (isG1 && isG2) {
                              h2_Wmass_inWindow_ptpair_gg->Fill(pt_pair, W_mass, w);
                              prof_W_inWindow_ptpair_gg->Fill(pt_pair, W_mass, w);
                              h_ptpair_inWindow_gg->Fill(pt_pair, w);
                              h2_Wmass_inWindow_ptpair_improved_gg->Fill(pt_pair_improved, W_mass, w);
                              prof_W_inWindow_ptpair_improved_gg->Fill(pt_pair_improved, W_mass, w);
                              h_ptpair_inWindow_improved_gg->Fill(pt_pair_improved, w);
                           }
                           else if (isG1 != isG2) {
                              h2_Wmass_inWindow_ptpair_qg->Fill(pt_pair, W_mass, w);
                              prof_W_inWindow_ptpair_qg->Fill(pt_pair, W_mass, w);
                              h_ptpair_inWindow_qg->Fill(pt_pair, w);
                              h2_Wmass_inWindow_ptpair_improved_qg->Fill(pt_pair_improved, W_mass, w);
                              prof_W_inWindow_ptpair_improved_qg->Fill(pt_pair_improved, W_mass, w);
                              h_ptpair_inWindow_improved_qg->Fill(pt_pair_improved, w);
                           }
                           else if (qqPairs(flav1, flav2)){
                              h2_Wmass_inWindow_ptpair_qq->Fill(pt_pair, W_mass, w);
                              prof_W_inWindow_ptpair_qq->Fill(pt_pair, W_mass, w);
                              h_ptpair_inWindow_qq->Fill(pt_pair, w);
                              h2_Wmass_inWindow_ptpair_improved_qq->Fill(pt_pair_improved, W_mass, w);
                              prof_W_inWindow_ptpair_improved_qq->Fill(pt_pair_improved, W_mass, w);
                              h_ptpair_inWindow_improved_qq->Fill(pt_pair_improved, w);
                           }
                           else {
                              h2_Wmass_inWindow_ptpair_qq_others->Fill(pt_pair, W_mass, w);
                              prof_W_inWindow_ptpair_qq_others->Fill(pt_pair, W_mass, w);
                              h_ptpair_inWindow_qq_others->Fill(pt_pair, w);
                              h2_Wmass_inWindow_ptpair_improved_qq_others->Fill(pt_pair_improved, W_mass, w);
                              prof_W_inWindow_ptpair_improved_qq_others->Fill(pt_pair_improved, W_mass, w);
                              h_ptpair_inWindow_improved_qq_others->Fill(pt_pair_improved, w);
                           }
                        }
                     }
                  }
               }

               // --- Diagonal window: Mass ≃ 2*ptpair + 4 (±10 GeV) and same Top_mass cut ---
               double center = 2.0 * pt_pair + 4.0;
               double halfWidth = 10.0;
               if (fabs(W_mass - center) < halfWidth && (Top_mass > 120. && Top_mass < 230.)) {
                   // Only fill when jets pass η cuts
                   if (!heavyFallback && secondLightJet >= 0 && fabs(lj1_p4.Eta()) < 1.3 && fabs(lj2_p4.Eta()) < 1.3) {
                       // Fill diagonal‐window histograms
                       h_Wmass_diagWindow->Fill(W_mass, w);
                       h_ptpair_diagWindow->Fill(pt_pair, w);
                       h2_Wmass_ptpair_diagWindow->Fill(pt_pair, W_mass, w);
                       prof_W_diagWindow_ptpair->Fill(pt_pair, W_mass, w);
                       prof_ptpair_diagWindow_mass->Fill(W_mass, pt_pair, w);
                   }
               }
               if ((W_mass > 40 && W_mass < 140.) && (Top_mass_improved > 120. && Top_mass_improved < 230.)){
                  h_Topmass_inWindow_improved->Fill(Top_mass_improved, w);
               }


               // 2 johtavaa kevyttä jettiä (exclude muon ja b jetit)
               // leptonID + isolaatio 
               // jettien IDt (b + light)
               // Loop over all jets in this event
            /*for (int i = 0; i != nJet; ++i) {
               p4jet.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
               double jet_pt = p4jet.Pt();
               // histograms
               h1_test->Fill(jet_pt, 1);
               // WMass 
               // topMass (2 hypoteesia) (spaNet? Mikko?)
            }*/
         } 
      } // jentry
         //cout << "Skipped " << _nbadevents_json << " events due to JSON ("
	      //<< (100.*_nbadevents_json/_nevents) << "%) \n";
      if (fmap) {
         fmap->cd();

         // --- inclusive ---
         p_partonPt_vs_partonPt->Write();
         p_genJetPt_vs_partonPt->Write();
         p_genJetPtOverPartonPt_vs_partonPt->Write();

         // --- inclusive FSR up/down ---
         p_partonPt_vs_partonPt_fsrUp->Write();
         p_partonPt_vs_partonPt_fsrDown->Write();
         p_genJetPt_vs_partonPt_fsrUp->Write();
         p_genJetPt_vs_partonPt_fsrDown->Write();
         p_genJetPtOverPartonPt_vs_partonPt_fsrUp->Write();
         p_genJetPtOverPartonPt_vs_partonPt_fsrDown->Write();

         p_genJetPtOverPartonPt_vs_partonPt_dlt0p5R->Write();
         p_genJetPtOverPartonPt_vs_partonPt_d0p5to2R->Write();
         p_genJetPtOverPartonPt_vs_partonPt_dgt2R->Write();

         // --- GenJet flavour splits (b / light) ---
         p_partonPt_vs_partonPt_b->Write();
         p_genJetPt_vs_partonPt_b->Write();
         p_genJetPtOverPartonPt_vs_partonPt_b->Write();

         p_partonPt_vs_partonPt_light->Write();
         p_genJetPt_vs_partonPt_light->Write();
         p_genJetPtOverPartonPt_vs_partonPt_light->Write();

         // --- GenJet flavour splits FSR up/down ---
         p_partonPt_vs_partonPt_b_fsrUp->Write();
         p_partonPt_vs_partonPt_b_fsrDown->Write();
         p_genJetPt_vs_partonPt_b_fsrUp->Write();
         p_genJetPt_vs_partonPt_b_fsrDown->Write();
         p_genJetPtOverPartonPt_vs_partonPt_b_fsrUp->Write();
         p_genJetPtOverPartonPt_vs_partonPt_b_fsrDown->Write();

         p_partonPt_vs_partonPt_light_fsrUp->Write();
         p_partonPt_vs_partonPt_light_fsrDown->Write();
         p_genJetPt_vs_partonPt_light_fsrUp->Write();
         p_genJetPt_vs_partonPt_light_fsrDown->Write();
         p_genJetPtOverPartonPt_vs_partonPt_light_fsrUp->Write();
         p_genJetPtOverPartonPt_vs_partonPt_light_fsrDown->Write();

         // --- selection-based splits (your analysis-selected jets) ---
         p_genJetPtOverPartonPt_vs_partonPt_selB->Write();
         p_genJetPtOverPartonPt_vs_partonPt_selLight->Write();

         // --- selection-based splits FSR up/down ---
         p_genJetPtOverPartonPt_vs_partonPt_selB_fsrUp->Write();
         p_genJetPtOverPartonPt_vs_partonPt_selB_fsrDown->Write();
         p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrUp->Write();
         p_genJetPtOverPartonPt_vs_partonPt_selLight_fsrDown->Write();

         fmap->Close();
      }   
      // Write run vs BX map and close file
      if (frunbx) {
         frunbx->cd();
         if (h2runbx) h2runbx->Write();
         frunbx->Close();
      }
      std::cout << "\n===== Events after run split =====\n";
for (const auto &kv : runCountAfterSplit) {
   std::cout << "run " << kv.first << " : " << kv.second << std::endl;
}
std::cout << "==================================\n";
   fout->cd();
   fout->Write();
   fout->Close();
   std::cout << "\n================ Summary of Jet Criteria =================" << std::endl;
   std::cout << "Total processed events: " << _nevents << std::endl;
   std::cout << "Events with >= 3 jets (pT > 30, |eta| < 2.4): " << _nPass_3jets
             << " (" << 100.0 * _nPass_3jets / _nevents << "%)" << std::endl;
   std::cout << "Events with >= 1 b-tagged jet: " << _nPass_1btag
             << " (" << 100.0 * _nPass_1btag / _nevents << "%)" << std::endl;
   std::cout << "Events with >= 3 jets and >= 1 b-tag: " << _nPass_3jets_1btag
             << " (" << 100.0 * _nPass_3jets_1btag / _nevents << "%)" << std::endl;
   std::cout << "=========================================================" << std::endl;
   exit(0);
} //WMassRun3 Loop

// Code originally from jetphys/HistosFill.C
void WMassRun3::PrintInfo(string info, bool printcout) //(.h!)
{
  //*ferr << info << endl << flush;
  if (printcout) cout << info << endl << flush;
}

// Code originally from jetphys/HistosFill.C
bool WMassRun3::LoadJSON(string json)
{
  PrintInfo(string("Processing LoadJSON() with ") + json + " ...",true);
  ifstream file(json, ios::in);
  if (!file.is_open()) { assert(false); return false; }
  char c;
  string s, s2, s3;
  char s1[256];
  int rn(0), ls1(0), ls2(0), nrun(0), nls(0);
  file.get(c);
  if (c!='{') return false;
  while (file >> s and sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    if (_gh_debug) PrintInfo(Form("\"%d\": ",rn),true);

    while (file.get(c) and c==' ') {};
    if (_gh_debug) { PrintInfo(Form("%c",c),true); assert(c=='['); }
    ++nrun;

    bool endrun = false;
    while (!endrun and file >> s >> s2 and (sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3 or (file >> s3 and sscanf((s+s2+s3).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3))) {

      s2 = s1;
      if (s2=="]") { file >> s3; s2 += s3; }

      if (_gh_debug) PrintInfo(Form("[%d,%d,'%s']",ls1,ls2,s1),true);

      for (int ls = ls1; ls != ls2+1; ++ls) {
        _json[rn][ls] = 1;
        ++nls;
      }

      endrun = (s2=="]," || s2=="]}");
      if (_gh_debug and !endrun and s2!=",") { PrintInfo(string("s1: ")+s2,true); assert(s2==","); }
    } // while ls
    if (_gh_debug) PrintInfo("",true);

    if (s2=="]}") continue;
    else if (_gh_debug and s2!="],") PrintInfo(string("s2: ")+s2,true); assert(s2=="],");
  } // while run
  if (s2!="]}") { PrintInfo(string("s3: ")+s2,true); return false; }

  PrintInfo(string("Called LoadJSON() with ") + json + ":",true);
  PrintInfo(Form("Loaded %d good runs and %d good lumi sections",nrun,nls),true);
  return true;
} // LoadJSON