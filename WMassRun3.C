#define WMassRun3_cxx
#include "WMassRun3.h"
#include <TH2.h>
#include <TH3.h>
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
#include <algorithm>
#include <TLegend.h>
#include <TColor.h>
#include <TStopwatch.h>
#include "tdrstyle_mod22.C"

bool _gh_debug = false;




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
  FactorizedJetCorrector *jec = new FactorizedJetCorrector(v);

  return jec;
} // getJFC 


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
   fChain->SetBranchStatus("Jet_jetId", 1);
   //fChain->SetBranchStatus("",1);
   //fChain->SetBranchStatus("",1);

   TDirectory *curdir = gDirectory;
   // Create the output file based on a condition
   TFile* fout;

   fout = new TFile("output.root", "RECREATE");

   TH1::SetDefaultSumw2();

   TH1D *h1_test = new TH1D("h1_test", ";pt;N", 100,0,400);
   TH1D* h_Wmass = new TH1D("h_Wmass", "; W Mass (GeV); Events", 80, 0, 400);
   TH1D* h_Top1mass = new TH1D("h_Top1mass", "; Top Mass b1 (GeV); Events", 80, 0, 400);
   TH1D* h_Top2mass = new TH1D("h_Top2mass", "; Top Mass b2 (GeV); Events", 80, 0, 400);
   TH1D* h_Topmass = new TH1D("h_Topmass", "; Top Mass (GeV); Events", 80, 0, 400);
   TH1D* h_Wmass_inWindow = new TH1D("h_Wmass_inWindow", "; W Mass in Window (GeV); Events", 80, 0, 400);
   TH1D* h_Topmass_inWindow = new TH1D("h_Topmass_inWindow", "; Top Mass in Window (GeV); Events", 80, 0, 400);

   TH1D* h_mbl1 = new TH1D("h_mbl1", ";; Events", 80, 0, 400);
   TH1D* h_mbl2 = new TH1D("h_mbl2", ";; Events", 80, 0, 400);
   TH1D* h_mbl = new TH1D("h_mbl", ";; Events", 80, 0, 400);
   TH1D* h_mbl_red = new TH1D("h_mbl_red", ";; Events", 100, 0, 5);
   TH1D* h_mbl_inWindow = new TH1D("h_mbl_inWindow", ";; Events", 80, 0, 400);
   TH1D* h_rbq = new TH1D("h_rbq", ";; Events", 100, 0, 5);
   TH1D* h_rbl = new TH1D("h_rbl", ";; Events", 100, 0, 5);


   curdir->cd();
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

   TLorentzVector p4jet, lhe, jet, jet2, jetn;
   TLorentzVector gamorig; // for QCD bkg
   TLorentzVector met, met1, metn, metu, metnu, rawmet, corrmet, rawgam;
   TLorentzVector jeti, corrjets, rawjet, rawjets, rcjet, rcjets, rcoffsets;
   TLorentzVector geni, genjet, genjet2;
   TLorentzVector fox; // for isQCD
   TLorentzVector lj1_p4, lj2_p4, W_p4, top1_p4, top2_p4, b1_p4, b2_p4, muon_p4, b1l_p4, b2l_p4, toph_p4, mbl_p4, bl_p4, bh_p4;

   int _ntot(0), _nevents(0), _nbadevents_json(0), _nbadevents_trigger(0);
   int _nbadevents_veto(0);


   
   // JSON valinta (Photon + jet code) (Bettina)

   // Load JSON files
   LoadJSON("Cert_Collisions2024_378981_386951_Golden.json");

    // Jetti korjaukset + tyyppi 1 Met korjaukset (JEC) //Winter24Prompt24 (Bettina)

    FactorizedJetCorrector *jec(0);
    jec = getFJC("", "Winter24Prompt24_RunG_V2_DATA_L2Relative_AK4PFPuppi", "Winter24Prompt24_RunG_V2_DATA_L2L3Residual_AK4PFPuppi");
    assert(jec);

   
    TFile *fjv(0);
    fjv = new TFile("jet_veto_maps/Winter24Prompt24/Winter24Prompt24_2024BCDEFGHI.root","READ");
    if (!fjv) cout << "Jetvetomap file not found for 2024G"  << endl << flush;
    assert(fjv);

    TH2D *h2jv = 0;
    h2jv = (TH2D*)fjv->Get("jetvetomap");
    if (!h2jv) cout << "Jetvetomap histo not found for 2024G" << endl << flush;
    assert(h2jv);

   std::vector<int> goodMuonIndices;
   goodMuonIndices.reserve(nMuon);

   std::vector<int> bjetIndices;
   bjetIndices.reserve(nJet);

   std::vector<int> lightJetIndices;
   lightJetIndices.reserve(nJet);

   std::vector<int> gluonJetIndices;
   gluonJetIndices.reserve(nJet);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%nlap==0) {
         cout << "." << flush;
      }
      if (jentry%nlap2==0 && jentry!=0) {
         double time = t.RealTime();
      if (time>0) cout << Form("\n\%1.0f ev/s\n",nlap2/time) << flush;
         t.Reset();
         t.Start();
      }

      // Trigger IsoMu24
      bool pass_trig = (HLT_IsoMu24);

         // Does the run/LS pass the latest JSON selection?
         if (_json[run][luminosityBlock]==0) {
            //_badjson.insert(pair<int, int>(run, lbn));
            ++_nbadevents_json;
            continue;
         }
         else 
            ++_nevents;

         // Select leading jets. Just exclude muon, don't apply JetID yet
         static const int nJetMax = 200;
         Float_t         Jet_resFactor[nJetMax]; // Custom addition
         Float_t         Jet_deltaJES[nJetMax]; // Custom addition
         Float_t         Jet_CF[nJetMax]; // Custom addition
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
            if (jec!=0) {

               double rawJetPt = Jet_pt[i] * (1.0 - Jet_rawFactor[i]);
               double rawJetMass = Jet_mass[i] * (1.0 - Jet_rawFactor[i]);
               jec->setJetPt(rawJetPt);
               jec->setJetEta(Jet_eta[i]);
               jec->setJetPhi(Jet_phi[i]);
               //double corr = jec->getCorrection();
               vector<float> v = jec->getSubCorrections();
               double corr = v.back();
               double res = (v.size()>1 ? v[v.size()-1]/v[v.size()-2] : 1.);
               //Jet_RES[i] = 1./res;
               Jet_deltaJES[i] = (1./corr) / (1.0 - Jet_rawFactor[i]);
               Jet_pt[i] = corr * rawJetPt;
               Jet_mass[i] = corr * rawJetMass;
               Jet_rawFactor[i] = (1.0 - 1.0/corr);
               Jet_resFactor[i] = (1.0 - 1.0/res);
            }
         } // for i in nJet 
      
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
   

         //JetVetoMap //Summer22EE_23Sep2023
         bool pass_jetveto = true;
         if (true) { // jet veto
            int i1 = h2jv->GetXaxis()->FindBin(jet.Eta());
            int j1 = h2jv->GetYaxis()->FindBin(jet.Phi());
            if (h2jv->GetBinContent(i1,j1)>0) {
               ++_nbadevents_veto;
               pass_jetveto = false;
            }
         } // jet veto
         // Met leikkaus?
         rawmet.SetPtEtaPhiM(RawPuppiMET_pt, 0, RawPuppiMET_phi, 0);

         bool pass_basic = (pass_trig && pass_filt && pass_jetveto);


        
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
                  if (passMuonIsoId && passMuonId && (mu_pt > 30.) && (fabs(mu_eta) < 2.4)){
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

               // Example threshold for "looes" btagThreshold
               double btagThreshold = 0.4319; //0.0849; //0.8482; 

               for (int j = 0; j < nJet; j++) {
                  // Exclude if matched to a selected muon
                  //if (Jet_muonIdx1[j] == goodMuonIndices[0]) continue;
                  //if (Jet_muonIdx2[j] == goodMuonIndices[0]) continue;
                  if (j == Muon_jetIdx[goodMuonIndices[0]]) continue;
                 
                  // If this jet passes the b-tag threshold => b-jet
                  if ((Jet_btagUParTAK4B[j] > btagThreshold) && (Jet_pt[j]  > 30.) && (fabs(Jet_eta[j]) < 2.4)) {
                     bjetIndices.push_back(j);
                  } 
                  else if (lightJetIndices.size() < 2 && (Jet_pt[j]  > 30.) && (fabs(Jet_eta[j]) < 2.4)){
                     //bool passesPt = (Jet_pt[j]  > 15.);
                     //bool passesEta  = (fabs(Jet_eta[j]) < 2.4); //myös blle^ > 30GeV jos pt alle 30GeV mutta b -> ei light!
                     // otherwise, treat it as a "light jet" 
                     lightJetIndices.push_back(j); //
                  }
                  else if ((Jet_pt[j]  > 30.)){
                     gluonJetIndices.push_back(j);
                  }
               }
               // 3) Now pick the leading two from each set (if available)
               if (bjetIndices.size() != 2) {
               // Not enough b-jets => skip event or do something else
               continue;
               }
               if (lightJetIndices.size() < 2) {
               // Not enough light jets => skip or do something else
               continue;
               }

               // The leading two b-jets
               int leadBJet   = bjetIndices[0];
               int secondBJet = bjetIndices[1];

               // The leading two light jets
               int leadLightJet   = lightJetIndices[0];
               int secondLightJet = lightJetIndices[1];

               // -------------------------------------------------
               // 4) Lastly, check Jet ID for the chosen jets
               // -------------------------------------------------
               // In NanoAOD: bit 2 = tight, bit 3 = tightLeptonVeto
               bool passJetID = true;

               //if ((Jet_jetId[j] < 4)) continue;

               // Macro-like function to test a single jet’s ID
               auto checkJetID = [&](int j) {
                  // If bit 2 is not set => fails tight
                  return ((Jet_jetId[j] & 4) != 0);
               };

               // Check the two b-jets
               if (!checkJetID(leadBJet))   passJetID = false;
               if (!checkJetID(secondBJet)) passJetID = false;

               // Check the two light jets
               if (!checkJetID(leadLightJet))   passJetID = false;
               if (!checkJetID(secondLightJet)) passJetID = false;

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

               double b1_pt = Jet_pt[leadBJet];
               double b2_pt = Jet_pt[secondBJet];
               double lj1_pt = Jet_pt[leadLightJet];
               double lj2_pt = Jet_pt[secondLightJet];

               // Build the 4-vectors for the two selected light jets
               lj1_p4.SetPtEtaPhiM(
                  Jet_pt[leadLightJet],
                  Jet_eta[leadLightJet],
                  Jet_phi[leadLightJet],
                  Jet_mass[leadLightJet]
               );

               lj2_p4.SetPtEtaPhiM(
                  Jet_pt[secondLightJet],
                  Jet_eta[secondLightJet],
                  Jet_phi[secondLightJet],
                  Jet_mass[secondLightJet]
               );

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



               top1_p4 = b1_p4 + lj1_p4 + lj2_p4;
               double Top1_mass = top1_p4.M();  // invariant mass in GeV

               top2_p4 = b2_p4 + lj1_p4 + lj2_p4;
               double Top2_mass = top2_p4.M();  // invariant mass in GeV

               // Sum the 4-vectors to get the W candidate
               W_p4 = lj1_p4 + lj2_p4;
               double W_mass = W_p4.M();  // invariant mass in GeV

               b1l_p4 = b1_p4 + muon_p4;
               double mbl1 =  b1l_p4.M();

               b2l_p4 = b2_p4 + muon_p4;
               double mbl2 =  b2l_p4.M();

               bh_p4 = (fabs(Top1_mass-172.5)< fabs(Top2_mass-172.5) ? b1_p4 : b2_p4);
               bl_p4 = (fabs(Top1_mass-172.5)< fabs(Top2_mass-172.5) ? b2_p4 : b1_p4);

               toph_p4 = bh_p4 + lj1_p4 + lj2_p4;
               double Top_mass = toph_p4.M();  // invariant mass in GeV

               mbl_p4 = bl_p4 + muon_p4;
               double mbl = mbl_p4.M();  // invariant mass in GeV

               double mbl_red = mbl/Top_mass;

               // Fill the histogram
               h_Wmass->Fill(W_mass);
               // histograms: top1 = b1 + lj1 + lj2, top2 = b2 + lj1+ lj2
               h_Top2mass->Fill(Top2_mass);
               h_Top1mass->Fill(Top1_mass);
               h_Topmass->Fill(Top_mass);
                // massb&lep = mbl
               h_mbl1->Fill(mbl1);
               h_mbl2->Fill(mbl2);
               h_mbl->Fill(mbl);
               h_mbl_red->Fill(mbl_red);

               double rbq = (b1_p4.Pt()+b2_p4.Pt())/(lj1_p4.Pt()+lj2_p4.Pt());
               double rbl = (b1_p4.Pt()+b2_p4.Pt())/(2*muon_p4.Pt());
               // pt: b1+b2/(lj1+lj2) = rbq
               h_rbq->Fill(rbq);
               // rbl: b1+b2/(2*l1)
               h_rbl->Fill(rbl);

               double top_mass = (fabs(Top1_mass-172.5)< fabs(Top2_mass-172.5) ? Top1_mass : Top2_mass);
               //if (fabs(W_mass-80.4) < 20. && fabs(top_mass-172.5) < 40.){
               if ((W_mass > 60 && W_mass < 110) && (top_mass > 120 && top_mass < 230)){
                  h_Wmass_inWindow->Fill(W_mass);
                  h_Topmass_inWindow->Fill(top_mass);
                  h_mbl_inWindow->Fill(mbl);
               }

               // 2 johtavaa kevyttä jettiä (exclude muon ja b jetit)
               // leptonID + isolaatio 
               // jettien IDt (b + light)
               // Loop over all jets in this event
            for (int i = 0; i != nJet; ++i) {
               p4jet.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
               double jet_pt = p4jet.Pt();
               // histograms
               h1_test->Fill(jet_pt, 1);
               // WMass 
               // topMass (2 hypoteesia) (spaNet? Mikko?)
            }
         } 
      } // jentry
         //cout << "Skipped " << _nbadevents_json << " events due to JSON ("
	      //<< (100.*_nbadevents_json/_nevents) << "%) \n";

   fout->Write();
   fout->Close();
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