#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"
#include "PhysicsTools/KinFitter/interface/KinFitEfficiency.hh"

#include "TLorentzVector.h"

#include "TH1.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Modified version of the ExampleEtEtaPhi2CFit.C macro from CMS AN 2005/025.
// Jet error parametrization from CMS AN 2005/005.
//
// To run this macro in a root session do:
// root [0] gSystem->Load("libPhysicsToolsKinFitter.so");
// root [1] .x PhysicsTools/KinFitter/test/testKinFitFWLite.C+
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void testKinFitFWLite()
{
  // Read the signal files
  /*TFile *genSignal = TFile::Open("signalSelectedEvents.root");
  TTree *genSignalTree = (TTree *) genSignal->Get("mutau_tree");*/
  TFile *recoSignal = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_Skim.root");
  TTree *recoSignalTree = (TTree *) recoSignal->Get("mutau_tree");
  double signalCrossSection = 48.61 * 0.1;

  // Open Background Files
  /*TFile *genDYJets = TFile::Open("dyJetsSelectedEvents.root");
  TTree *genDYJetsTree = (TTree *) genDYJets->Get("mutau_tree");
  TFile *genTTLeptonic = TFile::Open("ttLeptonicSelectedEvents.root");
  TTree *genTTLeptonicTree = (TTree *) genTTLeptonic->Get("mutau_tree");
  TFile *genTTSemiLeptonic = TFile::Open("ttSemiLeptonicSelectedEvents.root");
  TTree *genTTSemiLeptonicTree = (TTree *) genTTSemiLeptonic->Get("mutau_tree");*/
  TFile *recoDYJets = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/DYJetsToLL_Skim.root");
  TTree *recoDYJetsTree = (TTree *) recoDYJets->Get("mutau_tree");
  TFile *recoTTLeptonic = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/TTTo2L2Nu_Skim.root");
  TTree *recoTTLeptonicTree = (TTree *) recoTTLeptonic->Get("mutau_tree");
  TFile *recoTTSemiLeptonic = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/TTToSemiLeptonic_Skim.root");
  TTree *recoTTSemiLeptonicTree = (TTree *) recoTTSemiLeptonic->Get("mutau_tree");
  double dyJetsCrossSection = 5343.0;
  double ttLeptonicCrossSection = 88.29;
  double ttSemiLeptonicCrossSection = 365.35;
  
  std::cout << "Creating KinFit Histograms...\n";
  auto signalKinFitOutputMod = KinFitOutputModule(recoSignalTree, "RecoSignalKinFitHistograms.root");
  signalKinFitOutputMod.run();
  auto dyJetsKinFitOutputMod = KinFitOutputModule(recoDYJetsTree, "RecoDYJetsKinFitHistograms.root");
  dyJetsKinFitOutputMod.run();
  auto ttLeptonicKinFitOutputMod = KinFitOutputModule(recoTTLeptonicTree, "RecoTTLeptonicKinFitHistograms.root");
  ttLeptonicKinFitOutputMod.run();
  auto ttSemiLeptonicKinFitOutputMod = KinFitOutputModule(recoTTSemiLeptonicTree, "RecoTTSemiLeptonicKinFitHistograms.root");
  ttSemiLeptonicKinFitOutputMod.run();

  std::cout << "Calculating S/B and S/sqrt(S+B) ratios...\n";
  double lumi = 59.74;
  auto kinFitEfficiency = KinFitEfficiency(signalKinFitOutputMod, signalCrossSection, lumi);
  kinFitEfficiency.addBackground(dyJetsKinFitOutputMod, dyJetsCrossSection, lumi);
  kinFitEfficiency.addBackground(ttLeptonicKinFitOutputMod, ttLeptonicCrossSection, lumi);
  kinFitEfficiency.addBackground(ttSemiLeptonicKinFitOutputMod, ttSemiLeptonicCrossSection, lumi);
  std::cout << "Reco events:\n";
  kinFitEfficiency.run(45., 5., 5.);
  kinFitEfficiency.run(45., 3., 3.);
}