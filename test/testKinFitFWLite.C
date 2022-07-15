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
  TFile *genSignal = TFile::Open("signalSelectedEvents.root");
  TFile *recoSignal = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_Skim.root");
  TTree *genSignalTree = (TTree *) genSignal->Get("mutau_tree");
  TTree *recoSignalTree = (TTree *) recoSignal->Get("mutau_tree")
  double signalCrossSection = 48.61 * 0.1;

  // Open Background Files
  TFile *genDYJets = TFile::Open("dyJetsSelectedEvents.root");
  TFile *recoDYJets = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/DYJetsToLL_Skim.root");
  TTree *genDYJetsTree = (TTree *) genDYJets->Get("mutau_tree");
  TTree *recoDYJetsTree = (TTree *) recoDYJets->Get("mutau_tree");
  double dyJetsCrossSection = 5343.0;
  TFile *genTTLeptonic = TFile::Open("ttLeptonicSelectedEvents.root");
  TFile *recoTTLeptonic = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/TTTo2L2Nu_Skim.root");
  TTree *genTTLeptonicTree = (TTree *) genTTLeptonic->Get("mutau_tree");
  TTree *recoTTLeptonicTree = (TTree *) recoTTLeptonic->Get("mutau_tree");
  double ttLeptonicCrossSection = 88.29;
  TFile *genTTSemiLeptonic = TFile::Open("ttSemiLeptonicSelectedEvents.root");
  TFile *recoTTSemiLeptonic = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/TTToSemiLeptonic_Skim.root");
  TTree *genTTSemiLeptonicTree = (TTree *) genTTSemiLeptonic->Get("mutau_tree");
  TTree *recoTTSemiLeptonicTree = (TTree *) recoTTSemiLeptonic->Get("mutau_tree");
  double ttSemiLeptonicCrossSection = 365.35;
  
  std::cout << "Creating KinFit Histograms...\n";
  auto genSignalKinFitOutputMod = KinFitOutputModule(genSignalTree, "GenSignalKinFitHistograms.root");
  auto recoSignalKinFitOutputMod = KinFitOutputModule(recoSignalTree, "RecoSignalKinFitHistograms.root");
  auto genDYJetsKinFitOutputMod = KinFitOutputModule(genDYJetsTree, "GenDYJetsKinFitHistograms.root");
  auto recoDYJetsKinFitOutputMod = KinFitOutputModule(recoDYJetsTree, "RecoDYJetsKinFitHistograms.root");
  auto genTTLeptonicKinFitOutputMod = KinFitOutputModule(genTTLeptonicTree, "GenTTLeptonicKinFitHistograms.root");
  auto genTTLeptonicKinFitOutputMod = KinFitOutputModule(recoTTLeptonicTree, "RecoTTLeptonicKinFitHistograms.root");
  auto genTTSemiLeptonicKinFitOutputMod = KinFitOutputModule(genTTSemiLeptonicTree, "GenTTSemiLeptonicKinFitHistograms.root");
  auto genTTSemiLeptonicKinFitOutputMod = KinFitOutputModule(recoTTSemiLeptonicTree, "RecoTTSemiLeptonicKinFitHistograms.root");

  std::cout << "Calculating S/B and S/sqrt(S+B) ratios...\n";
  double lumi = 59.74;
  auto genKinFitEfficiency = KinFitEfficiency(genSignalKinFitOutputMod, signalCrossSection, lumi);
  genKinFitEfficiency.addBackground(genDYJetsKinFitOutputMod, dyJetsCrossSection, lumi);
  genKinFitEfficiency.addBackground(genTTLeptonicKinFitOutputMod, ttLeptonicCrossSection, lumi);
  genKinFitEfficiency.addBackground(genTTSemiLeptonicKinFitOutputMod, ttSemiLeptonicCrossSection, lumi);
  genKinFitEfficiency.run(45., 5., 5.);
  genKinFitEfficiency.run(45., 3., 3.);

  auto recoKinFitEfficiency = KinFitEfficiency(recoSignalKinFitOutputMod, signalCrossSection, lumi);
  recoKinFitEfficiency.addBackground(recoDYJetsKinFitOutputMod, dyJetsCrossSection, lumi);
  recoKinFitEfficiency.addBackground(recoTTLeptonicKinFitOutputMod, ttLeptonicCrossSection, lumi);
  recoKinFitEfficiency.addBackground(recoTTSemiLeptonicKinFitOutputMod, ttSemiLeptonicCrossSection, lumi);
  recoKinFitEfficiency.run(45., 5., 5.);
  recoKinFitEfficiency.run(45., 3., 3.);
}