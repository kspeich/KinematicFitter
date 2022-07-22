#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"
#include "PhysicsTools/KinFitter/interface/SVFitKinFitOutputModule.hh"
#include "PhysicsTools/KinFitter/interface/KinFitEfficiency.hh"

#include "TLorentzVector.h"

#include "TH1.h"

#include <iostream>

void testKinFitFWLite()
{
  // Read the signal files
  TFile *recoSignal = TFile::Open("InputFiles/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_Skim.root");
  TTree *recoSignalTree = (TTree *) recoSignal->Get("mutau_tree");
  TFile *recoSVFitSignal = TFile::Open("InputFiles/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_Skim_svfitted.root");
  TTree *recoSVFitSignalTree = (TTree *) recoSVFitSignal->Get("mutau_tree");
  double signalCrossSection = 48.61 * 0.1;

  // Open Background Files
  TFile *recoDYJets = TFile::Open("InputFiles/DYJetsToLL_Skim.root");
  TTree *recoDYJetsTree = (TTree *) recoDYJets->Get("mutau_tree");
  TFile *recoTTLeptonic = TFile::Open("InputFiles/TTTo2L2Nu_Skim.root");
  TTree *recoTTLeptonicTree = (TTree *) recoTTLeptonic->Get("mutau_tree");
  TFile *recoTTSemiLeptonic = TFile::Open("InputFiles/TTToSemiLeptonic_Skim.root");
  TTree *recoTTSemiLeptonicTree = (TTree *) recoTTSemiLeptonic->Get("mutau_tree");
  TFile *recoSVFitDYJets = TFile::Open("InputFiles/DYJetsToLL_Skim_svfitted.root");
  TTree *recoSVFitDYJetsTree = (TTree *) recoSVFitDYJets->Get("mutau_tree");
  TFile *recoSVFitTTLeptonic = TFile::Open("InputFiles/TTTo2L2Nu_Skim_svfitted.root");
  TTree *recoSVFitTTLeptonicTree = (TTree *) recoSVFitTTLeptonic->Get("mutau_tree");
  TFile *recoSVFitTTSemiLeptonic = TFile::Open("InputFiles/TTToSemiLeptonic_Skim_svfitted.root");
  TTree *recoSVFitTTSemiLeptonicTree = (TTree *) recoSVFitTTSemiLeptonic->Get("mutau_tree");
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

  auto signalSVFitKinFitOutputMod = SVFitKinFitOutputModule(recoSVFitSignalTree, "RecoSVFitSignalKinFitHistograms.root");
  signalSVFitKinFitOutputMod.run();
  auto dyJetsSVFitKinFitOutputMod = SVFitKinFitOutputModule(recoSVFitDYJetsTree, "RecoSVFitDYJetsKinFitHistograms.root");
  dyJetsSVFitKinFitOutputMod.run();
  auto ttLeptonicSVFitKinFitOutputMod = SVFitKinFitOutputModule(recoSVFitTTLeptonicTree, "RecoSVFitTTLeptonicKinFitHistograms.root");
  ttLeptonicSVFitKinFitOutputMod.run();
  auto ttSemiLeptonicSVFitKinFitOutputMod = SVFitKinFitOutputModule(recoSVFitTTSemiLeptonicTree, "RecoSVFitTTSemiLeptonicKinFitHistograms.root");
  ttSemiLeptonicSVFitKinFitOutputMod.run();


  std::cout << "Calculating S/B and S/sqrt(S+B) ratios...\n";
  double lumi = 59.74;
  std::cout << "Reco values: \n";
  auto kinFitEfficiency = KinFitEfficiency(signalKinFitOutputMod, signalCrossSection, lumi);
  kinFitEfficiency.addBackground(dyJetsKinFitOutputMod, dyJetsCrossSection, lumi);
  kinFitEfficiency.addBackground(ttLeptonicKinFitOutputMod, ttLeptonicCrossSection, lumi);
  kinFitEfficiency.addBackground(ttSemiLeptonicKinFitOutputMod, ttSemiLeptonicCrossSection, lumi);
  kinFitEfficiency.run(45., 5., 5.);
  kinFitEfficiency.run(45., 3., 3.);

  std::cout << "SVFit Reco values: \n";
  auto svFitKinFitEfficiency = KinFitEfficiency(signalSVFitKinFitOutputMod, signalCrossSection, lumi);
  svFitKinFitEfficiency.addBackground(dyJetsSVFitKinFitOutputMod, dyJetsCrossSection, lumi);
  svFitKinFitEfficiency.addBackground(ttLeptonicSVFitKinFitOutputMod, ttLeptonicCrossSection, lumi);
  svFitKinFitEfficiency.addBackground(ttSemiLeptonicSVFitKinFitOutputMod, ttSemiLeptonicCrossSection, lumi);
  svFitKinFitEfficiency.run(45., 5., 5.);
  svFitKinFitEfficiency.run(45., 3., 3.);
}