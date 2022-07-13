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
  // Read the file
  // TFile *data = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_Skim.root");    // Reco File
  TFile *data = TFile::Open("selectedEvents.root");         // Gen file that has been processed by EventSelection.py
  TTree *tree = (TTree *) data->Get("mutau_tree");
  double signalCrossSection = 48.61 * 0.1;

  // Open Background Files
  TFile *dyJets = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/DYJetsToLL_Skim.root");
  TTree *dyJetsTree = (TTree *) dyJets->Get("mutau_tree");
  double dyJetsCrossSection = 5343.0;
  TFile *ttLeptonic = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/TTTo2L2Nu_Skim.root");
  TTree *ttLeptonicTree = (TTree *) ttLeptonic->Get("mutau_tree");
  double ttLeptonicCrossSection = 88.29;
  TFile *ttSemiLeptonic = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/TTToSemiLeptonic_Skim.root");
  TTree *ttSemiLeptonicTree = (TTree *) ttSemiLeptonic->Get("mutau_tree");
  double ttSemiLeptonicCrossSection = 365.35;
  
  std::cout << "Creating KinFit Histograms...\n";
  auto kinFitOutputMod = KinFitOutputModule(tree);
  kinFitOutputMod.run();

  /*std::cout << "Calculating S/B and S/sqrt(S+B) ratios...\n";
  auto kinFitEfficiency = KinFitEfficiency(tree, signalCrossSection, 45., 5., 5., 59.74);
  std::cout << "Adding DYJets Background\n";
  kinFitEfficiency.addBackground(dyJetsTree, dyJetsCrossSection);
  std::cout << "Adding TTLeptonic Background\n";
  kinFitEfficiency.addBackground(ttLeptonicTree, ttLeptonicCrossSection);
  std::cout << "Adding TTSemiLeptonic Background\n";
  kinFitEfficiency.addBackground(ttSemiLeptonicTree, ttSemiLeptonicCrossSection);
  std::cout << "Backgrounds added\n";
  kinFitEfficiency.run();*/
}