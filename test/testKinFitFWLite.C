#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"

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
  //TFile *data = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_Skim.root");
  //TTree *tree = (TTree *) data->Get("mutau_tree");

  TFile *data = TFile::Open("root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18NanoAODv7/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/AE82D64D-12A8-B547-9A00-339213036849.root");
  TTree *tree = (TTree *) data->Get("Events");

  auto kinFitOutputMod = KinFitOutputModule(tree, true);
  kinFitOutputMod.run();
}