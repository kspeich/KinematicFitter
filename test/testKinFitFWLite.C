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
  // TFile *data = TFile::Open("/afs/cern.ch/work/s/skkwan/public/forKodai/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_Skim.root");
  TFile *data = TFile::Open("selectedEvents.root");
  TTree *tree = (TTree *) data->Get("mutau_tree");

  auto kinFitOutputMod = KinFitOutputModule(tree);
  kinFitOutputMod.run();
}