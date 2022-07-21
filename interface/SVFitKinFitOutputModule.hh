#ifndef SVFITKINFITOUTPUTMODULE_HH
#define SVFITKINFITOUTPUTMODULE_HH

#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/Particle.hh"
#include "PhysicsTools/KinFitter/interface/Particles.hh"

#include "TLorentzVector.h"

#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <vector>

class SVFitKinFitOutputModule : public KinFitOutputModule
{
public:
  SVFitKinFitOutputModule(TTree* iTree, std::string iOutputFile);
  ~SVFitKinFitOutputModule() {};
  
  void run() override;

private:
  Particles fitEvent(Particles event) override;
  void runFitter() override;
};

#endif