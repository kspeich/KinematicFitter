#ifndef KINFITEFFICIENCY_HH
#define KINFITEFFICIENCY_HH

#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"
#include "PhysicsTools/KinFitter/interface/Particle.hh"
#include "PhysicsTools/KinFitter/interface/Particles.hh"

#include "TLorentzVector.h"

#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <vector>

class KinFitEfficiency
{
public:
  KinFitEfficiency(KinFitOutputModule signal, double iSignalCrossSection, double iPeakMass, double iLowerWidth, double iUpperWidth, double iLuminosity);
  void addBackground(KinFitOutputModule background, double crossSection);

  void run();

private:
  std::vector<Particles> signalEvents;
  std::vector<Particles> fittedSignalEvents;
  double signalCrossSection;
  std::vector<std::vector<Particles>> backgroundEvents;
  std::vector<std::vector<Particles>> fittedBackgroundEvents;
  std::vector<double> backgroundCrossSections;

  double peakMass;
  double lowerWidth;
  double upperWidth;

  double luminosity; 

  void calculateRatio(bool fit);
  std::pair<double, double> countPassedEvents(std::vector<Particles> events, double weight);
};

#endif