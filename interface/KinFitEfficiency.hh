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
  KinFitEfficiency(KinFitOutputModule signal, double iSignalCrossSection, double iSignalLuminosity);
  void addBackground(KinFitOutputModule background, double crossSection, double luminosity);

  void run(double peakMass, double lowerWidth, double upperWidth);

private:
  std::vector<Particles> signalEvents;
  std::vector<Particles> fittedSignalEvents;
  double signalCrossSection;
  double signalLuminosity;
  std::vector<std::vector<Particles>> backgroundEvents;
  std::vector<std::vector<Particles>> fittedBackgroundEvents;
  std::vector<double> backgroundCrossSections;
  std::vector<double> backgroundLuminosities;


  void calculateRatio(bool fit, double peakMass, double lowerWidth, double upperWidth);
  std::pair<double, double> countPassedEvents(std::vector<Particles> events, double weight, double peakMass, double lowerWidth, double upperWidth);
};

#endif