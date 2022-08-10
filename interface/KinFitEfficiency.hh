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
  KinFitEfficiency(KinFitOutputModule signal, double iSignalCrossSection, double iSignalLuminosity, std::vector<KinFitOutputModule> backgrounds, std::vector<double> backgroundCrossSections, std::vector<double> backgroundLuminosities);
  void addBackground(KinFitOutputModule background, double crossSection, double luminosity);
  void addBackgrounds(std::vector<KinFitOutputModule> backgrounds, std::vector<double> backgroundCrossSections, std::vector<double> backgroundLuminosities);

  void run(double pseudoscalarMass, double pseudoscalarLowerWidth, double pseudoscalarUpperWidth, double higgsMass, double higgsLowerWidth, double higgsUpperWidth);

private:
  std::vector<Particles> signalEvents;
  std::vector<Particles> fittedSignalEvents;
  double signalCrossSection;
  double signalLuminosity;
  std::vector<std::vector<Particles>> backgroundEvents;
  std::vector<std::vector<Particles>> fittedBackgroundEvents;
  std::vector<double> backgroundCrossSections;
  std::vector<double> backgroundLuminosities;


  void calculateRatio(bool fit, double pseudoscalarMass, double pseudoscalarLowerWidth, double pseudoscalarUpperWidth, double higgsMass, double higgsLowerWidth, double higgsUpperWidth);
  std::vector<double> countPassedEvents(std::vector<Particles> events, double weight, double pseudoscalarMass, double pseudoscalarLowerWidth, double pseudoscalarUpperWidth, double higgsMass, double higgsLowerWidth, double higgsUpperWidth);
};

#endif