#ifndef KINFITOUTPUTMODULE_HH
#define KINFITOUTPUTMODULE_HH

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

class KinFitOutputModule
{
public:
  KinFitOutputModule(TTree* iTree, std::string iOutputFile);
  virtual ~KinFitOutputModule() {};
  
  std::vector<Particles> getUnfittedEvents() {return unfittedEvents;};
  std::vector<Particles> getFittedEvents() {return fittedEvents;};

  virtual void run();

protected:
  TTree* tree;

  std::vector<Particles> unfittedEvents;
  std::vector<Particles> fittedEvents;
  std::vector<TH1F*> histograms;

  std::string outputFile;

  Double_t ErrEt(TLorentzVector particleVec);
  Double_t ErrEta(TLorentzVector particleVec);
  Double_t ErrPhi(TLorentzVector particleVec);
  
  Float_t calculatePt(Float_t Et, Float_t eta, Float_t phi, Float_t m);

  void print(TKinFitter *fitter);

  TFitParticleEtEtaPhi* convertParticle(Particle particle);

  virtual Particles fitEvent(Particles event);
  virtual void runFitter();
  void fillHistograms(Particles event, TH1F* hEt, TH1F* hEta, TH1F* hPhi, TH1F* hTauTauInvMass, TH1F* hBBInvMass);
  void fillHistogramOneParticle(Particle particle, std::string kinematic, TH1F* h);
  void fillKinematicHistogramsByLeg(Particles event, std::string kinematic, TH1F* hLeadingB, TH1F* hNTLB, TH1F* hMTau, TH1F* hHTau);

  void makeHistograms();
  void drawHistograms();
};

#endif