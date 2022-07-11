#ifndef KINFITOUTPUTMODULE_HH
#define KINFITOUTPUTMODULE_HH

#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

#include "TLorentzVector.h"

#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <vector>

class KinFitOutputModule
{
public:
  KinFitOutputModule(TTree* iTree);
  ~KinFitOutputModule() {};
  
  std::vector<std::vector<TLorentzVector>> getUnfittedEvents() {return unfittedEvents;};
  std::vector<std::vector<TLorentzVector>> getFittedEvents() {return fittedEvents;};

  void run();

private:
  TTree* tree;
  std::vector<std::vector<TLorentzVector>> unfittedEvents;
  std::vector<std::vector<TLorentzVector>> fittedEvents;
  std::vector<TH1F*> histograms;

  Double_t ErrEt(Float_t Et, Float_t Eta);
  Double_t ErrEta(Float_t Et, Float_t Eta);
  Double_t ErrPhi(Float_t Et, Float_t Eta);
  
  Float_t calculatePt(Float_t Et, Float_t eta, Float_t phi, Float_t m);

  void print(TKinFitter *fitter);

  std::vector<TLorentzVector> fitEvent(std::vector<TLorentzVector> particleVectors);
  void runFitter();
  void fillHistograms(std::vector<TLorentzVector> particleVectors, TH1F* hEt, TH1F* hEta, TH1F* hPhi, TH1F* hTauTauInvMass, TH1F* hBBInvMass);

  void makeHistograms();
  void drawHistograms();
};

#endif