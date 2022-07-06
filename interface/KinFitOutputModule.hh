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
  KinFitOutputModule(TTree* iTree) {tree = iTree;};
  ~KinFitOutputModule() {};
  
  void run();

private:
  TTree* tree;
  std::vector<TH1F*> histograms;

  Double_t ErrEt(Float_t Et, Float_t Eta);
  Double_t ErrEta(Float_t Et, Float_t Eta);
  Double_t ErrPhi(Float_t Et, Float_t Eta);
  
  Float_t CalculatePt(Float_t Et, Float_t eta, Float_t phi, Float_t m);

  void print(TKinFitter *fitter);

  std::vector<const TMatrixD*> fitEvent(std::vector<TLorentzVector> particleVectors);

  void makeHistograms();
  void drawHistograms();
};

#endif