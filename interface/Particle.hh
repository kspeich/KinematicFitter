#ifndef PARTICLE_HH
#define PARTICLE_HH

#include "TLorentzVector.h"

#include <iostream>
#include <vector>

class Particle
{
public:
  Particle(TLorentzVector iFourVector, int iPdgId);

  TLorentzVector getFourVector() {return fourVector;};
  int getPdgId() {return pdgId;};

  double Pt() {return fourVector.Pt();};
  double Et() {return fourVector.Et();};
  double Eta() {return fourVector.Eta();};
  double Phi() {return fourVector.Phi();};
  double M() {return fourVector.M();};

private:
  TLorentzVector fourVector;
  int pdgId;
};

#endif