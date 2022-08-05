#ifndef PARTICLE_HH
#define PARTICLE_HH

#include "TLorentzVector.h"

#include <iostream>
#include <vector>

class Particle
{
public:
  Particle(TLorentzVector iFourVector, int iPdgId);
  Particle(TLorentzVector iFourVector, int iPdgId, std::vector<std::string> iTags);
  ~Particle() {};

  TLorentzVector getFourVector() {return fourVector;};
  int getPdgId() {return pdgId;};
  std::vector<std::string> getTags() {return tags;};

  void addTag(std::string tag) {tags.push_back(tag);};
  bool containsTag(std::string tag);

  double Pt() {return fourVector.Pt();};
  double Et() {return fourVector.Et();};
  double Eta() {return fourVector.Eta();};
  double Phi() {return fourVector.Phi();};
  double M() {return fourVector.M();};
  double Theta() {return fourVector.Theta();};
  double E() {return fourVector.E();};
  double P() {return fourVector.P();};

private:
  TLorentzVector fourVector;
  int pdgId;
  std::vector<std::string> tags;
};

#endif