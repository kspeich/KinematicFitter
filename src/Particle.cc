#include "PhysicsTools/KinFitter/interface/Particle.hh"

#include <iostream>
#include <vector>

Particle::Particle(TLorentzVector iFourVector, int iPdgId) :
  fourVector(iFourVector),
  pdgId(iPdgId)
{
}

bool Particle::containsTag(std::string tag)
{
  bool containsTag = false;

  for (auto particleTag : tags)
  {
    if (tag == particleTag)
    {
      containsTag = true;
    }
  }

  return containsTag;
}