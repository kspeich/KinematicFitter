#include "PhysicsTools/KinFitter/interface/Particle.hh"

#include <iostream>
#include <vector>

Particle::Particle(TLorentzVector iFourVector, int iPdgId) :
  fourVector(iFourVector),
  pdgId(iPdgId)
{
}