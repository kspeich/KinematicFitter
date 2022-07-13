#include "PhysicsTools/KinFitter/interface/Particles.hh"

#include <iostream>
#include <vector>
#include <cmath>

void Particles::addParticle(TLorentzVector particleVec, int pdgId)
{
  auto particle = Particle(particleVec, pdgId);
  particles.push_back(particle);
}

double Particles::getInvariantMass()
{
  TLorentzVector sum;

  for (auto particle : particles)
  {
    sum += particle.getFourVector();
  }
  
  return sum.M();
}

Particles Particles::getParticles(int pdgId)
{
  Particles specificParticles;

  for (auto particle : particles)
  {
    if (abs(particle.getPdgId()) == abs(pdgId))
    {
      specificParticles.addParticle(particle);
    }
  }

  return specificParticles;
}

std::vector<TLorentzVector> Particles::getParticleVectors()
{
  std::vector<TLorentzVector> particleVectors;

  for (auto particle : particles)
  {
    particleVectors.push_back(particle.getFourVector());
  }

  return particleVectors;
}

std::vector<int> Particles::getPdgIds()
{
  std::vector<int> pdgIds;

  for (auto particle : particles)
  {
    pdgIds.push_back(particle.getPdgId());
  }

  return pdgIds;
}