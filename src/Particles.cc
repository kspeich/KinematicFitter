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

Particles Particles::getParticles(int pdgId, std::string tag)
{
  Particles tagParticles;

  for (auto particle : getParticles(pdgId).getParticles())
  {
    if (particle.containsTag(tag))
    {
      tagParticles.addParticle(particle);
    }
  }

  return tagParticles;
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

Particle Particles::getLeading(int pdgId)
{
  auto relevantParticles = getParticles(pdgId).getParticles();
  sort(relevantParticles.begin(), relevantParticles.end(), [](Particle &p1, Particle &p2){return p1.Pt() > p2.Pt();});
  return relevantParticles[0];
}

Particle Particles::getNextToLeading(int pdgId)
{
  auto relevantParticles = getParticles(pdgId).getParticles();
  sort(relevantParticles.begin(), relevantParticles.end(), [](Particle &p1, Particle &p2){return p1.Pt() > p2.Pt();});
  return relevantParticles[1];
}