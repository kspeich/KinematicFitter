#ifndef PARTICLES_HH
#define PARTICLES_HH

#include "TLorentzVector.h"

#include "PhysicsTools/KinFitter/interface/Particle.hh"

#include <iostream>
#include <vector>

class Particles
{
public:
  Particles() {};
  Particles(std::vector<Particle> iParticles) : particles(iParticles) {};
  virtual ~Particles() {};

  Particle operator[](const int i) {return particles[i];};

  void addParticle(Particle particle) {particles.push_back(particle);};
  void addParticle(TLorentzVector particle, int pdgId);

  std::vector<Particle> getParticles() {return particles;};
  Particles getParticles(int pdgId);
  std::vector<TLorentzVector> getParticleVectors();
  std::vector<TLorentzVector> getParticleVectors(int pdgId) {return getParticles(pdgId).getParticleVectors();};

  std::vector<int> getPdgIds();

  int getNumParticles() {return particles.size();};
  int getNumParticles(int pdgId) {return getParticles(pdgId).getNumParticles();};
  
  double getInvariantMass();
  double getInvariantMass(int pdgId) {return getParticles(pdgId).getInvariantMass();};

private:
  std::vector<Particle> particles;
};

#endif