#include "PhysicsTools/KinFitter/interface/KinFitEfficiency.hh"

#include <iostream>
#include <vector>
#include <cmath>

KinFitEfficiency::KinFitEfficiency(KinFitOutputModule signal, double iSignalCrossSection, double iSignalLuminosity) :
  signalCrossSection(iSignalCrossSection),
  signalLuminosity(iSignalLuminosity)
{
  signalEvents = signal.getUnfittedEvents();
  fittedSignalEvents = signal.getFittedEvents();
}

KinFitEfficiency::KinFitEfficiency(KinFitOutputModule signal, double iSignalCrossSection, double iSignalLuminosity, std::vector<KinFitOutputModule> backgrounds, std::vector<double> backgroundCrossSections, std::vector<double> backgroundLuminosities) :
  signalCrossSection(iSignalCrossSection),
  signalLuminosity(iSignalLuminosity)
{
  signalEvents = signal.getUnfittedEvents();
  fittedSignalEvents = signal.getFittedEvents();
  addBackgrounds(backgrounds, backgroundCrossSections, backgroundLuminosities);
}

void KinFitEfficiency::addBackground(KinFitOutputModule background, double crossSection, double luminosity)
{
  backgroundEvents.push_back(background.getUnfittedEvents());
  fittedBackgroundEvents.push_back(background.getFittedEvents());
  backgroundCrossSections.push_back(crossSection);
  backgroundLuminosities.push_back(luminosity);
}

void KinFitEfficiency::addBackgrounds(std::vector<KinFitOutputModule> backgrounds, std::vector<double> backgroundCrossSections, std::vector<double> backgroundLuminosities)
{
  if (backgrounds.size() == backgroundCrossSections.size() && backgrounds.size() == backgroundLuminosities.size())
  {
    for (unsigned long int i = 0; i < backgrounds.size(); i++)
    {
      addBackground(backgrounds[i], backgroundCrossSections[i], backgroundLuminosities[i]);
    }
  }
}

void KinFitEfficiency::run(double pseudoscalarMass, double pseudoscalarLowerWidth, double pseudoscalarUpperWidth, double higgsMass, double higgsLowerWidth, double higgsUpperWidth)
{
  std::cout << "\nPseudoscalar Mass: " << pseudoscalarMass << " GeV; Reconstruction from " << pseudoscalarMass - pseudoscalarLowerWidth << " GeV to " << pseudoscalarMass + pseudoscalarUpperWidth << " GeV\n";
  std::cout << "Higgs Mass: " << higgsMass << " GeV; Reconstruction from " << higgsMass - higgsLowerWidth << " GeV to " << higgsMass + higgsUpperWidth << " GeV\n";
  std::cout << "\nBefore kinematic fit: \n";
  calculateRatio(false, pseudoscalarMass, pseudoscalarLowerWidth, pseudoscalarUpperWidth, higgsMass, higgsLowerWidth, higgsUpperWidth);
  std::cout << "\nAfter kinematic fit: \n";
  calculateRatio(true, pseudoscalarMass, pseudoscalarLowerWidth, pseudoscalarUpperWidth, higgsMass, higgsLowerWidth, higgsUpperWidth);
}

void KinFitEfficiency::calculateRatio(bool fit, double pseudoscalarMass, double pseudoscalarLowerWidth, double pseudoscalarUpperWidth, double higgsMass, double higgsLowerWidth, double higgsUpperWidth)
{
  double signalTauTau = 0;
  double signalBB = 0;
  double backgroundTauTau = 0;
  double backgroundBB = 0;

  double signalHiggs = 0;
  double backgroundHiggs = 0;

  // Pre-fit values
  std::vector<Particles> signal = signalEvents;
  std::vector<std::vector<Particles>> backgrounds = backgroundEvents;

  if (fit)    // Post-fit
  {
    signal = fittedSignalEvents;
    backgrounds = fittedBackgroundEvents;
  }

  auto signalWeight = signalCrossSection * signalLuminosity / signal.size();
  auto signalCounts = countPassedEvents(signal, signalWeight, pseudoscalarMass, pseudoscalarLowerWidth, pseudoscalarUpperWidth, higgsMass, higgsLowerWidth, higgsUpperWidth);
  
  signalTauTau += signalCounts[0];
  signalBB += signalCounts[1];
  signalHiggs += signalCounts[2];

  for (unsigned long int i = 0; i < backgrounds.size(); i++)
  {
    auto backgroundWeight = backgroundCrossSections[i] * backgroundLuminosities[i] / backgrounds[i].size();
    auto backgroundCounts = countPassedEvents(backgrounds[i], backgroundWeight, pseudoscalarMass, pseudoscalarLowerWidth, pseudoscalarUpperWidth, higgsMass, higgsLowerWidth, higgsUpperWidth);

    backgroundTauTau += backgroundCounts[0];
    backgroundBB += backgroundCounts[1];
    backgroundHiggs += backgroundCounts[2];
  }

  std::cout << "Tau Tau Signal Count: " << signalTauTau << "\t BB Signal Count: " << signalBB << std::endl;
  std::cout << "Tau Tau Background Count: " << backgroundTauTau << "\t BB Background Count: " << backgroundBB << std::endl;
  std::cout << "Tau Tau S/B ratio: " << (signalTauTau / backgroundTauTau) << std::endl;
  std::cout << "BB S/B ratio: " << (signalBB / backgroundBB) << std::endl;
  std::cout << "Tau Tau S/sqrt(S+B) ratio: " << (signalTauTau / pow(signalTauTau + backgroundTauTau, 0.5)) << std::endl;
  std::cout << "BB S/sqrt(S+B) ratio: " << (signalBB / pow(signalBB + backgroundBB, 0.5)) << std::endl;
  std::cout << "Higgs Signal Count: " << signalHiggs << "\nHiggs Background Count: " << backgroundHiggs << std::endl;
  std::cout << "Higgs S/B ratio: " << (signalHiggs / backgroundHiggs) << std::endl;
  std::cout << "Higgs S/sqrt(S+B) ratio: " << (signalHiggs / pow(signalHiggs + backgroundHiggs, 0.5)) << std::endl;
}

std::vector<double> KinFitEfficiency::countPassedEvents(std::vector<Particles> events, double weight, double pseudoscalarMass, double pseudoscalarLowerWidth, double pseudoscalarUpperWidth, double higgsMass, double higgsLowerWidth, double higgsUpperWidth)
{
  double tautauCount = 0;
  double bbCount = 0;
  double higgsCount = 0;

  for(auto event : events)
  {
    double tautauInvariantMass = event.getInvariantMass(15);
    double bbInvariantMass = event.getInvariantMass(5);
    double allParticleInvariantMass = event.getInvariantMass();

    if ((tautauInvariantMass <= (pseudoscalarMass + pseudoscalarUpperWidth)) && (tautauInvariantMass >= (pseudoscalarMass - pseudoscalarLowerWidth)))
    {
      tautauCount += weight;
    }
    if ((bbInvariantMass <= (pseudoscalarMass + pseudoscalarUpperWidth)) && (bbInvariantMass >= (pseudoscalarMass - pseudoscalarLowerWidth)))
    {
      bbCount += weight;
    }
    if ((allParticleInvariantMass <= (higgsMass + higgsUpperWidth)) && (allParticleInvariantMass >= (higgsMass - higgsLowerWidth)))
    {
      higgsCount += weight;
    }
  }

  std::vector<double> counts = {tautauCount, bbCount, higgsCount};
  return counts;
}