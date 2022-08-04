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

void KinFitEfficiency::run(double peakMass, double lowerWidth, double upperWidth)
{
  std::cout << "\nPseudoscalar Mass: " << peakMass << " GeV; Reconstruction from " << peakMass - lowerWidth << " GeV to " << peakMass + upperWidth << " GeV\n";
  std::cout << "\nBefore kinematic fit: \n";
  calculateRatio(false, peakMass, lowerWidth, upperWidth);
  std::cout << "\nAfter kinematic fit: \n";
  calculateRatio(true, peakMass, lowerWidth, upperWidth);
}

void KinFitEfficiency::calculateRatio(bool fit, double peakMass, double lowerWidth, double upperWidth)
{
  double signalTauTau = 0;
  double signalBB = 0;
  double backgroundTauTau = 0;
  double backgroundBB = 0;

  // Pre-fit values
  std::vector<Particles> signal = signalEvents;
  std::vector<std::vector<Particles>> backgrounds = backgroundEvents;

  if (fit)    // Post-fit
  {
    signal = fittedSignalEvents;
    backgrounds = fittedBackgroundEvents;
  }

  auto signalWeight = signalCrossSection * signalLuminosity / signal.size();
  auto signalCounts = countPassedEvents(signal, signalWeight, peakMass, lowerWidth, upperWidth);
  
  signalTauTau += signalCounts.first;
  signalBB += signalCounts.second;

  for (unsigned long int i = 0; i < backgrounds.size(); i++)
  {
    auto backgroundWeight = backgroundCrossSections[i] * backgroundLuminosities[i] / backgrounds[i].size();
    auto backgroundCounts = countPassedEvents(backgrounds[i], backgroundWeight, peakMass, lowerWidth, upperWidth);

    backgroundTauTau += backgroundCounts.first;
    backgroundBB += backgroundCounts.second;
  }

  std::cout << "Tau Tau Signal Count: " << signalTauTau << "\t BB Signal Count: " << signalBB << std::endl;
  std::cout << "Tau Tau Background Count: " << backgroundTauTau << "\t BB Background Count: " << backgroundBB << std::endl;
  std::cout << "Tau Tau S/B ratio: " << (signalTauTau / backgroundTauTau) << std::endl;
  std::cout << "BB S/B ratio: " << (signalBB / backgroundBB) << std::endl;
  std::cout << "Tau Tau S/sqrt(S+B) ratio: " << (signalTauTau / pow(signalTauTau + backgroundTauTau, 0.5)) << std::endl;
  std::cout << "BB S/sqrt(S+B) ratio: " << (signalBB / pow(signalBB + backgroundBB, 0.5)) << std::endl;
}

std::pair<double, double> KinFitEfficiency::countPassedEvents(std::vector<Particles> events, double weight, double peakMass, double lowerWidth, double upperWidth)
{
  double tautauCount = 0;
  double bbCount = 0;

  for(auto event : events)
  {
    double tautauInvariantMass = event.getInvariantMass(15);
    double bbInvariantMass = event.getInvariantMass(5);

    if ((tautauInvariantMass <= (peakMass + upperWidth)) && (tautauInvariantMass >= (peakMass - lowerWidth)))
    {
      tautauCount += weight;
    }
    if ((bbInvariantMass <= (peakMass + upperWidth)) && (bbInvariantMass >= (peakMass - lowerWidth)))
    {
      bbCount += weight;
    }
  }

  return std::make_pair(tautauCount, bbCount);
}