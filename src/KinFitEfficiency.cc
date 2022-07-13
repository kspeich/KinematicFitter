#include "PhysicsTools/KinFitter/interface/KinFitEfficiency.hh"

#include <iostream>
#include <vector>
#include <cmath>

KinFitEfficiency::KinFitEfficiency(TTree* signal, double iSignalCrossSection, double iPeakMass, double iLowerWidth, double iUpperWidth, double iLuminosity) :
  signalCrossSection(iSignalCrossSection),
  peakMass(iPeakMass),
  lowerWidth(iLowerWidth),
  upperWidth(iUpperWidth),
  luminosity(iLuminosity)
{
  auto signalKinFitOutputMod = KinFitOutputModule(signal);
  signalEvents = signalKinFitOutputMod.getUnfittedEvents();
  fittedSignalEvents = signalKinFitOutputMod.getFittedEvents();
}

void KinFitEfficiency::addBackground(TTree* background, double crossSection)
{
  auto backgroundKinFitOutputMod = KinFitOutputModule(background);
  backgroundEvents.push_back(backgroundKinFitOutputMod.getUnfittedEvents());
  fittedBackgroundEvents.push_back(backgroundKinFitOutputMod.getFittedEvents());
  backgroundCrossSections.push_back(crossSection);
}

void KinFitEfficiency::run()
{
  std::cout << "\nPseudoscalar Mass: " << peakMass << " GeV; Reconstruction from " << peakMass - lowerWidth << " GeV to " << peakMass + upperWidth << " GeV\n";
  std::cout << "\nBefore kinematic fit: \n";
  calculateRatio(false);
  std::cout << "\nAfter kinematic fit: \n";
  calculateRatio(true);
}

void KinFitEfficiency::calculateRatio(bool fit)
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

  auto signalWeight = signalCrossSection * luminosity / signal.size();
  auto signalCounts = countPassedEvents(signal, signalWeight);
  
  signalTauTau += signalCounts.first;
  signalBB += signalCounts.second;

  for (unsigned long int i = 0; i < backgrounds.size(); i++)
  {
    auto backgroundWeight = backgroundCrossSections[i] * luminosity / backgrounds[i].size();
    auto backgroundCounts = countPassedEvents(backgrounds[i], backgroundWeight);

    backgroundTauTau += backgroundCounts.first;
    backgroundBB += backgroundCounts.second;
  }

  // std::cout << "Tau Tau Signal: " << signalTauTau << "\t\tBB Signal: " << signalBB << std::endl;
  // std::cout << "Tau Tau Background: " << backgroundTauTau << "\t\tBB Background: " << backgroundBB << std::endl;

  std::cout << "Tau Tau S/B ratio: " << (signalTauTau / backgroundTauTau) << std::endl;
  std::cout << "BB S/B ratio: " << (signalBB / backgroundBB) << std::endl;
  std::cout << "Tau Tau S/sqrt(S+B) ratio: " << (signalTauTau / pow(signalTauTau + backgroundTauTau, 0.5)) << std::endl;
  std::cout << "BB S/sqrt(S+B) ratio: " << (signalBB / pow(signalBB + backgroundBB, 0.5)) << std::endl;
}

std::pair<double, double> KinFitEfficiency::countPassedEvents(std::vector<Particles> events, double weight)
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