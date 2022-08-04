#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include <iostream>
#include <string>
#include <algorithm>

void overlayPlots(std::string canvasName, std::string xAxis, std::string yAxis, std::vector<TFile*> files, std::vector<std::string> histNames, std::vector<std::string> legends, std::vector<Color_t> colors);
void plotKinematics(TFile* file, std::string process, bool svFit = false);
void plotInvariantMasses(TFile* signalKinFitHistograms, TFile* dyJetsKinFitHistograms, TFile* ttLeptonicKinFitHistograms, TFile* ttSemiLeptonicKinFitHistograms, bool svFit = false);
void compareInvariantMasses(TFile* signalKinFitHistograms, TFile* svFitSignalKinFitHistograms, TFile* dyJetsKinFitHistograms, TFile* svFitDYJetsKinFitHistograms, TFile* ttLeptonicKinFitHistograms, TFile* svFitTTLeptonicKinFitHistograms, TFile* ttSemiLeptonicKinFitHistograms, TFile* svFitTTSemiLeptonicKinFitHistograms, bool kinFit, int numBJets = 0);
void compareLegs(TFile* file, std::string process, bool svFit = false, bool kinFit = false);

void plotOutput()
{
  gStyle->SetOptStat(0);

  TFile* recoSignal = TFile::Open("Histograms/RecoSignalKinFitHistograms.root");
  TFile* recoDYJets = TFile::Open("Histograms/RecoDYJetsKinFitHistograms.root");
  TFile* recoTTLeptonic = TFile::Open("Histograms/RecoTTLeptonicKinFitHistograms.root");
  TFile* recoTTSemiLeptonic = TFile::Open("Histograms/RecoTTSemiLeptonicKinFitHistograms.root");

  TFile* svFitSignal = TFile::Open("Histograms/RecoSVFitSignalKinFitHistograms.root");
  TFile* svFitDYJets = TFile::Open("Histograms/RecoSVFitDYJetsKinFitHistograms.root");
  TFile* svFitTTLeptonic = TFile::Open("Histograms/RecoSVFitTTLeptonicKinFitHistograms.root");
  TFile* svFitTTSemiLeptonic = TFile::Open("Histograms/RecoSVFitTTSemiLeptonicKinFitHistograms.root");

  plotKinematics(recoSignal, "H->aa->bbtautau");
  plotInvariantMasses(recoSignal, recoDYJets, recoTTLeptonic, recoTTSemiLeptonic);
  plotKinematics(svFitSignal, "H->aa->bbtautau", true);
  plotInvariantMasses(svFitSignal, svFitDYJets, svFitTTLeptonic, svFitTTSemiLeptonic, true);
  compareInvariantMasses(recoSignal, svFitSignal, recoDYJets, svFitDYJets, recoTTLeptonic, svFitTTLeptonic, recoTTSemiLeptonic, svFitTTSemiLeptonic, false);
  compareInvariantMasses(recoSignal, svFitSignal, recoDYJets, svFitDYJets, recoTTLeptonic, svFitTTLeptonic, recoTTSemiLeptonic, svFitTTSemiLeptonic, true);
  compareInvariantMasses(recoSignal, svFitSignal, recoDYJets, svFitDYJets, recoTTLeptonic, svFitTTLeptonic, recoTTSemiLeptonic, svFitTTSemiLeptonic, false, 1);
  compareInvariantMasses(recoSignal, svFitSignal, recoDYJets, svFitDYJets, recoTTLeptonic, svFitTTLeptonic, recoTTSemiLeptonic, svFitTTSemiLeptonic, true, 1);
  compareInvariantMasses(recoSignal, svFitSignal, recoDYJets, svFitDYJets, recoTTLeptonic, svFitTTLeptonic, recoTTSemiLeptonic, svFitTTSemiLeptonic, false, 2);
  compareInvariantMasses(recoSignal, svFitSignal, recoDYJets, svFitDYJets, recoTTLeptonic, svFitTTLeptonic, recoTTSemiLeptonic, svFitTTSemiLeptonic, true, 2);
  //compareLegs(recoSignal, "H->aa->bbtautau");
  //compareLegs(recoSignal, "H->aa->bbtautau", false, true);
  //compareLegs(svFitSignal, "H->aa->bbtautau", true);
  //compareLegs(svFitSignal, "H->aa->bbtautau", true, true);
}

void overlayPlots(std::string canvasName, std::string xAxis, std::string yAxis, std::vector<TFile*> files, std::vector<std::string> histNames, std::vector<std::string> legends, std::vector<Color_t> colors)
{
  if ((files.size() != histNames.size()) || (files.size() != legends.size()) || (histNames.size() != legends.size()))
  {
    return;
  }

  TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1400, 1000);
  
  std::vector<TH1*> histograms;
  for (unsigned long int i = 0; i < histNames.size(); i++)
  {
    histograms.push_back(dynamic_cast<TH1*>(files[i]->Get(histNames[i].c_str())));
    histograms[i]->SetLineColor(colors[i]);
  }
  
  histograms[0]->GetXaxis()->SetTitle(xAxis.c_str());
  histograms[0]->GetYaxis()->SetTitle(yAxis.c_str());
  histograms[0]->SetTitle(canvasName.c_str());

  histograms[0]->Draw("HIST");
  for (unsigned long int i = 1; i < histograms.size(); i++)
  {
    histograms[i]->Draw("SAMEHIST");
    histograms[i]->Scale(histograms[0]->Integral() / histograms[i]->Integral());  // Scale each histogram to the first 
  }

  // Set the max to the maximum of the histogram with the largest maximum
  double max = 0;
  for (auto histogram : histograms)
  {
    if (histogram->GetMaximum() > max)
    {
      max = histogram->GetMaximum();
    }
  }
  histograms[0]->SetMaximum(1.05 * max);

  auto legend = new TLegend(0.7, 0.8, 0.9, 0.9);
  for (unsigned long int i = 0; i < histograms.size(); i++)
  {
    legend->AddEntry(histograms[i], legends[i].c_str(), "L");
  }
  legend->Draw();
  
  canvas->Update();
  canvas->Write();
  canvas->SaveAs(("Plots/" + canvasName + ".png").c_str());
}

void plotKinematics(TFile* file, std::string process, bool svFit)
{
  std::string suffix = "(Reco)";
  if (svFit)
  {
    suffix = "(SVFitted Reco)";
  }

  overlayPlots(process + " Final State Transverse Energy " + suffix, "Transverse Energy (GeV)", "Number of Events", {file, file}, {"Unfitted Transverse Energy;1", "Fitted Transverse Energy;1"}, {"Before Kinematic Fit", "Kinematic Fit"}, {kBlue, kRed});
  overlayPlots(process + " Final State Pseudorapidity " + suffix, "Pseudorapidity", "Number of Events", {file, file}, {"Unfitted Eta;1", "Fitted Eta;1"}, {"Before Kinematic Fit", "Kinematic Fit"}, {kBlue, kRed});
  overlayPlots(process + " Final State Phi " + suffix, "Phi", "Number of Events", {file, file}, {"Unfitted Phi;1", "Fitted Phi;1"}, {"Before Kinematic Fit", "Kinematic Fit"}, {kBlue, kRed});
}

void plotInvariantMasses(TFile* signalKinFitHistograms, TFile* dyJetsKinFitHistograms, TFile* ttLeptonicKinFitHistograms, TFile* ttSemiLeptonicKinFitHistograms, bool svFit)
{
  std::string suffix = "(Reco)";
  if (svFit)
  {
    suffix = "(SVFitted Reco)";
  }

  overlayPlots("H->aa->bbtautau Tau Tau Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", {signalKinFitHistograms, signalKinFitHistograms}, {"Unfitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1"}, {"Before Kinematic Fit", "After Kinematic Fit"}, {kBlue, kRed});
  overlayPlots("H->aa->bbtautau BB Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", {signalKinFitHistograms, signalKinFitHistograms}, {"Unfitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1"}, {"Before Kinematic Fit", "After Kinematic Fit"}, {kBlue, kRed});
  overlayPlots("DY Jets Tau Tau Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", {dyJetsKinFitHistograms, dyJetsKinFitHistograms}, {"Unfitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1"}, {"Before Kinematic Fit", "After Kinematic Fit"}, {kBlue, kRed});
  overlayPlots("DY Jets BB Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", {dyJetsKinFitHistograms, dyJetsKinFitHistograms}, {"Unfitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1"}, {"Before Kinematic Fit", "After Kinematic Fit"}, {kBlue, kRed});
  overlayPlots("TT Leptonic Tau Tau Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", {ttLeptonicKinFitHistograms, ttLeptonicKinFitHistograms}, {"Unfitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1"}, {"Before Kinematic Fit", "After Kinematic Fit"}, {kBlue, kRed});
  overlayPlots("TT Leptonic BB Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", {ttLeptonicKinFitHistograms, ttLeptonicKinFitHistograms}, {"Unfitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1"}, {"Before Kinematic Fit", "After Kinematic Fit"}, {kBlue, kRed});
  overlayPlots("TT Semi-Leptonic Tau Tau Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", {ttSemiLeptonicKinFitHistograms, ttSemiLeptonicKinFitHistograms}, {"Unfitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1"}, {"Before Kinematic Fit", "After Kinematic Fit"}, {kBlue, kRed});
  overlayPlots("TT Semi-Leptonic BB Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", {ttSemiLeptonicKinFitHistograms, ttSemiLeptonicKinFitHistograms}, {"Unfitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1"}, {"Before Kinematic Fit", "After Kinematic Fit"}, {kBlue, kRed});


  std::vector<TFile*> files = {signalKinFitHistograms, dyJetsKinFitHistograms, ttLeptonicKinFitHistograms, ttSemiLeptonicKinFitHistograms};
  std::vector<std::string> legends = {"H->aa->bbtautau", "DY Jets", "TT Leptonic", "TT Semi-Leptonic"};
  std::vector<Color_t> colors = {kBlue, kRed, kGreen, kGray};
  overlayPlots("Unfitted Tau Tau Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", files, {"Unfitted Tau Tau Invariant Mass;1", "Unfitted Tau Tau Invariant Mass;1", "Unfitted Tau Tau Invariant Mass;1", "Unfitted Tau Tau Invariant Mass;1"}, legends, colors);
  overlayPlots("Unfitted BB Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", files, {"Unfitted BB Invariant Mass;1", "Unfitted BB Invariant Mass;1", "Unfitted BB Invariant Mass;1", "Unfitted BB Invariant Mass;1"}, legends, colors);
  overlayPlots("Fitted Tau Tau Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", files, {"Fitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1"}, legends, colors);
  overlayPlots("Fitted BB Invariant Mass " + suffix, "Invariant Mass (GeV)", "Number of Events", files, {"Fitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1"}, legends, colors);

}

void compareInvariantMasses(TFile* signalKinFitHistograms, TFile* svFitSignalKinFitHistograms, TFile* dyJetsKinFitHistograms, TFile* svFitDYJetsKinFitHistograms, TFile* ttLeptonicKinFitHistograms, TFile* svFitTTLeptonicKinFitHistograms, TFile* ttSemiLeptonicKinFitHistograms, TFile* svFitTTSemiLeptonicKinFitHistograms, bool kinFit, int numBJets)
{
  std::string suffixPrefix = "";
  std::string histNamesModifier = "";
  std::string suffix = "(Unfitted Reco vs. Unfitted SVFit Reco)";

  if (numBJets == 1)
  {
    suffixPrefix = "for Events with 1 B-Jet ";
    histNamesModifier = " (1 b-jet)";
  }
  else if (numBJets == 2)
  {
    suffixPrefix = "for Events with 2 B-Jets ";
    histNamesModifier = " (2 b-jets)";
  }

  std::vector<std::string> tautauHistNames = {"Unfitted Tau Tau Invariant Mass" + histNamesModifier + ";1", "Unfitted Tau Tau Invariant Mass" + histNamesModifier + ";1"};
  std::vector<std::string> bbHistNames = {"Unfitted BB Invariant Mass" + histNamesModifier + ";1", "Unfitted BB Invariant Mass" + histNamesModifier + ";1"};

  if (kinFit)
  {
    suffix = "(Fitted Reco vs. Fitted SVFit Reco)";
    tautauHistNames = {"Fitted Tau Tau Invariant Mass" + histNamesModifier + ";1", "Fitted Tau Tau Invariant Mass" + histNamesModifier + ";1"};
    bbHistNames = {"Fitted BB Invariant Mass" + histNamesModifier + ";1", "Fitted BB Invariant Mass" + histNamesModifier + ";1"};
  }

  std::cout << tautauHistNames[0] << '\t' << tautauHistNames[1] << '\n' << bbHistNames[0] << '\t' << bbHistNames[1] << std::endl;

  overlayPlots("H->aa->bbtautau Tau Tau Invariant Mass " + suffixPrefix + suffix, "Invariant Mass (GeV)", "Number of Events", {signalKinFitHistograms, svFitSignalKinFitHistograms}, tautauHistNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("H->aa->bbtautau BB Invariant Mass " + suffixPrefix + suffix, "Invariant Mass (GeV)", "Number of Events", {signalKinFitHistograms, svFitSignalKinFitHistograms}, bbHistNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("DY Jets Tau Tau Invariant Mass " + suffixPrefix + suffix, "Invariant Mass (GeV)", "Number of Events", {dyJetsKinFitHistograms, svFitDYJetsKinFitHistograms}, tautauHistNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("DY Jets BB Invariant Mass " + suffixPrefix + suffix, "Invariant Mass (GeV)", "Number of Events", {dyJetsKinFitHistograms, svFitDYJetsKinFitHistograms}, bbHistNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("TT Leptonic Tau Tau Invariant Mass " + suffixPrefix + suffix, "Invariant Mass (GeV)", "Number of Events", {ttLeptonicKinFitHistograms, svFitTTLeptonicKinFitHistograms}, tautauHistNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("TT Leptonic BB Invariant Mass " + suffixPrefix + suffix, "Invariant Mass (GeV)", "Number of Events", {ttLeptonicKinFitHistograms, svFitTTLeptonicKinFitHistograms}, bbHistNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("TT Semi-Leptonic Tau Tau Invariant Mass " + suffixPrefix + suffix, "Invariant Mass (GeV)", "Number of Events", {ttSemiLeptonicKinFitHistograms, svFitTTSemiLeptonicKinFitHistograms}, tautauHistNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("TT Semi-Leptonic BB Invariant Mass " + suffixPrefix + suffix, "Invariant Mass (GeV)", "Number of Events", {ttSemiLeptonicKinFitHistograms, svFitTTSemiLeptonicKinFitHistograms}, bbHistNames, {"Reco", "SVFit"}, {kBlue, kRed});
}

void compareLegs(TFile* file, std::string process, bool svFit, bool kinFit)
{
  std::string prefix = "Unfitted ";
  if (kinFit)
  {
    prefix = "Fitted ";
  }

  std::string suffix = "(" + prefix + "Reco)";
  std::vector<std::string> legends = {"Leading b-quark", "Next-To-Leading b-quark", "Muonic Tau", "Hadronic Tau"};
  std::vector<TFile*> files = {file, file, file, file};
  std::vector<Color_t> colors = {kBlue, kRed, kGreen, kGray};
  std::vector<std::string> etHistNames = {prefix + "Leading b-quark Transverse Energy;1", prefix + "Next-To-Leading b-quark Transverse Energy;1", prefix + "Muonic Tau Transverse Energy;1", prefix + "Hadronic Tau Transverse Energy;1"};
  std::vector<std::string> etaHistNames = {prefix + "Leading b-quark Pseudorapidity;1", prefix + "Next-To-Leading b-quark Pseudorapidity;1", prefix + "Muonic Tau Pseudorapidity;1", prefix + "Hadronic Tau Pseudorapidity;1"};
  std::vector<std::string> phiHistNames = {prefix + "Leading b-quark Phi;1", prefix + "Next-To-Leading b-quark Phi;1", prefix + "Muonic Tau Phi;1", prefix + "Hadronic Tau Phi;1"};
  std::vector<std::string> massHistNames = {prefix + "Leading b-quark Mass;1", prefix + "Next-To-Leading b-quark Mass;1", prefix + "Muonic Tau Mass;1", prefix + "Hadronic Tau Mass;1"};

  if (svFit)
  {
    suffix = "(" + prefix + "SVFit Reco)";
    legends = {"Leading b-quark", "Next-To-Leading b-quark", "Di-Tau"};
    files.pop_back();
    colors.pop_back();
    etHistNames = {prefix + "Leading b-quark Transverse Energy;1", prefix + "Next-To-Leading b-quark Transverse Energy;1", prefix + "Di-Tau Transverse Energy;1"};
    etaHistNames = {prefix + "Leading b-quark Pseudorapidity;1", prefix + "Next-To-Leading b-quark Pseudorapidity;1", prefix + "Di-Tau Pseudorapidity;1"};
    phiHistNames = {prefix + "Leading b-quark Phi;1", prefix + "Next-To-Leading b-quark Phi;1", prefix + "Di-Tau Phi;1"};
    massHistNames = {prefix + "Leading b-quark Mass;1", prefix + "Next-To-Leading b-quark Mass;1", prefix + "Di-Tau Mass;1"};
  }

  overlayPlots(process + " Transverse Energy by Particle " + suffix, "Transverse Energy (GeV)", "Number of Events", files, etHistNames, legends, colors);
  overlayPlots(process + " Pseudorapidity by Particle " + suffix, "Pseudorapidity", "Number of Events", files, etaHistNames, legends, colors);
  overlayPlots(process + " Phi by Particle " + suffix, "Phi", "Number of Events", files, phiHistNames, legends, colors);
  overlayPlots(process + " Mass by Particle " + suffix, "Mass (GeV)", "Number of Events", files, massHistNames, legends, colors);
}