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
void compareInvariantMasses(TFile* signalKinFitHistograms, TFile* svFitSignalKinFitHistograms, TFile* dyJetsKinFitHistograms, TFile* svFitDYJetsKinFitHistograms, TFile* ttLeptonicKinFitHistograms, TFile* svFitTTLeptonicKinFitHistograms, TFile* ttSemiLeptonicKinFitHistograms, TFile* svFitTTSemiLeptonicKinFitHistograms, bool kinFit);

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
  plotInvariantMasses(svFitSignal, svFitDYJets, svFitTTLeptonic, svFitTTSemiLeptonic);
  compareInvariantMasses(recoSignal, svFitSignal, recoDYJets, svFitDYJets, recoTTLeptonic, svFitTTLeptonic, recoTTSemiLeptonic, svFitTTSemiLeptonic, false);
  compareInvariantMasses(recoSignal, svFitSignal, recoDYJets, svFitDYJets, recoTTLeptonic, svFitTTLeptonic, recoTTSemiLeptonic, svFitTTSemiLeptonic, true);
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
}

void compareInvariantMasses(TFile* signalKinFitHistograms, TFile* svFitSignalKinFitHistograms, TFile* dyJetsKinFitHistograms, TFile* svFitDYJetsKinFitHistograms, TFile* ttLeptonicKinFitHistograms, TFile* svFitTTLeptonicKinFitHistograms, TFile* ttSemiLeptonicKinFitHistograms, TFile* svFitTTSemiLeptonicKinFitHistograms, bool kinFit)
{
  std::vector<std::string> histNames = {"Unfitted Tau Tau Invariant Mass;1", "Unfitted Tau Tau Invariant Mass;1"};
  if (kinFit)
  {
    histNames = {"Fitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1"};
  }

  overlayPlots("H->aa->bbtautau Tau Tau Invariant Mass (Reco vs. SVFit Reco)", "Invariant Mass (GeV)", "Number of Events", {signalKinFitHistograms, svFitSignalKinFitHistograms}, histNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("H->aa->bbtautau BB Invariant Mass (Reco vs. SVFit Reco)", "Invariant Mass (GeV)", "Number of Events", {signalKinFitHistograms, svFitSignalKinFitHistograms}, histNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("DY Jets Tau Tau Invariant Mass (Reco vs. SVFit Reco)", "Invariant Mass (GeV)", "Number of Events", {dyJetsKinFitHistograms, svFitDYJetsKinFitHistograms}, histNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("DY Jets BB Invariant Mass (Reco vs. SVFit Reco)", "Invariant Mass (GeV)", "Number of Events", {dyJetsKinFitHistograms, svFitDYJetsKinFitHistograms}, histNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("TT Leptonic Tau Tau Invariant Mass (Reco vs. SVFit Reco)", "Invariant Mass (GeV)", "Number of Events", {ttLeptonicKinFitHistograms, svFitTTLeptonicKinFitHistograms}, histNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("TT Leptonic BB Invariant Mass (Reco vs. SVFit Reco)", "Invariant Mass (GeV)", "Number of Events", {ttLeptonicKinFitHistograms, svFitTTLeptonicKinFitHistograms}, histNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("TT Semi-Leptonic Tau Tau Invariant Mass (Reco vs. SVFit Reco)", "Invariant Mass (GeV)", "Number of Events", {ttSemiLeptonicKinFitHistograms, svFitTTSemiLeptonicKinFitHistograms}, histNames, {"Reco", "SVFit"}, {kBlue, kRed});
  overlayPlots("TT Semi-Leptonic BB Invariant Mass (Reco vs. SVFit Reco)", "Invariant Mass (GeV)", "Number of Events", {ttSemiLeptonicKinFitHistograms, svFitTTSemiLeptonicKinFitHistograms}, histNames, {"Reco", "SVFit"}, {kBlue, kRed});
}