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

void plotOutput()
{
  gStyle->SetOptStat(0);
  TFile* kinFitHistograms = TFile::Open("SignalKinFitHistograms.root");
  TFile* 
  
  //overlayPlots("H->aa->bbtautau Final State Transverse Energy (Gen)", "Transverse Energy (GeV)", "Number of Events", {kinFitHistograms, kinFitHistograms}, {"Unfitted Transverse Energy;1", "Fitted Transverse Energy;1"}, {"Original", "Kinematic Fit"}, {kBlue, kRed});
  //overlayPlots("H->aa->bbtautau Final State Pseudorapidity (Gen)", "Pseudorapidity", "Number of Events", {kinFitHistograms, kinFitHistograms}, {"Unfitted Eta;1", "Fitted Eta;1"}, {"Original", "Kinematic Fit"}, {kBlue, kRed});
  //overlayPlots("H->aa->bbtautau Final State Phi (Gen)", "Phi", "Number of Events", {kinFitHistograms, kinFitHistograms}, {"Unfitted Phi;1", "Fitted Phi;1"}, {"Original", "Kinematic Fit"}, {kBlue, kRed});
  //overlayPlots("H->aa->bbtautau TauTau Invariant Mass (Gen)", "Invariant Mass (GeV)", "Number of Events", {kinFitHistograms, kinFitHistograms}, {"Unfitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1"}, {"Original", "Kinematic Fit"}, {kBlue, kRed});
  //overlayPlots("H->aa->bbtautau BB Invariant Mass (Gen)", "Invariant Mass (GeV)", "Number of Events", {kinFitHistograms, kinFitHistograms}, {"Unfitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1"}, {"Original", "Kinematic Fit"}, {kBlue, kRed});

  overlayPlots("Unfitted BB Invariant Mass (Gen)", "Invairant Mass (GeV)", "Number of Events", {kinFitHistograms, })
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
  canvas->SaveAs((canvasName + ".png").c_str());
}