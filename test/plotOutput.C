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

void displayPlot(std::string canvasName, std::string xAxis, std::string yAxis, TFile* file1, std::string histName1);
void overlayTwoPlots(std::string canvasName, std::string xAxis, std::string yAxis, TFile* file1, TFile* file2, std::string histName1, std::string histName2, std::string legend1, std::string legend2);

void plotOutput()
{
  gStyle->SetOptStat(0);
  TFile* kinFitHistograms = TFile::Open("KinFitHistograms.root");

  overlayTwoPlots("H->aa->bbtautau Final State Transverse Energy (Gen)", "Transverse Energy (GeV)", "Number of Events", kinFitHistograms, kinFitHistograms, "Unfitted Transverse Energy;1", "Fitted Transverse Energy;1", "Original", "Kinematic Fit");
  overlayTwoPlots("H->aa->bbtautau Final State Pseudorapidity (Gen)", "Pseudorapidity", "Number of Events", kinFitHistograms, kinFitHistograms, "Unfitted Eta;1", "Fitted Eta;1", "Original", "Kinematic Fit");
  overlayTwoPlots("H->aa->bbtautau Final State Phi (Gen)", "Phi", "Number of Events", kinFitHistograms, kinFitHistograms, "Unfitted Phi;1", "Fitted Phi;1", "Original", "Kinematic Fit");
  overlayTwoPlots("H->aa->bbtautau TauTau Invariant Mass (Gen)", "Invariant Mass", "Number of Events", kinFitHistograms, kinFitHistograms, "Unfitted Tau Tau Invariant Mass;1", "Fitted Tau Tau Invariant Mass;1", "Original", "Kinematic Fit");
  overlayTwoPlots("H->aa->bbtautau BB Invariant Mass (Gen)", "Invariant Mass", "Number of Events", kinFitHistograms, kinFitHistograms, "Unfitted BB Invariant Mass;1", "Fitted BB Invariant Mass;1", "Original", "Kinematic Fit");
}

void displayPlot(std::string canvasName, std::string xAxis, std::string yAxis, TFile* file1, std::string histName1)
{
  TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1400, 1000);
  
  auto *hist1 = dynamic_cast<TH1*>(file1->Get(histName1.c_str()));

  hist1->SetLineColor(kBlue);
  
  hist1->GetXaxis()->SetTitle(xAxis.c_str());
  hist1->GetYaxis()->SetTitle(yAxis.c_str());
  hist1->SetTitle(canvasName.c_str());

  hist1->Draw("HIST");

  // Scale hist2 to hist1 and set the max to whichever has larger values
  //hist1->SetMaximum(1.05 * std::max(hist1->GetMaximum(), hist2->GetMaximum()));

  /*auto legend = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend->AddEntry(hist1, legend1.c_str(), "L");
  legend->Draw();*/
  
  canvas->Update();
  canvas->Write();
  canvas->SaveAs((canvasName + ".png").c_str());
}

void overlayTwoPlots(std::string canvasName, std::string xAxis, std::string yAxis, TFile* file1, TFile* file2, std::string histName1, std::string histName2, std::string legend1, std::string legend2)
{
  TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1400, 1000);
  
  auto *hist1 = dynamic_cast<TH1*>(file1->Get(histName1.c_str()));
  auto *hist2 = dynamic_cast<TH1*>(file2->Get(histName2.c_str()));

  hist1->SetLineColor(kBlue);
  hist2->SetLineColor(kRed);
  
  hist1->GetXaxis()->SetTitle(xAxis.c_str());
  hist1->GetYaxis()->SetTitle(yAxis.c_str());
  hist1->SetTitle(canvasName.c_str());

  hist1->Draw("HIST");
  hist2->Draw("SAMEHIST");

  // Scale hist2 to hist1 and set the max to whichever has larger values
  hist2->Scale(hist1->Integral() / hist2->Integral());
  hist1->SetMaximum(1.05 * std::max(hist1->GetMaximum(), hist2->GetMaximum()));

  auto legend = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend->AddEntry(hist1, legend1.c_str(), "L");
  legend->AddEntry(hist2, legend2.c_str(), "L");
  legend->Draw();
  
  canvas->Update();
  canvas->Write();
  canvas->SaveAs((canvasName + ".png").c_str());
}