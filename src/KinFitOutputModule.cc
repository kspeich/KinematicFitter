#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"

#include <iostream>
#include <vector>
#include <cmath>

KinFitOutputModule::KinFitOutputModule(TTree* iTree, std::string iOutputFile) :
  tree(iTree),
  outputFile(iOutputFile)
{}

void KinFitOutputModule::run()
{
  runFitter();
  makeHistograms();
  drawHistograms();
}

Double_t KinFitOutputModule::ErrEt(TLorentzVector particleVec)
{  
  return 0.05 * particleVec.Et();
}

Double_t KinFitOutputModule::ErrEta(TLorentzVector particleVec)
{
  return 0.05 * particleVec.Eta();
}

Double_t KinFitOutputModule::ErrPhi(TLorentzVector particleVec)
{
  return 0.05 * particleVec.Phi();
}

Float_t KinFitOutputModule::calculatePt(Float_t Et, Float_t eta, Float_t phi, Float_t m)
{
  Float_t E = Et * cosh(eta);
  Float_t p = pow(E * E - m * m, 0.5);

  return (p * Et / E);
}

void KinFitOutputModule::print(TKinFitter *fitter)
{
  std::cout << "=============================================" << std ::endl;
  std::cout << "-> Number of measured Particles  : " << fitter->nbMeasParticles() << std::endl;
  std::cout << "-> Number of unmeasured particles: " << fitter->nbUnmeasParticles() << std::endl;
  std::cout << "-> Number of constraints         : " << fitter->nbConstraints() << std::endl;
  std::cout << "-> Number of degrees of freedom  : " << fitter->getNDF() << std::endl;
  std::cout << "-> Number of parameters A        : " << fitter->getNParA() << std::endl;
  std::cout << "-> Number of parameters B        : " << fitter->getNParB() << std::endl;
  std::cout << "-> Maximum number of iterations  : " << fitter->getMaxNumberIter() << std::endl;
  std::cout << "-> Maximum deltaS                : " << fitter->getMaxDeltaS() << std::endl;
  std::cout << "-> Maximum F                     : " << fitter->getMaxF() << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std ::endl;
  std::cout << "-> Status                        : " << fitter->getStatus() << std::endl;
  std::cout << "-> Number of iterations          : " << fitter->getNbIter() << std::endl;
  std::cout << "-> S                             : " << fitter->getS() << std::endl;
  std::cout << "-> F                             : " << fitter->getF() << std::endl;
  std::cout << "=============================================" << std ::endl;
}

TFitParticleEtEtaPhi* KinFitOutputModule::convertParticle(Particle particle)
{
  auto vec = particle.getFourVector();
  
  TMatrixD covMatrix(3,3);
  covMatrix.Zero();

  covMatrix(0,0) = ErrEt(vec);
  covMatrix(1,1) = ErrEta(vec);
  covMatrix(2,2) = ErrPhi(vec);

  TFitParticleEtEtaPhi *tFitParticle = new TFitParticleEtEtaPhi(&vec, &covMatrix);
  return tFitParticle;
}

Particles KinFitOutputModule::fitEvent(Particles event)
{
  auto mTauPart = event.getParticles(15)[0];
  auto hTauPart = event.getParticles(15)[1];
  auto bJet1Part = event.getParticles(5)[0];

  auto mTau = convertParticle(mTauPart);
  auto hTau = convertParticle(hTauPart);
  auto bJet1 = convertParticle(bJet1Part);

  // mTau and hTau must make an a pseudoscalar
  TFitConstraintM *mCons1 = new TFitConstraintM( "AMassConstraint1", "AMass-Constraint1", 0, 0 , 45.);
  mCons1->addParticles1(mTau, hTau);

  std::vector<TFitParticleEtEtaPhi*> particles = {mTau, hTau, bJet1};
  std::vector<TFitConstraintM*> constraints = {mCons1};

  if (event.getNumParticles(5) == 2)  // Do everything including the second b-jet
  {
    auto bJet2Part = event.getParticles(5)[1];
    auto bJet2 = convertParticle(bJet2Part);

    // bJet1 and bJet2 must make an a pseudoscalar: this only happens when there is a second b-jet
    TFitConstraintM *mCons2 = new TFitConstraintM( "AMassConstraint2", "AMass-Constraint2", 0, 0 , 45.);
    mCons2->addParticles1(bJet1, bJet2);

    // All four particles must make a Higgs
    TFitConstraintM *mCons3 = new TFitConstraintM( "HiggsMassConstraint", "HiggsMass-Constraint", 0, 0, 125.);
    mCons3->addParticles1(mTau, hTau, bJet1, bJet2);

    particles.push_back(bJet2);
    constraints.push_back(mCons2);
    constraints.push_back(mCons3);
  }

  // Make the fitter, add particles and constraints
  TKinFitter* fitter = new TKinFitter("fitter", "fitter");
  for (auto particle : particles)
  {
    fitter->addMeasParticle(particle);
  }
  for(auto constraint : constraints)
  {
    fitter->addConstraint(constraint);
  }

  // Set convergence criteria
  fitter->setMaxNbIter( 30 );
  fitter->setMaxDeltaS( 1e-2 );
  fitter->setMaxF( 1e-1 );
  fitter->setVerbosity(3);

  // Perform the fit
  // std::cout << "Performing kinematic fit..." << std::endl;
  // print(fitter);
  fitter->fit();
  // std::cout << "Done." << std::endl;
  // print(fitter);
  
  Particles params;
  
  for (unsigned long int i; i < particles.size(); i++)
  {
    auto fittedValues = particles[i]->getParCurr();
    int pdgId = event.getPdgIds()[i];

    TLorentzVector v;
    double et = fittedValues->GetMatrixArray()[0];
    double eta = fittedValues->GetMatrixArray()[1];
    double phi = fittedValues->GetMatrixArray()[2];
    double m = event.getParticles()[i].M();
    double pt = calculatePt(et, eta, phi, m);
    v.SetPtEtaPhiM(pt, eta, phi, m);

    params.addParticle(v, pdgId);
  }
  
  for (auto constraint : constraints)
  {
    delete constraint;
  }

  delete fitter;

  return params;
}

void KinFitOutputModule::runFitter()
{
  // Set the addresses of the branches to elsewhere
  Float_t pt1, eta1, phi1, m1, pt2, eta2, phi2, m2, pt3, eta3, phi3, m3, pt4, eta4, phi4, m4;
  tree->SetBranchAddress("pt_1", &pt1);
  tree->SetBranchAddress("eta_1", &eta1);
  tree->SetBranchAddress("phi_1", &phi1);
  tree->SetBranchAddress("m_1", &m1);
  tree->SetBranchAddress("pt_2", &pt2);
  tree->SetBranchAddress("eta_2", &eta2);
  tree->SetBranchAddress("phi_2", &phi2);
  tree->SetBranchAddress("m_2", &m2);
  tree->SetBranchAddress("bpt_deepflavour_1", &pt3);
  tree->SetBranchAddress("beta_deepflavour_1", &eta3);
  tree->SetBranchAddress("bphi_deepflavour_1", &phi3);
  tree->SetBranchAddress("bm_deepflavour_1", &m3);
  tree->SetBranchAddress("bpt_deepflavour_2", &pt4);
  tree->SetBranchAddress("beta_deepflavour_2", &eta4);
  tree->SetBranchAddress("bphi_deepflavour_2", &phi4);
  tree->SetBranchAddress("bm_deepflavour_2", &m4);

  // Loop through each event, perform necessary calculations, and fill the histograms
  int max = tree->GetEntries();
  if (max > 100000)
  {
    max = 100000;
  }
  for(int i = 0; i < max; i++)   // GetEntries() returns the # of entries in the branch
  {
    if ((i % 1000) == 0)
    {
      std::cout << "Event #: " << i << " / " << tree->GetEntries() << std::endl;
    }
    
    tree->GetEntry(i);

    TLorentzVector v1, v2, v3, v4;
    Particles particles;

    v1.SetPtEtaPhiM(pt1, eta1, phi1, m1);  // v1 contains the values of the (muonic) tau
    v2.SetPtEtaPhiM(pt2, eta2, phi2, m2);  // v2 contains the values of the (hadronic) tau
    v3.SetPtEtaPhiM(pt3, eta3, phi3, m3);  // v3 contains the values of the first b-jet

    // The PDGIDs here represent that the particles are taus and b's, but do not indicate whether they are anti-taus or anti-b's
    particles.addParticle(v1, 15);
    particles.addParticle(v2, 15);
    particles.addParticle(v3, 5);

    if (pt4 != -9999 && eta4 != -9999 && phi4 != -9999 && m4 != -9999)
    {
      v4.SetPtEtaPhiM(pt4, eta4, phi4, m4);  // v4 contains the values of the second b-jet ONLY IF IT EXISTS
      particles.addParticle(v4, 5);
    }

    unfittedEvents.push_back(particles);
    fittedEvents.push_back(fitEvent(particles));
  }
}

void KinFitOutputModule::fillHistograms(Particles event, TH1F* hEt, TH1F* hEta, TH1F* hPhi, TH1F* hTauTauInvMass, TH1F* hBBInvMass)
{
  for (auto particle : event.getParticles())
  {
    hEt->Fill(particle.Et());
    hEta->Fill(particle.Eta());
    hPhi->Fill(particle.Phi());
  }

  hTauTauInvMass->Fill(event.getInvariantMass(15));    // Filling the tau tau invariant mass is unrelated to whether or not a second b-jet exists

  if (event.getNumParticles(5) == 2)
  {
    hBBInvMass->Fill(event.getInvariantMass(5));     // Only fill the BB invariant mass if there are TWO b-jets
  }
}

void KinFitOutputModule::fillKinematicHistogramsByLeg(Particles event, int pdgId, std::string kinematic, TH1F* hLeading, TH1F* hNTL)
{
  auto particles = event.getParticles(pdgId);
  auto leading = event.getLeading(pdgId);

  if (kinematic == "et")
  {
    hLeading->Fill(leading.Et());
  }
  else if (kinematic == "eta")
  {
    hLeading->Fill(leading.Eta());
  }
  else if (kinematic == "phi")
  {
    hLeading->Fill(leading.Phi());
  }

  if (hNTL != nullptr && event.getNumParticles(pdgId) >= 2)
  {
    auto nextToLeading = event.getNextToLeading(pdgId);
    if (kinematic == "et")
    {
      hNTL->Fill(nextToLeading.Et());
    }
    else if (kinematic == "eta")
    {
      hNTL->Fill(nextToLeading.Eta());
    }
    else if (kinematic == "phi")
    {
      hNTL->Fill(nextToLeading.Phi());
    }
  }
}

void KinFitOutputModule::makeHistograms()
{
  // Initialize the unfitted histograms
  auto *hEt = new TH1F("Unfitted Transverse Energy", "Unfitted Transverse Energy", 100, 0, 200);
  auto *hEta = new TH1F("Unfitted Eta", "Unfitted Eta", 100, -10, 10);
  auto *hPhi = new TH1F("Unfitted Phi", "Unfitted Phi", 100, -4, 4);
  auto *hTauTauInvMass = new TH1F("Unfitted Tau Tau Invariant Mass", "Unfitted Tau Tau Invariant Mass", 100, 0, 200);
  auto *hBBInvMass = new TH1F("Unfitted BB Invariant Mass", "Unfitted BB Invariant Mass", 100, 0, 200);

  // Initialize the fitted histograms
  auto *hEtFit = new TH1F("Fitted Transverse Energy", "Fitted Transverse Energy", 100, 0, 200);
  auto *hEtaFit = new TH1F("Fitted Eta", "Fitted Eta", 100, -10, 10);
  auto *hPhiFit = new TH1F("Fitted Phi", "Fitted Phi", 100, -4, 4);
  auto *hTauTauInvMassFit = new TH1F("Fitted Tau Tau Invariant Mass", "Fitted Tau Tau Invariant Mass", 100, 0, 200);
  auto *hBBInvMassFit = new TH1F("Fitted BB Invariant Mass", "Fitted BB Invariant Mass", 100, 0, 200);

  // Initialize the debug/test histograms
  auto *hEtB = new TH1F("Unfitted Leading b-quark Transverse Energy", "Unfitted Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hEtNTLB = new TH1F("Unfitted Next-To-Leading b-quark Transverse Energy", "Unfitted Next-To-Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hEtTau = new TH1F("Unfitted Leading Tau Transverse Energy", "Unfitted Leading Tau Transverse Energy", 100, 0, 200);
  auto *hEtNTLTau = new TH1F("Unfitted Next-To-Leading Tau Transverse Energy", "Unfitted Next-To-Leading Tau Transverse Energy", 100, 0, 200);
  auto *hEtBFit = new TH1F("Fitted Leading b-quark Transverse Energy", "Fitted Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hEtNTLBFit = new TH1F("Fitted Next-To-Leading b-quark Transverse Energy", "Fitted Next-To-Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hEtTauFit = new TH1F("Fitted Leading Tau Transverse Energy", "Fitted Leading Tau Transverse Energy", 100, 0, 200);
  auto *hEtNTLTauFit = new TH1F("Fitted Next-To-Leading Tau Transverse Energy", "Fitted Next-To-Leading Tau Transverse Energy", 100, 0, 200);
  auto *hEtaBFit = new TH1F("Fitted Leading b-quark Pseudorapidity", "Fitted Leading b-quark Pseudorapidity", 100, -10, 10);
  auto *hEtaNTLBFit = new TH1F("Fitted Next-To-Leading b-quark Pseudorapidity", "Fitted Next-To-Leading b-quark Pseudorapidity", 100, -10, 10);
  auto *hEtaTauFit = new TH1F("Fitted Leading Tau Pseudorapidity", "Fitted Leading Tau Pseudorapidity", 100, -10, 10);
  auto *hEtaNTLTauFit = new TH1F("Fitted Next-To-Leading Tau Pseudorapidity", "Fitted Next-To-Leading Tau Pseudorapidity", 100, -10, 10);
  auto *hPhiBFit = new TH1F("Fitted Leading b-quark Phi", "Fitted Leading b-quark Phi", 100, -4, 4);
  auto *hPhiNTLBFit = new TH1F("Fitted Next-To-Leading b-quark Phi", "Fitted Next-To-Leading b-quark Phi", 100, -4, 4);
  auto *hPhiTauFit = new TH1F("Fitted Leading Tau Phi", "Fitted Leading Tau Phi", 100, -4, 4);
  auto *hPhiNTLTauFit = new TH1F("Fitted Next-To-Leading Tau Phi", "Fitted Next-To-Leading Tau Phi", 100, -4, 4);

  // Fill the histograms
  for (auto unfittedEvent : unfittedEvents)
  {
    fillHistograms(unfittedEvent, hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass);
    fillKinematicHistogramsByLeg(unfittedEvent, 5, "et", hEtB, hEtNTLB);
    fillKinematicHistogramsByLeg(unfittedEvent, 15, "et", hEtTau, hEtNTLTau);
  }
  for (auto fittedEvent : fittedEvents)
  {
    fillHistograms(fittedEvent, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit);
    fillKinematicHistogramsByLeg(fittedEvent, 5, "et", hEtBFit, hEtNTLBFit);
    fillKinematicHistogramsByLeg(fittedEvent, 5, "eta", hEtaBFit, hEtaNTLBFit);
    fillKinematicHistogramsByLeg(fittedEvent, 5, "phi", hPhiBFit, hPhiNTLBFit);
    fillKinematicHistogramsByLeg(fittedEvent, 15, "et", hEtTauFit, hEtNTLTauFit);
    fillKinematicHistogramsByLeg(fittedEvent, 15, "eta", hEtaTauFit, hEtaNTLTauFit);
    fillKinematicHistogramsByLeg(fittedEvent, 15, "phi", hPhiTauFit, hPhiNTLTauFit);
  }

  // Add
  histograms = {hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit, hEtB, hEtNTLB, hEtTau, hEtNTLTau, hEtBFit, hEtNTLBFit, hEtTauFit, hEtNTLTauFit, hEtaBFit, hEtaNTLBFit, hEtaTauFit, hEtaNTLTauFit, hPhiBFit, hPhiNTLBFit, hPhiTauFit, hPhiNTLTauFit};
}

void KinFitOutputModule::drawHistograms()
{
  TFile* outFile = new TFile(("Histograms/" + outputFile).c_str(), "RECREATE");
  
  for (auto histogram : histograms)
  {
    histogram->Write();
  }

  std::cout << "Histograms written to Histograms/" << outputFile << '\n';
  outFile->Close();
  delete outFile;
}