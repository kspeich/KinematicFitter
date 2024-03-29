#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"

#include <iostream>
#include <vector>
#include <cmath>

KinFitOutputModule::KinFitOutputModule(TTree* iTree, std::string iOutputFile) :
  tree(iTree),
  outputFile(iOutputFile),
  maxEvents(100000)
{}

void KinFitOutputModule::run()
{
  runFitter();
  makeHistograms();
  drawHistograms();
}

Double_t KinFitOutputModule::ErrEt(Particle particle)
{  
  if (particle.getPdgId() == 5)                                // B-jet
  {
    return (0.2 * particle.Et());
  }
  else if (particle.getPdgId() == 15)                          // Tau
  {
    if (particle.containsTag("hadronic"))                   // Hadronic Tau
    {
      return 1;
    }
    else if (particle.containsTag("muonic"))                // Muonic Tau
    {
      double ptUncertainty = 0.05 * particle.Pt();
      return (particle.E() * ptUncertainty / particle.P());
    }
    else if (particle.containsTag("ditau"))                 // Di-tau object
    {
      double ptUncertainty = 10;
      return (particle.E() * ptUncertainty / particle.P());
    }
  }
  
  return 0;          // This should error out the kinematic fitter and the cause the event to not be fitted
}

Double_t KinFitOutputModule::ErrEta(Particle particle)
{
  if (particle.getPdgId() == 5)                                // B-jet
  {
    return 0.04;
  }
  else if (particle.getPdgId() == 15)                          // Tau
  {
    if (particle.containsTag("hadronic"))                   // Hadronic Tau
    {
      return 0.02;
    }
    else if (particle.containsTag("muonic"))                // Muonic Tau
    {
      double spatialUncertainty = 0.001;      // Spacial uncertainty in eta is roughly 1000 µm => 0.001 m
      double radius = 7.5;                    // Radius of the muon detector is 7.5m

      double thetaUncertainty = abs(atan(spatialUncertainty / radius));

      double etaDerivative = -1.0 / sin(particle.Theta());                // d(eta)/d(theta) = -csc(theta)

      return thetaUncertainty * etaDerivative;                            // deltaEta = d(eta)/d(theta) * deltaTheta (roughly)
    }
    else if (particle.containsTag("ditau"))                 // Di-tau object
    {
      return 0.4;
    }
  }

  return 0;          // This should error out the kinematic fitter and the cause the event to not be fitted
}

Double_t KinFitOutputModule::ErrPhi(Particle particle)
{
  if (particle.getPdgId() == 5)                                // B-jet
  {
    return 0.04;
  }
  else if (particle.getPdgId() == 15)                          // Tau
  {
    if (particle.containsTag("hadronic"))                   // Hadronic Tau
    {
      return 0.02;
    }
    else if (particle.containsTag("muonic"))                // Muonic Tau
    {
      double spatialUncertainty = 0.0003;     // Spacial uncertainty in eta is roughly 300 µm => 0.0003 m
      double radius = 7.5;                    // Radius of the muon detector is 7.5m

      return (atan(spatialUncertainty / radius));
    }
    else if (particle.containsTag("ditau"))                 // Di-tau object
    {
      return 0.8;
    }
  }
  
  return 0;          // This should error out the kinematic fitter and the cause the event to not be fitted
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

  covMatrix(0,0) = ErrEt(particle);
  covMatrix(1,1) = ErrEta(particle);
  covMatrix(2,2) = ErrPhi(particle);

  TFitParticleEtEtaPhi *tFitParticle = new TFitParticleEtEtaPhi(&vec, &covMatrix);
  return tFitParticle;
}

Particles KinFitOutputModule::fitEvent(Particles event)
{
  auto mTauPart = event.getParticles(15, "muonic")[0];
  auto hTauPart = event.getParticles(15, "hadronic")[0];
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
  
  for (unsigned long int i = 0; i < particles.size(); i++)
  {
    auto fittedValues = particles[i]->getParCurr();
    int pdgId = event.getPdgIds()[i];
    auto tags = event[i].getTags();

    TLorentzVector v;
    double et = fittedValues->GetMatrixArray()[0];
    double eta = fittedValues->GetMatrixArray()[1];
    double phi = fittedValues->GetMatrixArray()[2];
    double m = event.getParticles()[i].M();
    double pt = calculatePt(et, eta, phi, m);
    v.SetPtEtaPhiM(pt, eta, phi, m);

    params.addParticle(v, pdgId, tags);
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
  if (max > maxEvents)
  {
    max = maxEvents;
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

    v1.SetPtEtaPhiM(pt1, eta1, phi1, m1);  // v1 contains the values of the muonic tau
    v2.SetPtEtaPhiM(pt2, eta2, phi2, m2);  // v2 contains the values of the hadronic tau
    v3.SetPtEtaPhiM(pt3, eta3, phi3, m3);  // v3 contains the values of the first b-jet

    // The PDGIDs here represent that the particles are taus and b's, but do not indicate whether they are anti-taus or anti-b's
    particles.addParticle(v1, 15, {"muonic"});
    particles.addParticle(v2, 15, {"hadronic"});
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

void KinFitOutputModule::fillHistograms(Particles event, TH1F* hEt, TH1F* hEta, TH1F* hPhi, TH1F* hTauTauInvMass, TH1F* hBBInvMass, TH1F* hAllInvMass)
{
  fillKinematicHistograms(event, hEt, hEta, hPhi);
  fillInvariantMassHistograms(event, hTauTauInvMass, hBBInvMass, hAllInvMass);
}

void KinFitOutputModule::fillKinematicHistograms(Particles event, TH1F* hEt, TH1F* hEta, TH1F* hPhi)
{
  for (auto particle : event.getParticles())
  {
    hEt->Fill(particle.Et());
    hEta->Fill(particle.Eta());
    hPhi->Fill(particle.Phi());
  }
}

void KinFitOutputModule::fillInvariantMassHistograms(Particles event, TH1F* hTauTauInvMass, TH1F* hBBInvMass, TH1F* hAllInvMass)
{
  hTauTauInvMass->Fill(event.getInvariantMass(15));    // Filling the tau tau invariant mass is unrelated to whether or not a second b-jet exists

  if (event.getNumParticles(5) == 2)
  {
    hBBInvMass->Fill(event.getInvariantMass(5));     // Only fill the BB invariant mass if there are TWO b-jets
  }

  hAllInvMass->Fill(event.getInvariantMass());
}

void KinFitOutputModule::fillInvariantMassHistogramsByNumBJets(Particles event, TH1F* hTauTauInvMass1Jet, TH1F* hBBInvMass1Jet, TH1F* hTauTauInvMass2Jet, TH1F* hBBInvMass2Jet, TH1F* hAllInvMass1Jet, TH1F* hAllInvMass2Jet)
{
  // Fills either the histograms for 1 b-jet or the histograms for 2 b-jets
  
  if (event.getNumParticles(5) == 1)
  {
    fillInvariantMassHistograms(event, hTauTauInvMass1Jet, hBBInvMass1Jet, hAllInvMass1Jet);
  }
  else if (event.getNumParticles(5) == 2)
  {
    fillInvariantMassHistograms(event, hTauTauInvMass2Jet, hBBInvMass2Jet, hAllInvMass2Jet);
  }
}

void KinFitOutputModule::fillHistogramOneParticle(Particle particle, std::string kinematic, TH1F* h)
{
  if (kinematic == "et")
  {
    h->Fill(particle.Et());
  }
  else if (kinematic == "eta")
  {
    h->Fill(particle.Eta());
  }
  else if (kinematic == "phi")
  {
    h->Fill(particle.Phi());
  }
  else if (kinematic == "m")
  {
    h->Fill(particle.M());
  }
}

void KinFitOutputModule::fillKinematicHistogramsByLeg(Particles event, std::string kinematic, TH1F* hLeadingB, TH1F* hNTLB, TH1F* hMTau, TH1F* hHTau, TH1F* hDTau)
{
  auto bquarks = event.getParticles(5);
  auto leadingB = event.getLeading(5);
  fillHistogramOneParticle(leadingB, kinematic, hLeadingB);

  if (bquarks.getNumParticles() >= 2)
  {
    auto nextToLeadingB = event.getNextToLeading(5);
    fillHistogramOneParticle(nextToLeadingB, kinematic, hNTLB);
  }

  auto mTaus = event.getParticles(15, "muonic");
  auto hTaus = event.getParticles(15, "hadronic");
  auto dTaus = event.getParticles(15, "ditau");

  if (mTaus.getNumParticles() >= 1)
  {
    auto mTau = mTaus[0];
    fillHistogramOneParticle(mTau, kinematic, hMTau);
  }
  if (hTaus.getNumParticles() >= 1)
  {
    auto hTau = hTaus[0];
    fillHistogramOneParticle(hTau, kinematic, hHTau);
  }
  if (dTaus.getNumParticles() >= 1)
  {
    auto dTau = dTaus[0];
    fillHistogramOneParticle(dTau, kinematic, hDTau);
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
  auto *hAllInvMass = new TH1F("Unfitted All-Particle Invariant Mass", "Unfitted All-Particle Invariant Mass", 200, 0, 400);

  // Initialize the fitted histograms
  auto *hEtFit = new TH1F("Fitted Transverse Energy", "Fitted Transverse Energy", 100, 0, 200);
  auto *hEtaFit = new TH1F("Fitted Eta", "Fitted Eta", 100, -10, 10);
  auto *hPhiFit = new TH1F("Fitted Phi", "Fitted Phi", 100, -4, 4);
  auto *hTauTauInvMassFit = new TH1F("Fitted Tau Tau Invariant Mass", "Fitted Tau Tau Invariant Mass", 100, 0, 200);
  auto *hBBInvMassFit = new TH1F("Fitted BB Invariant Mass", "Fitted BB Invariant Mass", 100, 0, 200);
  auto *hAllInvMassFit = new TH1F("Fitted All-Particle Invariant Mass", "Fitted All-Particle Invariant Mass", 200, 0, 400);

  // Initialize the debug/test histograms
  auto *hEtB = new TH1F("Unfitted Leading b-quark Transverse Energy", "Unfitted Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hEtNTLB = new TH1F("Unfitted Next-To-Leading b-quark Transverse Energy", "Unfitted Next-To-Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hEtMTau = new TH1F("Unfitted Muonic Tau Transverse Energy", "Unfitted Muonic Tau Transverse Energy", 100, 0, 200);
  auto *hEtHTau = new TH1F("Unfitted Hadronic Tau Transverse Energy", "Unfitted Hadronic Tau Transverse Energy", 100, 0, 200);
  auto *hEtDTau = new TH1F("Unfitted Di-Tau Transverse Energy", "Unfitted Di-Tau Transverse Energy", 100, 0, 200);

  auto *hEtBFit = new TH1F("Fitted Leading b-quark Transverse Energy", "Fitted Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hEtNTLBFit = new TH1F("Fitted Next-To-Leading b-quark Transverse Energy", "Fitted Next-To-Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hEtMTauFit = new TH1F("Fitted Muonic Tau Transverse Energy", "Fitted Muonic Tau Transverse Energy", 100, 0, 200);
  auto *hEtHTauFit = new TH1F("Fitted Hadronic Tau Transverse Energy", "Fitted Hadronic Tau Transverse Energy", 100, 0, 200);
  auto *hEtDTauFit = new TH1F("Fitted Di-Tau Transverse Energy", "Fitted Di-Tau Transverse Energy", 100, 0, 200);

  auto *hEtaB = new TH1F("Unfitted Leading b-quark Pseudorapidity", "Unfitted Leading b-quark Pseudorapidity", 100, -10, 10);
  auto *hEtaNTLB = new TH1F("Unfitted Next-To-Leading b-quark Pseudorapidity", "Unfitted Next-To-Leading b-quark Pseudorapidity", 100, -10, 10);
  auto *hEtaMTau = new TH1F("Unfitted Muonic Tau Pseudorapidity", "Unfitted Muonic Tau Pseudorapidity", 100, -10, 10);
  auto *hEtaHTau = new TH1F("Unfitted Hadronic Tau Pseudorapidity", "Unfitted Hadronic Tau Pseudorapidity", 100, -10, 10);
  auto *hEtaDTau = new TH1F("Unfitted Di-Tau Pseudorapidity", "Unfitted Di-Tau Pseudorapidity", 100, -10, 10);

  auto *hEtaBFit = new TH1F("Fitted Leading b-quark Pseudorapidity", "Fitted Leading b-quark Pseudorapidity", 100, -10, 10);
  auto *hEtaNTLBFit = new TH1F("Fitted Next-To-Leading b-quark Pseudorapidity", "Fitted Next-To-Leading b-quark Pseudorapidity", 100, -10, 10);
  auto *hEtaMTauFit = new TH1F("Fitted Muonic Tau Pseudorapidity", "Fitted Muonic Tau Pseudorapidity", 100, -10, 10);
  auto *hEtaHTauFit = new TH1F("Fitted Hadronic Tau Pseudorapidity", "Fitted Hadronic Tau Pseudorapidity", 100, -10, 10);
  auto *hEtaDTauFit = new TH1F("Fitted Di-Tau Pseudorapidity", "Fitted Di-Tau Pseudorapidity", 100, -10, 10);

  auto *hPhiB = new TH1F("Unfitted Leading b-quark Phi", "Unfitted Leading b-quark Phi", 100, -4, 4);
  auto *hPhiNTLB = new TH1F("Unfitted Next-To-Leading b-quark Phi", "Unfitted Next-To-Leading b-quark Phi", 100, -4, 4);
  auto *hPhiMTau = new TH1F("Unfitted Muonic Tau Phi", "Unfitted Muonic Tau Phi", 100, -4, 4);
  auto *hPhiHTau = new TH1F("Unfitted Hadronic Tau Phi", "Unfitted Hadronic Tau Phi", 100, -4, 4);
  auto *hPhiDTau = new TH1F("Unfitted Di-Tau Phi", "Unfitted Di-Tau Phi", 100, -4, 4);

  auto *hPhiBFit = new TH1F("Fitted Leading b-quark Phi", "Fitted Leading b-quark Phi", 100, -4, 4);
  auto *hPhiNTLBFit = new TH1F("Fitted Next-To-Leading b-quark Phi", "Fitted Next-To-Leading b-quark Phi", 100, -4, 4);
  auto *hPhiMTauFit = new TH1F("Fitted Muonic Tau Phi", "Fitted Muonic Tau Phi", 100, -4, 4);
  auto *hPhiHTauFit = new TH1F("Fitted Hadronic Tau Phi", "Fitted Hadronic Tau Phi", 100, -4, 4);
  auto *hPhiDTauFit = new TH1F("Fitted Di-Tau Phi", "Fitted Di-Tau Phi", 100, -4, 4);

  auto *hMB = new TH1F("Unfitted Leading b-quark Mass", "Unfitted Leading b-quark Mass", 100, 0, 200);
  auto *hMNTLB = new TH1F("Unfitted Next-To-Leading b-quark Mass", "Unfitted Next-To-Leading b-quark Mass", 100, 0, 200);
  auto *hMMTau = new TH1F("Unfitted Muonic Tau Mass", "Unfitted Muonic Tau Mass", 100, 0, 200);
  auto *hMHTau = new TH1F("Unfitted Hadronic Tau Mass", "Unfitted Hadronic Tau Mass", 100, 0, 200);
  auto *hMDTau = new TH1F("Unfitted Di-Tau Mass", "Unfitted Di-Tau Mass", 100, 0, 200);

  auto *hMBFit = new TH1F("Fitted Leading b-quark Mass", "Fitted Leading b-quark Mass", 100, 0, 200);
  auto *hMNTLBFit = new TH1F("Fitted Next-To-Leading b-quark Mass", "Fitted Next-To-Leading b-quark Mass", 100, 0, 200);
  auto *hMMTauFit = new TH1F("Fitted Muonic Tau Mass", "Fitted Muonic Tau Mass", 100, 0, 200);
  auto *hMHTauFit = new TH1F("Fitted Hadronic Tau Mass", "Fitted Hadronic Tau Mass", 100, 0, 200);
  auto *hMDTauFit = new TH1F("Fitted Di-Tau Mass", "Fitted Di-Tau Mass", 100, 0, 200);

  auto *hTauTauInvMass1BJet = new TH1F("Unfitted Tau Tau Invariant Mass (1 b-jet)", "Unfitted Tau Tau Invariant Mass (1 b-jet)", 100, 0, 200);
  auto *hTauTauInvMass2BJet = new TH1F("Unfitted Tau Tau Invariant Mass (2 b-jets)", "Unfitted Tau Tau Invariant Mass (2 b-jets)", 100, 0, 200);
  auto *hBBInvMass1BJet = new TH1F("Unfitted BB Invariant Mass (1 b-jet)", "Unfitted BB Invariant Mass (1 b-jet)", 100, 0, 200);
  auto *hBBInvMass2BJet = new TH1F("Unfitted BB Invariant Mass (2 b-jets)", "Unfitted BB Invariant Mass (2 b-jets)", 100, 0, 200);
  auto *hAllInvMass1BJet = new TH1F("Unfitted All-Particle Invariant Mass (1 b-jet)", "Unfitted All-Particle Invariant Mass (1 b-jet)", 200, 0, 400);
  auto *hAllInvMass2BJet = new TH1F("Unfitted All-Particle Invariant Mass (2 b-jets)", "Unfitted All-Particle Invariant Mass (2 b-jets)", 200, 0, 400);

  auto *hTauTauInvMass1BJetFit = new TH1F("Fitted Tau Tau Invariant Mass (1 b-jet)", "Fitted Tau Tau Invariant Mass (1 b-jet)", 100, 0, 200);
  auto *hTauTauInvMass2BJetFit = new TH1F("Fitted Tau Tau Invariant Mass (2 b-jets)", "Fitted Tau Tau Invariant Mass (2 b-jets)", 100, 0, 200);
  auto *hBBInvMass1BJetFit = new TH1F("Fitted BB Invariant Mass (1 b-jet)", "Fitted BB Invariant Mass (1 b-jet)", 100, 0, 200);
  auto *hBBInvMass2BJetFit = new TH1F("Fitted BB Invariant Mass (2 b-jets)", "Fitted BB Invariant Mass (2 b-jets)", 100, 0, 200);
  auto *hAllInvMass1BJetFit = new TH1F("Fitted All-Particle Invariant Mass (1 b-jet)", "Fitted All-Particle Invariant Mass (1 b-jet)", 200, 0, 400);
  auto *hAllInvMass2BJetFit = new TH1F("Fitted All-Particle Invariant Mass (2 b-jets)", "Fitted All-Particle Invariant Mass (2 b-jets)", 200, 0, 400);

  // Fill the histograms
  for (auto unfittedEvent : unfittedEvents)
  {
    fillHistograms(unfittedEvent, hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass, hAllInvMass);
    fillKinematicHistogramsByLeg(unfittedEvent, "et", hEtB, hEtNTLB, hEtMTau, hEtHTau, hEtDTau);
    fillKinematicHistogramsByLeg(unfittedEvent, "eta", hEtaB, hEtaNTLB, hEtaMTau, hEtaHTau, hEtaDTau);
    fillKinematicHistogramsByLeg(unfittedEvent, "phi", hPhiB, hPhiNTLB, hPhiMTau, hPhiHTau, hPhiDTau);
    fillKinematicHistogramsByLeg(unfittedEvent, "m", hMB, hMNTLB, hMMTau, hMHTau, hMDTau);
    fillInvariantMassHistogramsByNumBJets(unfittedEvent, hTauTauInvMass1BJet, hBBInvMass1BJet, hTauTauInvMass2BJet, hBBInvMass2BJet, hAllInvMass1BJet, hAllInvMass2BJet);
  }
  for (auto fittedEvent : fittedEvents)
  {
    fillHistograms(fittedEvent, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit, hAllInvMassFit);
    fillKinematicHistogramsByLeg(fittedEvent, "et", hEtBFit, hEtNTLBFit, hEtMTauFit, hEtHTauFit, hEtDTauFit);
    fillKinematicHistogramsByLeg(fittedEvent, "eta", hEtaBFit, hEtaNTLBFit, hEtaMTauFit, hEtaHTauFit, hEtaDTauFit);
    fillKinematicHistogramsByLeg(fittedEvent, "phi", hPhiBFit, hPhiNTLBFit, hPhiMTauFit, hPhiHTauFit, hPhiDTauFit);
    fillKinematicHistogramsByLeg(fittedEvent, "m", hMBFit, hMNTLBFit, hMMTauFit, hMHTauFit, hMDTauFit);
    fillInvariantMassHistogramsByNumBJets(fittedEvent, hTauTauInvMass1BJetFit, hBBInvMass1BJetFit, hTauTauInvMass2BJetFit, hBBInvMass2BJetFit, hAllInvMass1BJetFit, hAllInvMass2BJetFit);
  }

  // Add
  histograms = {hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass, hAllInvMass, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit, hAllInvMassFit, hEtB, hEtNTLB, hEtMTau, hEtHTau, hEtDTau, hEtBFit, hEtNTLBFit, hEtMTauFit, hEtHTauFit, hEtDTauFit, hEtaB, hEtaNTLB, hEtaMTau, hEtaHTau, hEtaDTau, hEtaBFit, hEtaNTLBFit, hEtaMTauFit, hEtaHTauFit, hEtaDTauFit, hPhiB, hPhiNTLB, hPhiMTau, hPhiHTau, hPhiDTau, hPhiBFit, hPhiNTLBFit, hPhiMTauFit, hPhiHTauFit, hPhiDTauFit, hMB, hMNTLB, hMMTau, hMHTau, hMDTau, hMBFit, hMNTLBFit, hMMTauFit, hMHTauFit, hMDTauFit, hTauTauInvMass1BJet, hBBInvMass1BJet, hTauTauInvMass2BJet, hBBInvMass2BJet, hAllInvMass1BJet, hAllInvMass2BJet, hTauTauInvMass1BJetFit, hBBInvMass1BJetFit, hTauTauInvMass2BJetFit, hBBInvMass2BJetFit, hAllInvMass1BJetFit, hAllInvMass2BJetFit};
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