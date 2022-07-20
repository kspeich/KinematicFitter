#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"

#include <iostream>
#include <vector>
#include <cmath>

KinFitOutputModule::KinFitOutputModule(TTree* tree, std::string iOutputFile) :
  outputFile(iOutputFile)
{
  runFitter(tree);
}

void KinFitOutputModule::run()
{
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

Particles KinFitOutputModule::fitEvent(Particles event)
{
  auto mTauVec = event.getParticleVectors(15)[0];
  auto hTauVec = event.getParticleVectors(15)[1];
  auto bJet1Vec = event.getParticleVectors(5)[0];

  TMatrixD m1(3,3);
  TMatrixD m2(3,3);
  TMatrixD m3(3,3);
  m1.Zero();
  m2.Zero();
  m3.Zero();

  m1(0,0) = ErrEt (mTauVec); // et
  m1(1,1) = ErrEta(mTauVec); // eta
  m1(2,2) = ErrPhi(mTauVec); // phi
  m2(0,0) = ErrEt (hTauVec); // et
  m2(1,1) = ErrEta(hTauVec); // eta
  m2(2,2) = ErrPhi(hTauVec); // phi
  m3(0,0) = ErrEt (bJet1Vec); // et
  m3(1,1) = ErrEta(bJet1Vec); // eta
  m3(2,2) = ErrPhi(bJet1Vec); // phi

  TFitParticleEtEtaPhi *mTau = new TFitParticleEtEtaPhi( "mTau", "mTau", &mTauVec, &m1 );
  TFitParticleEtEtaPhi *hTau = new TFitParticleEtEtaPhi( "hTau", "hTau", &hTauVec, &m2 );
  TFitParticleEtEtaPhi *bJet1 = new TFitParticleEtEtaPhi( "bJet1", "bJet1", &bJet1Vec, &m3 );

  // mTau and hTau must make an a pseudoscalar
  TFitConstraintM *mCons1 = new TFitConstraintM( "AMassConstraint1", "AMass-Constraint1", 0, 0 , 45.);
  mCons1->addParticles1( mTau, hTau );

  std::vector<TFitParticleEtEtaPhi*> particles = {mTau, hTau, bJet1};
  std::vector<TFitConstraintM*> constraints = {mCons1};

  if (event.getNumParticles(5) == 2)  // Do everything including the second b-jet
  {
    auto bJet2Vec = event.getParticleVectors(5)[1];
    
    TMatrixD m4(3,3);
    m4.Zero();
    
    m4(0,0) = ErrEt (bJet2Vec); // et
    m4(1,1) = ErrEta(bJet2Vec); // eta
    m4(2,2) = ErrPhi(bJet2Vec); // phi

    TFitParticleEtEtaPhi *bJet2 = new TFitParticleEtEtaPhi( "bJet2", "bJet2", &bJet2Vec, &m4 );

    // bJet1 and bJet2 must make an a pseudoscalar: this only happens when there is a second b-jet
    TFitConstraintM *mCons2 = new TFitConstraintM( "AMassConstraint2", "AMass-Constraint2", 0, 0 , 45.);
    mCons2->addParticles1( bJet1, bJet2 );

    // All four particles must make a Higgs
    TFitConstraintM *mCons3 = new TFitConstraintM( "HiggsMassConstraint", "HiggsMass-Constraint", 0, 0, 125.);
    mCons3->addParticles1( mTau, hTau, bJet1, bJet2 );

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
  
  Particles params = {};
  
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

    auto fittedParticle = Particle(v, pdgId);
    params.addParticle(fittedParticle);
  }
  
  for (auto constraint : constraints)
  {
    delete constraint;
  }

  delete fitter;

  return params;
}

void KinFitOutputModule::runFitter(TTree* tree)
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

  // Fill the histograms
  for (auto unfittedEvent : unfittedEvents)
  {
    fillHistograms(unfittedEvent, hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass);
  }
  for (auto fittedEvent : fittedEvents)
  {
    fillHistograms(fittedEvent, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit);
  }

  // Add
  histograms = {hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit};
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