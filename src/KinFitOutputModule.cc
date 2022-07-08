#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"

#include <iostream>
#include <vector>
#include <cmath>

void KinFitOutputModule::run()
{
  makeHistograms();
  drawHistograms();
}

Double_t KinFitOutputModule::ErrEt(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 5.6;
    b = 1.25;
    c = 0.033;
  }
  else{
    a = 4.8;
    b = 0.89;
    c = 0.043;
  }
  InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
  return InvPerr2;
}

Double_t KinFitOutputModule::ErrEta(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 1.215;
    b = 0.037;
    c = 7.941 * 0.0001;
  }
  else{
    a = 1.773;
    b = 0.034;
    c = 3.56 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

Double_t KinFitOutputModule::ErrPhi(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 6.65;
    b = 0.04;
    c = 8.49 * 0.00001;
  }
  else{
    a = 2.908;
    b = 0.021;
    c = 2.59 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
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

std::vector<TLorentzVector> KinFitOutputModule::fitEvent(std::vector<TLorentzVector> particleVectors)
{
  auto mTauVec = particleVectors[0];
  auto hTauVec = particleVectors[1];
  auto bJet1Vec = particleVectors[2];

  TMatrixD m1(3,3);
  TMatrixD m2(3,3);
  TMatrixD m3(3,3);
  m1.Zero();
  m2.Zero();
  m3.Zero();

  // In this example the covariant matrix depends on the transverse energy and eta of the jets
  m1(0,0) = ErrEt (mTauVec.Et(), mTauVec.Eta()); // et
  m1(1,1) = ErrEta(mTauVec.Et(), mTauVec.Eta()); // eta
  m1(2,2) = ErrPhi(mTauVec.Et(), mTauVec.Eta()); // phi
  m2(0,0) = ErrEt (hTauVec.Et(), hTauVec.Eta()); // et
  m2(1,1) = ErrEta(hTauVec.Et(), hTauVec.Eta()); // eta
  m2(2,2) = ErrPhi(hTauVec.Et(), hTauVec.Eta()); // phi
  m3(0,0) = ErrEt (bJet1Vec.Et(), bJet1Vec.Eta()); // et
  m3(1,1) = ErrEta(bJet1Vec.Et(), bJet1Vec.Eta()); // eta
  m3(2,2) = ErrPhi(bJet1Vec.Et(), bJet1Vec.Eta()); // phi

  TFitParticleEtEtaPhi *mTau = new TFitParticleEtEtaPhi( "mTau", "mTau", &mTauVec, &m1 );
  TFitParticleEtEtaPhi *hTau = new TFitParticleEtEtaPhi( "hTau", "hTau", &hTauVec, &m2 );
  TFitParticleEtEtaPhi *bJet1 = new TFitParticleEtEtaPhi( "bJet1", "bJet1", &bJet1Vec, &m3 );

  // mTau and hTau must make an a pseudoscalar
  TFitConstraintM *mCons1 = new TFitConstraintM( "AMassConstraint1", "AMass-Constraint1", 0, 0 , 45.);
  mCons1->addParticles1( mTau, hTau );

  std::vector<TFitParticleEtEtaPhi*> particles = {mTau, hTau, bJet1};
  std::vector<TFitConstraintM*> constraints = {mCons1};

  if (particleVectors.size() == 4)  // Do everything including the second b-jet
  {
    auto bJet2Vec = particleVectors[3];
    
    TMatrixD m4(3,3);
    m4.Zero();
    
    m4(0,0) = ErrEt (bJet2Vec.Et(), bJet2Vec.Eta()); // et
    m4(1,1) = ErrEta(bJet2Vec.Et(), bJet2Vec.Eta()); // eta
    m4(2,2) = ErrPhi(bJet2Vec.Et(), bJet2Vec.Eta()); // phi

    TFitParticleEtEtaPhi *bJet2 = new TFitParticleEtEtaPhi( "bJet2", "bJet2", &bJet2Vec, &m4 );

    // bJet1 and bJet2 must make an a pseudoscalar: this only happens when there is a second b-jet
    TFitConstraintM *mCons2 = new TFitConstraintM( "AMassConstraint2", "AMass-Constraint2", 0, 0 , 45.);
    mCons2->addParticles1( bJet1, bJet2 );

    particles.push_back(bJet2);
    constraints.push_back(mCons2);
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
  
  std::vector<TLorentzVector> params = {};
  
  for (unsigned long int i; i < particles.size(); i++)
  {
    auto fittedValues = particles[i]->getParCurr();

    TLorentzVector v;

    Float_t et = fittedValues->GetMatrixArray()[0];
    Float_t eta = fittedValues->GetMatrixArray()[1];
    Float_t phi = fittedValues->GetMatrixArray()[2];
    Float_t m = particleVectors[i].M();
    Float_t pt = calculatePt(et, eta, phi, m);

    v.SetPtEtaPhiM(pt, eta, phi, m);
    params.push_back(v);
  }
  
  for (auto constraint : constraints)
  {
    delete constraint;
  }

  delete fitter;

  return params;
}

void KinFitOutputModule::fillHistograms(std::vector<TLorentzVector> particleVectors, TH1F* hEt, TH1F* hEta, TH1F* hPhi, TH1F* hTauTauInvMass, TH1F* hBBInvMass)
{
  // particleVectors is organized as (m)tau (h)tau bjet (bjet)
  for (auto vec : particleVectors)
  {
    hEt->Fill(vec.Et());
    hEta->Fill(vec.Eta());
    hPhi->Fill(vec.Phi());
  }

  hTauTauInvMass->Fill((particleVectors[0] + particleVectors[1]).M());    // Filling the tau tau invariant mass is unrelated to whether or not a second b-jet exists

  if (particleVectors.size() == 4)
  {
    hBBInvMass->Fill((particleVectors[2] + particleVectors[3]).M());     // Only fill the BB invariant mass if there are TWO b-jets
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

  // Tests
  auto *hLeadingBPt = new TH1F("Unfitted Leading b-quark Transverse Momentum", "Unfitted Leading b-quark Transverse Momentum", 100, 0, 200);
  auto *hNextToLeadingBPt = new TH1F("Unfitted Next-To-Leading b-quark Transverse Momentum", "Unfitted Next-To-Leading b-quark Transverse Momentum", 100, 0, 200);
  auto *hAPt = new TH1F("Unfitted b-quark Mother Pseudoscalar Transverse Momentum", "Unfitted b-quark Mother Pseudoscalar Transverse Momentum", 100, 0, 200);
  auto *hLeadingBEt = new TH1F("Unfitted Leading b-quark Transverse Energy", "Unfitted Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hNextToLeadingBEt = new TH1F("Unfitted Next-To-Leading b-quark Transverse Energy", "Unfitted Next-To-Leading b-quark Transverse Energy", 100, 0, 200);
  auto *hAEt = new TH1F("Unfitted b-quark Mother Pseudoscalar Transverse Energy", "Unfitted b-quark Mother Pseudoscalar Transverse Energy", 100, 0, 200);
  auto *hLeadingBEta = new TH1F("Unfitted Leading b-quark Pseudorapidity", "Unfitted Leading b-quark Pseudorapidity", 100, -10, 10);
  auto *hNextToLeadingBEta = new TH1F("Unfitted Next-To-Leading b-quark Pseudorapidity", "Unfitted Next-To-Leading b-quark Pseudorapidity", 100, -10, 10);
  auto *hAEta = new TH1F("Unfitted b-quark Mother Pseudoscalar Pseudorapidity", "Unfitted b-quark Mother Pseudoscalar Pseudorapidity", 100, -10, 10);
  auto *hLeadingBPhi = new TH1F("Unfitted Leading b-quark Phi", "Unfitted Leading b-quark Phi", 100, -4, 4);
  auto *hNextToLeadingBPhi = new TH1F("Unfitted Next-To-Leading b-quark Phi", "Unfitted Next-To-Leading b-quark Phi", 100, -4, 4);
  auto *hAPhi = new TH1F("Unfitted b-quark Mother Pseudoscalar Phi", "Unfitted b-quark Mother Pseudoscalar Pseudorapidity", 100, -4, 4);

  // Set the addresses of the branches to elsewhere
  Float_t pt1, eta1, phi1, m1, pt2, eta2, phi2, m2, pt3, eta3, phi3, m3, pt4, eta4, phi4, m4, pt5, eta5, phi5, m5;
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
  tree->SetBranchAddress("pt_atobb", &pt5);
  tree->SetBranchAddress("eta_atobb", &eta5);
  tree->SetBranchAddress("phi_atobb", &phi5);
  tree->SetBranchAddress("m_atobb", &m5);

  // Loop through each event, perform necessary calculations, and fill the histograms
  for(int i = 0; i < tree->GetEntries(); i++)   // GetEntries() returns the # of entries in the branch
  {
    // std::cout << "Event #: " << i << std::endl;

    tree->GetEntry(i);

    TLorentzVector v1, v2, v3, v4, motherA, v1Fit, v2Fit, v3Fit, v4Fit;

    v1.SetPtEtaPhiM(pt1, eta1, phi1, m1);  // v1 contains the values of the (muonic) tau
    v2.SetPtEtaPhiM(pt2, eta2, phi2, m2);  // v2 contains the values of the (hadronic) tau
    v3.SetPtEtaPhiM(pt3, eta3, phi3, m3);  // v3 contains the values of the first b-jet

    std::vector<TLorentzVector> particleVectors = {v1, v2, v3};

    if (pt4 != -9999 && eta4 != -9999 && phi4 != -9999 && m4 != -9999)
    {
      v4.SetPtEtaPhiM(pt4, eta4, phi4, m4);  // v4 contains the values of the second b-jet ONLY IF IT EXISTS
      particleVectors.push_back(v4);

      motherA.SetPtEtaPhiM(pt5, eta5, phi5, m5);

      TLorentzVector leadingB, nextToLeadingB;
      if (v3.Pt() >= v4.Pt())
      {
        leadingB = v3;
        nextToLeadingB = v4;
      }
      else
      {
        leadingB = v4;
        nextToLeadingB = v3;
      }
      hLeadingBPt->Fill(leadingB.Pt());
      hNextToLeadingBPt->Fill(nextToLeadingB.Pt());
      hLeadingBEt->Fill(leadingB.Et());
      hNextToLeadingBEt->Fill(nextToLeadingB.Et());
      hLeadingBEta->Fill(leadingB.Eta());
      hNextToLeadingBEta->Fill(nextToLeadingB.Eta());
      hLeadingBPhi->Fill(leadingB.Phi());
      hNextToLeadingBPhi->Fill(nextToLeadingB.Phi());

      hAPt->Fill(motherA.Pt());
      hAEt->Fill(motherA.Et());
      hAEta->Fill(motherA.Eta());
      hAPhi->Fill(motherA.Phi());
    }

    fillHistograms(particleVectors, hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass);

    auto fittedParticleVectors = fitEvent(particleVectors);
    fillHistograms(fittedParticleVectors, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit);
  }

  histograms = {hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit, hLeadingBPt, hNextToLeadingBPt, hAPt, hLeadingBEt, hNextToLeadingBEt, hAEt, hLeadingBEta, hNextToLeadingBEta, hAEta, hLeadingBPhi, hNextToLeadingBPhi, hAPhi};
}

void KinFitOutputModule::drawHistograms()
{
  TFile* outputFile = new TFile("KinFitHistograms.root", "RECREATE");
  
  for (auto histogram : histograms)
  {
    histogram->Write();
  }

  std::cout << "Histograms written to KinFitHistograms.root\n";
  outputFile->Close();
  delete outputFile;
}