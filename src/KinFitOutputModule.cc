#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"

#include <iostream>
#include <vector>

KinFitOutputModule::KinFitOutputModule(TTree* iTree, bool iGenEvents) :
  tree(iTree),
  genEvents(iGenEvents)
{}

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

std::vector<const TMatrixD*> KinFitOutputModule::fitEvent(std::vector<TLorentzVector> particleVectors)
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
  TFitConstraintM *mCons1 = new TFitConstraintM( "WMassConstraint", "WMass-Constraint", 0, 0 , 45.);
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
    TFitConstraintM *mCons2 = new TFitConstraintM( "TopMassConstraint", "TopMass-Constraint", 0, 0 , 45.);
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
  std::cout << "Performing kinematic fit..." << std::endl;
  print(fitter);
  fitter->fit();
  std::cout << "Done." << std::endl;
  print(fitter);
  
  std::vector<const TMatrixD*> params = {};
  
  for (auto particle : particles)
  {
    params.push_back(particle->getParCurr());
    //delete particle;
  }
  
  for (auto constraint : constraints)
  {
    delete constraint;
  }

  delete fitter;

  return params;
}

int KinFitOutputModule::findParent(std::vector<std::vector<Float_t>> rawData, int rowNum, int motherPdgId)
{
  // row = {instance, event, pdgId, motherInstance, pt, eta, phi, m}
  // so event is the 1st (2nd) element in the row motherInstance is the 3rd (4th) element in the row
  auto event = rawData[rowNum][1];
  auto motherInstance = rawData[rowNum][3];

  for (unsigned long int motherRow = 0; motherRow < rawData.size(); motherRow++)
  {
    if (rawData[motherRow][1] == event && rawData[motherRow][3] == motherInstance && rawData[motherRow][4] == motherPdgId)
    {
      return motherRow;
    }
  }

  return -9999;
}

void KinFitOutputModule::makeHistograms()
{
  // Initialize the unfitted histograms
  auto *hEt = new TH1F("Unfitted Transverse Energy", "Unfitted Transverse Energy", 50, 0, 200);
  auto *hEta = new TH1F("Unfitted Eta", "Unfitted Eta", 50, -10, 10);
  auto *hPhi = new TH1F("Unfitted Phi", "Unfitted Phi", 50, -4, 4);
  auto *hTauTauInvMass = new TH1F("Unfitted Tau Tau Invariant Mass", "Unfitted Tau Tau Invariant Mass", 50, 0, 200);
  auto *hBBInvMass = new TH1F("Unfitted BB Invariant Mass", "Unfitted BB Invariant Mass", 50, 0, 200);

  // Initialize the fitted histograms
  auto *hEtFit = new TH1F("Fitted Transverse Energy", "Fitted Transverse Energy", 50, 0, 200);
  auto *hEtaFit = new TH1F("Fitted Eta", "Fitted Eta", 50, -10, 10);
  auto *hPhiFit = new TH1F("Fitted Phi", "Fitted Phi", 50, -4, 4);
  auto *hTauTauInvMassFit = new TH1F("Fitted Tau Tau Invariant mass", "Fitted Tau Tau Invariant Mass", 50, 0, 200);
  auto *hBBInvMassFit = new TH1F("Fitted BB Invariant Mass", "Fitted BB Invariant Mass", 50, 0, 200);

  if (!genEvents)
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
    tree->SetBranchAddress("bpt_deepflavour_2", &pt4);
    tree->SetBranchAddress("beta_deepflavour_2", &eta4);
    tree->SetBranchAddress("bphi_deepflavour_2", &phi4);

    // Loop through each event, perform necessary calculations, and fill the histograms
    for(int i = 0; i < tree->GetEntries(); i++)   // GetEntries() returns the # of entries in the branch
    {
      std::cout << "Event #: " << i << std::endl;

      tree->GetEntry(i);

      TLorentzVector v1, v2, v3, v4, v1Fit, v2Fit, v3Fit, v4Fit;

      v1.SetPtEtaPhiM(pt1, eta1, phi1, m1);  // v1 contains the values of the muonic tau
      v2.SetPtEtaPhiM(pt2, eta2, phi2, m2);  // v2 contains the values of the hadronic tau
      v3.SetPtEtaPhiM(pt3, eta3, phi3, m3);  // v3 contains the values of the first b-jet

      std::vector<TLorentzVector> particleVectors = {v1, v2, v3};

      if (pt4 != -9999 && eta4 != -9999 && phi4 != -9999 && m4 != -9999)
      {
        v4.SetPtEtaPhiM(pt4, eta4, phi4, m4);  // v4 contains the values of the second b-jet ONLY IF IT EXISTS
        particleVectors.push_back(v4);

        hBBInvMass->Fill((v3 + v4).M());     // Only fill the BB invariant mass if there are TWO b-jets
      }
    
      hTauTauInvMass->Fill((v1 + v2).M());  // Filling the tau tau invariant mass is unrelated to whether or not a second b-jet exists

      for (auto vec : particleVectors)       // Don't fill v4 into any histogram if it does not exist
      {
        hEt->Fill(vec.Et());
        hEta->Fill(vec.Eta());
        hPhi->Fill(vec.Phi());
      }

      auto fittedValues = fitEvent(particleVectors);

      for (const TMatrixD* vecFit : fittedValues)
      {
        hEtFit->Fill(vecFit->GetMatrixArray()[0]);
        hEtaFit->Fill(vecFit->GetMatrixArray()[1]);
        hPhiFit->Fill(vecFit->GetMatrixArray()[2]);
      }
    }
  }
  
  else          // Running on generated samples => must consider the pdgIds of the mother particles
  {
    // Set the addresses of the branches to elsewhere
    Int_t instance, event, pdgId, motherInstance;
    Float_t pt, eta, phi, m;
    tree->SetBranchAddress("Instance", &instance);
    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("GenPart_pdgId", &pdgId);
    tree->SetBranchAddress("GenPart_genPartIdxMother", &motherInstance);
    tree->SetBranchAddress("GenPart_pt", &pt);
    tree->SetBranchAddress("GenPart_eta", &eta);
    tree->SetBranchAddress("GenPart_phi", &phi);
    tree->SetBranchAddress("GenPart_mass", &m);
    
    // Collect the data from the trees/branches into a readable format
    std::vector<std::vector<Float_t>> rawData;

    for(int i = 0; i < tree->GetEntries(); i++)   // GetEntries() returns the # of entries in the branch
    {
      tree->GetEntry(i);

      if (abs(pdgId) == 5 || abs(pdgId) == 15 || pdgId == 25 || pdgId == 36)    // Only include b, tau, a, and H
      {
        std::vector<Float_t> row = {static_cast<Float_t>(instance), static_cast<Float_t>(event), static_cast<Float_t>(pdgId), static_cast<Float_t>(motherInstance), pt, eta, phi, m};
        rawData.push_back(row);
      }
    }

    // Keep the useful data
    std::vector<std::tuple<Float_t, Float_t, TLorentzVector>> importantData;              // event #, pdgId, and TLorentzVector
    std::vector<Float_t> eventNums;                                                       // just the event #

    for (unsigned long int i = 0; i < rawData.size(); i++)
    {
      if (abs(rawData[i][2]) == 5 || abs(rawData[i][2] == 15))                            // b and tau
      {
        if (findParent(rawData, i, 36) != -9999)                                          // has pseudoscalar as parent
        {
          if (findParent(rawData, findParent(rawData, i, 36), 25) != -9999)               // which has a higgs as a parent
          {
            TLorentzVector vec; 
            vec.SetPtEtaPhiM(rawData[i][4], rawData[i][5], rawData[i][6], rawData[i][7]);

            std::tuple<Float_t, Float_t, TLorentzVector> particleData (rawData[i][1], rawData[i][2], vec);
            importantData.push_back(particleData);
            eventNums.push_back(rawData[i][1]);
          }
        }
      }
    }

    // Now organize the useful data by event
    std::map<Float_t, std::vector<std::pair<Float_t, TLorentzVector>>> events;
    for (unsigned long int i = 0; i < importantData.size(); i++)
    {
      auto [eventNum, pdgId, vec] = importantData[i];

      if (events.count(eventNum) == 0)            // std::map::count() counts the # of times a key exists in the map
      {
        std::vector<std::pair<Float_t, TLorentzVector>> v;
        events.insert({eventNum, v});
      }

      events[eventNum].push_back(std::make_pair(pdgId, vec));
    }

    // Showtime
    for (auto eventNum : eventNums)
    {
      std::vector<TLorentzVector> particles, bJets, taus;
      for (auto particle : events[eventNum])
      {
        particles.push_back(particle.second);
        if (abs(particle.first) == 5)
        {
          bJets.push_back(particle.second);
        }
        else if (abs(particle.first) == 15)
        {
          taus.push_back(particle.second);
        }
      }

      for (auto particle : particles)
      {
        hEt->Fill(particle.Et()); 
        hEta->Fill(particle.Eta());
        hPhi->Fill(particle.Phi());
      }

      TLorentzVector bSum, tauSum;
      for (auto bJet : bJets)
      {
        bSum += bJet;
      }
      for (auto tau : taus)
      {
        tauSum += tau;
      }
      
      hBBInvMass->Fill(bSum.M());
      hTauTauInvMass->Fill(tauSum.M());

      // Run the fitter
      auto fittedValues = fitEvent(particles);
      for (const TMatrixD* vecFit : fittedValues)
      {
        hEtFit->Fill(vecFit->GetMatrixArray()[0]);
        hEtaFit->Fill(vecFit->GetMatrixArray()[1]);
        hPhiFit->Fill(vecFit->GetMatrixArray()[2]);
      }
    }
  }  

  histograms = {hEt, hEta, hPhi, hTauTauInvMass, hBBInvMass, hEtFit, hEtaFit, hPhiFit, hTauTauInvMassFit, hBBInvMassFit};
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