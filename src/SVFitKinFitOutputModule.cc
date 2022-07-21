#include "PhysicsTools/KinFitter/interface/SVFitKinFitOutputModule.hh"
#include "PhysicsTools/KinFitter/interface/KinFitOutputModule.hh"

SVFitKinFitOutputModule::SVFitKinFitOutputModule(TTree* iTree, std::string iOutputFile) : 
  KinFitOutputModule(iTree, iOutputFile)
{}

void SVFitKinFitOutputModule::run()
{
  runFitter();
  makeHistograms();
  drawHistograms();
}

Particles SVFitKinFitOutputModule::fitEvent(Particles event)
{
  auto diTauVec = event.getParticleVectors(15)[0];
  auto bJet1Vec = event.getParticleVectors(5)[0];

  TMatrixD m1(3,3);
  TMatrixD m2(3,3);
  m1.Zero();
  m2.Zero();

  m1(0,0) = ErrEt (diTauVec); // et
  m1(1,1) = ErrEta(diTauVec); // eta
  m1(2,2) = ErrPhi(diTauVec); // phi
  m2(0,0) = ErrEt (bJet1Vec); // et
  m2(1,1) = ErrEta(bJet1Vec); // eta
  m2(2,2) = ErrPhi(bJet1Vec); // phi

  TFitParticleEtEtaPhi *diTau = new TFitParticleEtEtaPhi("diTau", "diTau", &diTauVec, &m1);
  TFitParticleEtEtaPhi *bJet1 = new TFitParticleEtEtaPhi("bJet1", "bJet1", &bJet1Vec, &m2);

  // diTau must make an a pseudoscalar
  TFitConstraintM *mCons1 = new TFitConstraintM("AMassConstraint1", "AMass-Constraint1", 0, 0, 45.);
  //mCons1->addParticles1(diTau);

  std::vector<TFitParticleEtEtaPhi*> particles = {diTau, bJet1};
  std::vector<TFitConstraintM*> constraints;

  if (event.getNumParticles(5) == 2)  // Do everything including the second b-jet
  {
    auto bJet2Vec = event.getParticleVectors(5)[1];
    
    TMatrixD m3(3,3);
    m3.Zero();
    
    m3(0,0) = ErrEt (bJet2Vec); // et
    m3(1,1) = ErrEta(bJet2Vec); // eta
    m3(2,2) = ErrPhi(bJet2Vec); // phi

    TFitParticleEtEtaPhi *bJet2 = new TFitParticleEtEtaPhi("bJet2", "bJet2", &bJet2Vec, &m3);

    // bJet1 and bJet2 must make an a pseudoscalar: this only happens when there is a second b-jet
    TFitConstraintM *mCons2 = new TFitConstraintM("AMassConstraint2", "AMass-Constraint2", 0, 0, 45.);
    mCons2->addParticles1(bJet1, bJet2);

    // All four particles must make a Higgs
    TFitConstraintM *mCons3 = new TFitConstraintM( "HiggsMassConstraint", "HiggsMass-Constraint", 0, 0, 125.);
    mCons3->addParticles1(diTau, bJet1, bJet2);

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

void SVFitKinFitOutputModule::runFitter()
{
  // Set the addresses of the branches to elsewhere
  Float_t pt1, eta1, phi1, m1, pt2, eta2, phi2, m2, m3, pt3, eta3, phi3;
  tree->SetBranchAddress("bpt_deepflavour_1", &pt1);
  tree->SetBranchAddress("beta_deepflavour_1", &eta1);
  tree->SetBranchAddress("bphi_deepflavour_1", &phi1);
  tree->SetBranchAddress("bm_deepflavour_1", &m1);
  tree->SetBranchAddress("bpt_deepflavour_2", &pt2);
  tree->SetBranchAddress("beta_deepflavour_2", &eta2);
  tree->SetBranchAddress("bphi_deepflavour_2", &phi2);
  tree->SetBranchAddress("bm_deepflavour_2", &m2);
  tree->SetBranchAddress("m_sv", &m3);
  tree->SetBranchAddress("pt_sv", &pt3);
  tree->SetBranchAddress("eta_sv", &eta3);
  tree->SetBranchAddress("phi_sv", &phi3);

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

    TLorentzVector bJet1Vec, bJet2Vec, diTauVec;
    Particles particles;

    bJet1Vec.SetPtEtaPhiM(pt1, eta1, phi1, m1);  // v1 contains the values of the (muonic) tau
    diTauVec.SetPtEtaPhiM(pt3, eta3, phi3, m3);  // v3 contains the values of the first b-jet

    // The PDGIDs here represent that the particles are taus and b's, but do not indicate whether they are anti-taus or anti-b's
    particles.addParticle(diTauVec, 15);
    particles.addParticle(bJet1Vec, 5);

    if (pt2 != -9999 && eta2 != -9999 && phi2 != -9999 && m2 != -9999)
    {
      bJet2Vec.SetPtEtaPhiM(pt2, eta2, phi2, m2);  // v4 contains the values of the second b-jet ONLY IF IT EXISTS
      particles.addParticle(bJet2Vec, 5);
    }

    unfittedEvents.push_back(particles);
    fittedEvents.push_back(fitEvent(particles));
  }
}