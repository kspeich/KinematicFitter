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
  auto diTauPart = event.getParticles(15, "ditau")[0];
  auto bJet1Part = event.getParticles(5)[0];
  auto diTau = convertParticle(diTauPart);
  auto bJet1 = convertParticle(bJet1Part);

  std::vector<TFitParticleEtEtaPhi*> particles = {diTau, bJet1};
  std::vector<TFitConstraintM*> constraints;

  if (event.getNumParticles(5) == 2)  // Do everything including the second b-jet
  {
    auto bJet2Part = event.getParticles(5)[1];
    auto bJet2 = convertParticle(bJet2Part);

    // bJet1 and bJet2 must make an a pseudoscalar: this only happens when there is a second b-jet
    TFitConstraintM *mCons1 = new TFitConstraintM("AMassConstraint2", "AMass-Constraint2", 0, 0, 45.);
    mCons1->addParticles1(bJet1, bJet2);

    // All four particles must make a Higgs
    TFitConstraintM *mCons2 = new TFitConstraintM( "HiggsMassConstraint", "HiggsMass-Constraint", 0, 0, 125.);
    mCons2->addParticles1(diTau, bJet1, bJet2);

    particles.push_back(bJet2);
    constraints.push_back(mCons1);
    constraints.push_back(mCons2);
  }
  else    // Since we aren't fitting the ditau, and there is no constraint since there is no second b-jet to put a constraint on, just return the two objects
  {
    return event;
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
  params.addParticle(diTauPart);
  
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

    bJet1Vec.SetPtEtaPhiM(pt1, eta1, phi1, m1);
    diTauVec.SetPtEtaPhiM(pt3, eta3, phi3, m3);

    // The PDGIDs here represent that the particles are taus and b's, but do not indicate whether they are anti-taus or anti-b's
    particles.addParticle(diTauVec, 15, {"ditau"});
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