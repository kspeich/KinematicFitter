from ROOT import TFile, TTree, TH1D, TLorentzVector
from array import array

inFile = TFile.Open("SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8.root")
inTree = inFile.Get("Events")           # Get the Events TTree
nentries = inTree.GetEntries()          # Number of entries in the Events TTree
# print("nentries={0:d}".format(nentries))

# Create the output file and tree
outFile = TFile("selectedEvents.root", 'recreate')
outTree = TTree("mutau_tree", "mutau_tree")

# Create arrays
# 1, 2: taus; 3, 4: b quarks
pt1 = array('f', [0])
eta1 = array('f', [0])
phi1 = array('f', [0])
m1 = array('f', [0])
pt2 = array('f', [0])
eta2 = array('f', [0])
phi2 = array('f', [0])
m2 = array('f', [0])
pt3 = array('f', [0])
eta3 = array('f', [0])
phi3 = array('f', [0])
m3 = array('f', [0])
pt4 = array('f', [0])
eta4 = array('f', [0])
phi4 = array('f', [0])
m4 = array('f', [0])
pt5 = array('f', [0])
eta5 = array('f', [0])
phi5 = array('f', [0])
m5 = array('f', [0])

# Create TBranches with the arrays
outTree.Branch("pt_1", pt1, "pt_1/F")
outTree.Branch("eta_1", eta1, "eta_1/F")
outTree.Branch("phi_1", phi1, "phi_1/F")
outTree.Branch("m_1", m1, "m_1/F")
outTree.Branch("pt_2", pt2, "pt_2/F")
outTree.Branch("eta_2", eta2, "eta_2/F")
outTree.Branch("phi_2", phi2, "phi_2/F")
outTree.Branch("m_2", m2, "m_2/F")
outTree.Branch("bpt_deepflavour_1", pt3, "bpt_deepflavour_1/F")
outTree.Branch("beta_deepflavour_1", eta3, "beta_deepflavour_1/F")
outTree.Branch("bphi_deepflavour_1", phi3, "bphi_deepflavour_1/F")
outTree.Branch("bm_deepflavour_1", m3, "bm_deepflavour_1/F")
outTree.Branch("bpt_deepflavour_2", pt4, "bpt_deepflavour_2/F")
outTree.Branch("beta_deepflavour_2", eta4, "beta_deepflavour_2/F")
outTree.Branch("bphi_deepflavour_2", phi4, "bphi_deepflavour_2/F")
outTree.Branch("bm_deepflavour_2", m4, "bm_deepflavour_2/F")
outTree.Branch("pt_atobb", pt5, "pt_atobb/F")
outTree.Branch("eta_atobb", eta5, "eta_atobb/F")
outTree.Branch("phi_atobb", phi5, "phi_atobb/F")
outTree.Branch("m_atobb", m5, "m_atobb/F")

tauMass = 1.77686
bMass = 4.18
aMass = 45.

passed = 0

for count, event in enumerate(inTree):             # count is the index, e is the event data
    # print("Event " + str(count) + " out of " + str(nentries))

    nGenPart = event.nGenPart

    pdgIdList = []
    motherList = []
    ptList = []
    etaList = []
    phiList = []
    massList = []

    for i in range(nGenPart):
        pdgIdList.append(event.GenPart_pdgId[i])
        motherList.append(event.GenPart_genPartIdxMother[i])
        ptList.append(event.GenPart_pt[i])
        etaList.append(event.GenPart_eta[i])
        phiList.append(event.GenPart_phi[i])
        massList.append(event.GenPart_mass[i])
    
    # Lists containing indices
    higgsList = []              # All Higgs
    aList = []                  # All pseudoscalars that are from a Higgs
    abList = []                 # All pseudoscalars that are from a Higgs and decay to bb
    tauList = []                # All taus that are from a pseudoscalar from Higgs
    bList = []                  # All b quarks that are from a pseudoscalar from Higgs

    # Event selection

    # First selection: Higgs
    for i in range(nGenPart):                
        if (pdgIdList[i] == 25):
            higgsList.append(i)
    
    # Second selection: Higgs must decay to two pseudoscalars
    for higgs in higgsList:
        for i in range(nGenPart):
            if (motherList[i] == higgs and pdgIdList[i] == 36):
                aList.append(i)

    # Third selection: pseudoscalars decay into b quarks and taus
    for a in aList:
        for i in range(nGenPart):
            if (motherList[i] == a and abs(pdgIdList[i]) == 15):
                tauList.append(i)
            elif (motherList[i] == a and abs(pdgIdList[i]) == 5):
                bList.append(i)
                abList.append(a)
    
    # Fourth selection: 2 b quarks and 2 taus
    if (len(tauList) == 2 and len(bList) == 2):
        passed += 1

        # Add elements to each list
        pt1[0] = ptList[tauList[0]]
        eta1[0] = etaList[tauList[0]]
        phi1[0] = phiList[tauList[0]]
        m1[0] = tauMass
        pt2[0] = ptList[tauList[1]]
        eta2[0] = etaList[tauList[1]]
        phi2[0] = phiList[tauList[1]]
        m2[0] = tauMass
        pt3[0] = ptList[bList[0]]
        eta3[0] = etaList[bList[0]]
        phi3[0] = phiList[bList[0]]
        m3[0] = bMass
        pt4[0] = ptList[tauList[1]]
        eta4[0] = etaList[tauList[1]]
        phi4[0] = phiList[tauList[1]]
        m4[0] = bMass
        pt5[0] = ptList[abList[0]]
        eta5[0] = etaList[abList[0]]
        phi5[0] = phiList[abList[0]]
        m5[0] = aMass

        outTree.Fill()

print(str(passed) + "/" + str(nentries) + " events passed")

outFile.Write()
outFile.Close()