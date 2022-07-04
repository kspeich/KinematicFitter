from ROOT import TFile, TTree, TH1D, TLorentzVector
from array import array

inFile = TFile.Open("root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18NanoAODv7/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/AE82D64D-12A8-B547-9A00-339213036849.root")
inTree = inFile.Get("Events")           # Get the Events TTree
nentries = inTree.GetEntries()          # Number of entries in the Events TTree
# print("nentries={0:d}".format(nentries))

# Create the output file and tree
outFile = TFile("selectedEvents.root", 'recreate')
outTree = TTree("mutau_tree", "mutau_tree")

# Create arrays => these will be used to create the TBranches later
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

for count, event in enumerate(inTree):             # count is the index, e is the event data
    pdgIdList = event.GenPart_pdgId
    motherList = event.GenPart_genPartIdxMother
    ptList = event.GenPart_pt
    etaList = event.GenPart_eta
    phiList = event.GenPart_phi
    massList = event.GenPart_mass

    nHiggs = 0                  # number of Higgs

    # Lists containing indices
    higgsList = []              # All Higgs
    aList = []                  # All pseudoscalars that are from a Higgs
    tauList = []                # All taus that are from a pseudoscalar from Higgs
    bList = []                  # All b quarks that are from a pseudoscalar from Higgs

    # Event selection

    # First selection: Higgs
    for pdgId in pdgIdList:                
        if (abs(pdgId) == 25):
            nHiggs += 1
            higgsList.append(pdgIdList.index(pdgId))

    # Second selection: Higgs must decay to two pseudoscalars
    for higgs in higgsList:
        for mother in motherList:
            if (mother == higgs and pdgIdList[motherList.index(mother)] == 36):
                aList.append(motherList.index(mother))

    # Third selection: pseudoscalars decay into b quarks and taus
    for a in aList:
        for mother in motherList:
            if (mother == a and abs(pdgIdList[motherList.index(mother)]) == 15):
                tauList.append(motherList.index(mother))
            elif (mother == a and abs(pdgIdList[motherList.index(mother)]) == 5):
                bList.append(motherList.index(mother))

    # Fourth selection: 2 b quarks and 2 taus
    if (len(tauList) == 2 and len(bList) == 2):
        # Add elements to each list
        pt1[0] = ptList[tauList[0]]
        eta1[0] = etaList[tauList[0]]
        phi1[0] = phiList[tauList[0]]
        m1[0] = massList[tauList[0]]
        pt2[0] = ptList[tauList[1]]
        eta2[0] = etaList[tauList[1]]
        phi2[0] = phiList[tauList[1]]
        m2[0] = massList[tauList[1]]
        pt3[0] = ptList[bList[0]]
        eta3[0] = etaList[bList[0]]
        phi3[0] = phiList[bList[0]]
        m3[0] = massList[bList[0]]
        pt4[0] = ptList[tauList[1]]
        eta4[0] = etaList[tauList[1]]
        phi4[0] = phiList[tauList[1]]
        m4[0] = massList[tauList[1]]
        
        outTree.Fill()

outFile.Write()
outFile.Close()        
