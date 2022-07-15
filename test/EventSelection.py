from ROOT import TFile, TTree, TH1D, TLorentzVector
from array import array

inSignal = "SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_FilterTauTauTrigger_TuneCP5_13TeV_madgraph_pythia8.root"
inDYJets = "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root"
inTTLeptonic = "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"
inTTSemiLeptonic = "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root"

outSignal = "signalSelectedEvents.root"
outDYJets = "dyJetsSelectedEvents.root"
outTTLeptonic = "ttLeptonicSelectedEvents.root"
outTTSemiLeptonic = "ttSemiLeptonicSelectedEvents.root"

inFile = TFile.Open(inDYJets)
inTree = inFile.Get("Events")           # Get the Events TTree
nentries = inTree.GetEntries()          # Number of entries in the Events TTree

# Create the output file and tree
outFile = TFile(outDYJets, 'recreate')
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

tauMass = 1.77686
bMass = 4.18
aMass = 45.

passed = 0

for count, event in enumerate(inTree):             # count is the index, e is the event data
    if (count % 1000 == 0):
        print("Event " + str(count) + " / " + str(nentries))

    nGenPart = event.nGenPart

    pdgIdList = []
    motherList = []
    ptList = []
    etaList = []
    phiList = []
    massList = []
    statusFlagsList = []

    for i in range(nGenPart):
        pdgIdList.append(event.GenPart_pdgId[i])
        motherList.append(event.GenPart_genPartIdxMother[i])
        ptList.append(event.GenPart_pt[i])
        etaList.append(event.GenPart_eta[i])
        phiList.append(event.GenPart_phi[i])
        massList.append(event.GenPart_mass[i])
        statusFlagsList.append(bin(event.GenPart_statusFlags[i]))           # Convert the integer to a binary string
    
    # Lists containing indices
    higgsList = []              # All Higgs
    aList = []                  # All pseudoscalars that are from a Higgs
    tauList = []                # All taus that are from a pseudoscalar from Higgs
    bList = []                  # All b quarks that are from a pseudoscalar from Higgs

    # Status Flags
    isFirstCopy = 12
    fromHardProcess = 8

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

    # Third selection: pseudoscalars decay into b quarks and taus, with the status flag isFirstCopy (12) and fromHardProcess(8)
    for a in aList:
        for i in range(nGenPart):
            if (len(statusFlagsList[i]) >= isFirstCopy + 2):
                if (statusFlagsList[i][-1 * (isFirstCopy + 1)] == "1" and statusFlagsList[i][-1 * (fromHardProcess + 1)] == "1"):
                    if (motherList[i] == a and abs(pdgIdList[i]) == 15):
                        tauList.append(i)
                    elif (motherList[i] == a and abs(pdgIdList[i]) == 5):
                        bList.append(i)
    
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

        outTree.Fill()

print(str(passed) + "/" + str(nentries) + " events passed")

outFile.Write()
outFile.Close()