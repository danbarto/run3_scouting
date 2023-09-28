
import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8ConcurrentGeneratorFilter",
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13600.0),
                         PythiaParameters = cms.PSet(
                             pythia8CommonSettingsBlock,
                             pythia8CP5SettingsBlock,
                             processParameters = cms.vstring(
                                 "HiggsSM:gg2H = on",
                                 "25:m0 =125",
                                 "25:addChannel = 1 0.5 102 4900101 -4900101",
                                 "25:addChannel = 1 0.5 102 4900102 -4900102",
                                 "25:0:onMode=0",
                                 "25:1:onMode=0",
                                 "25:2:onMode=0",
                                 "25:3:onMode=0",
                                 "25:4:onMode=0",
                                 "25:5:onMode=0",
                                 "25:6:onMode=0",
                                 "25:7:onMode=0",
                                 "25:8:onMode=0",
                                 "25:9:onMode=0",
                                 "25:10:onMode=0",
                                 "25:11:onMode=0",
                                 "25:12:onMode=0",
                                 "25:13:onMode=0",
                                 "HiddenValley:Ngauge = 3",
                                 "HiddenValley:nFlav = 2",
                                 "HiddenValley:fragment = on",
                                 "HiddenValley:FSR = on",
                                 "HiddenValley:alphaOrder = 1",
                                 "HiddenValley:setLambda = on",
                                 "HiddenValley:Lambda = 30",
                                 "HiddenValley:pTminFSR = 33.0",
                                 "HiddenValley:spinFv=0",
                                 "HiddenValley:probVector=0.75",
                                 "HiddenValley:separateFlav = on",
                                 "HiddenValley:probKeepEta1 = 1.0",
                                 "4900101:m0 = 31.551517243857766",
                                 "4900102:m0 = 31.715482756142233",
                                 "4900111:m0 = 9.9",
                                 "4900211:m0 = 9.9",
                                 "4900221:m0 = 30",
                                 "4900113:m0 = 30",
                                 "4900213:m0 = 30",
                                 "4900113:addChannel = 1 1.00 91 4900211 -4900211",
                                 "4900213:addChannel = 1 1.00 91 4900211 4900111",
                                 "999999:all = GeneralResonance void 1 0 0 4.8 0.001 0. 0. 0.",
                                 "999999:addChannel = 1 0.0541 91 1 -1",
                                 "999999:addChannel = 1 0.217 91 2 -2",
                                 "999999:addChannel = 1 0.0541 91 3 -3",
                                 "999999:addChannel = 1 0.211 91 4 -4",
                                 "999999:addChannel = 1 0.162 91 11 -11",
                                 "999999:addChannel = 1 0.162 91 13 -13",
                                 "999999:addChannel = 1 0.139 91 15 -15",
                                 "999999:tau0 = 1.1e-05",
                                 "4900221:addChannel = 1 0.5 91 4900211 4900211 4900111",
                                 "4900221:addChannel = 1 0.5 91 -4900211 -4900211 4900111",
                                 "4900221:tau0 = 0.0",
                                 "4900111:addChannel = 1 1.0 91 999999 999999",
                                 "4900111:tau0 = 26.099999999999998",
                                 "4900211:onMode = 0",
                             ),
                             parameterSets = cms.vstring('pythia8CommonSettings',
                                                         'pythia8CP5Settings',
                                                         'processParameters',
                                                         )
                         )
                         )

# -- Require at least one muon in the final state. Muon from taus and HF decays are not considered.
MuFilter = cms.EDFilter("PythiaFilter",
    Status = cms.untracked.int32(0),
    MotherID = cms.untracked.int32(999999),
    MinPt = cms.untracked.double(2),
    ParticleID = cms.untracked.int32(13),
    MaxEta = cms.untracked.double(10),
    MinEta = cms.untracked.double(-10)
)
ProductionFilterSequence = cms.Sequence(generator*MuFilter)
