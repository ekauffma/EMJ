import FWCore.ParameterSet.Config as cms

HcalDepthEnergyFractionNtuplizer = cms.EDAnalyzer(
    'HcalDepthEnergyFractionNtuplizer',
    PackedCandidateToken = cms.InputTag("packedPFCandidates"),
    HcalDepthEFToken = cms.InputTag("packedPFCandidates", "hcalDepthEnergyFractions")
)
