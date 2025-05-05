import FWCore.ParameterSet.Config as cms

EMJMCNtuplizer = cms.EDAnalyzer(
    'EMJMCNtuplizer',
    PFCandidateToken = cms.InputTag("particleFlow"),
    PFJetToken = cms.InputTag("ak4PFJetsPuppi"),
    GenJetToken = cms.InputTag("ak4GenJets"),
    GenParticleToken = cms.InputTag("genParticles"),
    hbRecHitsToken = cms.InputTag("reducedHcalRecHits", "hbhereco")
)
