import FWCore.ParameterSet.Config as cms

JetInfoNtuplizer = cms.EDAnalyzer(
    'JetInfoNtuplizer',
    PFCandidateToken = cms.InputTag("particleFlow"),
    PFJetToken = cms.InputTag("ak4PFJetsPuppi"),
    GenJetToken = cms.InputTag("ak4GenJets"),
    GenParticleToken = cms.InputTag("genParticles"),
)
