import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('DrellYann',
                      verbosity = cms.bool(False),
                      isMC = cms.bool(True),
                      vertexSrc = cms.InputTag('selectedPrimaryVertices'),
                      muonSrc = cms.InputTag('slimmedMuons')
                      )
