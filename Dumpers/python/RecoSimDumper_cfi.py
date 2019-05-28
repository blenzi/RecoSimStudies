import FWCore.ParameterSet.Config as cms

recosimdumper = cms.EDAnalyzer("RecoSimDumper",

    genParticleCollection             = cms.InputTag("genParticles","","HLT"),
    caloParticleCollection            = cms.InputTag("mix","MergedCaloTruth","HLT"),
    ebRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection                = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection               = cms.InputTag("particleFlowEGamma","EBEEClusters","RECO"),
    ebSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RECO"), 
    eeSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RECO"), 
    
    useRechits                        = cms.bool(True),
    usePFRechits                      = cms.bool(True),
    usePFCluster                      = cms.bool(True),
    useSuperCluster                   = cms.bool(True),
    motherID                          = cms.int32(22)
    #motherID                          = cms.int32(0)
)
