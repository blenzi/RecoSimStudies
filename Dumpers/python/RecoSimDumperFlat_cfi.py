import FWCore.ParameterSet.Config as cms

import Geometry.CaloEventSetup.caloTowerConstituents_cfi
CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
   MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
   )

recosimdumperFlat = cms.EDAnalyzer("RecoSimDumperFlat",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll"),
    vertexCollection                = cms.InputTag("offlinePrimaryVertices"),
    genParticleCollection           = cms.InputTag("genParticles",""),
    caloParticleCollection          = cms.InputTag("mix","MergedCaloTruth"),
    ebRechitCollection              = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection              = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),

    doCompression                   = cms.bool(True),  #do the compression of floats
    nBits                           = cms.int32(23),   #nbits for float compression (<=23)

    saveGenParticles                = cms.bool(True),  #save genParticles information
    saveCaloParticles               = cms.bool(True),  #save caloParticles information
    saveSimhits                     = cms.bool(True), #save simHits information
    saveRechits                     = cms.bool(True), #save recHits information

    genID                           = cms.vint32(22,11, -11), #save only caloParticles with this pdgId
    #genID                          = cms.vdouble(0),  #save only caloParticles with this pdgId

)
