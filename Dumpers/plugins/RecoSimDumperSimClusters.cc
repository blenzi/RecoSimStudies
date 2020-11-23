#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// #include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

// #include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
// #include "FWCore/Framework/interface/LooperFactory.h"
// #include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

// #include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
// #include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// #include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
// #include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Math/interface/libminifloat.h"

// #include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
// #include "DataFormats/EgammaCandidates/interface/Conversion.h"
// #include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
// #include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
// #include "DataFormats/EgammaReco/interface/ElectronSeed.h"
// #include "DataFormats/EgammaReco/interface/SuperCluster.h"
// #include "DataFormats/EgammaCandidates/interface/Photon.h"
// #include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
// #include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
// #include "FWCore/Utilities/interface/isFinite.h"
// #include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
// #include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoSimStudies/Dumpers/plugins/RecoSimDumperSimClusters.h"

#include <iostream>
#include <vector>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <time.h>

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

//
// constructors and destructor
//
RecoSimDumperSimClusters::RecoSimDumperSimClusters(const edm::ParameterSet& iConfig)
{

   vtxToken_                      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   rhoToken_                      = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
   genToken_                      = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   caloPartToken_                 = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   ebRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   doCompression_                 = iConfig.getParameter<bool>("doCompression");
   nBits_                         = iConfig.getParameter<int>("nBits");
   saveGenParticles_              = iConfig.getParameter<bool>("saveGenParticles");
   saveCaloParticles_             = iConfig.getParameter<bool>("saveCaloParticles");
   saveSimhits_             	  = iConfig.getParameter<bool>("saveSimhits");
   saveOnlyCaloSimhits_        	  = iConfig.getParameter<bool>("saveOnlyCaloSimhits");
   saveRechits_                   = iConfig.getParameter<bool>("saveRechits");
   genID_                         = iConfig.getParameter<std::vector<int>>("genID");

   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree");
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   tree->Branch("rho", &rho, "rho/F");
   if(saveGenParticles_){
      tree->Branch("genParticle_id","std::vector<int>",&genParticle_id);
      tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
      tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
      tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
      tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
   }
   if(saveCaloParticles_){
      tree->Branch("caloParticle_pdgId","std::vector<int>",&caloParticle_pdgId);
      tree->Branch("caloParticle_id","std::vector<int>",&caloParticle_id);
      tree->Branch("caloParticle_genEnergy","std::vector<float>",&caloParticle_genEnergy);
      tree->Branch("caloParticle_simEnergy","std::vector<float>",&caloParticle_simEnergy);
      tree->Branch("caloParticle_genPt","std::vector<float>",&caloParticle_genPt);
      tree->Branch("caloParticle_simPt","std::vector<float>",&caloParticle_simPt);
      tree->Branch("caloParticle_genEta","std::vector<float>",&caloParticle_genEta);
      tree->Branch("caloParticle_simEta","std::vector<float>",&caloParticle_simEta);
      tree->Branch("caloParticle_genPhi","std::vector<float>",&caloParticle_genPhi);
      tree->Branch("caloParticle_simPhi","std::vector<float>",&caloParticle_simPhi);
      tree->Branch("caloParticle_simIeta","std::vector<int>",&caloParticle_simIeta);
      tree->Branch("caloParticle_simIphi","std::vector<int>",&caloParticle_simIphi);
      tree->Branch("caloParticle_simIz","std::vector<int>",&caloParticle_simIz);
      tree->Branch("caloParticle_nSimClusters","std::vector<int>",&caloParticle_nSimClusters);
      tree->Branch("simCluster_pdgId", "std::vector<int>", &simCluster_pdgId);
      tree->Branch("simCluster_energy", "std::vector<float>", &simCluster_energy);
      tree->Branch("simCluster_eta", "std::vector<float>", &simCluster_eta);
      tree->Branch("simCluster_phi", "std::vector<float>", &simCluster_phi);

   }
   if(saveSimhits_){
      tree->Branch("simHit_id","std::vector<uint32_t>",&simHit_id);
      tree->Branch("simHit_energy","std::vector<float>",&simHit_energy);
      tree->Branch("simHit_eta","std::vector<float>",&simHit_eta);
      tree->Branch("simHit_phi","std::vector<float>",&simHit_phi);
      tree->Branch("simHit_ieta","std::vector<int>",&simHit_ieta);
      tree->Branch("simHit_iphi","std::vector<int>",&simHit_iphi);
      tree->Branch("simHit_iz","std::vector<int>",&simHit_iz);
   }
   if(saveRechits_){
      tree->Branch("recHit_energy","std::vector<float>",&recHit_energy);
      tree->Branch("recHit_eta","std::vector<float>",&recHit_eta);
      tree->Branch("recHit_phi","std::vector<float>",&recHit_phi);
      tree->Branch("recHit_ieta","std::vector<int>",&recHit_ieta);
      tree->Branch("recHit_iphi","std::vector<int>",&recHit_iphi);
      tree->Branch("recHit_iz","std::vector<int>",&recHit_iz);
   }
}

RecoSimDumperSimClusters::~RecoSimDumperSimClusters()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void RecoSimDumperSimClusters::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{

   edm::Handle<double> rhos;
   ev.getByToken(rhoToken_,rhos);
   if (!rhos.isValid()) {
       std::cerr << "Analyze --> rhos not found" << std::endl;
       return;
   }

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl;
       return;
   }

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   ev.getByToken(genToken_,genParticles);
   if (!genParticles.isValid()) {
       std::cerr << "Analyze --> genParticles not found" << std::endl;
       return;
   }

   edm::Handle<std::vector<CaloParticle> > caloParticles;
   ev.getByToken(caloPartToken_,caloParticles);
   if (!caloParticles.isValid()) {
       std::cerr << "Analyze --> caloParticles not found" << std::endl;
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEB;
   if(saveRechits_) {
      ev.getByToken(ebRechitToken_, recHitsEB);
      if (!recHitsEB.isValid()) {
          std::cerr << "Analyze --> recHitsEB not found" << std::endl;
          return;
      }
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   if(saveRechits_) {
      ev.getByToken(eeRechitToken_, recHitsEE);
      if (!recHitsEE.isValid()) {
          std::cerr << "Analyze --> recHitsEE not found" << std::endl;
          return;
      }
   }

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();

   nVtx = vertices->size();
   rho = *(rhos.product());

   genParticle_id.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();

   caloParticle_id.clear();
   caloParticle_pdgId.clear();
   caloParticle_genEnergy.clear();
   caloParticle_simEnergy.clear();
   caloParticle_genPt.clear();
   caloParticle_simPt.clear();
   caloParticle_genEta.clear();
   caloParticle_simEta.clear();
   caloParticle_genPhi.clear();
   caloParticle_simPhi.clear();
   caloParticle_simIeta.clear();
   caloParticle_simIphi.clear();
   caloParticle_simIz.clear();

   caloParticle_nSimClusters.clear();
   simCluster_pdgId.clear();
   simCluster_energy.clear();
   simCluster_eta.clear();
   simCluster_phi.clear();

   simHit_id.clear();
   simHit_energy.clear();
   simHit_eta.clear();
   simHit_phi.clear();
   simHit_ieta.clear();
   simHit_iphi.clear();
   simHit_iz.clear();

   recHit_energy.clear();
   recHit_eta.clear();
   recHit_phi.clear();
   recHit_ieta.clear();
   recHit_iphi.clear();
   recHit_iz.clear();

   std::vector<GenParticle> genParts;
   for(const auto& iGen : *(genParticles.product()))
   {
       bool isGoodParticle = false;
       for(unsigned int id=0; id<genID_.size(); id++)
           if((iGen.pdgId()==genID_.at(id) || genID_.at(id)==0) && iGen.status()==1) isGoodParticle=true;

       if(!isGoodParticle) continue;
       genParticle_id.push_back(iGen.pdgId());
       genParticle_energy.push_back(iGen.energy());
       genParticle_pt.push_back(iGen.pt());
       genParticle_eta.push_back(iGen.eta());
       genParticle_phi.push_back(iGen.phi());

       genParts.push_back(iGen);
   }

   for (const auto& caloParticle : *(caloParticles.product()))
   {
       bool isGoodParticle = false;
       for(unsigned int id=0; id<genID_.size(); id++)
           if (caloParticle.pdgId()==genID_.at(id) || genID_.at(id)==0) isGoodParticle=true;

       if(!isGoodParticle) continue;

       caloParticle_pdgId.push_back(caloParticle.pdgId());
       caloParticle_nSimClusters.push_back(caloParticle.simClusters().size());
       for(unsigned int iSC = 0; iSC < caloParticle.simClusters().size() ; iSC++){
           const auto& simCluster = caloParticle.simClusters()[iSC];
           simCluster_pdgId.push_back(simCluster->pdgId());
           simCluster_energy.push_back(simCluster->energy());
           simCluster_eta.push_back(simCluster->eta());
           simCluster_phi.push_back(simCluster->phi());

           for (const auto& hit : simCluster->hits_and_energies())
           {
              DetId id(hit.first);
              if(saveOnlyCaloSimhits_ && id.subdetId()!=EcalBarrel && id.subdetId()!=EcalEndcap) continue;

              simHit_id.push_back(hit.first);
              simHit_energy.push_back(hit.second);
           }
       }
    }

   //fill tree for each event
   tree->Fill();
}

void RecoSimDumperSimClusters::beginJob()
{

}

void RecoSimDumperSimClusters::endJob()
{


}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
float RecoSimDumperSimClusters::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoSimDumperSimClusters);

