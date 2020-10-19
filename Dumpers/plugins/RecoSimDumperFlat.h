#ifndef RecoSimStudies_Dumpers_RecoSimDumperFlat_H
#define RecoSimStudies_Dumpers_RecoSimDumperFlat_H

// system include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
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
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"


#include "DataFormats/Math/interface/deltaR.h"
//#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

#include "TTree.h"

#include <vector>

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

class RecoSimDumperFlat : public edm::EDAnalyzer
{
      public:
         explicit RecoSimDumperFlat(const edm::ParameterSet&);
	 ~RecoSimDumperFlat();


      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
         virtual void endJob() ;

      // ----------additional functions-------------------
      float reduceFloat(float val, int bits);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesCaloPart(const CaloParticle* iCaloParticle, float simHitEnergy_cut);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      GlobalPoint calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES);

      // ----------collection tokens-------------------
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_;
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_;

      edm::Service<TFileService> iFile;
      const CaloSubdetectorGeometry* _ebGeom;
      const CaloSubdetectorGeometry* _eeGeom;
      const CaloSubdetectorGeometry* _esGeom;
      bool _esPlus;
      bool _esMinus;

      // ----------config inputs-------------------
      bool doCompression_;
      int nBits_;
      bool saveGenParticles_;
      bool saveCaloParticles_;
      bool saveSimhits_;
      bool saveRechits_;
      std::vector<int> genID_;

      // ----------DNN inputs-------------------
      std::vector<double> HLF_VectorVar_;
      std::vector<std::vector<double>> PL_VectorVar_;
      std::vector<double> x_mean_, x_std_, list_mean_, list_std_;

      // ----------histograms & trees & branches-------------------
      TTree* tree;

      long int eventId;
      int lumiId;
      int runId;
      int nVtx;
      float rho;

      std::vector<int> genParticle_id;
      std::vector<float> genParticle_energy;
      std::vector<float> genParticle_pt;
      std::vector<float> genParticle_eta;
      std::vector<float> genParticle_phi;

      std::vector<int> caloParticle_id;
      std::vector<float> caloParticle_genEnergy;
      std::vector<float> caloParticle_simEnergy;
      std::vector<float> caloParticle_genPt;
      std::vector<float> caloParticle_simPt;
      std::vector<float> caloParticle_genEta;
      std::vector<float> caloParticle_simEta;
      std::vector<float> caloParticle_genPhi;
      std::vector<float> caloParticle_simPhi;
      std::vector<int> caloParticle_simIeta;
      std::vector<int> caloParticle_simIphi;
      std::vector<int> caloParticle_simIz;

      std::vector<float> simHit_energy;
      std::vector<float> simHit_eta;
      std::vector<float> simHit_phi;
      std::vector<int> simHit_ieta;
      std::vector<int> simHit_iphi;
      std::vector<int> simHit_iz;
      std::vector<int> simHit_particleIndex1;
      std::vector<int> simHit_particleIndex2;
      std::vector<int> simHit_particleIndex3;
      std::vector<float> simHit_energyParticle1;
      std::vector<float> simHit_energyParticle2;
      std::vector<float> simHit_energyParticle3;

      std::vector<float> recHit_energy;
      std::vector<float> recHit_eta;
      std::vector<float> recHit_phi;
      std::vector<int> recHit_ieta;
      std::vector<int> recHit_iphi;
      std::vector<int> recHit_iz;

      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart;

};

#endif
