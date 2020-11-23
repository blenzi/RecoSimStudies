#ifndef RecoSimStudies_Dumpers_RecoSimDumperSimClusters_H
#define RecoSimStudies_Dumpers_RecoSimDumperSimClusters_H

// system include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/FWLite/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TTree.h"

#include <vector>


class RecoSimDumperSimClusters : public edm::EDAnalyzer
{
      public:
         explicit RecoSimDumperSimClusters(const edm::ParameterSet&);
	 ~RecoSimDumperSimClusters();


      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
         virtual void endJob() ;

      // ----------additional functions-------------------
      float reduceFloat(float val, int bits);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesCaloPart(const CaloParticle* iCaloParticle, float simHitEnergy_cut);
      GlobalPoint calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES);

      // ----------collection tokens-------------------
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_;
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_;

      edm::Service<TFileService> iFile;

      // ----------config inputs-------------------
      bool doCompression_;
      int nBits_;
      bool saveGenParticles_;
      bool saveCaloParticles_;
      bool saveSimhits_;
      bool saveOnlyCaloSimhits_;
      bool saveRechits_;
      std::vector<int> genID_;

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
      std::vector<int> caloParticle_pdgId;
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

      std::vector<int> caloParticle_nSimClusters;
      std::vector<int> simCluster_pdgId;
      std::vector<float> simCluster_energy;
      std::vector<float> simCluster_eta;
      std::vector<float> simCluster_phi;

      std::vector<uint32_t> simHit_id;
      std::vector<float> simHit_energy;
      std::vector<float> simHit_eta;
      std::vector<float> simHit_phi;
      std::vector<int> simHit_ieta;
      std::vector<int> simHit_iphi;
      std::vector<int> simHit_iz;

      std::vector<float> recHit_energy;
      std::vector<float> recHit_eta;
      std::vector<float> recHit_phi;
      std::vector<int> recHit_ieta;
      std::vector<int> recHit_iphi;
      std::vector<int> recHit_iz;

};

#endif
