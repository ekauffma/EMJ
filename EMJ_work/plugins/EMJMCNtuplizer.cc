#include <memory>
#include <iostream>
#include <string>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/CaloRecHitAuxSetter.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TTree.h"

typedef std::unordered_set<unsigned> PidSet;


class EMJMCNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit EMJMCNtuplizer(const edm::ParameterSet&);
  ~EMJMCNtuplizer() = default;

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};

  const edm::EDGetTokenT<std::vector<reco::PFCandidate>> tok_PFCandidate_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> tok_GenParticle_;
  const edm::EDGetTokenT<std::vector<reco::GenJet>> tok_GenJet_;
  const edm::EDGetTokenT<reco::PFJetCollection> tok_PFJet_;
  const edm::EDGetTokenT<GenEventInfoProduct> tok_genInfo_;
  const edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> tok_hcalRecHitsHBHE_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> tok_CaloGeometry_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> tok_Vertex_;

  const CaloSubdetectorGeometry *caloGeometry_HB;

  TTree* tree;
  edm::Service<TFileService> theFileService;
  
  int run;
  int lumi;
  int evt;

  float pthat;
  float weight;

  int nPV;

  int nPFCand;
  std::vector<double> rawHoverE;
  std::vector<double> HoverE;
  std::vector<double> hcalDepthEF1;
  std::vector<double> hcalDepthEF2;
  std::vector<double> hcalDepthEF3;
  std::vector<double> hcalDepthEF4;
  std::vector<double> hcalDepthEF5;
  std::vector<double> hcalDepthEF6;
  std::vector<double> hcalDepthEF7;
  
  int nJet;
  std::vector<double> Jet_Pt;
  std::vector<double> Jet_Eta;
  std::vector<double> Jet_Phi;
  std::vector<double> Jet_E;

  std::vector<int> Jet_nConstituent;
  std::vector<double> Constituent_Px;
  std::vector<double> Constituent_Py;
  std::vector<double> Constituent_Pz;
  std::vector<double> Constituent_E;
  std::vector<double> Constituent_rawHoverE;
  std::vector<double> Constituent_HoverE;
  std::vector<double> Constituent_hcalDepthEF1;
  std::vector<double> Constituent_hcalDepthEF2;
  std::vector<double> Constituent_hcalDepthEF3;
  std::vector<double> Constituent_hcalDepthEF4;
  std::vector<double> Constituent_hcalDepthEF5;
  std::vector<double> Constituent_hcalDepthEF6;
  std::vector<double> Constituent_hcalDepthEF7;
    
  std::vector<int> Jet_nJetRecHit;
  std::vector<int> JetRecHit_auxTDC3;
  std::vector<double> JetRecHit_energy;
  std::vector<double> JetRecHit_eta;
  std::vector<double> JetRecHit_phi;

  /*
  std::vector<int> Jet_nTrueJetRecHit;
  std::vector<int> TrueJetRecHit_auxTDC3;
  std::vector<double> TrueJetRecHit_energy;
  std::vector<double> TrueJetRecHit_eta;
  std::vector<double> TrueJetRecHit_phi;
  */  

  int nGenJet;
  std::vector<double> GenJet_Pt;
  std::vector<double> GenJet_Eta;
  std::vector<double> GenJet_Phi;
  std::vector<double> GenJet_E;
    
  std::vector<int> GenJet_nGenConstituent;
  std::vector<double> GenConstituent_Px;
  std::vector<double> GenConstituent_Py;
  std::vector<double> GenConstituent_Pz;
  std::vector<double> GenConstituent_E;
  std::vector<double> GenConstituent_vertex_x;
  std::vector<double> GenConstituent_vertex_y;
  std::vector<double> GenConstituent_vertex_z;
  std::vector<int> GenConstituent_status;
  std::vector<int> GenConstituent_pdgid;
  std::vector<int> GenConstituent_charge;
  std::vector<int> GenConstituent_idx;
    
  std::vector<int> GenJet_nDarkAncestor;
  std::vector<double> DarkAncestor_Px;
  std::vector<double> DarkAncestor_Py;
  std::vector<double> DarkAncestor_Pz;
  std::vector<double> DarkAncestor_E;
  std::vector<double> DarkAncestor_vertex_x;
  std::vector<double> DarkAncestor_vertex_y;
  std::vector<double> DarkAncestor_vertex_z;
  std::vector<int> DarkAncestor_status;
  std::vector<int> DarkAncestor_pdgid;
  std::vector<int> DarkAncestor_charge;
  std::vector<int> DarkAncestor_idx;
    
  int nGenParticle;
  std::vector<double> GenParticle_Px;
  std::vector<double> GenParticle_Py;
  std::vector<double> GenParticle_Pz;
  std::vector<double> GenParticle_E;
  std::vector<double> GenParticle_vertex_x;
  std::vector<double> GenParticle_vertex_y;
  std::vector<double> GenParticle_vertex_z;
  std::vector<int> GenParticle_status;
  std::vector<int> GenParticle_pdgid;
  std::vector<int> GenParticle_charge;
  std::vector<int> GenParticle_idx;
  std::vector<int> GenParticle_nParents;
  std::vector<int> GenParticle_parentIdx;
  std::vector<int> GenParticle_nDaughters;
  std::vector<int> GenParticle_daughterIdx;

  int nRecHit;
  std::vector<double> RecHit_time;
  std::vector<double> RecHit_E;
  std::vector<double> RecHit_eta;
  std::vector<double> RecHit_phi;
  std::vector<int> RecHit_auxTDC0;
  std::vector<int> RecHit_auxTDC1;
  std::vector<int> RecHit_auxTDC2;  
  std::vector<int> RecHit_auxTDC3;
  std::vector<int> RecHit_auxTDC4;
  std::vector<int> RecHit_auxTDC5;

};

EMJMCNtuplizer::EMJMCNtuplizer(const edm::ParameterSet& iConfig):
  tok_PFCandidate_( consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidateToken"))),
  tok_GenParticle_( consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("GenParticleToken"))),
  tok_GenJet_( consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("GenJetToken"))),
  tok_PFJet_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PFJetToken"))),
  tok_genInfo_( consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  tok_hcalRecHitsHBHE_(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>( iConfig.getParameter<edm::InputTag>("hbRecHitsToken") )),
  tok_CaloGeometry_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
  tok_Vertex_( consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexToken")))
{
  usesResource("TFileService");

  tree = theFileService->make<TTree>("Events","");

  tree->Branch("run", &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("evt", &evt);

  tree->Branch("pthat", &pthat);
  tree->Branch("weight", &weight);

  tree->Branch("nPV", &nPV);

  tree->Branch("nPFCand", &nPFCand);
  tree->Branch("PFCand_rawHoverE", &rawHoverE);
  tree->Branch("PFCand_HoverE", &HoverE);
  tree->Branch("PFCand_hcalDepthEF1", &hcalDepthEF1);
  tree->Branch("PFCand_hcalDepthEF2", &hcalDepthEF2);
  tree->Branch("PFCand_hcalDepthEF3", &hcalDepthEF3);
  tree->Branch("PFCand_hcalDepthEF4", &hcalDepthEF4);
  tree->Branch("PFCand_hcalDepthEF5", &hcalDepthEF5);
  tree->Branch("PFCand_hcalDepthEF6", &hcalDepthEF6);
  tree->Branch("PFCand_hcalDepthEF7", &hcalDepthEF7);

  tree->Branch("nJet", &nJet);
  tree->Branch("Jet_Pt", &Jet_Pt);
  tree->Branch("Jet_Eta", &Jet_Eta);
  tree->Branch("Jet_Phi", &Jet_Phi);
  tree->Branch("Jet_E", &Jet_E);

  tree->Branch("Jet_nConstituent", &Jet_nConstituent);
  tree->Branch("Constituent_Px", &Constituent_Px);
  tree->Branch("Constituent_Py", &Constituent_Py);
  tree->Branch("Constituent_Pz", &Constituent_Px);
  tree->Branch("Constituent_E", &Constituent_E);
  tree->Branch("Constituent_rawHoverE", &Constituent_rawHoverE);
  tree->Branch("Constituent_HoverE", &Constituent_HoverE);
  tree->Branch("Constituent_hcalDepthEF1", &Constituent_hcalDepthEF1);
  tree->Branch("Constituent_hcalDepthEF2", &Constituent_hcalDepthEF2);
  tree->Branch("Constituent_hcalDepthEF3", &Constituent_hcalDepthEF3);
  tree->Branch("Constituent_hcalDepthEF4", &Constituent_hcalDepthEF4);
  tree->Branch("Constituent_hcalDepthEF5", &Constituent_hcalDepthEF5);
  tree->Branch("Constituent_hcalDepthEF6", &Constituent_hcalDepthEF6);
  tree->Branch("Constituent_hcalDepthEF7", &Constituent_hcalDepthEF7);
    
  tree->Branch("Jet_nJetRecHit", &Jet_nJetRecHit);
  tree->Branch("JetRecHit_auxTDC3", &JetRecHit_auxTDC3);
  tree->Branch("JetRecHit_energy", &JetRecHit_energy);
  tree->Branch("JetRecHit_eta", &JetRecHit_eta);
  tree->Branch("JetRecHit_phi", &JetRecHit_phi);

  /*
  tree->Branch("Jet_nTrueJetRecHit", &Jet_nTrueJetRecHit);
  tree->Branch("JetRecHit_auxTDC3", &JetRecHit_auxTDC3);
  tree->Branch("JetRecHit_energy", &JetRecHit_energy);
  tree->Branch("JetRecHit_eta", &JetRecHit_eta);
  tree->Branch("JetRecHit_phi", &JetRecHit_phi);
  */

  tree->Branch("nGenJet", &nGenJet);
  tree->Branch("GenJet_Pt", &GenJet_Pt);
  tree->Branch("GenJet_Eta", &GenJet_Eta);
  tree->Branch("GenJet_Phi", &GenJet_Phi);
  tree->Branch("GenJet_E", &GenJet_E);
    
  tree->Branch("GenJet_nGenConstituent", &GenJet_nGenConstituent);
  tree->Branch("GenConstituent_Px", &GenConstituent_Px);
  tree->Branch("GenConstituent_Py", &GenConstituent_Py);
  tree->Branch("GenConstituent_Pz", &GenConstituent_Pz);
  tree->Branch("GenConstituent_E", &GenConstituent_E);
  tree->Branch("GenConstituent_vertex_x", &GenConstituent_vertex_x);
  tree->Branch("GenConstituent_vertex_y", &GenConstituent_vertex_y);
  tree->Branch("GenConstituent_vertex_z", &GenConstituent_vertex_z);
  tree->Branch("GenConstituent_status", &GenConstituent_status);
  tree->Branch("GenConstituent_pdgid", &GenConstituent_pdgid);
  tree->Branch("GenConstituent_charge", &GenConstituent_charge);
  tree->Branch("GenConstituent_idx", &GenConstituent_idx);
    
  tree->Branch("GenJet_nDarkAncestor", &GenJet_nDarkAncestor);
  tree->Branch("DarkAncestor_Px", &DarkAncestor_Px);
  tree->Branch("DarkAncestor_Py", &DarkAncestor_Py);
  tree->Branch("DarkAncestor_Pz", &DarkAncestor_Pz);
  tree->Branch("DarkAncestor_E", &DarkAncestor_E);
  tree->Branch("DarkAncestor_vertex_x", &DarkAncestor_vertex_x);
  tree->Branch("DarkAncestor_vertex_y", &DarkAncestor_vertex_y);
  tree->Branch("DarkAncestor_vertex_z", &DarkAncestor_vertex_z);
  tree->Branch("DarkAncestor_status", &DarkAncestor_status);
  tree->Branch("DarkAncestor_pdgid", &DarkAncestor_pdgid);
  tree->Branch("DarkAncestor_charge", &DarkAncestor_charge);
  tree->Branch("DarkAncestor_idx", &DarkAncestor_idx);

  tree->Branch("nGenParticle", &nGenParticle);
  tree->Branch("GenParticle_Px", &GenParticle_Px);
  tree->Branch("GenParticle_Py", &GenParticle_Py);
  tree->Branch("GenParticle_Pz", &GenParticle_Pz);
  tree->Branch("GenParticle_E", &GenParticle_E);
  tree->Branch("GenParticle_vertex_x", &GenParticle_vertex_x);
  tree->Branch("GenParticle_vertex_y", &GenParticle_vertex_y);
  tree->Branch("GenParticle_vertex_z", &GenParticle_vertex_z);
  tree->Branch("GenParticle_status", &GenParticle_status);
  tree->Branch("GenParticle_pdgid", &GenParticle_pdgid);
  tree->Branch("GenParticle_charge", &GenParticle_charge);
  tree->Branch("GenParticle_idx", &GenParticle_idx);
  tree->Branch("GenParticle_nParents", &GenParticle_nParents);
  tree->Branch("GenParticle_parentIdx", &GenParticle_parentIdx); 
  tree->Branch("GenParticle_nDaughters", &GenParticle_nDaughters);
  tree->Branch("GenParticle_daughterIdx", &GenParticle_daughterIdx); 

  tree->Branch("nRecHit", &nRecHit);
  tree->Branch("RecHit_time", &RecHit_time);
  tree->Branch("RecHit_E", &RecHit_E);
  tree->Branch("RecHit_eta", &RecHit_eta);
  tree->Branch("RecHit_phi", &RecHit_phi);
  tree->Branch("RecHit_auxTDC0", &RecHit_auxTDC0);
  tree->Branch("RecHit_auxTDC1", &RecHit_auxTDC1);
  tree->Branch("RecHit_auxTDC2", &RecHit_auxTDC2);
  tree->Branch("RecHit_auxTDC3", &RecHit_auxTDC3);
  tree->Branch("RecHit_auxTDC4", &RecHit_auxTDC4);
  tree->Branch("RecHit_auxTDC5", &RecHit_auxTDC5);
}

bool isDarkParticle(int pdgId) {
  return std::abs(pdgId)>490000;
}

template <typename T>
const T* getDarkAncestor(const T* part ) {
  if (!part) return nullptr;
  int pdgId = part->pdgId();
  if(isDarkParticle(pdgId)) return part;
  for(size_t i=0; i<part->numberOfMothers(); i++) {
    const reco::GenParticle* ancestor = getDarkAncestor(static_cast<const T*>(part->mother(i)));
      if (ancestor) return ancestor;
  }
  return nullptr;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = deltaPhi(phi1, phi2);
    return std::sqrt(deta * deta + dphi * dphi);
}

// Given the list of Gen particle pointers, find the where a given point is in the list
template <typename T_base, typename T_derived>
static unsigned get_index(const std::vector<const T_derived*>& list, const T_base* target, unsigned def) {
  auto it = std::find(list.begin(), list.end(), target);
  unsigned ret = it - list.begin();
   return (ret == list.size()) ? def : ret;
}

void EMJMCNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<reco::PFCandidate>> handle_PFCandidate;
  iEvent.getByToken(tok_PFCandidate_, handle_PFCandidate);
    
  edm::Handle<std::vector<reco::GenParticle>> handle_GenParticle;
  iEvent.getByToken(tok_GenParticle_, handle_GenParticle);

  pthat = 0.0;
  edm::Handle<GenEventInfoProduct> handle_genInfo;
  iEvent.getByToken(tok_genInfo_, handle_genInfo);
  if (handle_genInfo->hasBinningValues()) pthat = (float)handle_genInfo->binningValues()[0];
  weight = (float)handle_genInfo->weight();

  const edm::Handle<reco::PFJetCollection> handle_PFJet = iEvent.getHandle(tok_PFJet_);
    
  edm::Handle<std::vector<reco::GenJet>> handle_GenJet;
  iEvent.getByToken(tok_GenJet_, handle_GenJet);

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> handle_hcalRecHitsHBHE;
  iEvent.getByToken(tok_hcalRecHitsHBHE_, handle_hcalRecHitsHBHE);

  auto const& geoHandle = iSetup.getData(tok_CaloGeometry_);
  caloGeometry_HB = geoHandle.getSubdetectorGeometry(DetId::Hcal, HcalBarrel);

  edm::Handle<std::vector<reco::Vertex>> handle_Vertex;
  iEvent.getByToken(tok_Vertex_, handle_Vertex);

  nPV = handle_Vertex->size();

  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  evt = iEvent.id().event();

  nJet = handle_PFJet->size();
  for (
    reco::PFJetCollection::const_iterator it_jet = handle_PFJet->begin(); 
    it_jet != handle_PFJet->end(); 
    ++it_jet
  ) {

    const reco::PFJet* jet = &(*it_jet);
    Jet_Pt.push_back(jet->pt());
    Jet_Eta.push_back(jet->eta());
    Jet_Phi.push_back(jet->phi());
    Jet_E.push_back(jet->energy());

    // get constituent pfcandidates
    std::vector<reco::PFCandidatePtr> jet_const = jet->getPFConstituents();
    Jet_nConstituent.push_back(jet_const.size());
    // int Jet_nTrueJetRecHit_current = 0;
    for (
      std::vector<reco::PFCandidatePtr>::const_iterator it_cand = jet_const.begin(); 
      it_cand != jet_const.end(); 
      ++it_cand
    ) {
      Constituent_Px.push_back((*it_cand)->px());
      Constituent_Py.push_back((*it_cand)->py());
      Constituent_Pz.push_back((*it_cand)->pz());
      Constituent_E.push_back((*it_cand)->energy());
      Constituent_rawHoverE.push_back((*it_cand)->rawHcalEnergy()/(*it_cand)->rawEcalEnergy());
      Constituent_HoverE.push_back((*it_cand)->hcalEnergy()/(*it_cand)->ecalEnergy());
      Constituent_hcalDepthEF1.push_back((*it_cand)->hcalDepthEnergyFraction(1));
      Constituent_hcalDepthEF2.push_back((*it_cand)->hcalDepthEnergyFraction(2));
      Constituent_hcalDepthEF3.push_back((*it_cand)->hcalDepthEnergyFraction(3));
      Constituent_hcalDepthEF4.push_back((*it_cand)->hcalDepthEnergyFraction(4));
      Constituent_hcalDepthEF5.push_back((*it_cand)->hcalDepthEnergyFraction(5));
      Constituent_hcalDepthEF6.push_back((*it_cand)->hcalDepthEnergyFraction(6));
      Constituent_hcalDepthEF7.push_back((*it_cand)->hcalDepthEnergyFraction(7));

      /*
      for (const reco::PFCandidate::ElementInBlock &eib : (*it_cand)->elementsInBlocks()) {
        const reco::PFBlockRef &block = eib.first;
        for (const reco::PFBlockElement &element : block->elements()) {
          if (element.type() != reco::PFBlockElement::Type::HCAL) {
            continue;
          }
	  const reco::PFClusterRef &cluster_ref = element.clusterRef();
      	  if (!cluster_ref.isAvailable()) {
            continue;
          }
          for (const std::pair<DetId, float> &p : cluster_ref->hitsAndFractions()) {
            const det_Id = p.first;
            const fraction = p.second;
            std::cout << "Out:" << det_Id.rawId() << " " << fraction << std::endl;
          }
        }
      }
      */
    }
      
    int Jet_nJetRecHit_current = 0;
    for( uint ih = 0; ih < handle_hcalRecHitsHBHE->size(); ih++ ){
      const HBHERecHit *recHit = &(*handle_hcalRecHitsHBHE)[ih];
      const HcalDetId recHitId = recHit->detid();
      if( recHit->detid().subdetId() != HcalBarrel ) continue; // ensure recHit is in the HCAL barrel
      const auto recHitPos = caloGeometry_HB->getGeometry(recHitId)->getPosition();
      if (deltaR(recHitPos.eta(), recHitPos.phi(), jet->eta(), jet->phi())>0.4) continue; // ensure recHit is within jet radius
      const uint32_t auxTDC = recHit->auxTDC();
      if (auxTDC) {
        const unsigned six_bits_mask = 0x3f;
        ++Jet_nJetRecHit_current;
        JetRecHit_auxTDC3.push_back(CaloRecHitAuxSetter::getField(auxTDC, six_bits_mask, 3 * 6));
        JetRecHit_energy.push_back(recHit->energy());
        JetRecHit_eta.push_back(recHitPos.eta());
        JetRecHit_phi.push_back(recHitPos.phi());
        JetRecHit_energy.push_back(recHit->energy());
      }
    }
    Jet_nJetRecHit.push_back(Jet_nJetRecHit_current);

  }

  nPFCand = handle_PFCandidate->size();
  for(
    auto theObject = handle_PFCandidate->begin();
    theObject != handle_PFCandidate->end();
    ++theObject
  ) {
    rawHoverE.push_back(theObject->rawHcalEnergy()/theObject->rawEcalEnergy());
    HoverE.push_back(theObject->hcalEnergy()/theObject->ecalEnergy());
    hcalDepthEF1.push_back(theObject->hcalDepthEnergyFraction(1));
    hcalDepthEF2.push_back(theObject->hcalDepthEnergyFraction(2));
    hcalDepthEF3.push_back(theObject->hcalDepthEnergyFraction(3));
    hcalDepthEF4.push_back(theObject->hcalDepthEnergyFraction(4));
    hcalDepthEF5.push_back(theObject->hcalDepthEnergyFraction(5));
    hcalDepthEF6.push_back(theObject->hcalDepthEnergyFraction(6));
    hcalDepthEF7.push_back(theObject->hcalDepthEnergyFraction(7));
  }
    
  // create unfiltered list of pointers for index
  std::vector<const reco::GenParticle *> gen_list;
  for (
    std::vector<reco::GenParticle>::const_iterator it_gpart = handle_GenParticle->begin();
    it_gpart != handle_GenParticle->end();
    ++it_gpart
  ) {
    const reco::GenParticle* gpart = &(*it_gpart);
    gen_list.push_back(gpart);
  }

  nGenParticle = handle_GenParticle->size();
  for (
    std::vector<reco::GenParticle>::const_iterator it_gpart = handle_GenParticle->begin();
    it_gpart != handle_GenParticle->end();
    ++it_gpart
  ) {
    const reco::GenParticle* gpart = &(*it_gpart);
    GenParticle_Px.push_back(gpart->px());
    GenParticle_Py.push_back(gpart->py());
    GenParticle_Pz.push_back(gpart->pz());
    GenParticle_E.push_back(gpart->energy());
    GenParticle_vertex_x.push_back(gpart->vertex().x());
    GenParticle_vertex_y.push_back(gpart->vertex().y());
    GenParticle_vertex_z.push_back(gpart->vertex().z());
    GenParticle_status.push_back(gpart->status());
    GenParticle_pdgid.push_back(gpart->pdgId());
    GenParticle_charge.push_back(gpart->charge());
    size_t self_idx = get_index(gen_list, gpart, -1);
    GenParticle_idx.push_back(self_idx);
    GenParticle_nParents.push_back(gpart->numberOfMothers());
    for(size_t i = 0; i < gpart->numberOfMothers(); i++) {
      GenParticle_parentIdx.push_back(get_index(gen_list, gpart->mother(i), self_idx)); 
    }
    GenParticle_nDaughters.push_back(gpart->numberOfDaughters());
    for(size_t i = 0; i < gpart->numberOfDaughters(); i++) {
      GenParticle_daughterIdx.push_back(get_index(gen_list, gpart->daughter(i), self_idx));
    }
  }
    
  nGenJet = handle_GenJet->size();
  for (
    std::vector<reco::GenJet>::const_iterator it_gjet = handle_GenJet->begin();
    it_gjet != handle_GenJet->end();
    ++it_gjet
  ) {
    int nDarkAncestors = 0;

    const reco::GenJet* gjet = &(*it_gjet);
    GenJet_Pt.push_back(gjet->pt());
    GenJet_Eta.push_back(gjet->eta());
    GenJet_Phi.push_back(gjet->phi());
    GenJet_E.push_back(gjet->energy());
      
    std::vector<const reco::GenParticle*> genjet_const = gjet->getGenConstituents();
    GenJet_nGenConstituent.push_back(genjet_const.size());
    for (
      std::vector<const reco::GenParticle*>::const_iterator it_gconst = genjet_const.begin();
      it_gconst != genjet_const.end();
      ++it_gconst
    ) {
      GenConstituent_Px.push_back((*it_gconst)->px());
      GenConstituent_Py.push_back((*it_gconst)->py());
      GenConstituent_Pz.push_back((*it_gconst)->pz());
      GenConstituent_E.push_back((*it_gconst)->energy());
      GenConstituent_vertex_x.push_back((*it_gconst)->vertex().x());
      GenConstituent_vertex_y.push_back((*it_gconst)->vertex().y());
      GenConstituent_vertex_z.push_back((*it_gconst)->vertex().z());
      GenConstituent_status.push_back((*it_gconst)->status());
      GenConstituent_pdgid.push_back((*it_gconst)->pdgId());
      GenConstituent_charge.push_back((*it_gconst)->charge());
      GenConstituent_idx.push_back(get_index(gen_list, *it_gconst, -1));

      const reco::GenParticle* ancestor = getDarkAncestor(*it_gconst );
      if(ancestor) {
      	std::cout<<"ANCESTOR PDGID = "<<ancestor->pdgId()<<std::endl;
	nDarkAncestors++;
	DarkAncestor_Px.push_back(ancestor->px());
	DarkAncestor_Py.push_back(ancestor->py());
        DarkAncestor_Pz.push_back(ancestor->pz());
	DarkAncestor_E.push_back(ancestor->energy());
	DarkAncestor_vertex_x.push_back(ancestor->vertex().x());
        DarkAncestor_vertex_y.push_back(ancestor->vertex().y());
        DarkAncestor_vertex_z.push_back(ancestor->vertex().z());
	DarkAncestor_status.push_back(ancestor->status());
	DarkAncestor_pdgid.push_back(ancestor->pdgId());
	DarkAncestor_charge.push_back(ancestor->charge());
	DarkAncestor_idx.push_back(get_index(gen_list, ancestor, -1));
      }
      else {
      	std::cout<<"NO DARK ANCESTOR FOUND"<<std::endl;
      }
    }
    GenJet_nDarkAncestor.push_back(nDarkAncestors);
  }

  nRecHit = 0;
  for( uint ih = 0; ih < handle_hcalRecHitsHBHE->size(); ih++ ){
    const HBHERecHit *recHit = &(*handle_hcalRecHitsHBHE)[ih];
    const HcalDetId recHitId = recHit->detid();
    if( recHit->detid().subdetId() != HcalBarrel ) continue;
    nRecHit++;
    RecHit_time.push_back(recHit->time());
    RecHit_E.push_back(recHit->energy());
    const auto recHitPos = caloGeometry_HB->getGeometry(recHitId)->getPosition();
    RecHit_eta.push_back(recHitPos.eta());
    RecHit_phi.push_back(recHitPos.phi());
    const uint32_t auxTDC = recHit->auxTDC();
    if (auxTDC) {
      const unsigned six_bits_mask = 0x3f;
      RecHit_auxTDC0.push_back(CaloRecHitAuxSetter::getField(auxTDC, six_bits_mask, 0 * 6));
      RecHit_auxTDC1.push_back(CaloRecHitAuxSetter::getField(auxTDC, six_bits_mask, 1 * 6));
      RecHit_auxTDC2.push_back(CaloRecHitAuxSetter::getField(auxTDC, six_bits_mask, 2 * 6));
      RecHit_auxTDC3.push_back(CaloRecHitAuxSetter::getField(auxTDC, six_bits_mask, 3 * 6));
      RecHit_auxTDC4.push_back(CaloRecHitAuxSetter::getField(auxTDC, six_bits_mask, 4 * 6));
      RecHit_auxTDC5.push_back(CaloRecHitAuxSetter::getField(auxTDC, six_bits_mask, 5 * 6));
    }
  }

  tree->Fill();

  rawHoverE.clear();
  HoverE.clear();
  hcalDepthEF1.clear();
  hcalDepthEF2.clear();
  hcalDepthEF3.clear();
  hcalDepthEF4.clear();
  hcalDepthEF5.clear();
  hcalDepthEF6.clear();
  hcalDepthEF7.clear();

  Jet_Pt.clear();
  Jet_Eta.clear();
  Jet_Phi.clear();
  Jet_E.clear();
  Jet_nConstituent.clear();
  Constituent_Px.clear();
  Constituent_Py.clear();
  Constituent_Pz.clear();
  Constituent_E.clear();
  Constituent_rawHoverE.clear();
  Constituent_HoverE.clear();
  Constituent_hcalDepthEF1.clear();
  Constituent_hcalDepthEF2.clear();
  Constituent_hcalDepthEF3.clear();
  Constituent_hcalDepthEF4.clear();
  Constituent_hcalDepthEF5.clear();
  Constituent_hcalDepthEF6.clear();
  Constituent_hcalDepthEF7.clear();
    
  Jet_nJetRecHit.clear();
  JetRecHit_auxTDC3.clear();
  JetRecHit_energy.clear();
  JetRecHit_eta.clear();
  JetRecHit_phi.clear();
    
  GenJet_Pt.clear();
  GenJet_Eta.clear();
  GenJet_Phi.clear();
  GenJet_E.clear();
                                   
  GenJet_nGenConstituent.clear();
  GenConstituent_Px.clear();
  GenConstituent_Py.clear();
  GenConstituent_Pz.clear();
  GenConstituent_E.clear();
  GenConstituent_vertex_x.clear();
  GenConstituent_vertex_y.clear();
  GenConstituent_vertex_z.clear();
  GenConstituent_status.clear();
  GenConstituent_pdgid.clear();
  GenConstituent_charge.clear();
  GenConstituent_idx.clear();
                                   
  GenJet_nDarkAncestor.clear();
  DarkAncestor_Px.clear();
  DarkAncestor_Py.clear();
  DarkAncestor_Pz.clear();
  DarkAncestor_E.clear();
  DarkAncestor_vertex_x.clear();
  DarkAncestor_vertex_y.clear();
  DarkAncestor_vertex_z.clear();
  DarkAncestor_status.clear();
  DarkAncestor_pdgid.clear();
  DarkAncestor_charge.clear();
  DarkAncestor_idx.clear();
                                   
  GenParticle_Px.clear();
  GenParticle_Py.clear();
  GenParticle_Pz.clear();
  GenParticle_E.clear();
  GenParticle_vertex_x.clear();
  GenParticle_vertex_y.clear();
  GenParticle_vertex_z.clear();
  GenParticle_status.clear();
  GenParticle_pdgid.clear();
  GenParticle_charge.clear();
  GenParticle_idx.clear();
  GenParticle_nParents.clear();
  GenParticle_parentIdx.clear(); 
  GenParticle_nDaughters.clear();
  GenParticle_daughterIdx.clear();  

  RecHit_time.clear();
  RecHit_E.clear();
  RecHit_eta.clear();
  RecHit_phi.clear();
  RecHit_auxTDC0.clear();
  RecHit_auxTDC1.clear();
  RecHit_auxTDC2.clear();
  RecHit_auxTDC3.clear();
  RecHit_auxTDC4.clear();
  RecHit_auxTDC5.clear();
}

DEFINE_FWK_MODULE(EMJMCNtuplizer);
