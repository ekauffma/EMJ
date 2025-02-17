#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// #include "DataFormats/PatCandidates/interface/Jet.h"
// #include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

#include <string>


class JetInfoNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit JetInfoNtuplizer(const edm::ParameterSet&);
  ~JetInfoNtuplizer() = default;

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};

  edm::EDGetTokenT<std::vector<reco::PFCandidate>> tok_PFCandidate_;
  const edm::EDGetTokenT<reco::PFJetCollection> tok_PFJet_;
  const edm::EDGetTokenT<GenEventInfoProduct> tok_genInfo_;

  TTree* tree;
  edm::Service<TFileService> theFileService;
  
  int run;
  int lumi;
  int evt;

  float pthat;
  float weight;

  int nPFCand;
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

  std::vector<double> Jet_nConstituent;
  std::vector<double> Constituent_Px;
  std::vector<double> Constituent_Py;
  std::vector<double> Constituent_Pz;
  std::vector<double> Constituent_E;
  std::vector<double> Constituent_hcalDepthEF1;
  std::vector<double> Constituent_hcalDepthEF2;
  std::vector<double> Constituent_hcalDepthEF3;
  std::vector<double> Constituent_hcalDepthEF4;
  std::vector<double> Constituent_hcalDepthEF5;
  std::vector<double> Constituent_hcalDepthEF6;
  std::vector<double> Constituent_hcalDepthEF7;
};

JetInfoNtuplizer::JetInfoNtuplizer(const edm::ParameterSet& iConfig):
  tok_PFCandidate_( consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidateToken"))),
  tok_PFJet_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PFJetToken"))),
  tok_genInfo_( consumes<GenEventInfoProduct>(edm::InputTag("generator")))
{
  usesResource("TFileService");

  tree = theFileService->make<TTree>("PFCand","Particle Flow Candidate information");

  tree->Branch("run", &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("evt", &evt);

  tree->Branch("pthat", &pthat);
  tree->Branch("weight", &weight);

  tree->Branch("nPFCand", &nPFCand);
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
  tree->Branch("Constituent_hcalDepthEF1", &Constituent_hcalDepthEF1);
  tree->Branch("Constituent_hcalDepthEF2", &Constituent_hcalDepthEF2);
  tree->Branch("Constituent_hcalDepthEF3", &Constituent_hcalDepthEF3);
  tree->Branch("Constituent_hcalDepthEF4", &Constituent_hcalDepthEF4);
  tree->Branch("Constituent_hcalDepthEF5", &Constituent_hcalDepthEF5);
  tree->Branch("Constituent_hcalDepthEF6", &Constituent_hcalDepthEF6);
  tree->Branch("Constituent_hcalDepthEF7", &Constituent_hcalDepthEF7);
}

void JetInfoNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<reco::PFCandidate>> handle_PFCandidate;
  iEvent.getByToken(tok_PFCandidate_, handle_PFCandidate);

  pthat = 0.0;
  edm::Handle<GenEventInfoProduct> handle_genInfo;
  iEvent.getByToken(tok_genInfo_, handle_genInfo);
  if (handle_genInfo->hasBinningValues()) pthat = (float)handle_genInfo->binningValues()[0];
  weight = (float)handle_genInfo->weight();

  const edm::Handle<reco::PFJetCollection> handle_PFJet = iEvent.getHandle(tok_PFJet_);

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
    for (
      std::vector<reco::PFCandidatePtr>::const_iterator it_cand = jet_const.begin(); 
      it_cand != jet_const.end(); 
      ++it_cand
    ) {
      Constituent_Px.push_back((*it_cand)->px());
      Constituent_Py.push_back((*it_cand)->py());
      Constituent_Pz.push_back((*it_cand)->pz());
      Constituent_E.push_back((*it_cand)->energy());
      Constituent_hcalDepthEF1.push_back((*it_cand)->hcalDepthEnergyFraction(1));
      Constituent_hcalDepthEF2.push_back((*it_cand)->hcalDepthEnergyFraction(2));
      Constituent_hcalDepthEF3.push_back((*it_cand)->hcalDepthEnergyFraction(3));
      Constituent_hcalDepthEF4.push_back((*it_cand)->hcalDepthEnergyFraction(4));
      Constituent_hcalDepthEF5.push_back((*it_cand)->hcalDepthEnergyFraction(5));
      Constituent_hcalDepthEF6.push_back((*it_cand)->hcalDepthEnergyFraction(6));
      Constituent_hcalDepthEF7.push_back((*it_cand)->hcalDepthEnergyFraction(7));
    }
  }

  nPFCand = handle_PFCandidate->size();
  for(
    auto theObject = handle_PFCandidate->begin();
    theObject != handle_PFCandidate->end();
    ++theObject
  ) {
    hcalDepthEF1.push_back(theObject->hcalDepthEnergyFraction(1));
    hcalDepthEF2.push_back(theObject->hcalDepthEnergyFraction(2));
    hcalDepthEF3.push_back(theObject->hcalDepthEnergyFraction(3));
    hcalDepthEF4.push_back(theObject->hcalDepthEnergyFraction(4));
    hcalDepthEF5.push_back(theObject->hcalDepthEnergyFraction(5));
    hcalDepthEF6.push_back(theObject->hcalDepthEnergyFraction(6));
    hcalDepthEF7.push_back(theObject->hcalDepthEnergyFraction(7));
  }

  tree->Fill();

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
  Constituent_hcalDepthEF1.clear();
  Constituent_hcalDepthEF2.clear();
  Constituent_hcalDepthEF3.clear();
  Constituent_hcalDepthEF4.clear();
  Constituent_hcalDepthEF5.clear();
  Constituent_hcalDepthEF6.clear();
  Constituent_hcalDepthEF7.clear();
}

DEFINE_FWK_MODULE(JetInfoNtuplizer);
