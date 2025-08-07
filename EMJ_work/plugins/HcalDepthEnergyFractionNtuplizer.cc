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

//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/HcalDepthEnergyFractions.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"


class HcalDepthEnergyFractionNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit HcalDepthEnergyFractionNtuplizer(const edm::ParameterSet&);
  ~HcalDepthEnergyFractionNtuplizer() = default;

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};

  const edm::EDGetTokenT<pat::PackedCandidateCollection> tok_PackedCandidate_;
  const edm::EDGetTokenT<edm::ValueMap<pat::HcalDepthEnergyFractions>> tok_HcalDepthEF_;

  TTree* tree;
  edm::Service<TFileService> theFileService;
  
  int run;
  int lumi;
  int evt;

  int nPFCand;
  std::vector<double> hcalDepthEF1;
  std::vector<double> hcalDepthEF2;
  std::vector<double> hcalDepthEF3;
  std::vector<double> hcalDepthEF4;
  std::vector<double> hcalDepthEF5;
  std::vector<double> hcalDepthEF6;
  std::vector<double> hcalDepthEF7;
};

HcalDepthEnergyFractionNtuplizer::HcalDepthEnergyFractionNtuplizer(const edm::ParameterSet& iConfig):
  tok_PackedCandidate_( consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("PackedCandidateToken"))),
  tok_HcalDepthEF_( consumes<edm::ValueMap<pat::HcalDepthEnergyFractions>>(iConfig.getParameter<edm::InputTag>("HcalDepthEFToken")))
{
  usesResource("TFileService");

  tree = theFileService->make<TTree>("Events","");

  tree->Branch("run", &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("evt", &evt);

  tree->Branch("nPFCand", &nPFCand);
  tree->Branch("PFCand_hcalDepthEF1", &hcalDepthEF1);
  tree->Branch("PFCand_hcalDepthEF2", &hcalDepthEF2);
  tree->Branch("PFCand_hcalDepthEF3", &hcalDepthEF3);
  tree->Branch("PFCand_hcalDepthEF4", &hcalDepthEF4);
  tree->Branch("PFCand_hcalDepthEF5", &hcalDepthEF5);
  tree->Branch("PFCand_hcalDepthEF6", &hcalDepthEF6);
  tree->Branch("PFCand_hcalDepthEF7", &hcalDepthEF7);
}


void HcalDepthEnergyFractionNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<pat::PackedCandidateCollection> handle_PackedCandidate;
  iEvent.getByToken(tok_PackedCandidate_, handle_PackedCandidate);

  edm::Handle<edm::ValueMap<pat::HcalDepthEnergyFractions>> handle_HcalDepthEF;
  iEvent.getByToken(tok_HcalDepthEF_, handle_HcalDepthEF);

  std::cout << "PackedCandidate size: " << handle_PackedCandidate->size() << std::endl;
  std::cout << "HcalDepthEF isValid: " << handle_HcalDepthEF.isValid() << std::endl;

  auto ref = edm::Ref<pat::PackedCandidateCollection>(handle_PackedCandidate, 0);
  std::cout << "Ref is available? " << ref.isAvailable() << ", isNull? " << ref.isNull() << std::endl;

  std::cout << "Map size: " << handle_HcalDepthEF->size() << std::endl;

  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  evt = iEvent.id().event();

  nPFCand = handle_PackedCandidate->size();
  for(
    unsigned int i = 0;
    i < handle_PackedCandidate->size();
    ++i
  ) {
    edm::Ref<pat::PackedCandidateCollection> pfCandRef(handle_PackedCandidate, i);
    const auto& hcalEF = (*handle_HcalDepthEF)[pfCandRef];

    hcalDepthEF1.push_back(hcalEF.fraction(0));
    hcalDepthEF2.push_back(hcalEF.fraction(1));
    hcalDepthEF3.push_back(hcalEF.fraction(2));
    hcalDepthEF4.push_back(hcalEF.fraction(3));
    hcalDepthEF5.push_back(hcalEF.fraction(4));
    hcalDepthEF6.push_back(hcalEF.fraction(5));
    hcalDepthEF7.push_back(hcalEF.fraction(6));
    }
    
  tree->Fill();

  hcalDepthEF1.clear();
  hcalDepthEF2.clear();
  hcalDepthEF3.clear();
  hcalDepthEF4.clear();
  hcalDepthEF5.clear();
  hcalDepthEF6.clear();
  hcalDepthEF7.clear();
}

DEFINE_FWK_MODULE(HcalDepthEnergyFractionNtuplizer);
