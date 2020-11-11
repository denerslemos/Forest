
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/HeavyIonEvent/interface/HFFilterInfo.h" //this line is needed to access the HF Filters

#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <HepMC/PdfInfo.h>

#include "TTree.h"

//
// class declaration
//

class GOForest : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit GOForest(const edm::ParameterSet&);
//  ~GOForest();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  //to save output

  edm::Service<TFileService> fs;

  //access the collections


  edm::EDGetTokenT<reco::Centrality> CentralityTag_;
  edm::EDGetTokenT<int> CentralityBinTag_;

  edm::EDGetTokenT<reco::HFFilterInfo> HFfilters_;

  edm::EDGetTokenT<std::vector<reco::Vertex>> VertexTag_;

//  edm::EDGetTokenT<edm::GenHIEvent> HiMCTag_;
//  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfoToken_;
//  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
//  edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;

  //make output TTrees

  TTree *Ev_TTree;
  TTree *Skim_TTree;
  TTree *GO_TTree;

  //useful booleans

  bool doCentrality_;
  bool doVertex_;
  bool doMC_;
  bool doHiMC_;
  bool useHepMC_;
 // bool doVertex_;

  //variables for filters

  unsigned short int pprimaryvertexfilter, clustercompatibilityfilter, beamscramping;
  unsigned short int numMinHFTower2, numMinHFTower3, numMinHFTower4, numMinHFTower5;

  //centrality variables

  int hiBin;
  int hiNpix, hiNpixelTracks, hiNtracks, hiNtracksPtCut, hiNtracksEtaCut, hiNtracksEtaPtCut;
  int hiNpixPlus, hiNpixMinus, hiNpixelTracksPlus, hiNpixelTracksMinus;
  float hiHF, hiHFplus, hiHFminus, hiHFplusEta4, hiHFminusEta4, hiHFhit, hiHFhitPlus, hiHFhitMinus;
  float hiHFECut, hiHFECutPlus, hiHFECutMinus;
  float hiEB, hiET, hiEE, hiEEplus, hiEEminus;
  float hiZDC, hiZDCplus, hiZDCminus;

  //MC variables

  float fNpart;
  float fNcoll;
  float fNhard;
  float fPhi0;
  float fb;

  int fNcharged;
  int fNchargedMR;
  float fMeanPt;
  float fMeanPtMR;
  float fEtMR;
  int fNchargedPtCut;
  int fNchargedPtCutMR;

  int proc_id;
  float pthat;
  float weight;
  float alphaQCD;
  float alphaQED;
  float qScale;
  int   nMEPartons;
  int   nMEPartonsFiltered;
  std::pair<int, int> pdfID;
  std::pair<float, float> pdfX;
  std::pair<float, float> pdfXpdf;

  std::vector<float> ttbar_w; //weights for systematics
  
  std::vector<int> npus;    //number of pileup interactions
  std::vector<float> tnpus; //true number of interactions

  //vertex variables

  float vx,vy,vz;

  //event variables

  unsigned long long event;
  unsigned int run;
  unsigned int lumi;

};

GOForest::GOForest(const edm::ParameterSet& iConfig): 
  CentralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralitySrc"))),
  CentralityBinTag_(consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc"))),
  HFfilters_(consumes<reco::HFFilterInfo>(iConfig.getParameter<edm::InputTag>("HFfilters"))),
  VertexTag_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("Vertex"))){
}
/*
  puInfoToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"))),
  genInfoToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  HiMCTag_(consumes<edm::GenHIEvent>(iConfig.getParameter<edm::InputTag>("HiMC"))),
  generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
  doCentrality_(iConfig.getParameter<bool> ("doCentrality")),
  doMC_(iConfig.getParameter<bool> ("doMC")),
  doHiMC_(iConfig.getParameter<bool> ("doHiMC")),
  useHepMC_(iConfig.getParameter<bool> ("useHepMC")),
  doVertex_(iConfig.getParameter<bool>("doVertex"))
*/


void GOForest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //cleanup previous event
  npus.clear();
  tnpus.clear();
  ttbar_w.clear();

  event = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();

  Ev_TTree->Fill();

  using namespace std;
  using namespace edm;
  using namespace reco;
//  using namespace pat;

  //access the vertex

  edm::Handle<std::vector<reco::Vertex>> vertex;
  iEvent.getByToken(VertexTag_, vertex);
  vx = vertex->begin()->x();
  vy = vertex->begin()->y();
  vz = vertex->begin()->z();

  //save primary vertex boolean 
  //
  std::cout << "vertex->begin()->isFake(): " << vertex->begin()->isFake() << endl;
  std::cout << "vertex->begin()->position().Rho(): "<< vertex->begin()->position().Rho() << endl;
  std::cout << "vertex->begin()->tracksSize(): " << vertex->begin()->tracksSize() << endl;
  pprimaryvertexfilter=0;
  if(!vertex->begin()->isFake() && fabs(vertex->begin()->z()) <= 25. && vertex->begin()->position().Rho() <= 2. /* || vertex->begin()->tracksSize() < 2*/){pprimaryvertexfilter=1;}

  //save the min number of tower, important to apply HF filter

  edm::Handle<reco::HFFilterInfo> HFfilter;
  iEvent.getByToken(HFfilters_, HFfilter);

  numMinHFTower2 = HFfilter->numMinHFTowers2; 
  numMinHFTower3 = HFfilter->numMinHFTowers3;
  numMinHFTower4 = HFfilter->numMinHFTowers4;
  numMinHFTower5 = HFfilter->numMinHFTowers5;

  //reserved for cluster compatibility filter


  Skim_TTree->Fill();

  //reserved for Isabela Studies

  edm::Handle<int> cbin_;
  iEvent.getByToken(CentralityBinTag_,cbin_);
  hiBin = *cbin_;

  GO_TTree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void GOForest::beginJob() {
  // please remove this method if not needed

  fNpart = -1;
  fNcoll = -1;
  fNhard = -1;
  fPhi0 = -1;
  fb = -1;
  fNcharged = -1;
  fNchargedMR = -1;
  fMeanPt = -1;
  fMeanPtMR = -1;

  fEtMR = -1;
  fNchargedPtCut = -1;
  fNchargedPtCutMR = -1;

  proc_id =   -1;
  pthat   =   -1.;
  weight  =   -1.;
  alphaQCD =  -1.;
  alphaQED =  -1.;
  qScale   =  -1.;
//  npu      =   1;
  hiBin = -1;
  vx = -100;
  vy = -100;
  vz = -100;
  numMinHFTower2 = -1;
  numMinHFTower3 = -1;
  numMinHFTower4 = -1;
  numMinHFTower5 = -1;

  TFileDirectory Event_Info = fs->mkdir("EventInfo");
  Ev_TTree = Event_Info.make<TTree>("EventInfo","EventInfo");
  Ev_TTree->Branch("run",&run,"run/I");
  Ev_TTree->Branch("evt",&event,"evt/l");
  Ev_TTree->Branch("lumi",&lumi,"lumi/I");

  TFileDirectory skim_analysis = fs->mkdir("skimanalysis");
  Skim_TTree = skim_analysis.make<TTree>("HltTree","HltTree");
  Skim_TTree->Branch("numMinHFTower2",&numMinHFTower2,"numMinHFTower2/I");
  Skim_TTree->Branch("numMinHFTower3",&numMinHFTower3,"numMinHFTower3/I");
  Skim_TTree->Branch("numMinHFTower4",&numMinHFTower4,"numMinHFTower4/I");
  Skim_TTree->Branch("numMinHFTower5",&numMinHFTower5,"numMinHFTower5/I");
  Skim_TTree->Branch("pprimaryvertexfilter",&pprimaryvertexfilter,"pprimaryvertexfilter/O");

  TFileDirectory hiEvt_Analyzer = fs->mkdir("hiEvtAnalyzer");
  GO_TTree = hiEvt_Analyzer.make<TTree>("HiTree","HiTree");
  GO_TTree->Branch("hiBin",&hiBin,"hiBin/I");
  GO_TTree->Branch("vx",&vx,"vx/F");
  GO_TTree->Branch("vy",&vy,"vy/F");
  GO_TTree->Branch("vz",&vz,"vz/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void GOForest::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GOForest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GOForest);
