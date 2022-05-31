#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"	
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

#include "SHarper/TrigNtup/interface/EGRegTreeStruct.hh"

#include "TFile.h"
#include "TTree.h"

#include <string>
#include <vector>
#include <array>

class EGRegTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {

private:
  EGRegTreeStruct egRegTreeData_;
  TTree* egRegTree_;
  std::string treeName_;
  
  edm::ESGetToken<CaloTopology,CaloTopologyRecord> caloTopologyToken_;
  edm::ESGetToken<EcalChannelStatus,EcalChannelStatusRcd> channelStatusToken_;
  edm::EDGetTokenT<reco::VertexCollection>  verticesToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPartsToken_;
  edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_;
  edm::EDGetTokenT<std::vector<reco::PFCluster> > pfClusterToken_; 
  std::vector<edm::EDGetTokenT<reco::SuperClusterCollection>> scTokens_;
  std::vector<edm::EDGetTokenT<reco::SuperClusterCollection>> scAltTokens_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalHitsEBToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalHitsEEToken_;
  edm::EDGetTokenT<std::vector<reco::GsfElectron> > elesToken_;
  edm::EDGetTokenT<std::vector<reco::Photon> > phosToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puSumToken_;

  std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectron> > > eleAltTokens_;
  std::vector<edm::EDGetTokenT<std::vector<reco::Photon> > > phoAltTokens_;

  EGRegTreeMaker(const EGRegTreeMaker& rhs)=delete;
  EGRegTreeMaker& operator=(const EGRegTreeMaker& rhs)=delete;

  std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart; 
  
  bool fillFromSC_;
  bool fillFromMC_;
  bool fillFromSIM_;
  
public:
  explicit EGRegTreeMaker(const edm::ParameterSet& iPara);
  virtual ~EGRegTreeMaker();
  
private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& iSetup);
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endRun(edm::Run const& iRun, edm::EventSetup const&){}
  virtual void endJob();
  template<typename T>
  void setToken(edm::EDGetTokenT<T>& token,const edm::ParameterSet& iPara,const std::string& tagName){
    token = consumes<T>(iPara.getParameter<edm::InputTag>(tagName));
  }
  template<typename T>
  void setToken(std::vector<edm::EDGetTokenT<T> > & tokens,const edm::ParameterSet& iPara,const std::string& tagName){
    for(auto& tag: iPara.getParameter<std::vector<edm::InputTag> >(tagName)){
      tokens.push_back(consumes<T>(tag));
    }
  }
  static const reco::GenParticle* matchGenPart(float eta,float phi,const std::vector<reco::GenParticle>& genParts);
  static const CaloParticle* matchCaloPart(const reco::GenParticle& genPart,const std::vector<reco::GenParticle>& genParts,const std::vector<CaloParticle>& caloParts);
  const std::vector<reco::PFCluster> getBestMatchedPFCluster(const std::vector<CaloParticle>& caloParts,const std::vector<reco::PFCluster>& pfClusters);
  const reco::SuperCluster* matchSCToCalo(const CaloParticle* caloToMatch,int iCalo,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles, std::vector<reco::PFCluster>* bestPFClusters);
  const reco::SuperCluster*  matchSC(const reco::SuperCluster* scToMatch,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles);
  static const reco::SuperCluster*  matchSC(float eta,float phi,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles);
  std::vector<std::pair<DetId, float> >* getHitsAndEnergiesCaloPart(const CaloParticle* iCaloParticle, float simHitEnergy_cut);
  std::vector<double> getScores(const reco::CaloClusterPtr seed, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart);
  std::vector<double> getScores(const reco::PFCluster* pfCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart);
  bool isInWindow(double genEta,double genPhi,double seedEta,double seedPhi);
  double deltaPhi(double phi1, double phi2);
  std::array<double,3> dynamicWindow(double eta);
  

};



EGRegTreeMaker::EGRegTreeMaker(const edm::ParameterSet& iPara):
  egRegTree_(nullptr),
  treeName_("egRegTree"),
  caloTopologyToken_(esConsumes()),
  channelStatusToken_(esConsumes()),
  fillFromSC_(false),
  fillFromMC_(false),
  fillFromSIM_(false)
{
  if(iPara.exists("fillFromSC")){
     fillFromSC_ = iPara.getParameter<bool>("fillFromSC"); 
  } 
  if(iPara.exists("fillFromMC")){
     fillFromMC_ = iPara.getParameter<bool>("fillFromMC"); 
  } 
  if(iPara.exists("fillFromSIM")){
     fillFromSIM_ = iPara.getParameter<bool>("fillFromSIM"); 
  } 
  if(iPara.exists("treeName")){
    treeName_ = iPara.getParameter<std::string>("treeName");
  }

  setToken(verticesToken_,iPara,"verticesTag");
  setToken(rhoToken_,iPara,"rhoTag");
  setToken(genPartsToken_,iPara,"genPartsTag");
  setToken(caloPartToken_,iPara,"caloPartsTag");
  setToken(pfClusterToken_,iPara,"pfClusTag");
  setToken(scTokens_,iPara,"scTag");
  setToken(scAltTokens_,iPara,"scAltTag");
  setToken(ecalHitsEBToken_,iPara,"ecalHitsEBTag");
  setToken(ecalHitsEEToken_,iPara,"ecalHitsEETag");
  setToken(elesToken_,iPara,"elesTag");
  setToken(phosToken_,iPara,"phosTag");
  setToken(puSumToken_,iPara,"puSumTag");
  setToken(eleAltTokens_,iPara,"elesAltTag");
  setToken(phoAltTokens_,iPara,"phosAltTag");

}

EGRegTreeMaker::~EGRegTreeMaker()
{

}


void EGRegTreeMaker::beginJob()
{
  edm::Service<TFileService> fs;
  fs->file().cd();
  egRegTree_ = new TTree(treeName_.c_str(),"");
  egRegTreeData_.setNrEnergies(eleAltTokens_.size(),phoAltTokens_.size());
  egRegTreeData_.createBranches(egRegTree_);

  if((fillFromSC_&& fillFromMC_) || (fillFromSC_&& fillFromSIM_) || (fillFromMC_&& fillFromSIM_) || (fillFromSC_&& fillFromMC_ && fillFromSIM_) || (!fillFromSC_&& !fillFromMC_ && !fillFromSIM_)){
     std::cout << "WARNING: Filling options are exclusive! Setting to FillFromMC" << std::endl;
     fillFromSC_  = false;
     fillFromMC_  = true;
     fillFromSIM_ = false;
  }
  if(fillFromSC_) std::cout <<  "Running in \"fillFromSC\" mode..." << std::endl; 
  if(fillFromMC_) std::cout <<  "Running in \"fillFromMC\" mode..." << std::endl; 
  if(fillFromSIM_) std::cout << "Running in \"fillFromSIM\" mode..." << std::endl; 
} 

void EGRegTreeMaker::beginRun(const edm::Run& run,const edm::EventSetup& iSetup)
{ 
 
}

namespace {
  template<typename T> 
  edm::Handle<T> getHandle(const edm::Event& iEvent,const edm::EDGetTokenT<T>& token)
  {
    edm::Handle<T> handle;
    iEvent.getByToken(token,handle);
    return handle;
  }
  template<typename T> 
  std::vector<edm::Handle<T> > getHandle(const edm::Event& iEvent,const std::vector<edm::EDGetTokenT<T> >& tokens)
  {
    std::vector<edm::Handle<T> > handles;
    for(auto& token : tokens){
      edm::Handle<T> handle;
      iEvent.getByToken(token,handle);
      handles.emplace_back(std::move(handle));
    }
    return handles;
  }
}
  
namespace{

  bool hasBasicClusters(const reco::SuperCluster& sc){
    if(!sc.seed().isAvailable()) return false;
    else if(!sc.clusters().isAvailable()) return false;
    else return true;
  }
  template<typename T> 
  const T* matchBySCSeedId(unsigned int seedId,const std::vector<T>& objs){
    for(auto& obj : objs){
      if(obj.superCluster()->seed()->seed().rawId()==seedId) return &obj;
    }
    return nullptr;
  }
  const reco::GsfElectron* matchEle(unsigned int seedId,const std::vector<reco::GsfElectron>& eles){
    return matchBySCSeedId(seedId,eles);
  }
  const reco::Photon* matchPho(unsigned int seedId,const std::vector<reco::Photon>& phos){
    return matchBySCSeedId(seedId,phos);
  }
  
  template<typename T> std::vector<const T*> 
  matchToAltCollsBySCSeedId(const T* objToMatch,const std::vector<edm::Handle<std::vector<T> > >& objCollHandles){
    std::vector<const T*> matches;
    if(objToMatch){
      unsigned int seedId = objToMatch->superCluster()->seed()->seed().rawId();
      for(auto& objCollHandle : objCollHandles){
	if(objCollHandle.isValid()) matches.push_back(matchBySCSeedId(seedId,*objCollHandle));
	else matches.push_back(nullptr);
      }
    }
    else{
      matches.resize(objCollHandles.size(),nullptr);
    }
    return matches;
  }
}

void EGRegTreeMaker::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  auto scHandles = getHandle(iEvent,scTokens_);
  auto scAltHandles = getHandle(iEvent,scAltTokens_);
  auto ecalHitsEBHandle = getHandle(iEvent,ecalHitsEBToken_);
  auto ecalHitsEEHandle = getHandle(iEvent,ecalHitsEEToken_);
  auto pfClustersHandle = getHandle(iEvent,pfClusterToken_);
  auto genPartsHandle = getHandle(iEvent,genPartsToken_);
  auto caloPartsHandle = getHandle(iEvent,caloPartToken_);
  auto verticesHandle = getHandle(iEvent,verticesToken_);
  auto rhoHandle = getHandle(iEvent,rhoToken_);
  auto elesHandle = getHandle(iEvent,elesToken_);
  auto eleAltHandles = getHandle(iEvent,eleAltTokens_);
  auto phosHandle = getHandle(iEvent,phosToken_);
  auto phoAltHandles = getHandle(iEvent,phoAltTokens_);
  auto puSumHandle = getHandle(iEvent,puSumToken_);
  edm::ESHandle<CaloTopology> caloTopoHandle = iSetup.getHandle(caloTopologyToken_);
  edm::ESHandle<EcalChannelStatus> chanStatusHandle = iSetup.getHandle(channelStatusToken_);
  
  int nrVert = verticesHandle->size();
  float nrPUInt = -1;
  float nrPUIntTrue = -1;
  if(puSumHandle.isValid()){
    for(auto& puInfo : *puSumHandle){
      if(puInfo.getBunchCrossing()==0){
	nrPUInt = puInfo.getPU_NumInteractions();
	nrPUIntTrue = puInfo.getTrueNumInteractions();
	break;
      }
    }
  }

  if(fillFromSC_){
    //std::cout <<  "Running in \"fillFromSC\" mode..." << std::endl; 
    for(const auto& scHandle : scHandles){
      for(const auto& sc: *scHandle){
	const reco::GenParticle* genPart = genPartsHandle.isValid() ? matchGenPart(sc.eta(),sc.phi(),*genPartsHandle) : nullptr;
        const CaloParticle* caloPart = (genPartsHandle.isValid() && caloPartsHandle.isValid() && genPart!=nullptr) ? matchCaloPart(*genPart,*genPartsHandle,*caloPartsHandle) : nullptr;
	const reco::GsfElectron* ele = elesHandle.isValid() ? matchEle(sc.seed()->seed().rawId(),*elesHandle) : nullptr;
	const reco::Photon* pho = phosHandle.isValid() ? matchPho(sc.seed()->seed().rawId(),*phosHandle) : nullptr;
	const reco::SuperCluster* scAlt = matchSC(&sc,scAltHandles);
	
	const std::vector<const reco::GsfElectron*> altEles = matchToAltCollsBySCSeedId(ele,eleAltHandles);
	const std::vector<const reco::Photon*> altPhos = matchToAltCollsBySCSeedId(pho,phoAltHandles);

	if(hasBasicClusters(sc)){
	  egRegTreeData_.fill(iEvent,nrVert,*rhoHandle,nrPUInt,nrPUIntTrue,
			      *ecalHitsEBHandle,*ecalHitsEEHandle,
			      *caloTopoHandle,
			      *chanStatusHandle,
			      &sc,genPart,caloPart,ele,pho,scAlt,
			      altEles,altPhos);
	  egRegTree_->Fill();
	}
      }
    }
  }else if(fillFromSIM_){
    //std::cout <<  "Running in \"fillFromSIM\" mode..." << std::endl; 
    std::vector<reco::PFCluster> bestPFCluster = getBestMatchedPFCluster(*caloPartsHandle,*pfClustersHandle);

    int iCalo=0;
    for(const auto& caloPart : *caloPartsHandle){
      if((std::abs(caloPart.pdgId())==11 || caloPart.pdgId()==22)){
        const reco::SuperCluster* sc = matchSCToCalo(&caloPart,iCalo,scHandles,&bestPFCluster);
	const reco::SuperCluster* scAlt = matchSC(sc,scAltHandles);
	const reco::GsfElectron* ele = elesHandle.isValid() && sc ? matchEle(sc->seed()->seed().rawId(),*elesHandle) : nullptr;
	const reco::Photon* pho = phosHandle.isValid() && sc ? matchPho(sc->seed()->seed().rawId(),*phosHandle) : nullptr;
        const reco::GenParticle* genPart = genPartsHandle.isValid() ? &((*genPartsHandle)[caloPart.g4Tracks()[0].genpartIndex()-1]) : nullptr;
	
	const std::vector<const reco::GsfElectron*> altEles = matchToAltCollsBySCSeedId(ele,eleAltHandles);
	const std::vector<const reco::Photon*> altPhos = matchToAltCollsBySCSeedId(pho,phoAltHandles);


	egRegTreeData_.fill(iEvent,nrVert,*rhoHandle,nrPUInt,nrPUIntTrue,
			    *ecalHitsEBHandle,*ecalHitsEEHandle,
			    *caloTopoHandle,
			    *chanStatusHandle,
			    sc,genPart,&caloPart,ele,pho,scAlt,
			    altEles,altPhos);
	egRegTree_->Fill();
      }
      iCalo++;
    }
  }else if(fillFromMC_){
    //std::cout <<  "Running in \"fillFromMC\" mode..." << std::endl; 
    for(const auto& genPart : *genPartsHandle){
      if((std::abs(genPart.pdgId())==11 || genPart.pdgId()==22) && genPart.statusFlags().isPrompt() && genPart.statusFlags().isFirstCopy()){
	const reco::SuperCluster* sc = matchSC(genPart.eta(),genPart.phi(),scHandles);
	const reco::SuperCluster* scAlt = matchSC(sc,scAltHandles);
	const reco::GsfElectron* ele = elesHandle.isValid() && sc ? matchEle(sc->seed()->seed().rawId(),*elesHandle) : nullptr;
	const reco::Photon* pho = phosHandle.isValid() && sc ? matchPho(sc->seed()->seed().rawId(),*phosHandle) : nullptr;
        const CaloParticle* caloPart = caloPartsHandle.isValid() ? matchCaloPart(genPart,*genPartsHandle,*caloPartsHandle) : nullptr;
	
	const std::vector<const reco::GsfElectron*> altEles = matchToAltCollsBySCSeedId(ele,eleAltHandles);
	const std::vector<const reco::Photon*> altPhos = matchToAltCollsBySCSeedId(pho,phoAltHandles);


	egRegTreeData_.fill(iEvent,nrVert,*rhoHandle,nrPUInt,nrPUIntTrue,
			    *ecalHitsEBHandle,*ecalHitsEEHandle,
			    *caloTopoHandle,
			    *chanStatusHandle,
			    sc,&genPart,caloPart,ele,pho,scAlt,
			    altEles,altPhos);
	egRegTree_->Fill();
      }
    }
    
  } 
}

const reco::GenParticle*  EGRegTreeMaker::matchGenPart(float eta,float phi,const std::vector<reco::GenParticle>& genParts)
{
  const reco::GenParticle* bestMatch=nullptr;
  float bestDR2=0.3*0.3;
  for(const auto& genPart : genParts){
    if(std::abs(genPart.pdgId())==11){
      if(genPart.statusFlags().isPrompt() && genPart.statusFlags().isFirstCopy()){
	float dR2 = reco::deltaR2(genPart.eta(),genPart.phi(),eta,phi);
	if(dR2<bestDR2){
	  bestMatch = &genPart;
	  bestDR2 = dR2;
	}
	  
      }
    }
  }
  return bestMatch;
}

const CaloParticle*  EGRegTreeMaker::matchCaloPart(const reco::GenParticle& genPart,const std::vector<reco::GenParticle>& genParts,const std::vector<CaloParticle>& caloParts)
{
  const CaloParticle* bestMatch=nullptr;
  for(const auto& caloPart : caloParts){
    const auto genParticle = genParts[caloPart.g4Tracks()[0].genpartIndex()-1]; 
    float dR2 = reco::deltaR2(genParticle.eta(),genParticle.phi(),genPart.eta(),genPart.phi()); 
    float dE = fabs(genParticle.energy()-genPart.energy()); 
    if(genPart.statusFlags().isPrompt() && genPart.statusFlags().isFirstCopy() && dR2<1.e-6 && dE<1.e-6) bestMatch = &caloPart;
  }
  return bestMatch;
}

const std::vector<reco::PFCluster> EGRegTreeMaker::getBestMatchedPFCluster(const std::vector<CaloParticle>& caloParts,const std::vector<reco::PFCluster>& pfClusters)
{
  std::vector<reco::PFCluster> bestPFClusters;
  bestPFClusters.resize(caloParts.size());
 
  int iCalo=0;
  for(const auto& caloPart : caloParts){
    double maxScore = -1.;
    std::vector<std::pair<DetId, float> > hitsAndEnergies = *getHitsAndEnergiesCaloPart(&caloPart,-1.);
    for(const auto& pfCluster : pfClusters){
      std::vector<double> scores = getScores(&pfCluster,&hitsAndEnergies);   
      if(scores[0]>maxScore){
         bestPFClusters[iCalo] = pfCluster;
         maxScore = scores[0];
      }
    }
    iCalo++; 
  }
  return bestPFClusters;
}

const reco::SuperCluster* EGRegTreeMaker::matchSCToCalo(const CaloParticle* caloToMatch,int iCalo,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles, std::vector<reco::PFCluster>* bestPFClusters)
{
  const reco::SuperCluster* bestMatch=nullptr;
  std::vector<std::pair<DetId, float> > hitsAndEnergies = *getHitsAndEnergiesCaloPart(caloToMatch,-1.);
  double maxScore = -1.;
  for(const auto& scHandle : scHandles){
      if(scHandle.isValid()){
	for(const auto& sc: *scHandle){
          std::vector<double> scores = getScores(sc.seed(),&hitsAndEnergies);  
          bool sameCluster = reco::CaloCluster((*bestPFClusters)[iCalo])==*sc.seed();
          if(scores[0]>maxScore && scores[0]>1e-2 && sameCluster){
             bestMatch = &sc;
             maxScore = scores[0];
          }
        }
      }
  }

  if(bestMatch!=nullptr && !isInWindow(caloToMatch->eta(),caloToMatch->phi(),bestMatch->seed()->eta(),bestMatch->seed()->phi())) return nullptr;
  return bestMatch;  
}


//matched by seed det id which should be unique
const reco::SuperCluster* EGRegTreeMaker::matchSC(const reco::SuperCluster* scToMatch,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles)
{
  if(scToMatch){
    auto seedId = scToMatch->seed()->seed().rawId();
    for(const auto& scHandle : scHandles){
      if(scHandle.isValid()){
	for(const auto& sc: *scHandle){
	  if(sc.seed()->seed().rawId()==seedId){
	    return &sc;
	  }
	}
      }
    }
  }
  return nullptr;
}

const reco::SuperCluster*  EGRegTreeMaker::matchSC(float eta,float phi,const std::vector<edm::Handle<reco::SuperClusterCollection> >& scHandles)
{
  const reco::SuperCluster* bestMatch=nullptr;
  float bestDR2=0.3*0.3;
  for(const auto& scHandle : scHandles){
    for(const auto& sc: *scHandle){
      float dR2 = reco::deltaR2(sc.eta(),sc.phi(),eta,phi);
      if(dR2<bestDR2){
	bestMatch = &sc;
	bestDR2 = dR2;
      }
    }
  }
  return bestMatch;
}

std::vector<std::pair<DetId, float> >* EGRegTreeMaker::getHitsAndEnergiesCaloPart(const CaloParticle* iCaloParticle, float simHitEnergy_cut)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_final = new std::vector<std::pair<DetId, float> >;
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;
    
    const auto& simClusters = iCaloParticle->simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
        auto simCluster = simClusters[iSC];  
        auto hits_and_energies = simCluster->hits_and_energies();
        for(unsigned int i = 0; i < hits_and_energies.size(); i++){ 
            if(hits_and_energies[i].second < simHitEnergy_cut) continue; 
            HitsAndEnergies_tmp->push_back(std::make_pair(DetId(hits_and_energies[i].first),hits_and_energies[i].second));  
        }  
    }

    for(unsigned int i = 0; i < HitsAndEnergies_tmp->size(); i++){  
        if (HitsAndEnergies_map.find(HitsAndEnergies_tmp->at(i).first) == HitsAndEnergies_map.end()) {
            HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]=HitsAndEnergies_tmp->at(i).second;      
        }else{
            HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]=HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]+HitsAndEnergies_tmp->at(i).second; 
        }
    }

    for(auto const& hit : HitsAndEnergies_map) 
        HitsAndEnergies_final->push_back(std::make_pair(hit.first,hit.second));

    return HitsAndEnergies_final;
}

std::vector<double> EGRegTreeMaker::getScores(const reco::CaloClusterPtr seed, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart)
{
    std::vector<double> scores;
    scores.resize(2);

    double simFraction=0.; 
    double simFraction_noHitsFraction=0.;
    double simEnergy=0.;
    double simEnergy_shared=0.;
    double simEnergy_shared_noHitsFraction=0.;
   
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart)
        simEnergy+=hit_CaloPart.second;
   
    const std::vector<std::pair<DetId,float> > hitsAndFractions = seed->hitsAndFractions();
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){
        for(const std::pair<DetId, float>& hit_Cluster : hitsAndFractions){     
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               simEnergy_shared+=hit_CaloPart.second*hit_Cluster.second;
               simEnergy_shared_noHitsFraction+=hit_CaloPart.second;
            } 
        }
    }

    if(simEnergy>0.) simFraction_noHitsFraction = simEnergy_shared_noHitsFraction/simEnergy;
    else simFraction_noHitsFraction = -1.; 

    if(simEnergy>0.) simFraction = simEnergy_shared/simEnergy;
    else simFraction = -1.;
   
    scores[0] = simFraction;
    scores[1] = simFraction_noHitsFraction;
    
    for(unsigned iVar=0; iVar<scores.size(); iVar++)
        if(std::isnan(scores.at(iVar))) std::cout << "score = " << iVar << " ---> NAN " << std::endl; 

    return scores;
}

std::vector<double> EGRegTreeMaker::getScores(const reco::PFCluster* pfCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart)
{
    std::vector<double> scores;
    scores.resize(2);

    double simFraction=0.; 
    double simFraction_noHitsFraction=0.;
    double simEnergy=0.;
    double simEnergy_shared=0.;
    double simEnergy_shared_noHitsFraction=0.;
   
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart)
        simEnergy+=hit_CaloPart.second;
   
    const std::vector<std::pair<DetId,float> > hitsAndFractions = pfCluster->hitsAndFractions();
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){
        for(const std::pair<DetId, float>& hit_Cluster : hitsAndFractions){     
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               simEnergy_shared+=hit_CaloPart.second*hit_Cluster.second;
               simEnergy_shared_noHitsFraction+=hit_CaloPart.second;
            } 
        }
    }

    if(simEnergy>0.) simFraction_noHitsFraction = simEnergy_shared_noHitsFraction/simEnergy;
    else simFraction_noHitsFraction = -1.; 

    if(simEnergy>0.) simFraction = simEnergy_shared/simEnergy;
    else simFraction = -1.;
   
    scores[0] = simFraction;
    scores[1] = simFraction_noHitsFraction;
    
    for(unsigned iVar=0; iVar<scores.size(); iVar++)
        if(std::isnan(scores.at(iVar))) std::cout << "score = " << iVar << " ---> NAN " << std::endl; 

    return scores;
}

bool EGRegTreeMaker::isInWindow(double genEta,double genPhi,double seedEta,double seedPhi)
{
    //Iz computation
    int genIz = 0;
    if(abs(genEta)>1.479 && genEta>0.) genIz = 1;
    if(abs(genEta)>1.479 && genEta<0.) genIz = -1;
    
    int seedIz = 0;
    if(abs(seedEta)>1.479 && seedEta>0.) seedIz = 1;
    if(abs(seedEta)>1.479 && seedEta<0.) seedIz = -1;

    if(genIz != seedIz) return false;

    //DeltaEta
    double etaw = seedEta-genEta;
    if(genEta<0.) etaw = -etaw; 
    
    //DeltaPhi
    double phiw = deltaPhi(genPhi, seedPhi);

    //Seed dynamic windows
    std::array<double,3> deltas = dynamicWindow(seedEta);

    if(etaw>=deltas[0] && etaw<=deltas[1] && abs(phiw)<=deltas[2]) return true;
    else return false;
}

double EGRegTreeMaker::deltaPhi(double phi1, double phi2)
{
    double dphi = phi1 - phi2;
    if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
    return dphi;
}

std::array<double,3> EGRegTreeMaker::dynamicWindow(double eta)
{
    double aeta = abs(eta);

    double deta_up = 0.;
    if(aeta >= 0 && aeta < 0.1) deta_up = 0.075;
    if(aeta >= 0.1 && aeta < 1.3) deta_up = 0.0758929 -0.0178571* aeta + 0.0892857*(aeta*aeta); 
    else if(aeta >= 1.3 && aeta < 1.7) deta_up = 0.2;
    else if(aeta >=1.7 && aeta < 1.9) deta_up = 0.625 -0.25*aeta;
    else if(aeta >= 1.9) deta_up = 0.15;

    double deta_down = 0.;
    if(aeta < 2.1) deta_down = -0.075;
    else if(aeta >= 2.1 && aeta < 2.5) deta_down = -0.1875 *aeta + 0.31875;
    else if(aeta >=2.5) deta_down = -0.15;
    
    double dphi = 0.;        
    if(aeta < 1.9) dphi = 0.6;
    else if(aeta >= 1.9 && aeta < 2.7) dphi = 1.075 -0.25 * aeta;
    else if(aeta >= 2.7) dphi = 0.4;
     
    return std::array<double,3>{{deta_down,deta_up,dphi}};
}

void EGRegTreeMaker::endJob()
{ 

}


  


//define this as a plug-in
DEFINE_FWK_MODULE(EGRegTreeMaker);
