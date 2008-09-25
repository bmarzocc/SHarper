#include "SHarper/HEEPAnalyzer/interface/HEEPAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "TH1D.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


//delete these
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPDebug.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" 

HEEPAnalyzer::HEEPAnalyzer(const edm::ParameterSet& iPara):
  evtHelper_(),heepEvt_(),nrPass_(0),nrFail_(0)
{
  evtHelper_.setup(iPara);

}

void HEEPAnalyzer::beginJob(const edm::EventSetup& iSetup)
{
  edm::Service<TFileService> fs;
  massHist_ = fs->make<TH1D>("massHist","Di-Electron Mass;M_{ee} (GeV/c^{2});# Events / 5 GeV/c^{2}",290,50,1500); 
}

void HEEPAnalyzer::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup)
{ 
//heep::listAllProducts<EcalRecHitCollection>(iEvent,"HEEPAnalyzer");
 
  //make the heep event (see easy isnt it)
  evtHelper_.makeHeepEvent(iEvent,iSetup,heepEvt_);
  std::cout <<"made event "<<std::endl;
 
  //edm::Handle<reco::IsoDepositMap> ecalIsolHandle;
  //iEvent.getByLabel("eleIsoDepositEcalFromHitsFast",ecalIsolHandle);
  //const reco::IsoDepositMap& ecalIsol = *ecalIsolHandle.product();

  

  //do what ever you want
  //lets get the heep electrons and count the number that pass / fail cuts
  const std::vector<heep::Ele>& eles = heepEvt_.heepElectrons();
  for(size_t eleNr=0;eleNr<eles.size();eleNr++){
    if(eles[eleNr].cutCode()==0x0) nrPass_++;
    else nrFail_++;
  }

  //lets fill a mass hist
  for(size_t ele1Nr=0;ele1Nr<eles.size();ele1Nr++){
    for(size_t ele2Nr=ele1Nr+1;ele2Nr<eles.size();ele2Nr++){
      const heep::Ele& ele1 = eles[ele1Nr];
      const heep::Ele& ele2 = eles[ele2Nr];
      
      std::cout <<"ele1 ecal isol "<<ele1.isolEm()<<" ele 2 ecal isol "<<ele2.isolEm()<<std::endl;
      
      //const reco::IsoDeposit& iso = ecalIsol[ele1.gsfEle()];
      //std::cout <<"isol "<<iso.sumWithin(0.4)<<std::endl;

      float mass = (ele1.p4()+ele2.p4()).mag();
      massHist_->Fill(mass);
    }
  }
    
}


void HEEPAnalyzer::endJob()
{
  std::cout <<"Nr pass "<<nrPass_<<" nr fail "<<nrFail_<<" eff "<<static_cast<float>(nrPass_)/static_cast<float>(nrFail_)<<std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HEEPAnalyzer);
