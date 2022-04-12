
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

class L1EGammaFilter : public edm::stream::EDFilter<> {

private:
  edm::EDGetTokenT<BXVector<l1t::EGamma> > l1EGToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  int nrPass_;
  int nrTot_;
  int nrL1EGsRequired_;
  float l1EGEtCut_;
  
public:
  explicit L1EGammaFilter(const edm::ParameterSet& para);
  virtual ~L1EGammaFilter(){}
  
  virtual bool filter(edm::Event& event,const edm::EventSetup& setup);
  virtual void endJob(){};
};

L1EGammaFilter::L1EGammaFilter(const edm::ParameterSet& para):nrPass_(0),nrTot_(0)

{
  l1EGToken_ = consumes<BXVector<l1t::EGamma>>(para.getParameter<edm::InputTag>("l1EGTag"));
  genEventInfoToken_ = consumes<GenEventInfoProduct>(edm::InputTag("generator", ""));
  nrL1EGsRequired_= para.getParameter<int>("nrL1EGsRequired"); 
  l1EGEtCut_ = para.getParameter<double>("l1EGEtCut");
}

bool L1EGammaFilter::filter(edm::Event& event,const edm::EventSetup& setup)
{   
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  event.getByToken(genEventInfoToken_, genEvtInfo);
  
  edm::Handle<BXVector<l1t::EGamma> > l1EGammas;
  event.getByToken(l1EGToken_, l1EGammas);

  int nrEGPass = 0;
  for(auto egIt = l1EGammas->begin(0);egIt != l1EGammas->end(0);egIt++){
    if(egIt->et()>=nrL1EGsRequired_){
      nrEGPass++;
    }
  }
  return nrEGPass>=nrL1EGsRequired_;
}  

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1EGammaFilter);
