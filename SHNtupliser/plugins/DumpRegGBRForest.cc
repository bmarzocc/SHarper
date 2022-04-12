
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CondFormats/DataRecord/interface/GBRDWrapperRcd.h"
#include "CondFormats/GBRForest/interface/GBRForestD.h"

#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


class DumpRegGBRForest : public edm::one::EDAnalyzer<edm::one::SharedResources> { 


public:
  explicit DumpRegGBRForest(const edm::ParameterSet& iPara);
  virtual ~DumpRegGBRForest(){}
  
  DumpRegGBRForest(const DumpRegGBRForest& rhs)=delete;
  DumpRegGBRForest& operator=(const DumpRegGBRForest& rhs)=delete;

 private:
  template<typename T>
  void getToken(edm::EDGetTokenT<T>& token,const edm::ParameterSet& iPara,const std::string& paraName){
    token = consumes<T>(iPara.getParameter<edm::InputTag>(paraName));
  }
  
  virtual void beginJob(){}
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& iSetup){}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob(){}
  
private:
  edm::ESGetToken<GBRForestD, GBRDWrapperRcd> forestToken_;
  std::vector<std::string> names_;
  bool written_;
};


DumpRegGBRForest::DumpRegGBRForest(const edm::ParameterSet& iPara):
  names_(iPara.getParameter<std::vector<std::string>>("names")),
  written_(false)
{


}

void DumpRegGBRForest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(!written_){
    edm::Service<TFileService> fs;
    edm::ESHandle<GBRForestD> gbrHandle;
    for(const auto& name : names_){
      forestToken_ = esConsumes<GBRForestD, GBRDWrapperRcd>(edm::ESInputTag("", name.c_str()));
      gbrHandle = iSetup.getHandle(forestToken_); 
      fs->file().WriteObject(gbrHandle.product(),name.c_str());
    }
    written_=true;
  }
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DumpRegGBRForest);
