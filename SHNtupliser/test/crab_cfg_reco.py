# CRAB3 config template for flashgg
# More options available on the twiki :
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial

from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName       = 'FourGammasGunPt1-100_DeepSC_algoD_idealICs'
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName          = 'egRegTreeMakerRefinedAOD_DeepSC_algoD_idealICs.py'
config.JobType.priority          = 30
config.JobType.maxMemoryMB       = 3000
config.JobType.numCores          = 2

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
#config.Data.inputDataset         = '/FourElectronsGunPt1-100_pythia8_StdMix_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15/dpg_ecal-GEN-SIM-RAW-Reduced-f73d3671898879ae5e675fcbfe8fb906/USER'
config.Data.inputDataset         = '/FourGammasGunPt1-100_pythia8_StdMix_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15/dpg_ecal-GEN-SIM-RAW-Reduced-f73d3671898879ae5e675fcbfe8fb906/USER'
config.JobType.outputFiles       = ["egmRegTree.root"]
config.Data.inputDBS             = 'global'   
config.Data.inputDBS             = 'phys03'   
config.Data.splitting            = 'FileBased'
config.Data.unitsPerJob          = 1
config.Data.publication          = False
config.Data.ignoreLocality       = True
config.Data.outputDatasetTag     = 'RECO_UL18_pfRechitThres_DeepSC_algoD_idealICs'

# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase        =  '/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite         = 'T2_CH_CERN'
config.Site.whitelist           = ['T2_CH_CERN']


## config.Data.allowNonValidInputDataset=True
