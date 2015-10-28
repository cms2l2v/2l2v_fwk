from WMCore.Configuration import Configuration

import os
binfile = os.path.join(os.environ['CMSSW_BASE'], 'bin',
                       os.environ['SCRAM_ARCH'], 'runPhoJetAnalysis')

#print binfile

config = Configuration()

config.section_("General")
#config.General.requestName = 'GJets_HT-100to200'
#config.General.requestName = 'GJet_Pt-15to3000'
#config.General.requestName = 'QCD_Pt-15to3000'
config.General.requestName = 'DYJetsToLL_M-50'
#config.General.requestName = 'GJet_Pt-15To6000'
config.General.workArea = 'crab_phojet'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'phojet_cfg.py'
config.JobType.scriptExe = 'crabExe.sh'
config.JobType.inputFiles = [#'FrameworkJobReport.xml', #'x509_proxy',
                             binfile, 'phojet_cfg.py']
config.JobType.outputFiles = ['output.root']

config.section_("Data")
#config.Data.inputDataset = '/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = "/GJet_Pt-15to3000_Tune4C_13TeV_pythia8/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM"
#config.Data.inputDataset = "/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM"
#config.Data.inputDataset = "/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
config.Data.inputDataset = "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM"
#config.Data.inputDataset = "/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM"
#config.Data.inputDataset = "/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM"

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.publication = True
#config.Data.publishDBS = 'phys03'
#config.Data.publishDataName = 'phojet_analysis'

config.section_("Site")
config.Site.storageSite = 'T2_US_Purdue'
