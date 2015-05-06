import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.PileupJetIDParams_cfi import cutbased as pu_jetid

process = cms.Process("PhoJet")

def lfn_to_pfn(f):
    import subprocess
    proc = subprocess.Popen(["edmFileUtil -d %s" %f],
                            stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    pfn = out.strip()
    return pfn

inputFiles = cms.untracked.vstring("file:input.root")
try:
    import PSet
    print "loading from PSet..."
    fnames = list(PSet.process.source.fileNames)
    fnames = [ lfn_to_pfn(f) for f in fnames]    
    inputFiles =  cms.untracked.vstring(fnames)                
except:
    print "No PSet found, using local file..."
    inputFiles =  cms.untracked.vstring('file:input.root')  

print inputFiles

process.source = cms.Source("PoolSource", fileNames = inputFiles)

runProcess = cms.PSet(
    input = inputFiles, 
    outfile = cms.string("output.root"),
    debug = cms.bool(False),
    isMC = cms.bool(True),
    xsec = cms.double(9999.99),
    #mctruthmode = cms.int32(22), # for photon
    mctruthmode = cms.int32(0), # DY process 
    maxevents = cms.int32(-1), # set to -1 when running on grid. 
    pujetidparas = cms.PSet(pu_jetid)
)
