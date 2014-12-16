import FWCore.ParameterSet.Config as cms

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
    #inputFiles =  PSet.process.source.fileNames
    fnames = list(PSet.process.source.fileNames)

    #fnames = [ 'root://cms-xrd-global.cern.ch/%s' %f for f in fnames]
    # import sys 
    # for f in fnames:
    #     lfn_to_pfn(f)
    #     sys.exit()
    fnames = [ lfn_to_pfn(f) for f in fnames]    
    inputFiles =  cms.untracked.vstring(fnames)                
except:
    print "not able to import" 
    pass

print inputFiles

process.source = cms.Source("PoolSource", fileNames = inputFiles)

runProcess = cms.PSet(
    #input = cms.untracked.vstring("file:input.root"),
    input = inputFiles, 
    #outdir = cms.string("results"),
    output = cms.string("output.root"),
    debug = cms.bool(True),
    isMC = cms.bool(True),
    xsec = cms.double(9999.99),
    mctruthmode = cms.int32(22)
)
