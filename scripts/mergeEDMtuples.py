#!/usr/bin/env python

import os,sys
import json
import optparse
import commands
from UserCode.llvv_fwk.storeTools_cff import fillFromStore, addPrefixSuffixToFileList
import LaunchOnCondor




usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--inDir'      ,    dest='inDir'              , help='input directory containing the crab output'           , default='/store/cmst3/user/querten/2013_Jul_EDMtuples')
parser.add_option('-o', '--outDir'     ,    dest='outDir'             , help='output directory where the merged file will be moved' , default='/store/user/querten/2013_Jul_EDMtuples_merged')
parser.add_option('-s', '--sub'        ,    dest='queue'              , help='batch queue'                                          , default='2nd')
parser.add_option('-R', '--R'          ,    dest='requirementtoBatch' , help='requirement for batch queue'                          , default='pool>30000')
parser.add_option('-j', '--json'       ,    dest='samplesDB'          , help='samples json file'                                    , default='')
parser.add_option('-t', '--tag'        ,    dest='onlytag'            , help='process only samples matching this tag'               , default='all')
parser.add_option('-n', '--n'          ,    dest='fperjob'            , help='input files per job'                                  , default=-1,  type=int)
(opt, args) = parser.parse_args()


scriptFile=os.path.expandvars('${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapLocalAnalysisRun.sh')

jsonFile = open(opt.samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()


os.system("cmsMkdir " + opt.outDir) #make sure that the output dir exist

#run over sample
for proc in procList :

    #run over processes
    for desc in proc[1] :
         
        #run over items in process
        data = desc['data']
        for d in data :

            #tag veto
            if(opt.onlytag!='all') :
                itag=d['dtag']
                if(itag.find(opt.onlytag)<0) : continue

            inputdir = opt.inDir+"/"+d['dtag'];
            print str(d['dtag'])+" --> "+inputdir


            filenames=fillFromStore(inputdir,0,-1,False)
            nfiles=len(filenames)
            filenames=addPrefixSuffixToFileList("   '", filenames, "',")

            LaunchOnCondor.Jobs_RunHere        = 0
            LaunchOnCondor.Jobs_Queue          = opt.queue
#            LaunchOnCondor.Jobs_LSFRequirement = '"'+opt.requirementtoBatch+'"'
            LaunchOnCondor.ListToFile(filenames, "/tmp/InputFile_"+d['dtag']+".txt")                
            LaunchOnCondor.Jobs_FinalCmds = ['ls', 'cmsStageOut out.root ' + opt.outDir+"/"+d['dtag']+".root"]
            LaunchOnCondor.SendCMSMergeJob("FARM_Merge", "Merge_"+d['dtag'], filenames, "'out.root'", "'keep *'")
