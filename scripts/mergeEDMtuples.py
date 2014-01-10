#!/usr/bin/env python

import os,sys
import json
import optparse
import commands
from UserCode.llvv_fwk.storeTools_cff import fillFromStore, addPrefixSuffixToFileList, removeDuplicates
import LaunchOnCondor

"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal


usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--inDir'      ,    dest='inDir'              , help='input directory containing the crab output'           , default='/store/cmst3/user/querten/2013_Jul_EDMtuples')
parser.add_option('-o', '--outDir'     ,    dest='outDir'             , help='output directory where the merged file will be moved' , default='/store/user/querten/2013_Jul_EDMtuples_merged')
parser.add_option('-s', '--sub'        ,    dest='queue'              , help='batch queue'                                          , default='2nd')
parser.add_option('-R', '--R'          ,    dest='requirementtoBatch' , help='requirement for batch queue'                          , default='pool>30000')
parser.add_option('-j', '--json'       ,    dest='samplesDB'          , help='samples json file'                                    , default='')
parser.add_option('-t', '--tag'        ,    dest='onlytag'            , help='process only samples matching this tag'               , default='all')
parser.add_option('-n', '--n'          ,    dest='fperjob'            , help='input files per job'                                  , default=-1,  type=int)
parser.add_option('-D', '--duplicates' ,    dest='duplicates'         , help='clean the input directory for duplicates'             , default=False)
(opt, args) = parser.parse_args()


scriptFile=os.path.expandvars('${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapLocalAnalysisRun.sh')

jsonFile = open(opt.samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

if(opt.outDir.find('/store/')==0):
   os.system("cmsMkdir " + opt.outDir) #make sure that the output dir exist
elif(opt.outDir.find('/lustre/ncg.ingrid.pt/')==0):
   os.system("mkdir " + opt.outDir) 
elif(opt.outDir.find('srm://')==0):
   os.system('source /nfs/soft/grid/ui/setup/grid-env.sh; voms-proxy-init -voms cms -valid 192:00; mkdir -p /nfs/home/fynu/quertenmont/x509_user_proxy; cp $X509_USER_PROXY /nfs/home/fynu/quertenmont/x509_user_proxy/proxy')#all must be done in the same command to avoid environement problems

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

            if(opt.duplicates):
               dirToClean = inputdir
               if(dirToClean.find('/storage/data/cms/store/')): dirToClean = dirToClean.replace('/storage/data/cms/store/', '/storage_rw/data/cms/store/') #Hack for Louvain T2
               removeDuplicates(dirToClean);

            filenames=LaunchOnCondor.natural_sort(fillFromStore(inputdir,0,-1,False))
            nfiles=len(filenames)
            filenames=addPrefixSuffixToFileList("   '", filenames, "',")

            split=getByLabel(d,'split',1)
            NFilesToMerge = nfiles//split
            NFilesToMergeRemains = nfiles%split 
            startFile = 0
            endFile = 0 
            for segment in range(0,split) :
                startFile = endFile 
                endFile   = endFile + NFilesToMerge
                if(NFilesToMergeRemains>0):
                    endFile+=1
                    NFilesToMergeRemains-=1

                mergedFileName = d['dtag']
                if(split>1): mergedFileName+='_' + str(segment)
                mergedFileName+= '.root'
                mergedFilePath = opt.outDir + "/" + mergedFileName

                LaunchOnCondor.Jobs_RunHere        = 0
                LaunchOnCondor.Jobs_Queue          = opt.queue
                LaunchOnCondor.Jobs_LSFRequirement = '"'+opt.requirementtoBatch+'"'
                #LaunchOnCondor.ListToFile(filenames[startFile:endFile], "/tmp/InputFile_"+d['dtag']+".txt")                
                if(mergedFilePath.find('castor')>=0) :
                   LaunchOnCondor.Jobs_FinalCmds = ['pwd', 'ls -lth', 'rfcp '+mergedFileName+' ' + mergedFilePath]
                # lustre condition must be before /store/ condition (/lustre/ path contains /store/ string) 
                elif(mergedFilePath.find('/lustre/ncg.ingrid.pt/')==0):
                    LaunchOnCondor.Jobs_FinalCmds = ['pwd', 'ls -lth', 'cp '+mergedFileName+' ' + mergedFilePath]
                elif(mergedFilePath.find('/store/')==0):
                   LaunchOnCondor.Jobs_FinalCmds = ['pwd', 'ls -lth', 'cmsStageOut '+mergedFileName+' ' + mergedFilePath]
                elif(mergedFilePath.find('srm://')==0):
                   LaunchOnCondor.Jobs_FinalCmds = ['pwd', 'ls -lth', 'source /nfs/soft/grid/ui/setup/grid-env.sh', 'export X509_USER_PROXY=/nfs/home/fynu/quertenmont/x509_user_proxy/proxy','ls /$PWD/'+mergedFileName, 'for i in 5m 10m 15m 20m 25m 30m 45m 60m; do lcg-cp -v -D srmv2 -b file:/$PWD/'+mergedFileName+' '+mergedFilePath + ' && break || echo "Pause for $i" && sleep $i; done', 'mv ' + mergedFileName + ' /home/fynu/quertenmont/scratch/13_06_26_HTauTau/LoicFramework/CMSSW_5_3_11/src/UserCode/llvv_fwk/FARM_Merge/outputs/.', 'rm -rf ' + mergedFileName]
                else:
                   LaunchOnCondor.Jobs_FinalCmds = ['pwd', 'ls -lth', 'mv '+mergedFileName+' ' + os.getcwd()+"/FARM_Merge/outputs/"+mergedFileName, 'ls -lth '+opt.outDir]                
#                   LaunchOnCondor.Jobs_FinalCmds = ['pwd', 'ls -lth', 'mv '+mergedFileName+' ' + mergedFilePath, 'ls -lth '+opt.outDir]                
#                   mergedFileName = 'file:'+mergedFilePath
                LaunchOnCondor.SendCMSMergeJob("FARM_Merge", "Merge_"+d['dtag']+'_'+str(segment), filenames[startFile:endFile], "'"+mergedFileName+"'", "'keep *'")
