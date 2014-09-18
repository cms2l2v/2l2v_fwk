#!/usr/bin/env python
import os,sys
import json
import optparse
import commands
import LaunchOnCondor
import UserCode.llvv_fwk.storeTools_cff as storeTools

"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal


#configure
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-e', '--exe'        ,    dest='theExecutable'      , help='excecutable'                           , default='')
parser.add_option('-s', '--sub'        ,    dest='queue'              , help='batch queue'                           , default='')
parser.add_option('-R', '--R'          ,    dest='requirementtoBatch' , help='requirement for batch queue'           , default='pool>30000')
parser.add_option('-j', '--json'       ,    dest='samplesDB'          , help='samples json file'                     , default='')
parser.add_option('-d', '--dir'        ,    dest='indir'              , help='input directory or tag in json file'   , default='aoddir')
parser.add_option('-o', '--out'        ,    dest='outdir'             , help='output directory'                      , default='')
parser.add_option('-t', '--tag'        ,    dest='onlytag'            , help='process only samples matching this tag', default='all')
parser.add_option('-n', '--n'          ,    dest='fperjob'            , help='input files per job'                   , default=-1,  type=int)
parser.add_option('-p', '--pars'       ,    dest='params'             , help='extra parameters for the job'          , default='')
parser.add_option('-c', '--cfg'        ,    dest='cfg_file'           , help='base configuration file template'      , default='')
parser.add_option('-r', "--report"     ,    dest='report'             , help='If the report should be sent via email', default=False, action="store_true")
(opt, args) = parser.parse_args()
scriptFile=os.path.expandvars('${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapLocalAnalysisRun.sh')

split=1
segment=0
                                        
#open the file which describes the sample
jsonFile = open(opt.samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

FarmDirectory                      = opt.outdir+"/FARM"
JobName                            = opt.theExecutable
LaunchOnCondor.Jobs_RunHere        = 1
LaunchOnCondor.Jobs_Queue          = opt.queue
LaunchOnCondor.Jobs_LSFRequirement = '"'+opt.requirementtoBatch+'"'
LaunchOnCondor.Jobs_EmailReport    = opt.report
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)


#define local site
hostname = os.getenv("HOSTNAME", "");
localTier = ""
if(hostname.find("cern.ch")>0) : localTier = "CERN"



kInitDone = False;
initialCommand = '';
#run over sample
for proc in procList :

    #run over processes
    for desc in proc[1] :

        #run over items in process
        isdata=getByLabel(desc,'isdata',False)
        mctruthmode=getByLabel(desc,'mctruthmode',0)

        data = desc['data']
        for d in data :
            origdtag = getByLabel(d,'dtag','')
            dtag = origdtag
            xsec = getByLabel(d,'xsec',-1)
            br = getByLabel(d,'br',[])
            suffix = str(getByLabel(d,'suffix' ,""))
            if opt.onlytag!='all' and dtag.find(opt.onlytag)<0 : continue
            if mctruthmode!=0 : dtag+='_filt'+str(mctruthmode)      
                                
            if(xsec>0 and not isdata) :
                for ibr in br :  xsec = xsec*ibr
            split=getByLabel(d,'split',1)

            FileList = [];
            miniAODSamples = getByLabel(d,'miniAOD','')
            if(("/MINIAOD" in getByLabel(d,'dset','')) or len(getByLabel(d,'miniAOD',''))>0):
               listSites = commands.getstatusoutput('das_client.py --query="site dataset='+getByLabel(d,'dset','') + '" --limit=0')[1]

               list = []
               if(localTier in listSites and "CERN" in localTier):
                  list = commands.getstatusoutput('das_client.py --query="file dataset='+getByLabel(d,'dset','') + '" --limit=0')[1].split()
                  for i in range(0,len(list)): list[i] = "root://eoscms//eos/cms"+list[i]
               elif(len(getByLabel(d,'miniAOD',''))>0):
                  list = storeTools.fillFromStore(getByLabel(d,'miniAOD',''),0,-1,True);                  
               elif("/MINIAODSIM" in getByLabel(d,'dset','')):

                  if(not kInitDone):
                     print "You are going to run on a sample over grid using the AAA protocol, it is therefore needed to initialize your grid certificate"
                     os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init -voms cms -valid 192:00 --out ~/x509_user_proxy/proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain
                     initialCommand = 'export X509_USER_PROXY=~/x509_user_proxy/proxy;voms-proxy-init --noregen;'
                     kInitDone = True


                  print("Use das_client.py to list files from : " + getByLabel(d,'dset','') )
                  list = commands.getstatusoutput('das_client.py --query="file dataset='+getByLabel(d,'dset','') + '" --limit=0')[1].split()
                  for i in range(0,len(list)): list[i] = "root://cms-xrd-global.cern.ch/"+list[i]
               else:
                  list = storeTools.fillFromStore(getByLabel(d,'miniAOD',''),0,-1,True);

               ngroup = len(list)/split
               groupList = ''
               i=0;
               while(i <len(list) ):
                  groupList += '"'+list[i]+'",\\n';
                  if(i>0 and i%ngroup==0):
                     FileList.append(groupList)
                     groupList=''
                  i = i+1;                                      
            else:
 	       for segment in range(0,split) :
                  if(split==1): 
                     eventsFile=opt.indir + '/' + origdtag + '.root'
                  else:
                     eventsFile=opt.indir + '/' + origdtag + '_' + str(segment) + '.root'

                  if(eventsFile.find('/store/')==0)  : eventsFile = commands.getstatusoutput('cmsPfn ' + eventsFile)[1]
                  FileList.append('"'+eventsFile+'"')

            for s in range(0,len(FileList)):
                #create the cfg file
                eventsFile = FileList[s]
                eventsFile = eventsFile.replace('?svcClass=default', '')
            	sedcmd = 'sed \'s%"@input"%' + eventsFile +'%;s%@outdir%' + opt.outdir +'%;s%@isMC%' + str(not isdata) + '%;s%@mctruthmode%'+str(mctruthmode)+'%;s%@xsec%'+str(xsec)+'%;'
                sedcmd += 's%@cprime%'+str(getByLabel(d,'cprime',-1))+'%;'
                sedcmd += 's%@brnew%' +str(getByLabel(d,'brnew' ,-1))+'%;'
                sedcmd += 's%@suffix%' +suffix+'%;'
            	if(opt.params.find('@useMVA')<0) :          opt.params = '@useMVA=False ' + opt.params
                if(opt.params.find('@weightsFile')<0) :     opt.params = '@weightsFile= ' + opt.params
                if(opt.params.find('@evStart')<0) :         opt.params = '@evStart=0 '    + opt.params
                if(opt.params.find('@evEnd')<0) :           opt.params = '@evEnd=-1 '     + opt.params
            	if(opt.params.find('@saveSummaryTree')<0) : opt.params = '@saveSummaryTree=False ' + opt.params
            	if(opt.params.find('@runSystematics')<0) :  opt.params = '@runSystematics=False '  + opt.params
                if(opt.params.find('@jacknife')<0) :        opt.params = '@jacknife=-1 ' + opt.params
                if(opt.params.find('@jacks')<0) :           opt.params = '@jacks=-1 '    + opt.params
            	if(len(opt.params)>0) :
                    extracfgs = opt.params.split(' ')
                    for icfg in extracfgs :
                        varopt=icfg.split('=')
                        if(len(varopt)<2) : continue
                        sedcmd += 's%' + varopt[0] + '%' + varopt[1] + '%;'
            	sedcmd += '\''
		if(len(FileList)==1): 
                    cfgfile=opt.outdir +'/'+ dtag + suffix + '_cfg.py'
		else:
                    cfgfile=opt.outdir +'/'+ dtag + suffix + '_' + str(s) + '_cfg.py'
                os.system('cat ' + opt.cfg_file + ' | ' + sedcmd + ' > ' + cfgfile)

                #run the job
                if len(opt.queue)==0 :
                    os.system(opt.theExecutable + ' ' + cfgfile)
                else :
                    #old version
                    #localParams='-exe=%s -cfg=%s'%(opt.theExecutable,cfgfile)
                    #batchCommand='submit2batch.sh -q%s -R\"%s\" -J%s%d %s %s'%(opt.queue,opt.requirementtoBatch,d['dtag'],segment,scriptFile,localParams)
                    #os.system(batchCommand)
                    LaunchOnCondor.SendCluster_Push(["BASH", initialCommand + str(opt.theExecutable + ' ' + cfgfile)])

LaunchOnCondor.SendCluster_Submit()
