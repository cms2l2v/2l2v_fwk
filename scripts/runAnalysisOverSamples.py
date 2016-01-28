#!/usr/bin/env python
import os,sys,time
import json
import optparse
import commands
import LaunchOnCondor
import UserCode.llvv_fwk.storeTools_cff as storeTools

PROXYDIR = "~/x509_user_proxy"
DatasetFileDB = "DAS"  #DEFAULT: will use das_client.py command line interface
#DatasetFileDB = "DBS" #OPTION:  will use curl to parse https GET request on DBSserver

isLocalSample = False
nonLocalSamples = []

"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(proc,key,defaultVal=None) :
    try :
        return proc[key]
    except KeyError:
        return defaultVal


def getByLabelFromKeyword(proc,keyword,key,defaultVal=None) :
    try:
       print "found keyword option %s" % str(proc[keyword][key])
       return proc[keyword][key]
    except:
       try :
           return proc[key]
       except KeyError:
           return defaultVal


initialCommand = '';
def initProxy():
   global initialCommand
   validCertificate = True
   if(validCertificate and (not os.path.isfile(os.path.expanduser(PROXYDIR+'/x509_proxy')))):validCertificate = False
   if(validCertificate and (time.time() - os.path.getmtime(os.path.expanduser(PROXYDIR+'/x509_proxy')))>600): validCertificate = False
   # --voms cms, otherwise it does not work normally
   if(validCertificate and int(commands.getstatusoutput('(export X509_USER_PROXY='+PROXYDIR+'/x509_proxy;voms-proxy-init --voms cms --noregen;voms-proxy-info -all) | grep timeleft | tail -n 1')[1].split(':')[2])<8 ):validCertificate = False

   if(not validCertificate):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      if(not os.path.isfile(os.path.expanduser('~/.globus/mysecret.txt'))):
         if(hostname.find("iihe.ac.be")!=-1): os.system('mkdir -p '+PROXYDIR+'; voms-proxy-init --voms cms:/cms/becms  -valid 192:00 --out '+PROXYDIR+'/x509_proxy')
         else:                                os.system('mkdir -p '+PROXYDIR+'; voms-proxy-init --voms cms             -valid 192:00 --out '+PROXYDIR+'/x509_proxy')
      else:
         if(hostname.find("iihe.ac.be")!=-1): os.system('mkdir -p '+PROXYDIR+'; voms-proxy-init --voms cms:/cms/becms  -valid 192:00 --out '+PROXYDIR+'/x509_proxy -pwstdin < /home/fynu/quertenmont/.globus/mysecret.txt')
         else:                                os.system('mkdir -p '+PROXYDIR+'; voms-proxy-init --voms cms             -valid 192:00 --out '+PROXYDIR+'/x509_proxy -pwstdin < /home/fynu/quertenmont/.globus/mysecret.txt') 
   initialCommand = 'export X509_USER_PROXY='+PROXYDIR+'/x509_proxy;voms-proxy-init --voms cms --noregen; ' #no voms here, otherwise I (LQ) have issues

def getFileList(procData,DefaultNFilesPerJob):
   global nonLocalSamples
   global isLocalSample
   isLocalSample = False

   FileList = [];
   miniAODSamples = getByLabel(procData,'miniAOD','')
   isMINIAODDataset = ("/MINIAOD" in getByLabel(procData,'dset','')) or  ("amagitte" in getByLabel(procData,'dset',''))
   if(isMINIAODDataset or len(getByLabel(procData,'miniAOD',''))>0):
      instance = ""
      if(len(getByLabel(procData,'dbsURL',''))>0): instance =  "instance=prod/"+ getByLabel(procData,'dbsURL','')
      listSites = commands.getstatusoutput('das_client.py --query="site dataset='+getByLabel(procData,'dset','') + ' ' + instance + ' | grep site.name,site.dataset_fraction " --limit=0')[1]
      IsOnLocalTier=False
      MaxFraction=0;  FractionOnLocal=-1;
      for site in listSites.split('\n'):
         if(localTier==""):continue;
         try:
            MaxFraction = max(MaxFraction, float(site.split()[1].replace('%','')) )
         except:
            MaxFraction = max(MaxFraction, 0.0);
         if(localTier in site):
            FractionOnLocal = float(site.split()[1].replace('%',''));

      if(FractionOnLocal == MaxFraction):
            IsOnLocalTier=True            
            print ("Sample is found to be on the local grid tier %s (%f%%) for %s") %(localTier, FractionOnLocal, getByLabel(procData,'dset',''))

      isLocalSample = IsOnLocalTier

      if(localTier != "" and not IsOnLocalTier):
         nonLocalSamples += [getByLabel(procData,'dset','')]

      list = []
      if(IsOnLocalTier or isMINIAODDataset):
         list = []
         if(DatasetFileDB=="DAS"):
            list = commands.getstatusoutput('das_client.py --query="file dataset='+getByLabel(procData,'dset','') + ' ' + instance + '" --limit=0')[1].split()
         elif(DatasetFileDB=="DBS"):
            curlCommand="curl -ks --key $X509_USER_PROXY --cert $X509_USER_PROXY -X GET "
            dbsPath="https://cmsweb.cern.ch/dbs/prod/global/DBSReader"
            sedTheList=' | sed \"s#logical_file_name#\\nlogical_file_name#g\" | sed \"s#logical_file_name\': \'##g\" | sed \"s#\'}, {u\'##g\" | sed \"s#\'}]##g\" | grep store '
            list = commands.getstatusoutput(initialCommand + curlCommand+'"'+dbsPath+'/files?dataset='+getByLabel(procData,'dset','')+'"'+sedTheList)[1].split()

         list = [x for x in list if ".root" in x] #make sure that we only consider root files
         for i in range(0,len(list)):              
            if IsOnLocalTier:
               if  (hostname.find("iihe.ac.be")!=-1): list[i] = "dcap://maite.iihe.ac.be:/pnfs/iihe/cms/ph/sc4"+list[i]
               elif(hostname.find("ucl.ac.be" )!=-1): list[i] = "/storage/data/cms"+list[i]
               else:                                  list[i] = "root://eoscms//eos/cms"+list[i]            
            else:
               list[i] = "root://cms-xrd-global.cern.ch/"+list[i] #works worldwide
              #list[i] = "root://xrootd-cms.infn.it/"+list[i]    #optimal for EU side
              #list[i] = "root://cmsxrootd.fnal.gov/"+list[i]    #optimal for US side

      elif(len(getByLabel(procData,'miniAOD',''))>0):
         print "Processing private local sample: " + getByLabel(procData,'miniAOD','')
         list = storeTools.fillFromStore(getByLabel(procData,'miniAOD',''),0,-1,True);                  
      else:
         print "Processing an unknown type of sample (assuming it's a private local sample): " + getByLabel(procData,'miniAOD','')
         list = storeTools.fillFromStore(getByLabel(procData,'miniAOD',''),0,-1,True);

      list = storeTools.keepOnlyFilesFromGoodRun(list, os.path.expandvars(getByLabel(procData,'lumiMask','')))       
      split=getByLabel(procData,'split',-1)
      if(split>0):
         NFilesPerJob = max(1,len(list)/split)
      else:
         NFilesPerJob = DefaultNFilesPerJob

      for g in range(0, len(list), NFilesPerJob):
         groupList = ''
         for f in list[g:g+NFilesPerJob]:
            groupList += '"'+f+'",\\n';
         FileList.append(groupList)

   else:
      print "Processing a non EDM/miniAOD sample in : " + opt.indir + '/' + origdtag + '_' + str(segment) + '.root'
      for segment in range(0,split) :
         eventsFile=opt.indir + '/' + origdtag + '_' + str(segment) + '.root'
         if(eventsFile.find('/store/')==0)  : eventsFile = commands.getstatusoutput('cmsPfn ' + eventsFile)[1]
         FileList.append('"'+eventsFile+'"')
   return FileList

#configure
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-e', '--exe'        ,    dest='theExecutable'      , help='executable'                                , default='')
parser.add_option('-s', '--sub'        ,    dest='queue'              , help='batch queue OR "crab" to use crab3'        , default='8nh')
parser.add_option('-R', '--R'          ,    dest='requirementtoBatch' , help='requirement for batch queue'               , default='pool>30000')
parser.add_option('-j', '--json'       ,    dest='samplesDB'          , help='samples json file'                         , default='')
parser.add_option('-d', '--dir'        ,    dest='indir'              , help='input directory or tag in json file'       , default='aoddir')
parser.add_option('-o', '--out'        ,    dest='outdir'             , help='output directory'                          , default='')
parser.add_option('-t', '--tag'        ,    dest='onlytag'            , help='process only samples matching this tag'    , default='all')
parser.add_option('-k', '--key'        ,    dest='onlykeyword'        , help='process only samples matching this keyword', default='')
parser.add_option('-K', '--skipkey'    ,    dest='skipkeyword'        , help='skip process matching this keyword'        , default='')
parser.add_option('-p', '--pars'       ,    dest='params'             , help='extra parameters for the job'              , default='')
parser.add_option('-c', '--cfg'        ,    dest='cfg_file'           , help='base configuration file template'          , default='')
parser.add_option('-r', "--report"     ,    dest='report'             , help='If the report should be sent via email'    , default=False, action="store_true")
parser.add_option('-D', "--db"         ,    dest='db'                 , help='DB to get file list for a given dset'      , default=DatasetFileDB)
parser.add_option('-F', "--resubmit"   ,    dest='resubmit'           , help='resubmit jobs that failed'                 , default=False, action="store_true")
parser.add_option('-S', "--NFile"      ,    dest='NFile'              , help='default #Files per job (for autosplit)'    , default=6)
parser.add_option('-f', "--localnfiles",    dest='localnfiles'        , help='number of parallel jobs to run locally'    , default=8)
parser.add_option('-l', "--lfn"        ,    dest='crablfn'            , help='user defined directory for CRAB runs'      , default='')

(opt, args) = parser.parse_args()
scriptFile=os.path.expandvars('${CMSSW_BASE}/bin/${SCRAM_ARCH}/wrapLocalAnalysisRun.sh')
DatasetFileDB                      = opt.db

#define local site
localTier = ""
hostname = commands.getstatusoutput("hostname -f")[1]
if(hostname.find("ucl.ac.be")!=-1):localTier = "T2_BE_UCL"
if(hostname.find("iihe.ac.be")!=-1):localTier = "T2_BE_IIHE"
if(hostname.find("cern.ch")!=-1)  :localTier = "T2_CH_CERN"

FarmDirectory                      = opt.outdir+"/FARM"
PROXYDIR                           = FarmDirectory+"/inputs/"
initProxy()

JobName                            = opt.theExecutable
LaunchOnCondor.Jobs_RunHere        = 0
LaunchOnCondor.Jobs_Queue          = opt.queue
LaunchOnCondor.Jobs_LSFRequirement = '"'+opt.requirementtoBatch+'"'
LaunchOnCondor.Jobs_EmailReport    = opt.report
LaunchOnCondor.Jobs_LocalNJobs     = opt.localnfiles
LaunchOnCondor.Jobs_CRABLFN        = opt.crablfn
LaunchOnCondor.Jobs_ProxyDir       = FarmDirectory+"/inputs/" 

#open the file which describes the sample
jsonFile = open(opt.samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()
#run over sample
for procBlock in procList :
    #run over processes
    for proc in procBlock[1] :
        #run over items in process
        if(getByLabel(proc,'interpollation', '')!=''):continue #skip interpollated processes
        if(getByLabel(proc,'mixing'        , '')!=''):continue #skip mixed processes
        keywords = getByLabel(proc,'keys',[])
        if(opt.onlykeyword!='' and len(keywords)>0 and opt.onlykeyword not in keywords): continue #skip processes not maching the keyword
        if(opt.skipkeyword!='' and len(keywords)>0):
           skipThisProc = False
           for skip in opt.skipkeyword.split(','):
              if(skip in keywords): skipThisProc = True
           if(skipThisProc):continue #skip processes matching one of the skipkeyword

        isdata=getByLabelFromKeyword(proc,opt.onlykeyword,'isdata',False)
        isdatadriven=getByLabelFromKeyword(proc,opt.onlykeyword,'isdatadriven',False)       
        mctruthmode=getByLabelFromKeyword(proc,opt.onlykeyword,'mctruthmode',0)
        procSuffix=getByLabelFromKeyword(proc,opt.onlykeyword,'suffix' ,"")
        data = proc['data']

        for procData in data :
            LaunchOnCondor.Jobs_InitCmds = ['ulimit -c 0;'] 
 
            origdtag = getByLabel(procData,'dtag','')
            if(origdtag=='') : continue
            dtag = origdtag         

            xsec = getByLabel(procData,'xsec',-1)
            br = getByLabel(procData,'br',[])
            suffix = str(getByLabel(procData,'suffix' ,procSuffix))
            split=getByLabel(procData,'split',-1)
            if opt.onlytag!='all' and dtag.find(opt.onlytag)<0 : continue
            filt=''
            if mctruthmode!=0 : filt='_filt'+str(mctruthmode)      
            if(xsec>0 and not isdata):
                for ibr in br :  xsec = xsec*ibr


            if(opt.resubmit==False):
               FileList = ['"'+getByLabel(procData,'dset','UnknownDataset')+'"']
               LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + '_' + dtag)
               if(LaunchOnCondor.subTool!='crab'):FileList = getFileList(procData, int(opt.NFile) )
               if(not isLocalSample):LaunchOnCondor.Jobs_InitCmds      += [initialCommand]

               for s in range(0,len(FileList)):
                   #create the cfg file
                   eventsFile = FileList[s]
                   eventsFile = eventsFile.replace('?svcClass=default', '')
                   prodfilepath=opt.outdir +'/'+ dtag + suffix + '_' + str(s) + filt
               	   sedcmd = 'sed \''
                   sedcmd += 's%"@dtag"%"' + dtag +'"%;'
                   sedcmd += 's%"@input"%' + eventsFile+'%;'
            	   sedcmd += 's%@outfile%' + prodfilepath+'.root%;'
            	   sedcmd += 's%@isMC%' + str(not (isdata or isdatadriven) )+'%;'
            	   sedcmd += 's%@mctruthmode%'+str(mctruthmode)+'%;'
            	   sedcmd += 's%@xsec%'+str(xsec)+'%;'
                   sedcmd += 's%@cprime%'+str(getByLabel(procData,'cprime',-1))+'%;'
                   sedcmd += 's%@brnew%' +str(getByLabel(procData,'brnew' ,-1))+'%;'
                   sedcmd += 's%@suffix%' +suffix+'%;'
                   sedcmd += 's%@lumiMask%"' + os.path.expandvars(getByLabel(procData,'lumiMask',''))+'"%;'
              	   if(opt.params.find('@useMVA')<0) :          opt.params = '@useMVA=False ' + opt.params
                   if(opt.params.find('@weightsFile')<0) :     opt.params = '@weightsFile= ' + opt.params
                   if(opt.params.find('@evStart')<0) :         opt.params = '@evStart=0 '    + opt.params
                   if(opt.params.find('@evEnd')<0) :           opt.params = '@evEnd=-1 '     + opt.params
            	   if(opt.params.find('@saveSummaryTree')<0) : opt.params = '@saveSummaryTree=False ' + opt.params
            	   if(opt.params.find('@runSystematics')<0) :  opt.params = '@runSystematics=False '  + opt.params
                   if(opt.params.find('@jacknife')<0) :        opt.params = '@jacknife=-1 ' + opt.params
                   if(opt.params.find('@jacks')<0) :           opt.params = '@jacks=-1 '    + opt.params
                   if(opt.params.find('@trig')<0) :            opt.params = '@trig=False ' + opt.params
            	   if(len(opt.params)>0) :
                      extracfgs = opt.params.split(' ')
                      for icfg in extracfgs :
                        varopt=icfg.split('=')
                        if(len(varopt)<2) : continue
                        sedcmd += 's%' + varopt[0] + '%' + varopt[1] + '%;'
            	   sedcmd += '\''
                   cfgfile=prodfilepath + '_cfg.py'
                   os.system('cat ' + opt.cfg_file + ' | ' + sedcmd + ' > ' + cfgfile)

                   #run the job
                   if len(opt.queue)==0 :
                       os.system(opt.theExecutable + ' ' + cfgfile)
                   else:
                       if(LaunchOnCondor.subTool=='crab'):
                          LaunchOnCondor.Jobs_CRABDataset  = FileList[0]
                          LaunchOnCondor.Jobs_CRABcfgFile  = cfgfile
                          LaunchOnCondor.Jobs_CRABexe      = opt.theExecutable
                          if(commands.getstatusoutput("whoami")[1]=='vischia'):
                              LaunchOnCondor.Jobs_CRABStorageSite = 'T2_PT_NCG_Lisbon'
                          else:
                              LaunchOnCondor.Jobs_CRABStorageSite = 'T2_BE_UCL'
                          LaunchOnCondor.Jobs_CRABname     = dtag
                          LaunchOnCondor.Jobs_CRABInDBS    = getByLabel(procData,'dbsURL','global')
                          if(split>0):
                              LaunchOnCondor.Jobs_CRABUnitPerJob = 100 / split 
                          else:
                              LaunchOnCondor.Jobs_CRABUnitPerJob = int(opt.NFile)
                       LaunchOnCondor.SendCluster_Push(["BASH", str(opt.theExecutable + ' ' + cfgfile)])

               LaunchOnCondor.SendCluster_Submit()

            else:
               configList = commands.getstatusoutput('ls ' + opt.outdir +'/'+ dtag + suffix + '*_cfg.py')[1].split('\n')
               failedList = []
               for cfgfile in configList:
                  if( not os.path.isfile( cfgfile.replace('_cfg.py','.root'))):
                     failedList+= [cfgfile]

               if(len(failedList)>0):
                  LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + '_' + dtag)
                  for cfgfile in failedList:                  
                     LaunchOnCondor.SendCluster_Push(["BASH", str(opt.theExecutable + ' ' + cfgfile)])
                  LaunchOnCondor.SendCluster_Submit()


if(LaunchOnCondor.subTool=='criminal'):
    LaunchOnCondor.SendCluster_CriminalSubmit()


if(len(nonLocalSamples)>0):
   print "Some samples that you want to process are not currently available on the local tier ("+localTier+") you use."
   print "Consider transfering them via phedex to run your analysis faster:"
   for nonLocalSample in nonLocalSamples:
      print nonLocalSample

