#!/usr/bin/env python
import os,sys,time
import json
import optparse
import commands
import LaunchOnCondor
import UserCode.llvv_fwk.storeTools_cff as storeTools
import sqlite3
import pwd

PROXYDIR = "~/x509_user_proxy"
DatasetFileDB = "DAS"  #DEFAULT: will use das_client.py command line interface
#DatasetFileDB = "DBS" #OPTION:  will use curl to parse https GET request on DBSserver

cachedQueryDB = sqlite3.connect(os.path.expandvars('${CMSSW_BASE}/src/UserCode/llvv_fwk/data/das_query_cache.db') )
cachedQueryDBcursor = cachedQueryDB.cursor()
cachedQueryDBcursor.execute("""CREATE TABLE IF NOT EXISTS queries(id INTEGER PRIMARY KEY AUTOINCREMENT UNIQUE, date TEXT, query TEXT, result TEXT)""")
cachedQueryDB.commit()
#cachedQueryDBcursor.close()
#cachedQueryDB.close()


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

def DASQuery(query):
   cachedQueryDB = sqlite3.connect(os.path.expandvars('${CMSSW_BASE}/src/UserCode/llvv_fwk/data/das_query_cache.db') )
   cachedQueryDBcursor = cachedQueryDB.cursor()

   cachedQueryDBcursor.execute("""SELECT result FROM queries WHERE query=?""", (query,) )
   fetched = cachedQueryDBcursor.fetchall()
   if(len(fetched)>1):print("WARNING:  More than one query result found in cached DAS query DB, return first one")
   if(len(fetched)>=1):
      cachedQueryDB.commit()
      cachedQueryDB.close()
      return fetched[0][0]

   #get the result from DAS and cache it for future usage (only if there was no error with DAS)
   outputs = commands.getstatusoutput('/cvmfs/cms.cern.ch/common/das_client --query="' + query + '" --limit=0')
   result = outputs[1]
   if(outputs[0]==0):cachedQueryDBcursor.execute("""INSERT INTO queries(date, query, result) VALUES(?, ?, ?)""", (time.strftime('%Y-%m-%d %H:%M:%S'), query, result))
   cachedQueryDB.commit()  #commit changes to the DB

   return result 

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
         os.system('mkdir -p '+PROXYDIR+'; voms-proxy-init --voms cms             -valid 192:00 --out '+PROXYDIR+'/x509_proxy')
      else:
         os.system('mkdir -p '+PROXYDIR+'; voms-proxy-init --voms cms             -valid 192:00 --out '+PROXYDIR+'/x509_proxy -pwstdin < /home/fynu/quertenmont/.globus/mysecret.txt') 
   initialCommand = 'export X509_USER_PROXY='+PROXYDIR+'/x509_proxy;voms-proxy-init --voms cms --noregen; ' #no voms here, otherwise I (LQ) have issues
   os.system('cp '+PROXYDIR+'/x509_proxy /tmp/x509up_u%s' % pwd.getpwuid( os.getuid() ).pw_uid)

def getFileList(procData,DefaultNFilesPerJob):
   global nonLocalSamples
   global isLocalSample
   isLocalSample = False

   FileList = [];
   nonMiniAODSamples = []
   miniAODSamples = getByLabel(procData,'miniAOD',[])
   dsetSamples = getByLabel(procData,'dset',[])
   for s in dsetSamples: 
      if("/MINIAOD" in s): miniAODSamples+=[s]
      else: nonMiniAODSamples+=[s]

   for sample in miniAODSamples:
      
      instance = ""
      if(len(getByLabel(procData,'dbsURL',''))>0): instance =  "instance=prod/"+ getByLabel(procData,'dbsURL','')
      listSites = DASQuery('site dataset='+sample + ' ' + instance + ' | grep site.name,site.replica_fraction')
      IsOnLocalTier=False
      MaxFraction=0;  FractionOnLocal=-1;
      for site in listSites.split('\n'):
         if(localTier==""):continue;
         try:
            MaxFraction = max(MaxFraction, float(site.split()[2].replace('%','').replace('"','')) )
         except:
            MaxFraction = max(MaxFraction, 0.0);
         if(localTier in site):
            FractionOnLocal = float(site.split()[2].replace('%','').replace('"',''));

      if(FractionOnLocal == MaxFraction):
            IsOnLocalTier=True            
            print ("Sample is found to be on the local grid tier %s (%f%%) for %s") %(localTier, FractionOnLocal, sample)

      isLocalSample = IsOnLocalTier

      if(localTier != "" and not IsOnLocalTier):
         nonLocalSamples += [sample]

      list = [] 
      if(IsOnLocalTier or "/MINIAOD" in sample):
         list = []
         if(DatasetFileDB=="DAS"):
            list = DASQuery('file dataset='+sample + ' ' + instance).split()
         elif(DatasetFileDB=="DBS"):
            curlCommand="curl -ks --key $X509_USER_PROXY --cert $X509_USER_PROXY -X GET "
            dbsPath="https://cmsweb.cern.ch/dbs/prod/global/DBSReader"
            sedTheList=' | sed \"s#logical_file_name#\\nlogical_file_name#g\" | sed \"s#logical_file_name\': \'##g\" | sed \"s#\'}, {u\'##g\" | sed \"s#\'}]##g\" | grep store '
            list = commands.getstatusoutput(initialCommand + curlCommand+'"'+dbsPath+'/files?dataset='+sample+'"'+sedTheList)[1].split()

         list = [x for x in list if ".root" in x] #make sure that we only consider root files
         for i in range(0,len(list)):              
            if IsOnLocalTier:
               if  (hostname.find("iihe.ac.be")!=-1): list[i] = "dcap://maite.iihe.ac.be/pnfs/iihe/cms/ph/sc4"+list[i]
               elif(hostname.find("ucl.ac.be" )!=-1): list[i] = "/storage/data/cms"+list[i]
               else:                                  list[i] = "root://eoscms//eos/cms"+list[i]            
            else:
               list[i] = "root://cms-xrd-global.cern.ch/"+list[i] #works worldwide
              #list[i] = "root://xrootd-cms.infn.it/"+list[i]    #optimal for EU side
              #list[i] = "root://cmsxrootd.fnal.gov/"+list[i]    #optimal for US side

      else:
         print "Processing private local sample: " + sample 
         list = storeTools.fillFromStore(sample,0,-1,True);                  

      list = storeTools.keepOnlyFilesFromGoodRun(list, os.path.expandvars(getByLabel(procData,'lumiMask','')))       
      split=getByLabel(procData,'split',-1)
      if(split>0):
         NFilesPerJob = max(1,len(list)/split)
      else:
         NFilesPerJob = DefaultNFilesPerJob
         if(hostname.find("iihe.ac.be")!=-1): NFilesPerJob = DefaultNFilesPerJob #at iihe we don't want to have more than the defaultNFilesPerJob
         #elif((len(list)/NFilesPerJob)>100):NFilesPerJob=len(list)/100  #make sure the number of jobs isn't too big
         else: NFilesPerJob = DefaultNFilesPerJob #Hot fix to test removing the limit on the NFilesPerJob
            
      for g in range(0, len(list), NFilesPerJob):
         groupList = ''
         for f in list[g:g+NFilesPerJob]:
            groupList += '"'+f+'",\\n';
         FileList.append(groupList)

   for sample in nonMiniAODSamples:
      split=getByLabel(procData,'split',-1)
      for segment in range(0,split) :
         print "Processing a non EDM/miniAOD sample in : " + opt.indir + '/' + origdtag + '_' + str(segment) + '.root'
         eventsFile=opt.indir + '/' + origdtag + '_' + str(segment) + '.root'
         if(eventsFile.find('/store/')==0)  : eventsFile = commands.getstatusoutput('cmsPfn ' + eventsFile)[1]
         FileList.append('"'+eventsFile+'"')
   return FileList


def CacheInputs(FileList):
   NewList = ""
   CopyCommand = ""
   DeleteCommand = ""
   List = FileList.replace('"','').replace(',','').split('\\n')
   for l in List:
      if(len(l)<2):continue
      if("IIHE" in localTier): CopyCommand += "dccp " + l + " " + l[l.rfind('/')+1:] + ";\n"
      else: CopyCommand += "cp " + l + " " + l[l.rfind('/')+1:] + ";\n"
      DeleteCommand += "rm -f " + l[l.rfind('/')+1:] + ";\n"
      NewList += '"' +  l[l.rfind('/')+1:] + '",\\n'
   return (NewList, CopyCommand, DeleteCommand)

   


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
if(commands.getstatusoutput("hostname -f")[1].find("iihe.ac.be")!=-1): parser.add_option('-S', "--NFile"      ,    dest='NFile'              , help='default #Files per job (for autosplit)'    , default=6)
else: parser.add_option('-S', "--NFile"      ,    dest='NFile'              , help='default #Files per job (for autosplit)'    , default=8)
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
PROXYDIR                           = FarmDirectory+"/inputs"
initProxy()
doCacheInputs                      = False

if("IIHE" in localTier):
    print "kill big-submission and sleep"
    doCacheInputs =True  #False
    LaunchOnCondor.KillProcess("big-submission")
    LaunchOnCondor.KillProcess("sleep")


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
        if(getByLabel(proc,'nosample'      , '')!=''):continue #skip processes which do not have a real sample
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
        resonance=getByLabelFromKeyword(proc,opt.onlykeyword,'resonance',1);
        data = proc['data']

        for procData in data : 
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
               FileList = getByLabel(procData,'dset',['"UnknownDataset"'])
               LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + '_' + dtag+filt)
               if(LaunchOnCondor.subTool!='crab'):FileList = getFileList(procData, int(opt.NFile) )

               for s in range(0,len(FileList)):
                   LaunchOnCondor.Jobs_FinalCmds= []
                   LaunchOnCondor.Jobs_InitCmds = ['ulimit -c 0;'] 
                   if(not isLocalSample):LaunchOnCondor.Jobs_InitCmds      += [initialCommand]

                   #create the cfg file
                   eventsFile = FileList[s]
                   eventsFile = eventsFile.replace('?svcClass=default', '')
                   if(doCacheInputs and isLocalSample):
                      result = CacheInputs(eventsFile)
                      eventsFile = result[0]
                      if("IIHE" in localTier): LaunchOnCondor.Jobs_InitCmds.append('if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi;\n')
                      LaunchOnCondor.Jobs_InitCmds.append(result[1])
                      LaunchOnCondor.Jobs_FinalCmds.append(result[2])

                   prodfilepath=opt.outdir +'/'+ dtag + suffix + '_' + str(s) + filt
               	   sedcmd = 'sed \''
                   sedcmd += 's%"@dtag"%"' + dtag +'"%;'
                   sedcmd += 's%"@input"%' + eventsFile+'%;'
            	   sedcmd += 's%@outfile%' + prodfilepath+'.root%;'
            	   sedcmd += 's%@isMC%' + str(not (isdata or isdatadriven) )+'%;'
            	   sedcmd += 's%@mctruthmode%'+str(mctruthmode)+'%;'
                   sedcmd += 's%@resonance%'+str(resonance)+'%;'
            	   sedcmd += 's%@xsec%'+str(xsec)+'%;'
                   sedcmd += 's%@cprime%'+str(getByLabel(procData,'cprime',-1))+'%;'
                   sedcmd += 's%@brnew%' +str(getByLabel(procData,'brnew' ,-1))+'%;'
                   sedcmd += 's%@suffix%' +suffix+'%;'
                   sedcmd += 's%@lumiMask%"' + os.path.expandvars(getByLabel(procData,'lumiMask',''))+'"%;'
              	   if(opt.params.find('@useMVA')<0) :          opt.params = '@useMVA=False ' + opt.params
                   if(opt.params.find('@weightsFile')<0) :     opt.params = '@weightsFile= ' + opt.params
                   if(opt.params.find('@puWeightsFile')<0) :     opt.params = '@puWeightsFile= ' + opt.params
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
                          LaunchOnCondor.Jobs_CRABname     = dtag + '_' + str(s)
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

#if("IIHE" in localTier):
#    os.system("big-submission "+FarmDirectory+"/inputs/big.cmd") #now done in the submit.sh script
