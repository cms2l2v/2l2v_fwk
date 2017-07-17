#!/usr/bin/env python
import math
import os,sys
#import json
import getopt
import commands
import ROOT
from ROOT import TFile, TGraph, TCanvas, TF1, TH1
sys.path.append('../../../scripts/')
import LaunchOnCondor

#default value

signalSuffixVec = []
OUTName         = []
LandSArgOptions = []
BIN             = []
MODEL           = []

LaunchOnCondor.Jobs_Queue='8nh'
FarmDirectory  = "FARM"
JobName        = "computeLimits"
CMSSW_BASE=os.environ.get('CMSSW_BASE')
CWD=os.getcwd()
phase=-1

###################################################
##   VALUES TO BE EDITED BY THE USE ARE BELLOW   ##
###################################################

MODELS=["SM"] #,"RsGrav","BulkGrav","Rad"] #Models which can be used are: "RsGrav", "BulkGrav", "Rad", "SM"
based_key="2l2v_datadriven_" #mcbased_" #to run limits on MC use: 2l2v_mcbased_, to use data driven obj use: 2l2v_datadriven_
jsonPath='$CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples_full2016_GGH_VBF.json' #samples_full2016_GGH_WithoutBckg.json' #samples_full2016_GGH_WithoutWZ.json' #samples_full2016_GGH_WithoutWZandZVV.json' #samples_full2016_GGH_WithoutIrriducible.json' #samples_full2016_GGH.json'
inUrl='$CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plotter_2017_06_08_forLimits_Alessio.root' #plotter_2017_05_05_forLimits.root'
BESTDISCOVERYOPTIM=True #Set to True for best discovery optimization, Set to False for best limit optimization
ASYMTOTICLIMIT=True #Set to True to compute asymptotic limits (faster) instead of toy based hybrid-new limits
BINS = ["eq0jets,geq1jets,vbf"] # list individual analysis bins to consider as well as combined bins (separated with a coma but without space)

MASS = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000]
SUBMASS = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000] 

LandSArgCommonOptions="  --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001  --subNRB"

for model in MODELS:
   for shape in ["mt_shapes"]:# --histoVBF met_shapes"]:  #here run all the shapes you want to test.  '"mt_shapes --histoVBF met_shapes"' is a very particular case since we change the shape for VBF
      for bin in BINS:
         if(model=="SM" or model=="MELA"):
            for CP in [100.0,10.0,5.0]:
               if(CP!=100.0 and not ',' in bin):continue  #only do subchannel for SM
               for BRN in [0.0]:
    
                   suffix = "_cp%4.2f_brn%4.2f" % (CP, BRN)
 
                   #Run limit for ShapeBased GG+VBF
                   #signalSuffixVec += [ suffix ]
                   #OUTName         += ["SB13TeV_SM"]
                   #LandSArgOptions += [" --histo " + shape + " --histoVBF " + shape + " --systpostfix _13TeV --shape "]  #as this combine ggF and VBF, we need to scale VBF to the SM ratio expectation
                   #BIN             += [bin]
                   #MODEL           += [model]

                   signalSuffixVec += [ suffix ]
                   OUTName         += ["SB13TeV_SM_GGF"]
                   LandSArgOptions += [" --histo " + shape + " --histoVBF " + shape + "  --systpostfix _13TeV --shape --skipQQH "]
                   BIN             += [bin]
                   MODEL          += [model]

                   #signalSuffixVec += [ suffix ]
                   #OUTName         += ["SB13TeV_SM_VBF"]
                   #LandSArgOptions += [" --histo " + shape + " --histoVBF " + shape + "  --systpostfix _13TeV --shape --skipGGH "]
                   #BIN             += [bin]
                   #MODEL           += [model]

         elif(model=="RsGrav"):
                   signalSuffixVec += [ "" ]
                   OUTName         += ["SB13TeV_RsGrav"]
                   LandSArgOptions += [" --histo " + shape + " --systpostfix _13TeV --shape"]
                   BIN             += [bin]
                   MODEL           += [model]

         elif(model=="BulkGrav"):
                   signalSuffixVec += [ "" ]
                   OUTName         += ["SB13TeV_BulkGrav"]
                   LandSArgOptions += [" --histo " + shape + " --systpostfix _13TeV --shape"]
                   BIN             += [bin]
                   MODEL           += [model]

         elif(model=="Rad"):
                   signalSuffixVec += [ "" ]
                   OUTName         += ["SB13TeV_Rad"]
                   LandSArgOptions += [" --histo " + shape + " --systpostfix _13TeV --shape"]
                   BIN             += [bin]
                   MODEL           += [model]

###################################################
##   MAIN  CODE : MODIFY AT YOUR OWN RISK        ##
###################################################

def help() :
   print ' -p phase (no default value is assigned)'
   print '\t 1 --> compute limit for all possible selection cuts  (in view of optimization)'
   print '\t 2 --> wrap-up and order the results either by best significance or best limit'
   print '\t 3 --> choose the best cuts for the selection and produce limits and control plots for this selection'
   print '\t 4 --> produce final limit plot (brazilian flag plot)'
   print '\nNote: CMSSW_BASE must be set when launching optimize.py (current values is: ' + CMSSW_BASE + ')\n' 
   
#parse the options
try:
   # retrive command line options
   shortopts  = "p:f:m:i:s:j:o:h?"
   opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
   # print help information and exit:
   print "ERROR: unknown options in argument %s" % sys.argv[1:]
   help()
   sys.exit(1)

for o,a in opts:
   if o in("-?", "-h"):
      help()
      sys.exit(1)
   elif o in('-i'): inUrl = a
   elif o in('-p'): phase = int(a)
   elif o in('-o'): CWD=a
   elif o in('-j'): jsonUrl=a
      
if(phase<0 or len(CMSSW_BASE)==0):
   help()
   sys.exit(1)


#auxiliary function
def findCutIndex(cutsH, Gcut, m):
   for i in range(1, cutsH.GetXaxis().GetNbins()+1):
      passAllCuts=True
      for y in range(1, cutsH.GetYaxis().GetNbins()+1):
         if( (cutsH.GetYaxis().GetBinLabel(y).find("<")>=0 and cutsH.GetBinContent(i,y)>(Gcut[y-1].Eval(m,0,"")+0.000001)) or (cutsH.GetYaxis().GetBinLabel(y).find(">")>=0 and cutsH.GetBinContent(i,y)<(Gcut[y-1].Eval(m,0,"")-0.000001)) or  (cutsH.GetYaxis().GetBinLabel(y).find("<")==-1 and cutsH.GetYaxis().GetBinLabel(y).find(">")==-1 and cutsH.GetBinContent(i,y)<(Gcut[y-1].Eval(m,0,"")-0.000001)) ):
            passAllCuts=False;
            break;
      if(passAllCuts==True): return i;   
   return cutsH.GetXaxis().GetNbins()+1;

def findSideMassPoint(mass):
   global MASS
   LMass=0
   RMass=9999
   for m in MASS:
      if(m<=mass and m>=LMass):LMass=m
      if(m>=mass and m<=RMass):RMass=m
   return [LMass,RMass]

#######################

#Loop over all configurations
iConf = -1
cp = 0
for signalSuffix in signalSuffixVec : 
   iConf+=1;
   if(phase<=3 and ',' in BIN[iConf]):continue #only need individual bin for these phases

   jsonUrl = jsonPath + " --key " + based_key + MODEL[iConf]
   LandSArg = LandSArgCommonOptions + ' ' + LandSArgOptions[iConf];
   if(signalSuffix != ""):LandSArg+=' --signalSufix \"' + signalSuffix +'\" '
   binSuffix = ""
   if(',' not in BIN[iConf]):binSuffix="_"+ BIN[iConf]   

   DataCardsDir='cards_'+OUTName[iConf]+signalSuffix+binSuffix

   #prepare the output
   OUT = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+binSuffix+'/'
   os.system('mkdir -p ' + OUT)

   #get the cuts
   file = ROOT.TFile(inUrl)
   cutsH = file.Get('ggH(200)_BOnly/all_optim_cut') 

   #Cp
   if "100.0" in signalSuffix:
   	cp = 100.0
   elif "10.0" in signalSuffix:
	cp = 10.0
   elif "5.0" in signalSuffix:
        cp = 5.0     
 
   ###################################################
   ##   OPTIMIZATION LOOP                           ##
   ###################################################

   if( phase == 1 ):
      print '# RUN LIMITS FOR ALL POSSIBLE CUTS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])

      FILE = open(OUT+"/LIST.txt","w")
      i = 1
      while (i<cutsH.GetXaxis().GetNbins()+1):                   
          shapeCutMin_ = 0
          shapeCutMax_ = 9999
          SCRIPT = open(OUT+'script_'+str(i)+'_'+str(shapeCutMin_)+'_'+str(shapeCutMax_)+'.sh',"w")
          SCRIPT.writelines('cd ' + CMSSW_BASE + '/src;\n')
          SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
          SCRIPT.writelines("eval `scram r -sh`;\n")
          SCRIPT.writelines('cd -;\n')     
          for j in range(0, 1): #always run 1points per jobs
             for m in MASS:
                cardsdir = 'H'+ str(m) + '_' + OUTName[iConf] + '_' + str(i);
                SCRIPT.writelines('mkdir -p ' + cardsdir+';\ncd ' + cardsdir+';\n')
                SCRIPT.writelines("computeLimit --m " + str(m) + " --in " + inUrl + " --syst " + " --index " + str(i)     + " --json " + jsonUrl + " --shapeMin " + str(shapeCutMin_) + " --shapeMax " + str(shapeCutMax_) + " " + LandSArg + " --bins " + BIN[iConf] + " ;\n")
                SCRIPT.writelines("sh combineCards.sh;\n")
                SCRIPT.writelines("combine -M Asymptotic -m " +  str(m) + " --run expected card_combined.dat > LIMIT.log;\n") #limit computation
                SCRIPT.writelines("combine -M ProfileLikelihood  -m " +  str(m) + " --significance -t -1 --expectSignal=1 card_combined.dat  > SIGN.log;\n") #apriori significance computation
                SCRIPT.writelines('tail -n 100 LIMIT.log > ' +OUT+str(m)+'_'+str(i)+'_'+str(shapeCutMin_)+'_'+str(shapeCutMax_)+'.log;\n')
                SCRIPT.writelines('tail -n 100 SIGN.log >> ' +OUT+str(m)+'_'+str(i)+'_'+str(shapeCutMin_)+'_'+str(shapeCutMax_)+'.log;\n')
                SCRIPT.writelines('cat LIMIT.log; cat SIGN.log\n')
                SCRIPT.writelines('cd ..;\n')
                #SCRIPT.writelines('mv ' + cardsdir + ' ' + OUT + '/.\n')
                SCRIPT.writelines('rm -rd ' + cardsdir+';\n')            
          SCRIPT.close()
          LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_'+str(i)+'_'+str(shapeCutMin_)+'_'+str(shapeCutMax_)+'.sh'])
          i = i+1#increment the cut index
      FILE.close()
      LaunchOnCondor.SendCluster_Submit()
      
   ###################################################
   ##   WRAPPING UP RESULTS                         ##
   ###################################################
   elif(phase == 2):
      print '# SCANNING ALL SETS OF CUTS  for ' + DataCardsDir + ' (you may want to go for a coffee...)#\n'
      fileName = OUT + "/OPTIM"+signalSuffix
      FILE = open(fileName+".txt","w")
      for m in MASS:
         print 'Starting mass ' + str(m)
         FILE.writelines("------------------------------------------------------------------------------------\n")
         BestLimit = []
         fileList = commands.getstatusoutput("find " + OUT +" -name " + str(m)+"_*.log")[1].split();           
         for f in fileList:
            try:
               value = -1.0;
               if(BESTDISCOVERYOPTIM==True):
                  exp = commands.getstatusoutput("cat " + f + " | grep \"Significance:\"")[1]; 
                  if(len(exp)<=0):continue
                  value = exp.split()[1]
               else:
                  exp = commands.getstatusoutput("cat " + f + " | grep \"Expected 50.0%\"")[1];  
                  if(len(exp)<=0):continue
                  value = exp.split()[4]

               if(value=='matches'):continue             
               if(float(value)<=0.0):continue
               if(BESTDISCOVERYOPTIM and float(value)>1000.0):continue #too high for a significance --> something is going wrong here
               f = f.replace(".log","")
               fields = f.split('_')
               N = len(fields)
               index = fields[N-3] 
               Cuts = ''
               for c in range(1, cutsH.GetYaxis().GetNbins()+1): 
                  Cuts += str(cutsH.GetBinContent(int(index),c)).rjust(7) + " ("+str(cutsH.GetYaxis().GetBinLabel(c))+")   "
               BestLimit.append("mH="+str(m)+ " --> Limit/Significance=" + ('%010.6f' % float(value)) + "  Index: " + str(index)   + "  Cuts: " + Cuts + "   CutsOnShape: " + str(fields[N-2]).rjust(5) + " " + str(fields[N-1]).rjust(5))
            except:
               print "File %s does not contain a valid limit" % f

         #sort the limits for this mass
         BestLimit.sort(reverse=BESTDISCOVERYOPTIM)
         for s in BestLimit:
            FILE.writelines(s+"\n")

      #all done
      FILE.close()
      print("file "+fileName+".txt is written: it contains all selection points ordered by exp limit")

   ###################################################
   ##   CHOSE BEST SELECTION CUTS and SAVE IT       ##
   ###################################################

   elif(phase == 3 ):
      LaunchOnCondor.Jobs_RunHere        = 1
      print '# CHOOSE BEST SELECTION CUTS  for ' + DataCardsDir + '#\n'
      Gcut  = []
      for c in range(1, cutsH.GetYaxis().GetNbins()+1):
         Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
         Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
         Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

      fileName = OUT+"/OPTIM"+signalSuffix
      fileName+=".txt"
         
      mi=0
      for m in MASS:
         #if you want to display more than 3 options edit -m3 field
         cut_lines=commands.getstatusoutput("cat " + fileName + " | grep 'mH="+str(m)+"' -m20")[1].split('\n')
         if(len(cut_lines)<1):continue  #make sure that we get some lines
         if(len(cut_lines[0])<1):continue # make sure the line is not empty
         print 'mH='+str(m)+'\tOption \tR \tLimits and Cuts' 
         ictr=1
         for c in cut_lines: 
            cplit = c.split()
            print '\t #'+ str(ictr) + '\t' + cplit[2] + '\tIndex=' +  cplit[4] + '\t' + str(cplit[6:len(cplit)])
            ictr+=1
         print "Which option you want to keep?"
         #opt = int(raw_input(">"))-1   # remove prompting
         opt = 0 #always take the best cut

         #save cut chosen
         cutString = ''
         for c in range(1, cutsH.GetYaxis().GetNbins()+1):
            cutString += str(cutsH.GetYaxis().GetBinLabel(c)) + cut_lines[opt].split()[6+(c-1)*2] + '\t'
            Gcut[c-1].SetPoint(mi, m, float(cut_lines[opt].split()[6+(c-1)*2]) );
         cutString += cut_lines[opt].split()[6+(cutsH.GetYaxis().GetNbins()-1)*2 + 3 ] + '<shape<' + cut_lines[opt].split()[6+(cutsH.GetYaxis().GetNbins()-1)*2 + 4] 
         print cutString
         Gcut[cutsH.GetYaxis().GetNbins()+0].SetPoint(mi, m, float(cut_lines[opt].split()[6+(cutsH.GetYaxis().GetNbins()-1)*2 + 3 ]) );
         Gcut[cutsH.GetYaxis().GetNbins()+1].SetPoint(mi, m, float(cut_lines[opt].split()[6+(cutsH.GetYaxis().GetNbins()-1)*2 + 4 ]) );
         mi+=1

      #run limits for the cuts chosen (for intermediate masses use spline interpolation)
      print("carefully check that the cuts you have choosen are well listed below:\n")
      for m in SUBMASS:
           index = findCutIndex(cutsH, Gcut, m);
           Cuts = ''
           for c in range(1, cutsH.GetYaxis().GetNbins()+1):
              Cuts += str(cutsH.GetBinContent(int(index),c)).rjust(7) + " ("+str(cutsH.GetYaxis().GetBinLabel(c))+")   "

           print "M=%04i : Index=% 5i --> Cuts: %s"  % (m, index, Cuts)

      while True:
           #ans = raw_input('Use this fit and compute final limits? (y or n)\n')
           ans='y'
           if(ans=='y' or ans == 'Y'): break;
           else:			       sys.exit(0);           
      print 'YES'

      listcuts = open(OUT+'cuts.txt',"w")
      #LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:
           index = findCutIndex(cutsH, Gcut, m);
           listcuts.writelines(str(m)+' '+str(index)+' ');
           for c in range(1, cutsH.GetYaxis().GetNbins()+3):
              listcuts.writelines(str(Gcut[c-1].Eval(m,0,""))+' ');
           listcuts.writelines('\n');
      listcuts.close();


   elif(phase == 4 ):
      LaunchOnCondor.Jobs_RunHere        = 0
      print '# FINAL COMBINED LIMITS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+binSuffix+OUTName[iConf])
      for m in SUBMASS:

           SideMasses = findSideMassPoint(m)
           indexString = ' '
           indexLString = ' '
           indexRString = ' '
           for bin in BIN[iConf].split(',') :
               Gcut  = []
               for c in range(1, cutsH.GetYaxis().GetNbins()+1):
                 Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
               Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

               INbinSuffix = "_" + bin 
               IN = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+INbinSuffix+'/'
               try:
                  listcuts = open(IN+'cuts.txt',"r")
                  mi=0
                  for line in listcuts :
                     vals=line.split(' ')
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
   #                     Gcut[c-1].SetPoint(mi, float(vals[0]), float(vals[c+1]));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);
                  listcuts.close();          
               except:
                  mi=0
                  for mtmp in SUBMASS:
                     for c in range(1, cutsH.GetYaxis().GetNbins()+3):
                        #FIXME FORCE INDEX TO BE 17 (Met>125GeV)
                        Gcut[c-1].SetPoint(mi, 17, float(125));
                     mi+=1
                  for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);

               #add comma to index string if it is not empty
               if(indexString!=' '):
                  indexString+=','
                  if(not (SideMasses[0]==SideMasses[1])):
                     indexLString+=','
                     indexRString+=','

               #find the cut index for the current mass point
               indexString += str(findCutIndex(cutsH, Gcut, m));
               if(not (SideMasses[0]==SideMasses[1])):
                  indexLString = str(findCutIndex(cutsH, Gcut, SideMasses[0]));
                  indexRString = str(findCutIndex(cutsH, Gcut, SideMasses[1]));

           #print indexString

           cutStr = " "
           SideMassesArgs = ""
           if(not (SideMasses[0]==SideMasses[1])):
               SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + indexLString +  " --indexR " + indexRString + " "

           SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
           SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
           SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc6_amd64_gcc491")+";\n")
           SCRIPT.writelines("eval `scram r -sh`;\n")
           SCRIPT.writelines('cd -;\n')

           cardsdir=DataCardsDir+"/"+('%04.0f' % float(m));
           SCRIPT.writelines('if [[ $HOSTNAME =~ "iihe" ]]; then cd $TMPDIR; fi\n')
	   SCRIPT.writelines('mkdir '+DataCardsDir+'\n')
	   SCRIPT.writelines('mkdir '+cardsdir+'; cd '+cardsdir+'\n')
           SCRIPT.writelines('mkdir -p out;\ncd out;\n')
	   #SCRIPT.writelines("computeLimit --m " + str(m) + " --in " + inUrl + " " + " --index " + indexString + " --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +" ;\n")
           SCRIPT.writelines("computeLimit --m " + str(m) + " --in " + inUrl + " " + "--syst --index " + indexString + " --bins " + BIN[iConf] + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +" ;\n")
           SCRIPT.writelines("sh combineCards.sh;\n"); 
	   SCRIPT.writelines("text2workspace.py card_combined.dat -o workspace.root -P UserCode.llvv_fwk.HiggsWidth:higgswidth --PO verbose --PO \'is2l2nu\' --PO m=\'" + str(m) + "\' --PO w=\'" + str(cp) + "\' \n")
           #compute pvalue
           SCRIPT.writelines("combine -M ProfileLikelihood --signif --pvalue -m " +  str(m) + "  workspace.root > COMB.log;\n")
	   SCRIPT.writelines("combine -M MaxLikelihoodFit workspace.root  \n")
	   SCRIPT.writelines("python " + CMSSW_BASE + "/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py  mlfit.root -g Nuisance_CrossCheck.root \n")

	   ##fvbf=0 for GGH and 1 for VBF
           ### THIS IS FOR Asymptotic fit
           if(ASYMTOTICLIMIT==True):
              SCRIPT.writelines("combine -M Asymptotic -m " +  str(m) + " workspace.root >>  COMB.log;\n") 

           ### THIS is for toy (hybridNew) fit
           else:
              SCRIPT.writelines("combine -M Asymptotic -m " +  str(m) + " workspace.root > COMB.log;\n") #first run assymptotic limit to get quickly the range of interest
              SCRIPT.writelines("rm higgsCombineTest.Asymptotic*.root;\n")
              SCRIPT.writelines("RMIN=`cat COMB.log | grep 'Expected  2.5%' | awk '{print $5;}'`;\n") #get the low edge 2sigma band from the assymptotic --> will be used to know where to put points
              SCRIPT.writelines("RMAX=`cat COMB.log | grep 'Expected 97.5%' | awk '{print $5;}'`;\n") #get the high edge 2sigma band from the assymptotic --> will be used to know where to put points
              SCRIPT.writelines('echo "expected limit from the assymptotic in the 2sigma range [$RMIN, $RMAX]";\n')
              SCRIPT.writelines('RMIN=$(echo "$RMIN*0.5" | bc);\n'); #DIVIDE RMIN   BY 2 to make sure we are considering large space enough
              SCRIPT.writelines('RMAX=$(echo "$RMAX*3.0" | bc);\n'); #MULTIPLY RMAX BY 3 to make sure we are considering large space enough
              SCRIPT.writelines('echo "for the hybridNew, consider r to be in the range [$RMIN, $RMAX]";\n')
              SCRIPT.writelines("makeGridUsingCrab.py card_combined.dat $RMIN $RMAX -n 40 -m "+str(m)+" -o grid ;\n")
              SCRIPT.writelines("rm grid.root;\n")
              SCRIPT.writelines("sh grid.sh 1 16 &> /dev/null;\n")
              SCRIPT.writelines("rm higgsCombinegrid.HybridNew.*;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+";\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.025;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.160;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.500;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.840;\n")
              SCRIPT.writelines("combine workspace.root -M HybridNew --grid=grid.root -m "+str(m)+" --expectedFromGrid 0.975;\n")
              SCRIPT.writelines("hadd -f higgsCombineTest.HybridNewMerged.mH"+str(m)+".root  higgsCombineTest.HybridNew.mH"+str(m)+"*.root;\n")
              SCRIPT.writelines("rm higgsCombineTest.HybridNew.mH"+str(m)+"*.root;\n")

           SCRIPT.writelines('mkdir -p ' + CWD+'/'+cardsdir+';\n')
           SCRIPT.writelines('mv * ' + CWD+'/'+cardsdir+'/.;\n')
           SCRIPT.writelines('cd ..;\n\n') 
           SCRIPT.close()
           #os.system('sh ' + OUT+'script_mass_'+str(m)+'.sh ')  #uncomment this line to launch interactively (this may take a lot of time)
           LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_mass_'+str(m)+'.sh'])
      LaunchOnCondor.SendCluster_Submit()


   ######################################################################

   elif(phase == 5 ):
      print '# FINAL PLOT for ' + DataCardsDir + '#\n'
      os.system("hadd -f "+DataCardsDir+"/PValueTree.root "+DataCardsDir+"/*/higgsCombineTest.ProfileLikelihood.*.root > /dev/null")

      #THIS IS FOR ASYMPTOTIC
      if(ASYMTOTICLIMIT==True):
         os.system("hadd -f "+DataCardsDir+"/LimitTree.root "+DataCardsDir+"/*/higgsCombineTest.Asymptotic.*.root > /dev/null")
      #THIS IS FOR HYBRIDNEW
      else:
         os.system("hadd -f "+DataCardsDir+"/LimitTree.root "+DataCardsDir+"/*/higgsCombineTest.HybridNewMerged.*.root > /dev/null")

      os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Stength_\",\""+DataCardsDir+"/LimitTree.root\",\"\", false, true, 13 , 35914.143 )'")

   ######################################################################


if(phase>5):
      help()

