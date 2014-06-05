#!/usr/bin/env python

import os,sys
#import json
import getopt
import commands
import ROOT
from ROOT import TFile, TGraph, TCanvas, TF1, TH1
sys.path.append('../../../scripts/')
import LaunchOnCondor

#default values
shapeName='svfit_shapes '
inUrl='$CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau/plotter.root'
CWD=os.getcwd()
phase=-1
jsonUrl='$CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau/jsonForLucia.json'
CMSSW_BASE=os.environ.get('CMSSW_BASE')


MASS = [15,30,70]
SUBMASS = [15,30,70]

LandSArgCommonOptions=" --signalRescale 1000 "
signalSuffixVec = []
OUTName         = []
LandSArgOptions = []

#signalSuffixVec += [""];
#OUTName         += ["SB8TeV_comb"]
#LandSArgOptions += [" --bins _elmu,_elha,_muha,_haha       --systpostfix _8TeV"]# --shape "]

signalSuffixVec += [""];
OUTName         += ["SB8TeV_elmu"]
LandSArgOptions += [" --bins _OSelmu       --systpostfix _8TeV"]# --shape "]

signalSuffixVec += [""];
OUTName         += ["SB8TeV_elha"]
LandSArgOptions += [" --bins _OSelha       --systpostfix _8TeV"]# --shape "]

signalSuffixVec += [""];
OUTName         += ["SB8TeV_muha"]
LandSArgOptions += [" --bins _OSmuha       --systpostfix _8TeV"]# --shape "]

signalSuffixVec += [""];
OUTName         += ["SB8TeV_haha"]
LandSArgOptions += [" --bins _OShaha       --systpostfix _8TeV"]# --shape "]

FarmDirectory                      = "FARM"
JobName                            = "computeLimits"
LaunchOnCondor.Jobs_RunHere        = 1
#LaunchOnCondor.Jobs_Queue          = opt.queue
#LaunchOnCondor.Jobs_LSFRequirement = '"'+opt.requirementtoBatch+'"'



cutList='' 
def help() :
   print '\n\033[92m optimize.py \033[0m \n'
   print ' -p phase (no default value is assigned)'
   print '\t 1 --> submit landS jobs for all selection point'
   print '\t 2 --> check the logs to find the optimal selection point'
   print '\t      from the ouptut of the logs you can search for the optimal points yourself ;)'
   print '\t      and edit phase3 of this script with your optimal points (note: change this to be given as input)'
   print '\t 3 --> you are prompted to choose the best cuts for the selection: the limits will be comptued for the list of cuts'
   print '\t       if -f LIST.txt is given the LIST of cuts passed will be used instead'
   print '\t 4 --> once all the final limit jobs have been run, use this phase to build the brazilian flag plot'
   print '\t 0 --> cut and count based analysis'
   print '\t 1 --> shape based analysis'
   print ' -s shapename (default='+shapeName+')'
   print ' -i inputfile (default='+inUrl+')'
   print ' -o output (default='+CWD+')'
   print ' -j jsonfile (default='+jsonUrl+')'
   print '\nUsage example: \033[92m python optimize.py -m 0 -i ~/work/plotter.root -o ~/work/ -p 1 \033[0m'
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
   elif o in('-s'): shapeName=a
   elif o in('-f'): cutList=a
      
if(phase<0 or len(CMSSW_BASE)==0):
   help()
   sys.exit(1)


#auxiliary function
def findCutIndex(cutsH, Gcut, m):
   for i in range(1, cutsH.GetXaxis().GetNbins()+1):
      for y in range(1, cutsH.GetYaxis().GetNbins()+1):
         if(cutsH.GetBinContent(i,y)<Gcut[y-1].Eval(m,0,"")-0.000001):continue;
      return i;   
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
for signalSuffix in signalSuffixVec :
   iConf+=1;

   LandSArg = LandSArgCommonOptions + ' ' + LandSArgOptions[iConf];
   if(signalSuffix != ""):LandSArg+=' --signalSufix \"' + signalSuffix +'\" '
   DataCardsDir='cards_'+OUTName[iConf]+signalSuffix


   #prepare the output
   OUT = CWD+'/JOBS/'+OUTName[iConf]+signalSuffix+'/'
   os.system('mkdir -p ' + OUT)

   #get the cuts
   file = ROOT.TFile(inUrl)
   cutsH = file.Get('data/optim_cut') 

   ######################################################################

   if( phase == 1 ):
      print '# RUN LIMITS FOR ALL POSSIBLE CUTS  for ' + DataCardsDir + '#\n'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+OUTName[iConf])

      FILE = open(OUT+"/LIST.txt","w")
      i = 1
      while (i<cutsH.GetXaxis().GetNbins()+1):                   
#            FILE.writelines("index="+str(i).rjust(5) + " --> cut>" + str(cuts1.GetBinContent(i)).rjust(5) + " " + str(svfitmin_).rjust(5) + '<mt<'+str(svfitmax_).rjust(5) + "\n")
             #create wrappper script for each set of cuts ans submit it
          svfitmin_ = 0
          svfitmax_ = 9999
          SCRIPT = open(OUT+'script_'+str(i)+'_'+str(svfitmin_)+'_'+str(svfitmax_)+'.sh',"w")
          SCRIPT.writelines('echo "TESTING SELECTION : ' + str(i).rjust(5) + ' --> met>' + str(cutsH.GetBinContent(i,1)).rjust(5) + ' ' + str(svfitmin_) + '<mt<'+str(svfitmax_)+'";\n')
          SCRIPT.writelines('cd ' + CMSSW_BASE + '/src;\n')
          SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc5_amd64_gcc434")+";\n")
          SCRIPT.writelines("eval `scram r -sh`;\n")
          SCRIPT.writelines('cd /tmp/;\n')
          for j in range(0, 10): #always run 100points per jobs
             #if(not i%100==0):continue         #just run 1 job over 100 for debugging
             for m in MASS:
                shapeBasedOpt=''
                cardsdir = 'H'+ str(m);
                SCRIPT.writelines('mkdir -p ' + cardsdir+';\ncd ' + cardsdir+';\n')
                SCRIPT.writelines("computeLimit --m " + str(m) + " --histo " + shapeName + " --in " + inUrl + " --syst " + shapeBasedOpt + " --index " + str(i)     + " --json " + jsonUrl +" --fast " + " --shapeMin " + str(svfitmin_) + " --shapeMax " + str(svfitmax_) + " " + LandSArg + " ;\n")
                SCRIPT.writelines("sh combineCards.sh;\n")
                SCRIPT.writelines("combine -M Asymptotic -m " +  str(m) + " --run expected card_combined.dat > COMB.log;\n")
                SCRIPT.writelines('tail -n 100 COMB.log > ' +OUT+str(m)+'_'+str(i)+'_'+str(svfitmin_)+'_'+str(svfitmax_)+'.log;\n')
                SCRIPT.writelines('cd ..;\n\n')
             i = i+1#increment the cut index
          SCRIPT.close()
          LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_'+str(i)+'_'+str(svfitmin_)+'_'+str(svfitmax_)+'.sh &> '+OUT+'script_'+str(i)+'_'+str(svfitmin_)+'_'+str(svfitmax_)+'.log'])
      FILE.close()
      LaunchOnCondor.SendCluster_Submit()

         
   ######################################################################
   elif(phase == 2):
      print '# SCANNING ALL SETS OF CUTS  for ' + DataCardsDir + ' (you may want to go for a coffee...)#\n'
      
      fileName = OUT + "/OPTIM" + signalSuffix
      FILE = open(fileName+".txt","w")
      for m in MASS:
         print 'Starting mass ' + str(m)
         FILE.writelines("------------------------------------------------------------------------------------\n")
         BestLimit = []
         fileList = commands.getstatusoutput("ls " + OUT + str(m)+"_*.log")[1].split();           
         for f in fileList:
            exp = commands.getstatusoutput("cat " + f + " | grep \"Expected 50.0%\"")[1];
            if(len(exp)<=0):continue
            median = exp.split()[4]
#            print 'median is ' + median
	    if(median=='matches'):continue
            if(float(median)<=0.0):continue
            f = f.replace(".log","")
            fields = f.split('_')
            N = len(fields)
            index = fields[N-3] 
            Cuts = ''
            for c in range(1, cutsH.GetYaxis().GetNbins()+1): 
               Cuts += str(cutsH.GetBinContent(int(index),c)).rjust(7) + " ("+str(cutsH.GetYaxis().GetBinLabel(c))+")   "
            BestLimit.append("mH="+str(m)+ " --> Limit=" + ('%010.6f' % float(median)) + "  Cuts: " + Cuts + "   CutsOnShape: " + str(fields[N-2]).rjust(5) + " " + str(fields[N-1]).rjust(5))

         #sort the limits for this mass
         BestLimit.sort()
         for s in BestLimit:
            FILE.writelines(s+"\n")

      #all done
      FILE.close()
      print("file "+fileName+".txt is written: it contains all selection points ordered by exp limit")

   ######################################################################

   elif(phase == 3 ):
      print '# FINAL LIMITS  for ' + DataCardsDir + '#\n'
      Gcut  = []
      for c in range(1, cutsH.GetYaxis().GetNbins()+1):
         Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #add a graph for each cut
      Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMin
      Gcut.extend([ROOT.TGraph(len(SUBMASS))]) #also add a graph for shapeMax

      if(cutList=='') :
         fileName = OUT+"/OPTIM"+signalSuffix
         fileName+=".txt"
         
         mi=0
         for m in MASS:
            #if you want to display more than 3 options edit -m3 field
            cut_lines=commands.getstatusoutput("cat " + fileName + " | grep 'mH="+str(m)+"' -m20")[1].split('\n')
            print 'mH='+str(m)+'\tOption \tR \tLimits and Cuts' 
            ictr=1
            for c in cut_lines: 
               cplit = c.split()
               #print '\t #'+ str(ictr) + '\t' + c
               print '\t #'+ str(ictr) + '\t' + cplit[2] + '\t' + str(cplit[4:len(cplit)])
               ictr+=1
            print "Which option you want to keep?"
            opt = int(raw_input(">"))-1

            #save cut chosen
            cutString = ''
            for c in range(1, cutsH.GetYaxis().GetNbins()+1):
               cutString += str(cutsH.GetYaxis().GetBinLabel(c)) + cut_lines[opt].split()[4+(c-1)*2] + '\t'
               Gcut[c-1].SetPoint(mi, m, float(cut_lines[opt].split()[4+(c-1)*2]) );
            cutString += cut_lines[opt].split()[4+(cutsH.GetYaxis().GetNbins()-1)*2 + 3 ] + '<shape<' + cut_lines[opt].split()[4+(cutsH.GetYaxis().GetNbins()-1)*2 + 4] 
            print cutString
            Gcut[cutsH.GetYaxis().GetNbins()+0].SetPoint(mi, m, float(cut_lines[opt].split()[4+(cutsH.GetYaxis().GetNbins()-1)*2 + 3 ]) );
            Gcut[cutsH.GetYaxis().GetNbins()+1].SetPoint(mi, m, float(cut_lines[opt].split()[4+(cutsH.GetYaxis().GetNbins()-1)*2 + 4 ]) );
            mi+=1

         #display cuts chosen
#         c1 = ROOT.TCanvas("c1", "c1",300*(cutsH.GetYaxis().GetNbins()+2),300);
#         ROOT.gROOT.SetStyle('Plain')
#         ROOT.gStyle.SetOptStat(False);
#         c1.Divide(cutsH.GetYaxis().GetNbins()+2);
         for c in range(1, cutsH.GetYaxis().GetNbins()+3):
#            c1.cd(c);
            Gcut[c-1].SetMarkerStyle(20);
#            Gcut[c-1].Draw("APC");
            Gcut[c-1].GetXaxis().SetTitle("m_{A} (GeV/c^{2})");
            if(c<=cutsH.GetYaxis().GetNbins()+1):
               Gcut[c-1].GetYaxis().SetTitle(str(cutsH.GetYaxis().GetBinLabel(c)));
               Gcut[c-1].Set(mi);
            else:
               Gcut[c-1].GetYaxis().SetTitle("Shape");
#         c1.cd(0);
#         c1.Update();
#         c1.SaveAs("OptimizedCuts.png")

         #run limits for the cuts chosen (for intermediate masses use spline interpolation)
         for m in SUBMASS:
              index = findCutIndex(cutsH, Gcut, m);

         while True:
              ans = raw_input('Use this fit and compute final limits? (y or n)\n')
              if(ans=='y' or ans == 'Y'): break;
              else:			    sys.exit(0);           
         print 'YES'

      else :
         mi=0
         f= open(cutList,'r')
         for line in f :
            vals=line.split(' ')
            for c in range(1, cutsH.GetYaxis().GetNbins()+3):
               Gcut[c-1].SetPoint(mi, float(vals[0]), float(vals[c]));
            mi+=1
         f.close()
         for c in range(1, cutsH.GetYaxis().GetNbins()+3): Gcut[c-1].Set(mi);

      list = open(OUT+'list.txt',"w")
      listcuts = open(OUT+'cuts.txt',"w")
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName + "_"+signalSuffix+OUTName[iConf])
      for m in SUBMASS:
           index = findCutIndex(cutsH, Gcut, m);
           SCRIPT = open(OUT+'/script_mass_'+str(m)+'.sh',"w")
           SCRIPT.writelines('cd ' + CMSSW_BASE + ';\n')
           SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc5_amd64_gcc434")+";\n")
           SCRIPT.writelines("eval `scram r -sh`;\n")
           SCRIPT.writelines('cd ' + CWD + ';\n')
           shapeBasedOpt=''

           SideMassesArgs = ' '
           SideMasses = findSideMassPoint(m)
           if(not (SideMasses[0]==SideMasses[1])):
              #print "Side Mass for mass " + str(m) + " are " + str(SideMasses[0]) + " and " + str(SideMasses[1])
              Lindex = findCutIndex(cutsH, Gcut, SideMasses[0]);
              Rindex = findCutIndex(cutsH, Gcut, SideMasses[1]);
              SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + str(Lindex) +  " --indexR " + str(Rindex) + " "

           cutStr = " --rebin 2 --shapeMin " + str(Gcut[cutsH.GetYaxis().GetNbins()].Eval(m,0,"")) +" --shapeMax " + str(Gcut[cutsH.GetYaxis().GetNbins()+1].Eval(m,0,""));
           if(OUTName[iConf].find("SB")>=0):cutStr = " --rebin 8 "; #not cuts on the shape itself
           cutStr = " "

           cardsdir=DataCardsDir+"/"+('%04.0f' % float(m));
           SCRIPT.writelines('mkdir -p ' + cardsdir+';\ncd ' + cardsdir+';\n')
           SCRIPT.writelines("computeLimit --m " + str(m) + " --histo " + shapeName + " --in " + inUrl + " " + " --syst " + shapeBasedOpt + " --index " + str(index) + " --json " + jsonUrl + " " + SideMassesArgs + " " + LandSArg + cutStr  +" ;\n")
           SCRIPT.writelines("sh combineCards.sh;\n")
           SCRIPT.writelines("combine -M Asymptotic -m " +  str(m) + "  card_combined.dat > COMB.log;\n") 
           SCRIPT.writelines("combine -M MaxLikelihoodFit -m " +  str(m) + " --saveNormalizations card_combined.dat;\n")
           SCRIPT.writelines("extractFitNormalization.py mlfit.root hzz2l2v_"+str(m)+"_?TeV.root > fit.txt;\n")
           SCRIPT.writelines('cd ..;\n\n') 
           SCRIPT.close()

	   LaunchOnCondor.SendCluster_Push(["BASH", 'sh ' + OUT+'script_mass_'+str(m)+'.sh'])

           list.writelines('H'+str(m)+'_'+str(index)+'\n');
           listcuts.writelines(str(m)+' ');
           for c in range(1, cutsH.GetYaxis().GetNbins()+3):
              listcuts.writelines(str(Gcut[c-1].Eval(m,0,""))+' ');
           listcuts.writelines('\n');
      list.close();
      listcuts.close();
      LaunchOnCondor.SendCluster_Submit()

   ######################################################################

   elif(phase == 4 ):
      print '# FINAL PLOT for ' + DataCardsDir + '#\n'
      os.system("hadd -f "+DataCardsDir+"/LimitTree.root "+DataCardsDir+"/*/higgsCombineTest.Asymptotic.*.root > /dev/null")
      if(LandSArg.find('skip')<0):
         os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Stength_\",\""+DataCardsDir+"/LimitTree.root\",\"\",  true, false, 8 , 19.7 )'")
      else:
         os.system("getXSec "+DataCardsDir+"/XSecs.txt "+DataCardsDir+"/*/Efficiency.tex")
         os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/Stength_\",\""+DataCardsDir+"/LimitTree.root\",\""+DataCardsDir+"/XSecs.txt\",  true, false, 8 , 19.7 )'")
         os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"/XSec_\",\""+DataCardsDir+"/LimitTree.root\",\""+DataCardsDir+"/XSecs.txt\",  false, false, 8 , 19.7 )'")
#         os.system("root -l -b -q plotLimit.C+'(\""+DataCardsDir+"\",\""+DataCardsDir+"/LimitTree.root\",\"\",  true, 8 , 19.6 )'")
   ######################################################################

   else:
      help()


