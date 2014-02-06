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
shapeName='svfit_shapes'
inUrl='$CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau/plotter.root'
CWD=os.getcwd()
phase=-1
jsonUrl='$CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau/jsonForLucia.json'
CMSSW_BASE=os.environ.get('CMSSW_BASE')


MASS = [90,100,110,120,125,130,140,150,160]
SUBMASS = [90,100,110,120,125,130,140,150,160]

LandSArgCommonOptions=" "
signalSuffixVec = []
OUTName         = []
LandSArgOptions = []

signalSuffixVec += [""];
OUTName         += ["SB8TeV_comb"]
LandSArgOptions += [" --bins _elmu,_elha,_muha,_haha       --systpostfix _8TeV --shape "]

signalSuffixVec += [""];
OUTName         += ["SB8TeV_elmu"]
LandSArgOptions += [" --bins _elmu       --systpostfix _8TeV --shape "]

signalSuffixVec += [""];
OUTName         += ["SB8TeV_elha"]
LandSArgOptions += [" --bins _elha       --systpostfix _8TeV --shape "]

signalSuffixVec += [""];
OUTName         += ["SB8TeV_muha"]
LandSArgOptions += [" --bins _muha       --systpostfix _8TeV --shape "]

signalSuffixVec += [""];
OUTName         += ["SB8TeV_haha"]
LandSArgOptions += [" --bins _haha       --systpostfix _8TeV --shape "]



FarmDirectory                      = "FARM"
JobName                            = "computeLimits"
LaunchOnCondor.Jobs_RunHere        = 1
#LaunchOnCondor.Jobs_Queue          = opt.queue
#LaunchOnCondor.Jobs_LSFRequirement = '"'+opt.requirementtoBatch+'"'
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)



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
def findCutIndex(cutsH, cut1, cut2, cut3):
   for i in range(1, cutsH.GetXaxis().GetNbins()):
      if(cutsH.GetBinContent(i,1)<cut1-1):continue;
#      if(cutsH.GetBinContent(i,2)<cut2-5):continue;
#      if(cutsH.GetBinContent(i,3)<cut3-5):continue;
      return i;   
   return cutsH.GetXaxis().GetNbins();

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
   cutsH = file.Get('VV/optim_cut') 

   ######################################################################

   if( phase == 1 ):
      print '# RUN LIMITS FOR ALL POSSIBLE CUTS  for ' + DataCardsDir + '#\n'
      commandToRun = []

      FILE = open(OUT+"/LIST.txt","w")
      for i in range(1, cutsH.GetXaxis().GetNbins()):
#         FILE.writelines("index="+str(i).rjust(5) + " --> cut>" + str(cuts1.GetBinContent(i)).rjust(5) + " " + str(mtmin_).rjust(5) + '<mt<'+str(mtmax_).rjust(5) + "\n")
          #create wrappper script for each set of cuts ans submit it
          mtmin_ = 0
          mtmax_ = 9999
          SCRIPT = open(OUT+'script_'+str(i)+'_'+str(mtmin_)+'_'+str(mtmax_)+'.sh',"w")
          SCRIPT.writelines('echo "TESTING SELECTION : ' + str(i).rjust(5) + ' --> met>' + str(cutsH.GetBinContent(i,1)).rjust(5) + ' ' + str(mtmin_) + '<mt<'+str(mtmax_)+'";\n')
          SCRIPT.writelines('cd ' + CMSSW_BASE + '/src;\n')
          SCRIPT.writelines("export SCRAM_ARCH="+os.getenv("SCRAM_ARCH","slc5_amd64_gcc434")+";\n")
          SCRIPT.writelines("eval `scram r -sh`;\n")
          SCRIPT.writelines('cd /tmp/;\n')
          for m in MASS:
             if(m!=125):continue
             shapeBasedOpt=''
             cardsdir = 'H'+ str(m);
             SCRIPT.writelines('mkdir -p ' + cardsdir+';\ncd ' + cardsdir+';\n')
             SCRIPT.writelines("computeLimit --m " + str(m) + " --histo " + shapeName + " --in " + inUrl + " --syst " + shapeBasedOpt + " --index " + str(i)     + " --json " + jsonUrl +" --fast " + " --rebin 2 --shapeMin " + str(mtmin_) + " --shapeMax " + str(mtmax_) + " " + LandSArg + " ;\n")
             SCRIPT.writelines("sh combineCards.sh;\n")
             SCRIPT.writelines("combine -M Asymptotic -m " +  str(m) + " --run expected card_combined.dat > COMB.log;\n")
             SCRIPT.writelines('tail -n 100 COMB.log > ' +OUT+str(m)+'_'+str(i)+'_'+str(mtmin_)+'_'+str(mtmax_)+'.log;\n')
             SCRIPT.writelines('cd ..;\n\n')
             SCRIPT.close()
             commandToRun.append("bsub -o /dev/null -G u_zh -q 8nh -J optim"+str(i)+'_'+str(mtmin_)+'_'+str(mtmax_) + " 'sh " + OUT+"script_"+str(i)+'_'+str(mtmin_)+'_'+str(mtmax_)+".sh &> "+OUT+"script_"+str(i)+'_'+str(mtmin_)+'_'+str(mtmax_)+".log'")
      FILE.close()

      for c in commandToRun:
         print(c)
         os.system(c);

         
   ######################################################################
   elif(phase == 2):
      print '# SCANNING ALL SETS OF CUTS  for ' + DataCardsDir + ' (you may want to go for a coffee...)#\n'
      
      fileName = OUT + "/OPTIM" + signalSuffix
      FILE = open(fileName+".txt","w")
      for m in MASS:
         if(m!=125):continue
         print 'Starting mass ' + str(m)
         FILE.writelines("------------------------------------------------------------------------------------\n")
         BestLimit = []
         fileList = commands.getstatusoutput("ls " + OUT + str(m)+"_*.log")[1].split();           
         for f in fileList:
            exp = commands.getstatusoutput("cat " + f + " | grep \"Expected 50.0%\"")[1];
            if(len(exp)<=0):continue
            median = exp.split()[4]
            if(float(median)<=0.0):continue
            f = f.replace(".log","")
            fields = f.split('_')
            N = len(fields)
            index = fields[N-3] 
            BestLimit.append("mH="+str(m)+ " --> " + ('%07.3f' % float(median)) + " " + str(index) + "  " + str(cutsH.GetBinContent(int(index),1)).rjust(5) + " " + str(cutsH.GetBinContent(int(index),2)).rjust(5) + "  " + str(fields[N-2]).rjust(5) + " " + str(fields[N-1]).rjust(5))

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
      Gmet  = ROOT.TGraph(len(SUBMASS));
      Gtmin = ROOT.TGraph(len(SUBMASS));
      Gtmax = ROOT.TGraph(len(SUBMASS));


      if(cutList=='') :
         fileName = OUT+"/OPTIM"+signalSuffix
         fileName+=".txt"
         
         mi=0
         for m in MASS:
            #if you want to display more than 3 options edit -m3 field
            cut_lines=commands.getstatusoutput("cat " + fileName + " | grep 'mH="+str(m)+"' -m20")[1].split('\n')
            print 'mH='+str(m)+'\tOption \tR \tmin MET\tMT range' 
            ictr=1
            for c in cut_lines:
               print '\t #'+ str(ictr) + '\t' + c.split()[2] + '\t' + c.split()[4] + '\t(' + c.split()[5] + '-' + c.split()[6]+ ')'
               ictr+=1
            print "Which option you want to keep?"
            opt = int(raw_input(">"))-1

            #save cut chosen
            metCut=float(cut_lines[opt].split()[4])
            mtMinCut=float(cut_lines[opt].split()[5])
            mtMaxCut=float(cut_lines[opt].split()[6])
            Gmet .SetPoint(mi, m, metCut);
            Gtmin.SetPoint(mi, m, mtMinCut);
            Gtmax.SetPoint(mi, m, mtMaxCut);
            mi+=1

         #display cuts chosen
         c1 = ROOT.TCanvas("c1", "c1",900,300);
         ROOT.gROOT.SetStyle('Plain')
         ROOT.gStyle.SetOptStat(False);

         c1 = ROOT.TCanvas("c1", "c1",900,300);
         c1.Divide(3);
         c1.cd(1);
         Gmet.SetMarkerStyle(20);
         Gmet.SetTitle("MET");
         Gmet.Draw("APC");
         Gmet.GetXaxis().SetTitle("m_{H} (GeV/c^{2})");
         Gmet.GetYaxis().SetTitle("met cut");

         c1.cd(2);
         Gtmin.SetMarkerStyle(20);
         Gtmin.SetTitle("MT min");
         Gtmin.Draw("APC");
         Gtmin.GetXaxis().SetTitle("m_{H} (GeV/c^{2})");
         Gtmin.GetYaxis().SetTitle("mt_{min} cut");

         c1.cd(3);
         Gtmax.SetMarkerStyle(20);
         Gtmax.SetTitle("MT max");
         Gtmax.Draw("APC");
         Gtmax.GetXaxis().SetTitle("m_{H} (GeV/c^{2})");
         Gtmax.GetYaxis().SetTitle("mt_{max} cut");
         c1.cd(0);
         c1.Update();
         c1.SaveAs("OptimizedCuts.png")

         Gmet.Set(mi);
         Gtmin.Set(mi);
         Gtmax.Set(mi);

         #run limits for the cuts chosen (for intermediate masses use spline interpolation)
         for m in SUBMASS:
              index = findCutIndex(cutsH, Gmet.Eval(m,0,""), Gtmin.Eval(m,0,""),  Gtmax.Eval(m,0,""));
      #        print("mH="+str(m).rjust(3)+ " met>"+str(cutsH.GetBinContent(index)).rjust(5) + " " + str(cuts2.GetBinContent(index)).rjust(5) + "<mt<"+str(cuts3.GetBinContent(index)).rjust(5) )

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
            Gmet .SetPoint(mi, float(vals[0]), float(vals[1]));
            Gtmin.SetPoint(mi, float(vals[0]), float(vals[2]));
            Gtmax.SetPoint(mi, float(vals[0]), 0);#float(vals[3]));
            mi+=1
         f.close()
         Gmet.Set(mi);
         Gtmin.Set(mi);
         Gtmax.Set(mi);



      list = open(OUT+'list.txt',"w")
      listcuts = open(OUT+'cuts.txt',"w")
      for m in SUBMASS:
           index = 1;#findCutIndex(cutsH, Gmet.Eval(m,0,""), Gtmin.Eval(m,0,""), Gtmax.Eval(m,0,""));
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
              Lindex = findCutIndex(cutsH, Gmet.Eval(SideMasses[0],0,""), Gtmin.Eval(SideMasses[0],0,""),  Gtmax.Eval(SideMasses[0],0,""));
              Rindex = findCutIndex(cutsH, Gmet.Eval(SideMasses[1],0,""), Gtmin.Eval(SideMasses[1],0,""),  Gtmax.Eval(SideMasses[1],0,""));
              #print "cutIndex for sideBand are " + str(Lindex) + " and " + str(Rindex) 
              SideMassesArgs += "--mL " + str(SideMasses[0]) + " --mR " + str(SideMasses[1]) + " --indexL " + str(Lindex) +  " --indexR " + str(Rindex) + " "


           cutStr = " --rebin 2 --shapeMin " + str(Gtmin.Eval(m,0,"")) +" --shapeMax " + str(Gtmax.Eval(m,0,""));
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

#           os.system("bsub -o /dev/null -G u_zh -q 8nh 'sh " + OUT+"script_mass_"+str(m)+".sh'")
           list.writelines('H'+str(m)+'_'+str(index)+'\n');
           listcuts.writelines(str(m)+' ' + str(Gmet.Eval(m,0,"")) + ' ' + str(Gtmin.Eval(m,0,""))+' '+ str(Gtmax.Eval(m,0,"")) +'\n');
      list.close();
      listcuts.close();

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

LaunchOnCondor.SendCluster_Submit()

