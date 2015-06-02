#!/usr/bin/python

import os
import sys
import commands
from ROOT import TH1F,TGraph,TFile

type=sys.argv[1]

if(type=='VBF') :
    cafdir='/store/caf/user/rgerosa/powheg_lhe_production/qqH_CPS_corrected'
else :
    cafdir='/store/caf/user/rgerosa/powheg_lhe_production/ggH_CPS'

    
baseUrl='http://tier2.ihepa.ufl.edu/~tcheng/HighMass/ggH_shape_HTO/Tables/'
fOut = TFile.Open(type+"_LineShapes.root","RECREATE")

#
def getShapeFromFile(fname,col,hname) :
    hmasses = commands.getstatusoutput("cat %s | awk '{print $1}'"%fname)[1].strip().split()
    shape   = commands.getstatusoutput("cat %s | awk '{print $%d}'"%(fname,col))[1].strip().split()
    shapeGr=TGraph()
    shapeGr.SetName(hname)
    for i in xrange(0,len(hmasses)): shapeGr.SetPoint(shapeGr.GetN(),float(hmasses[i]),float(shape[i]))
    return shapeGr


allFiles=commands.getstatusoutput("cmsLs %s | grep -ir hmass | awk '{print $5}'"%cafdir)[1].strip().split()


for lhe_file in allFiles:

    #get the LHE
    print lhe_file
    mass=lhe_file.split('RG_hmass_')[1].split('_')[0]
    print 'Saving distributions and weights for mass: %s'%mass

    outdir=fOut.mkdir('H'+mass)
    outdir.cd()

    #read the generated Higgs mass lineshape
    print '    GEN shapes from %s'%lhe_file
    commands.getstatusoutput("cmsStage -f %s /tmp/tmp.lhe.tar.gz"%lhe_file)
    commands.getstatusoutput("cd /tmp && gtar -xzvf /tmp/tmp.lhe.tar.gz")
    grep_output = commands.getstatusoutput("grep -ir 25\ \ \ \ \  /tmp/pwgevents.lhe  | awk '{print $11}'")
    hmasses=grep_output[1].strip().split()
    genH=TH1F('gen','gen',1000,0,4000)
    for m in hmasses: genH.Fill(float(m))
    norm=genH.Integral()
    if(norm>0):genH.Scale(1./norm)
    genH.Write()
    commands.getstatusoutput("rm ./tmp.lhe")
    
    #get the CPS shape from Tonguang's files
    sfname='mZZ_Higgs%s_8TeV_Lineshape+Interference.txt'%mass
    url=baseUrl+'/'+sfname
    print '    CPS shapes from %s'%url
    commands.getstatusoutput('curl -O %s'%url)
    errorOut=commands.getstatusoutput('grep -ir Not\ found %s'%sfname)[1]

    #histos to build: the first entry is the column in the txt file
    histos=[[3,'cps'],[4,'cps_up'],[5,'cps_down'],[6,'nominal'],[7,'up'],[8,'down']]
    #no interference for VBF just take the nominal CPS again for now
    if(type=='VBF'):
        histos[3][0]=3
        histos[4][0]=4
        histos[5][0]=5

    #get shape + interference histograms
    for h in histos :

        #if not found clone the baseline so that the weights will be 1 - it typically happens for mH<400 GeV
        shape=genH.Clone(h[1])        
        shape.SetName(h[1]+"_shape")
        if(len(errorOut)<=0): 
            shapeGr=getShapeFromFile(sfname,h[0],h[1])
            shape.Reset('ICE')
            for ibin in xrange(1,shape.GetXaxis().GetNbins()+1):
                x=shape.GetXaxis().GetBinCenter(ibin)
                y=max(shapeGr.Eval(x),0.)                
                shape.SetBinContent(ibin,y)
            norm=shape.Integral()
            if norm>0 : shape.Scale(1./norm)
            shapeGr.Write(h[1]+'_orig')
        h.append(shape)
        shape.Write()
        
        refShape=genH
        if(h[1].find('cps')<0): refShape=histos[0][2]
        ratio=shape.Clone(h[1]+'_ratio')
        ratio.Divide(refShape)
        shapeWgtGr=TGraph(ratio)
        shapeWgtGr.SetName(h[1])
        shapeWgtGr.Write()
    commands.getstatusoutput("rm %s"%sfname)
    
    print ' ... all done: <M_H> Gen:'+str(genH.GetMean())+' CPS:'+str(histos[0][2].GetMean())

fOut.Close()
    
