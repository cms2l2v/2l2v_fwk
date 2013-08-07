#!/usr/bin/env python

import os,sys
import getopt
from math import sqrt
import commands
import ROOT
from array import array

if(len(sys.argv)<2) :

    samples=[
        ["q2down", "(Q/2)^{2}", "/store/cmst3/user/psilva/VBFNLO/eejj/8TeV/ScaleDown"],
        ["q2up", "(2Q)^{2}","/store/cmst3/user/psilva/VBFNLO/mmjj/8TeV/ScaleUp"],
        ["nominal", "nominal","/store/cmst3/user/psilva/VBFNLO/mmjj/8TeV/"]
        ]
    histos=[]

    for s in samples:
        name=s[0]
        title=s[1]
        url=s[2]

        files=commands.getstatusoutput("cmsLs "+url+" | grep -ir events_ | grep -ir lhe | awk '{print $5}'" )[1].split('\n')

        print 'Running on %d files for %s'%(len(files),title)


        ptbins = [0,5,10,20,30,40,50,60,70,80,90,100,150,200,300,500,1000]
        dilptH=ROOT.TH1D(name,title+';p_{T}(ll) [GeV]; Events;',16,array('d',ptbins))
        dilptH.Sumw2()

        for f in files:

            os.system("cmsStage -f " + f + " tmp.lhe")
            os.system("cat tmp.lhe | grep -ir 11\ \ \ \  > tmp.txt")
            os.system("cat tmp.lhe | grep -ir 13\ \ \ \  >> tmp.txt")
            
            leptons = [line.strip() for line in open('tmp.txt')]
            for l in xrange(0,len(leptons),2) :
                val1 = filter(None, leptons[l].split(" ") )
                val2 = filter(None, leptons[l+1].split(" ") )
                llp = [float(val1[6])+float(val2[6]), float(val1[7])+float(val2[7]), float(val1[8])+float(val2[8])]
                pt=sqrt(llp[0]*llp[0]+llp[1]*llp[1])
                ibin=dilptH.FindBin(pt)
                wid=dilptH.GetBinWidth(ibin)
                dilptH.Fill(pt,1./wid)
                
            os.system("rm tmp.txt")
            os.system("rm tmp.lhe")

        histos.append(dilptH)

    fOut=ROOT.TFile('dilPt.root','RECREATE')
    for h in histos: h.Write()
    fOut.Close()

else :

    fIn=ROOT.TFile(sys.argv[1])
    q2down=fIn.Get('q2down')
    q2down.SetLineColor(2)
    q2down.SetDirectory(0)
    q2down.Scale(1./q2down.Integral())
    
    q2up=fIn.Get('q2up')
    q2up.SetLineColor(3)
    q2up.SetDirectory(0)
    q2up.Scale(1./q2up.Integral())

    nominal=fIn.Get('nominal')
    nominal.SetLineColor(1)
    nominal.SetDirectory(0)
    nominal.Scale(1./nominal.Integral())

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.1)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.05)

    c=ROOT.TCanvas("c","c",600,600)
    c.Divide(1,2)
    p=c.cd(1)
    p.SetLogx(1)
    nominal.Draw("hist")
    nominal.GetYaxis().SetTitleSize(0.06)
    nominal.GetXaxis().SetTitleSize(0.06)
    nominal.GetXaxis().SetTitleOffset(0.7)
    q2up.Draw("histsame")
    q2down.Draw("histsame")
    leg=ROOT.TLegend(0.5,0.8,0.95,0.9)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.06)
    leg.AddEntry(nominal,nominal.GetTitle(),'f')
    leg.AddEntry(q2up,q2up.GetTitle(),'f')
    leg.AddEntry(q2down,q2down.GetTitle(),'f')
    leg.SetNColumns(3)
    leg.Draw()
    pt=ROOT.TPaveText(0.1,0.95,0.9,0.99,'brNDC')
    pt.SetFillStyle(0)
    pt.SetTextAlign(12)
    pt.SetBorderSize(0)
    pt.SetTextSize(0.06)
    pt.AddText('VBFNLO, #sqrt{s}=8 TeV')
    pt.Draw()

    p=c.cd(2)
    p.SetLogx(1)
    ratioUp=q2up.Clone('q2upWgt')
    ratioUp.Divide(nominal)
    ratioUpGr=ROOT.TGraphErrors(ratioUp)
    ratioUpGr.SetLineColor(ratioUp.GetLineColor())
    ratioUpGr.SetMarkerColor(ratioUp.GetLineColor())
    ratioUpGr.SetMarkerStyle(20)
    #upFunc=ROOT.TF1('q2upWgt_func','max([0]*((1.+[1]+x)/(1.+[2]*x)),0.)',1,250)
    upFunc=ROOT.TF1('q2upWgt_func','[0]*((1.+[1]+x)/(1.+[2]*x))',0,1000)
    upFunc.SetParameter(0,1.0)
    upFunc.SetParLimits(0,-1,1)
    upFunc.SetParameter(1,0)
    upFunc.SetParLimits(1,-200,200)
    upFunc.SetParameter(2,0)
    upFunc.SetParLimits(2,-1,1)
    upFunc.SetLineColor(ratioUpGr.GetLineColor())
    ratioUpGr.Fit(upFunc,'EX0M+')
    
    ratioDown=q2down.Clone('q2downWgt')
    ratioDown.Divide(nominal)
    ratioDownGr=ROOT.TGraphErrors(ratioDown)
    ratioDownGr.SetLineColor(ratioDown.GetLineColor())
    ratioDownGr.SetMarkerColor(ratioDown.GetLineColor())
    downFunc=upFunc.Clone('q2downWgt_func')
    downFunc.SetLineColor(ratioDownGr.GetLineColor())
    ratioDownGr.SetMarkerStyle(20)
    ratioDownGr.Fit(downFunc,'R+')

    ratioUpGr.Draw("ap")
    ratioUpGr.GetYaxis().SetTitle("Ratio")
    ratioUpGr.GetYaxis().SetTitleSize(0.06)
    ratioUpGr.GetXaxis().SetTitleSize(0.06)
    ratioUpGr.GetYaxis().SetRangeUser(0.74,1.16)
    ratioUpGr.GetXaxis().SetTitleOffset(0.7)
    ratioUpGr.GetXaxis().SetTitle( q2up.GetXaxis().GetTitle() )
    upFunc.Draw("same")
    ratioDownGr.Draw("p")
    downFunc.Draw("same")

    
    fOut=ROOT.TFile.Open('vbfnlo_q2var.root','RECREATE')
    ratioUpGr.Write()
    upFunc.Write()
    ratioDownGr.Write()
    downFunc.Write()
    fOut.Close()

    raw_input()
    
    
