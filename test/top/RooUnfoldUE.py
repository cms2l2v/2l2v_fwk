#!/usr/bin/env python


## first look at unfolding...
## first dummy script...
## documentation: http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
##
## command to run: ./RooUnfoldUE.py -i $dir/test.root
## needs a roor file conatining the histogram, here test.root
##
##
##




import os
import optparse
import sys


from ROOT import *
from ROOT import ROOT
ROOT.gSystem.Load("${CMSSW_BASE}/src/RooUnfold-1.1.1/libRooUnfold");

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
# from ROOT import RooUnfoldSvd
# from ROOT import RooUnfoldTUnfold


def openTFile(path, option='r'):
    f =  TFile.Open(path,option)
    if not f.__nonzero__() or not f.IsOpen():
        raise NameError('File '+path+' not open')
    return f

def getHist(file, hist):
    h = file.Get(hist)
    if not h.__nonzero__():
        raise NameError('Histogram '+hist+' doesn\'t exist in '+str(file))
    h.Sumw2()
    return h

def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inputfile'  ,    dest='inputfile'       , help='Name of the input tree file.'          , default=None)
    parser.add_option('-g', '--gen',            dest='genHisto',         help='Name of the histogram at generator level',   default='t#bar{t}172.5/emu_ngench')
    parser.add_option('-r', '--rec',            dest='recHisto',         help='Name of the histogram at reco level',        default='t#bar{t}172.5/emu_nch')
    parser.add_option('-c', '--cor',            dest='corHisto',         help='Name of the correlation histogram',          default='t#bar{t}172.5/emu_nchvsngench')
    parser.add_option('-m', '--method'     ,    dest='method'          , help='Unfolding method. Default = bayes'     , default='bayes')
    parser.add_option('-b', '--batch'      ,    dest='batch'           , help='Use this flag to run root in batch mode.'            , action='store_true'        , default=False)
 
    (opt, args) = parser.parse_args()
   
    if opt.inputfile is None:
        parser.error('No input directory defined!')

    if opt.batch:
        sys.argv.append( '-b' )
        ROOT.gROOT.SetBatch()


    inputfile = opt.inputfile
    method = opt.method
    genHisto = opt.genHisto
    recHisto = opt.recHisto
    corHisto = opt.corHisto
    rebin = 1
    reg_param = 2

    print 'Retrieving histograms from file'
    f = openTFile(inputfile)
    h_gen = getHist(f,genHisto).Rebin(rebin)
    h_gen.SetDirectory(0)
    h_reco = getHist(f,recHisto).Rebin(rebin)
    h_reco.SetDirectory(0)
    h_corr = getHist(f,corHisto).Rebin2D(rebin,rebin)
    h_corr.SetDirectory(0)
    f.Close();
    
    print 'Defining response'
    response= RooUnfoldResponse (h_reco, h_gen,h_corr);
  
    ## three different methods exist to perform the unfolding
    ## see the documentation for more information about those
    print 'Unfolding with %s'%method
    if method == 'svd':
        unfold = RooUnfoldSvd(response, h_reco, reg_param);   
    if method == 'bayes':
        unfold = RooUnfoldBayes(response, h_reco, reg_param);
    if method == 'bin':
        unfold = RooUnfoldBinByBin(response, h_reco)
    h_unfold= unfold.Hreco();
    unfold.PrintTable (cout, h_reco);

    print 'Saving plots'
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadBottomMargin(0.13)
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetPadRightMargin(0.05)
    
    c_corr = TCanvas('c_corr','c_corr',600,600)
    h_corr.Draw('colz')
    c_corr.Print('RooUnfoldUE_%s.png'%(h_corr.GetName()))

    c = TCanvas('c','c',600,600)
    h_unfold.SetMarkerColor(4)
    h_unfold.SetLineColor(4)
    h_unfold.SetMarkerStyle(24)
    h_unfold.SetMaximum(1.25*max(h_unfold.GetMaximum(), h_gen.GetMaximum(), h_reco.GetMaximum()))
    h_unfold.SetMinimum(0.)
    h_unfold.Draw('e1');
    h_reco.SetMarkerStyle(20)
    h_reco.Draw("e1 same");
    h_gen.SetLineColor(4);
    h_gen.Draw("hist same");

    leg=TLegend(0.6,0.6,0.9,0.9,'','brNDC')
    leg.AddEntry(h_unfold,'unfolded','p')
    leg.AddEntry(h_reco,'reconstructed','p')
    leg.AddEntry(h_gen,'generated','f')
    leg.SetTextFont(42)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw()
    
    c.Print('RooUnfoldUE_%s_%s.png'%(method,h_reco.GetName()))
        

if __name__ == '__main__':
    main()
    print 'done'
