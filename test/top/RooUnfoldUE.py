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

def GetKeyNames( file, dir = "" ):
        file.cd(dir)
        return [key.GetName() for key in gDirectory.GetListOfKeys()]



def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inputfile'  ,    dest='inputfile'       , help='Name of the input tree file.'          , default=None)
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


    rebin = 1

    #f = openTFile('/afs/cern.ch/user/j/jueugste/cmssw/Top_HLTweights/CMGTools/CMSSW_5_3_3_patch3/src/UserCode/TopMass//test.root')
    f = openTFile(inputfile)

    ## get the variables (for each variable one needs a gen, reco, and correlation histogram)
    keys = GetKeyNames(f)

    varis = list(set([key.split('_')[-1] for key in keys]))
    for var in varis:

        h_gen = getHist(f,'h_gen_'+var).Rebin(rebin)
        h_reco = getHist(f,'h_reco_'+var).Rebin(rebin)
        h_corr = getHist(f,'h_corr_'+var).Rebin2D(rebin,rebin)
        print 'defining response'
        response= RooUnfoldResponse (h_gen, h_reco,h_corr);

        c_corr = TCanvas('c_corr','c_corr',600,600)
        h_corr.Draw('colz')
        c_corr.Print('correlation.png')


        ## loop over different reg paramaeters to optimize, especially for svd
        # for i in xrange(10):


        ## this is the regularisation parameter
        reg_param = 2
        
        ## three different methods exist to perform the unfolding
        ## see the documentation for more information about those
        
        if method == 'svd':
            unfold = RooUnfoldSvd(response, h_reco, reg_param);   
        if method == 'bayes':
            unfold = RooUnfoldBayes(response, h_reco, reg_param);
        if method == 'bin':
            unfold = RooUnfoldBinByBin(response, h_reco)
             
        h_unfold= unfold.Hreco();
        unfold.PrintTable (cout, h_reco);
        c = TCanvas('c','c',600,600)
        c.SetLogy(1)
        h_unfold.SetMarkerStyle(4)
        h_unfold.SetMarkerColor(2)
        h_unfold.SetLineColor(2)
        h_unfold.SetMaximum(1.25*max(h_unfold.GetMaximum(), h_gen.GetMaximum(), h_reco.GetMaximum()))
        h_unfold.SetMinimum(1.)
        h_unfold.Draw('e1');
        h_reco.SetMarkerStyle(20)
        h_reco.Draw("e1 same");
        h_gen.SetLineColor(8);
        h_gen.Draw("hist same");
        
        c.Print('test_RooUnfoldUE_'+var+'.png')
        

if __name__ == '__main__':
    main()
    print 'done'
