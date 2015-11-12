import ROOT
import pickle

"""
Takes a directory on eos (starting from /store/...) and returns a list of all files with 'prepend' prepended
"""
def getEOSlslist(directory, tag, mask=''):
    from subprocess import Popen, PIPE
    print 'looking into: '+directory+'...'

    #eos_cmd = '/afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select'
    #data = Popen([eos_cmd, 'ls', '/eos/cms/'+directory],stdout=PIPE)
    data = Popen(['ls', directory],stdout=PIPE)
    out,err = data.communicate()

    full_list = []

    ## if input file was single root file:
    if directory.endswith('.root'):
        if len(out.split('\n')[0]) > 0:
            return [prepend + directory]
        
    ## instead of only the file name append the string to open the file in ROOT
    for line in out.split('\n'):
        if len(line.split()) == 0: continue
        if tag not in line: continue
        full_list.append( directory + '/' + line)

    ## strip the list of files if required
    if mask != '':
        stripped_list = [x for x in full_list if mask in x]
        return stripped_list

    ## return 
    return full_list

"""
Loops over a list of samples and produces a cache file to normalize MC
"""
def produceNormalizationCache(samplesList,inDir,cache,xsecWgts,integLumi,puWgts):

    #open pileup files
    puDists={}
    for puEstimate in ['nom','up','down']:
        try:
            fIn=ROOT.TFile.Open('${CMSSW_BASE}/src/UserCode/llvv_fwk/data/PileupData_%s.root'%puEstimate)
            puDists[puEstimate]=fIn.Get('pileup').Clone(puEstimate)
            puDists[puEstimate].SetDirectory(0)
        except:
            pass        
    hputrue=None
    try:
        hputrue=puDists['nom'].Clone('hputrue')
        hputrue.SetDirectory(0)
        print 'Pileup distributions will be computed, this make take a while...'
    except:
        print 'No pileup distributions found under data/ pileup weights won\'t be stored'

    #loop over samples
    for tag,sample in samplesList: 

        if sample[1]==1 : 
            xsecWgts[tag]=None
            continue

        if tag in xsecWgts:
            print '[Warning] won\'t override current definition for',tag,'. Use --resetCache option to override'
            continue
        

        input_list=getEOSlslist(directory=inDir, tag=tag)            
        xsec=sample[0]            
        norigEvents=None
        if hputrue : hputrue.Reset('ICE')
        for f in input_list:
            fIn=ROOT.TFile.Open(f)
            if norigEvents is None:
                norigEvents=fIn.Get('all_initNNLO').Clone('xsecwgts')
                norigEvents.SetDirectory(0)
                norigEvents.Reset('ICE')
            norigEvents.Add(fIn.Get('all_initNNLO'))
            if hputrue:
                hputrue.SetDirectory(fIn)
                fIn.Get('minievents').Draw('putrue>>+hputrue','ttbar_w[0]','goff')
                hputrue.SetDirectory(0)
            fIn.Close()
        try:
            for xbin in xrange(1,norigEvents.GetNbinsX()+1):
                norigEvents.SetBinContent(xbin,xsec/norigEvents.GetBinContent(xbin))
                norigEvents.SetBinError(xbin,0.)
        except:
            print 'No normalization histogram for ',tag
        xsecWgts[tag]  = norigEvents
        integLumi[tag] = 1./norigEvents.GetBinContent(1) if norigEvents else 0.
        puGr=[]
        if hputrue and hputrue.Integral()>0:
            hputrue.Scale(1./hputrue.Integral())
            for puEst in ['nom','up','down']:
                hwgt=puDists[puEst].Clone('htemp')
                if hwgt.Integral()==0 : 
                    puGr.append(None)
                    continue
                hwgt.Scale(1./hwgt.Integral())
                hwgt.Divide(hputrue)
                puGr.append( ROOT.TGraph(hwgt) )
                puGr[-1].SetName(puEst+'_wgt')
                hwgt.Delete()
        puWgts[tag] = puGr

        if norigEvents:
            print '... %s cross section=%f pb weights sum(initial events)=%3.0f lumi=%3.2f/fb' % (tag,xsec,xsec/norigEvents.GetBinContent(1),integLumi[tag]/1000.)
            
        
    #dump to file    
    cachefile=open(cache,'w')
    pickle.dump(xsecWgts, cachefile, pickle.HIGHEST_PROTOCOL)
    pickle.dump(integLumi, cachefile, pickle.HIGHEST_PROTOCOL)
    pickle.dump(puWgts, cachefile, pickle.HIGHEST_PROTOCOL)
    cachefile.close()
    print 'Produced normalization cache and pileup weights @ %s'%cache

    return xsecWgts,integLumi,puWgts
