import FWCore.ParameterSet.Config as cms
import os,sys

process = cms.Process("MTOPdiff")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')

#CONFIGURE FROM COMMAND LINE
myPYTHIAueSettings="PythiaP11Settings"
myMaxEvents=5000
outputFile="histos.root"
initSeedn=0
if len(sys.argv)>3 :
	outputFile=sys.argv[3]
	outFileBase=os.path.basename( outputFile )
	myPYTHIAueSettings=outFileBase.split('_')[0]
	initSeed=int((outFileBase.split('_')[1]).split('.')[0])
if len(sys.argv)>5 :
	myMaxEvents=int(sys.argv[5])

process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
process.RandomNumberGeneratorService.generator.initialSeed = 123456789+initSeed

print '*********************************************'
print 'Starting new generation with %d events'%myMaxEvents
print 'UE settings are %s'%myPYTHIAueSettings
print 'Output will be stored at %s'%outputFile
print 'Random seed is %s'%str(process.RandomNumberGeneratorService.generator.initialSeed)
print '*********************************************'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(myMaxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.source = cms.Source("EmptySource")

process.genstepfilter.triggerConditions=cms.vstring("generation_step")

#PYTHIA CONFIGURATION
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    ExternalDecays = cms.PSet(
        Tauola = cms.untracked.PSet(
            UseTauolaPolarization = cms.bool(True),
            InputCards = cms.PSet(
                mdtau = cms.int32(0),
                pjak2 = cms.int32(0),
                pjak1 = cms.int32(0)
            )
        ),
        parameterSets = cms.vstring('Tauola')
    ),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(8000.0),
    crossSection = cms.untracked.double(131.7),
    UseExternalGenerators = cms.untracked.bool(True),
    PythiaParameters = cms.PSet(
        PythiaP11Settings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
					'MSTJ(22)=2     ! Decay those unstable particles', 
					'PARJ(71)=10 .  ! for which ctau  10 mm', 
					'MSTP(51) =      7 ! PDF set', 
					'MSTP(52) =      1 ! PDF set internal (=1) or pdflib (=2)', 
					'MSTP( 3) =      1 ! INT switch for choice of LambdaQCD', 
					'MSTU(112)=      5 ! INT n(flavors) for LambdaQCD', 
					'PARU(112)= 0.1600 ! INT LambdaQCD', 
					'PARP( 1) = 0.1600 ! ME/UE LambdaQCD', 
					'PARJ(81) = 0.2600 ! FSR LambdaQCD (inside resonance decays)', 
					'PARP(72) = 0.2600 ! IFSR LambdaQCD (outside resonance decays', 
					'PARP(61) = 0.2600 ! ISR LambdaQCD', 
					'MSTP(64) =      2 ! ISR alphaS type', 
					'PARP(64) = 1.0000 ! ISR renormalization scale prefactor', 
					'MSTP(67) =      2 ! ISR coherence option for 1st emission', 
					'MSTP(68) =      3 ! ISR phase space choice & ME corrections', 
					'PARP(67) = 1.0000 ! ISR Q2max factor', 
					'MSTP(72) =      2 ! IFSR scheme for non-decay FSR', 
					'PARP(71) = 1.0000 ! IFSR Q2max factor in non-s-channel procs', 
					'MSTP(70) =      0 ! ISR IR regularization scheme', 
					'PARP(62) = 1.5000 ! ISR IR cutoff', 
					'PARJ(82) = 1.0000 ! FSR IR cutoff', 
					'MSTP(33) =      0 ! "K" switch for K-factor on/off & type', 
					'MSTP(81) =     21 ! UE model', 
					'PARP(82) = 2.9300 ! UE IR cutoff at reference ecm', 
					'PARP(89) = 7000.0 ! UE IR cutoff reference ecm', 
					'PARP(90) = 0.2650 ! UE IR cutoff ecm scaling power', 
					'MSTP(82) =      3 ! UE hadron transverse mass distribution', 
					'MSTP(88) =      0 ! BR composite scheme', 
					'MSTP(89) =      0 ! BR color scheme', 
					'PARP(79) = 2.0000 ! BR composite x enhancement', 
					'PARP(80) = 0.0150 ! BR breakup suppression', 
					'MSTP(91) =      1 ! BR primordial kT distribution', 
					'PARP(91) = 1.0000 ! BR primordial kT width <|kT|>', 
					'PARP(93) =      10.0000 ! BR primordial kT UV cutoff', 
					'MSTP(95) =      8 ! FSI color (re-)connection model', 
					'PARP(78) = 0.0360 ! FSI color reconnection strength', 
					'PARP(77) = 1.0000 ! FSI color reco high-pT damping strength', 
					'MSTJ(11) =      5 ! HAD choice of fragmentation function(s)', 
					'PARJ( 1) = 0.0870 ! HAD diquark suppression', 
					'PARJ( 2) = 0.1900 ! HAD strangeness suppression', 
					'PARJ( 3) = 0.9500 ! HAD strange diquark suppression', 
					'PARJ( 4) = 0.0430 ! HAD vector diquark suppression', 
					'PARJ( 5) = 0.5000 ! HAD P(popcorn)', 
					'PARJ( 6) = 1.0000 ! HAD extra popcorn B(s)-M-B(s) supp', 
					'PARJ( 7) = 1.0000 ! HAD extra popcorn B-M(s)-B supp', 
					'PARJ(11) = 0.3500 ! HAD P(vector meson), u and d only', 
					'PARJ(12) = 0.4000 ! HAD P(vector meson), contains s', 
					'PARJ(13) = 0.5400 ! HAD P(vector meson), heavy quarks', 
					'PARJ(21) = 0.3300 ! HAD fragmentation pT', 
					'PARJ(25) = 0.6300 ! HAD eta0 suppression', 
					'PARJ(26) = 0.1200 ! HAD eta0prim suppression', 
					'PARJ(41) = 0.3500 ! HAD string parameter a(Meson)', 
					'PARJ(42) = 0.8000 ! HAD string parameter b', 
					'PARJ(45) = 0.5500 ! HAD string a(Baryon)-a(Meson)', 
					'PARJ(46) = 1.0000 ! HAD Lund(=0)-Bowler(=1) rQ (rc)', 
					'PARJ(47) = 1.0000 ! HAD Lund(=0)-Bowler(=1) rb'),
	PythiaP11noCRSettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution',
					    'MSTJ(22)=2     ! Decay those unstable particles',
					    'PARJ(71)=10 .  ! for which ctau  10 mm',
					    'MSTP(51) =      7 ! PDF set',
					    'MSTP(52) =      1 ! PDF set internal (=1) or pdflib (=2)',
					    'MSTP( 3) =      1 ! INT switch for choice of LambdaQCD',
					    'MSTU(112)=      5 ! INT n(flavors) for LambdaQCD',
					    'PARU(112)= 0.1600 ! INT LambdaQCD',
					    'PARP( 1) = 0.1600 ! ME/UE LambdaQCD',
					    'PARJ(81) = 0.2600 ! FSR LambdaQCD (inside resonance decays)',
					    'PARP(72) = 0.2600 ! IFSR LambdaQCD (outside resonance decays',
					    'PARP(61) = 0.2600 ! ISR LambdaQCD',
					    'MSTP(64) =      2 ! ISR alphaS type',
					    'PARP(64) = 1.0000 ! ISR renormalization scale prefactor',
					    'MSTP(67) =      2 ! ISR coherence option for 1st emission',
					    'MSTP(68) =      3 ! ISR phase space choice & ME corrections',
					    'PARP(67) = 1.0000 ! ISR Q2max factor',
					    'MSTP(72) =      2 ! IFSR scheme for non-decay FSR',
					    'PARP(71) = 1.0000 ! IFSR Q2max factor in non-s-channel procs',
					    'MSTP(70) =      0 ! ISR IR regularization scheme',
					    'PARP(62) = 1.5000 ! ISR IR cutoff',
					    'PARJ(82) = 1.0000 ! FSR IR cutoff',
					    'MSTP(33) =      0 ! "K" switch for K-factor on/off & type',
					    'MSTP(81) =     21 ! UE model',
					    'PARP(82) = 3.0500 ! UE IR cutoff at reference ecm',
					    'PARP(89) = 7000.0 ! UE IR cutoff reference ecm',
					    'PARP(90) = 0.2650 ! UE IR cutoff ecm scaling power',
					    'MSTP(82) =      3 ! UE hadron transverse mass distribution',
					    'MSTP(88) =      0 ! BR composite scheme',
					    'MSTP(89) =      0 ! BR color scheme',
					    'PARP(79) = 2.0000 ! BR composite x enhancement',
					    'PARP(80) = 0.0150 ! BR breakup suppression',
					    'MSTP(91) =      1 ! BR primordial kT distribution',
					    'PARP(91) = 1.0000 ! BR primordial kT width <|kT|>',
					    'PARP(93) =      10.0000 ! BR primordial kT UV cutoff',
					    'MSTP(95) =      0 ! FSI color (re-)connection model',
					    'PARP(78) = 0.0000 ! FSI color reconnection strength',
					    'PARP(77) = 0.0000 ! FSI color reco high-pT damping strength',
					    'MSTJ(11) =      5 ! HAD choice of fragmentation function(s)',
					    'PARJ( 1) = 0.0870 ! HAD diquark suppression',
					    'PARJ( 2) = 0.1900 ! HAD strangeness suppression',
					    'PARJ( 3) = 0.9500 ! HAD strange diquark suppression',
					    'PARJ( 4) = 0.0430 ! HAD vector diquark suppression',
					    'PARJ( 5) = 0.5000 ! HAD P(popcorn)',
					    'PARJ( 6) = 1.0000 ! HAD extra popcorn B(s)-M-B(s) supp',
					    'PARJ( 7) = 1.0000 ! HAD extra popcorn B-M(s)-B supp',
					    'PARJ(11) = 0.3500 ! HAD P(vector meson), u and d only',
					    'PARJ(12) = 0.4000 ! HAD P(vector meson), contains s',
					    'PARJ(13) = 0.5400 ! HAD P(vector meson), heavy quarks',
					    'PARJ(21) = 0.3300 ! HAD fragmentation pT',
					    'PARJ(25) = 0.6300 ! HAD eta0 suppression',
					    'PARJ(26) = 0.1200 ! HAD eta0prim suppression',
					    'PARJ(41) = 0.3500 ! HAD string parameter a(Meson)',
					    'PARJ(42) = 0.8000 ! HAD string parameter b',
					    'PARJ(45) = 0.5500 ! HAD string a(Baryon)-a(Meson)',
					    'PARJ(46) = 1.0000 ! HAD Lund(=0)-Bowler(=1) rQ (rc)',
					    'PARJ(47) = 1.0000 ! HAD Lund(=0)-Bowler(=1) rb'),
	PythiaZ2Settings  = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution',
					'MSTJ(22)=2     ! Decay those unstable particles',
					'PARJ(71)=10 .  ! for which ctau  10 mm',
					'MSTP(33)=0     ! no K factors in hard cross sections',
					'MSTP(2)=1      ! which order running alphaS',
					'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)',
					'MSTP(52)=2     ! work with LHAPDF',
					'PARP(82)=1.832 ! pt cutoff for multiparton interactions',
					'PARP(89)=1800. ! sqrts for which PARP82 is set',
					'PARP(90)=0.275 ! Multiple interactions: rescaling power',
					'MSTP(95)=6     ! CR (color reconnection parameters)',
					'PARP(77)=1.016 ! CR',
					'PARP(78)=0.538 ! CR',
					'PARP(80)=0.1   ! Prob. colored parton from BBR',
					'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter',
					'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter',
					'PARP(62)=1.025 ! ISR cutoff',
					'MSTP(91)=1     ! Gaussian primordial kT',
					'PARP(93)=10.0  ! primordial kT-max',
					'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default',
					'MSTP(82)=4     ! Defines the multi-parton model'),
	processSettings = cms.vstring('MSEL=0         ! User defined processes',
				      'MSUB(81)  = 1     ! qqbar to QQbar',
				      'MSUB(82)  = 1     ! gg to QQbar',
				      'MSTP(7)   = 6     ! flavor = top',
				      'PMAS(6,1) = 172.5  ! top quark mass'),
        parameterSets = cms.vstring(myPYTHIAueSettings,
				    'processSettings')
	)
)



process.demo = cms.EDAnalyzer('TopMassDifferenceAnalyzer')

#output
process.TFileService = cms.Service("TFileService",
				   fileName = cms.string(outputFile),
				   closeFileFast = cms.untracked.bool(True)
				   )

#the sequence to run
process.ProductionFilterSequence = cms.Sequence(process.generator)
process.generation_step = cms.Path(process.pgen)
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

process.p = cms.Path(process.demo)
