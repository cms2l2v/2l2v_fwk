from HiggsAnalysis.CombinedLimit.PhysicsModel import *

class XsecRatioModel(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar('brratio[1,0,10]');
        self.modelBuilder.doVar('commonxsec[1,0,10]')
        self.modelBuilder.out.var("brratio").setRange(0,10)
        self.modelBuilder.out.var("brratio").setConstant(False)
        self.modelBuilder.out.var("commonxsec").setRange(0,10)
        self.modelBuilder.out.var("commonxsec").setConstant(False)
        self.modelBuilder.doSet("POI",'brratio')

        # mutau = commonxsec*ratio
        # emu = commonxsec
        # mutau/emu = commonxsec*ratio/commonxsec = ratio
        self.modelBuilder.factory_('expr::Scaling_ltau("@0*@1", commonxsec, brratio)')
        self.modelBuilder.factory_('expr::Scaling_emu("@0", commonxsec)')

        self.processScaling = { 'signal_emu':'emu', 'signal_singlemu':'ltau', 'signal_singlee':'ltau'}

        self.modelBuilder.out.Print()
        
    def getYieldScale(self,bin,process):

        for prefix, model in self.processScaling.iteritems():
            if process.startswith(prefix):
                print "Process ",process,"will be scaled according to the model"
                return 'Scaling_'+model
            
        return 1


xsecRatioModel = XsecRatioModel()

