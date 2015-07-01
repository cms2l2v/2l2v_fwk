from HiggsAnalysis.CombinedLimit.PhysicsModel import *

class XsecRatioModel(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar('ratio[1,0,2]');
        self.modelBuilder.doVar('commonxsec[1,0,2]');
        self.modelBuilder.doSet('POI','ratio')
        # mutau = commonxsec*ratio
        # emu = commonxsec
        # mutau/emu = commonxsec*ratio/commonxsec = ratio
        self.modelBuilder.factory_('expr::Scaling_ltau("@0*@1", commonxsec, ratio)')
        self.modelBuilder.factory_('expr::Scaling_emu("@0", commonxsec)')

        self.processScaling = { 'emu_signal':'emu', 'mutau_signal':'ltau', 'etau_signal':'ltau'}

        self.modelBuilder.out.Print()
        
    def getYieldScale(self,bin,process):

        for prefix, model in self.processScaling.iteritems():
            if process.startswith(prefix):
                return 'Scaling_'+model
            
        return 1


xsecRatioModel = XsecRatioModel()

