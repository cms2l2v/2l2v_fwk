from HiggsAnalysis.CombinedLimit.PhysicsModel import *


### This is the base python class to study the Higgs width

class HiggsWidth(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.GGsmfixed = False
        self.is2l2nu = False
        self.poiMap = []
        self.pois = {}
        self.verbose = False
        self.xsec= 1.0
        self.xsec_vbf= 1.0
        self.m= 200
        self.w= 10.
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        if process == "ggH_sonl": return "ggH_s_func"
        elif process == "ggH_bonl": return "ggH_b_func"
        elif process == "ggH_sand": return "ggH_sbi_func"
        #if process == "qqH_sonl": return "qqH_s_func"
        #elif process == "qqH_bonl": return "qqH_b_func"
        #elif process == "qqH_sand": return "qqH_sbi_func"
        else:
            return 1

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po == "GGsmfixed":
                print "Will fix CMS_zz4l_GGsm to 1 and float mu"
                self.GGsmfixed = True
            if po == "is2l2nu":
                print "Will consider cards in 2l2nu style (separated S, B, S+B+I)"
                self.is2l2nu = True
            if po == "muAsPOI":
                print "Will combine all ZZ cards"
                self.isCombine = True
            if po.startswith("m="):
                self.m = float(po.replace("m=",""))
                print "mass = ",self.m
            if po.startswith("w="):
                self.w = float(po.replace("w=",""))
                print "width = ",self.w




    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.modelBuilder.out.var("r"):
            print "have r inside"
        else:
            self.modelBuilder.doVar("r[10,0,1000]")
        if self.modelBuilder.out.var("fvbf"):
            print "have fvbf inside"
        else:
            self.modelBuilder.doVar("fvbf[0,0,1]")
            #     if self.is2l2nu:
            #self.setXsec()
            #self.modelBuilder.factory_( "expr::CMS_zz2l2nu_mu(\"@0*(1-@1)*0.0673*0.2*2/1000./%f\", r,fvbf)" %(self.xsec))
            #self.modelBuilder.factory_( "expr::CMS_zz2l2nu_mu_vbf(\"@0*@1*0.0673*0.2*2/1000./%f\", r,fvbf)" %(self.xsec_vbf))
            #poi = "r"

        poi = "r"
        
        self.modelBuilder.factory_( "expr::ggH_s_func(\"@0-sqrt(@0)\", r)")
        self.modelBuilder.factory_(  "expr::ggH_b_func(\"1-sqrt(@0)\", r)")
        self.modelBuilder.factory_(  "expr::ggH_sbi_func(\"sqrt(@0)\", r)")

#self.modelBuilder.factory_( "expr::qqH_s_func(\"@0-sqrt(@0)\", CMS_zz2l2nu_mu_vbf)")
#        self.modelBuilder.factory_(  "expr::qqH_b_func(\"1-sqrt(@0)\", CMS_zz2l2nu_mu_vbf)")
#        self.modelBuilder.factory_(  "expr::qqH_sbi_func(\"sqrt(@0)\", CMS_zz2l2nu_mu_vbf)")


	self.modelBuilder.doSet("POI",poi)



higgswidth = HiggsWidth()
