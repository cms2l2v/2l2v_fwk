from HiggsAnalysis.CombinedLimit.PhysicsModel import *

 
### This is the base python class to study the Higgs width
 
class Higgswidth(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.GGsmfixed = False
        self.is2l2nu = False
        self.poiMap = []
        self.pois = {}
        self.verbose = False
        self.xsec= 1.0 
        self.m= 200. 
        self.w= 10. 
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False
 
    def getYieldScale(self,bin,process):
        if process == "ggH_sonl": return "ggH_s_func"
        elif process == "ggH_bonl": return "ggH_b_func"
        elif process == "ggH_sand": return "ggH_sbi_func"
        if process == "qqH_sonl": return "qqH_s_func"
        elif process == "qqH_bonl": return "qqH_b_func"
        elif process == "qqH_sand": return "qqH_sbi_func"
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
           
    def setXsec(self):
        marr = [200 , 300 ,400 ,500 ,600 ,700 ,800 ,900 ,1000 ,1500 ,2000 ,2500 ,3000];
        XSec100 = [0.000716943, 0.00960252, 0.0176147, 0.0172144, 0.0121598, 0.00997052, 0.00723416, 0.00698813, 0.0054009, 0.00189508, 0.000560771, 0.000188943, 6.75111e-05];
        XSec010 = [0.00310052, 0.0659359, 0.208345, 0.206409, 0.139248, 0.108503, 0.0774346, 0.0731695, 0.0557424, 0.0192516, 0.00551768, 0.00176027, 0.000600868];
        XSec001 = [0.00573029, 0.127286, 0.422868, 0.418509, 0.281284, 0.217461, 0.155051, 0.14615, 0.112064, 0.038958, 0.0110115, 0.00348143, 0.00119788];

        for index in range(len(marr)):
		if self.m == marr[index]:
	  		if self.w == 100:
	  			self.xsec = XSec100[index]
	  		elif self.w == 10:
	  			self.xsec = XSec010[index]
	  		elif self.w == 5:
	  			self.xsec = XSec001[index]
			else:
				raise RuntimeError, "Don't know this mass, width %f, %f" %(self.m, self.w)
	print "mass: ", self.m
	print "width: ", self.w
	print "xsec: ", self.xsec

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.modelBuilder.out.var("r"):
            print "have r inside"
        else:
            self.modelBuilder.doVar("r[1,0,10000]")
        if self.is2l2nu:
            #self.modelBuilder.doVar("CMS_zz4l_GGsm[1.,0.,50.]")
            #self.modelBuilder.doVar("CMS_zz4l_mu[1.,0.,1000.]")
            #self.modelBuilder.doVar("CMS_widthH_kbkg[1.,0.,2.]")
            self.setXsec()
            self.modelBuilder.factory_( "expr::CMS_zz4l_mu(\"@0*0.0673*0.2*2/1000./%f\", r)" %(self.xsec))
            poi = "r"            
        #if self.GGsmfixed:
            #self.modelBuilder.out.var("CMS_zz4l_GGsm")
            #self.modelBuilder.out.var("CMS_zz4l_GGsm").setVal(1)
            #self.modelBuilder.out.var("CMS_zz4l_GGsm").setConstant(True)
            #self.modelBuilder.out.var("CMS_zz4l_mu")
            #print "Fixing CMS_zz4l_GGsm"
          #  poi = "CMS_zz4l_mu"
        #else:
           #poi = "r"
       

        self.modelBuilder.factory_( "expr::ggH_s_func(\"@0-sqrt(@0)\", CMS_zz4l_mu)")
        self.modelBuilder.factory_(  "expr::ggH_b_func(\"1-sqrt(@0)\", CMS_zz4l_mu)")
        self.modelBuilder.factory_(  "expr::ggH_sbi_func(\"sqrt(@0)\", CMS_zz4l_mu)")

        self.modelBuilder.factory_( "expr::qqH_s_func(\"@0-sqrt(@0)\", CMS_zz4l_mu)")
        self.modelBuilder.factory_(  "expr::qqH_b_func(\"1-sqrt(@0)\", CMS_zz4l_mu)")
        self.modelBuilder.factory_(  "expr::qqH_sbi_func(\"sqrt(@0)\", CMS_zz4l_mu)")       


	self.modelBuilder.doSet("POI",poi)


       
higgswidth = Higgswidth()
 
