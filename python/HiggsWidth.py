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
        XSec100 = [0.000722821, 0.00974505, 0.0186783, 0.0170494, 0.0127224, 0.0123448, 0.00807319, 0.00694976, 0.00492894, 0.00189508, 0.000569369, 0.000196548, 5.8855e-05];
        XSec010 = [0.00313343, 0.0670063, 0.218893, 0.204644, 0.145873, 0.137759, 0.086202, 0.072429, 0.0517826, 0.0192517, 0.00547137, 0.00184722, 0.000516527];
        XSec001 = [0.00578392, 0.129529, 0.443913, 0.414981, 0.29471, 0.27903, 0.172004, 0.142393, 0.10418, 0.0389582, 0.0110436, 0.00363833, 0.00101596];

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
 
