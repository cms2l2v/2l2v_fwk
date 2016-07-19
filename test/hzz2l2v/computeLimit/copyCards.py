#!/usr/bin/env python
import urllib
import string
import os,sys,time
import glob
import commands
import json

SVNdirectory = "/nfs/scratch/fynu/quertenmont/15_09_21_HZZ2l2v/SVN/HIG-16-001/"

SUBMASS = [200, 250, 300, 350, 400,450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1250, 1500]
for CP in [1.0, 0.6, 0.3, 0.1]:
   for M in SUBMASS:
      localDir = "cards_SB13TeV_cp%4.2f_brn%4.2f/%04i" % (CP, 0.0, M)
      svnDir   = SVNdirectory + "/%i/cp%03i_brnew%02i/" % (M,CP*100.0, 0.0*100)
      os.system("mkdir -p " + svnDir)
      print("cp -rd " + localDir+"/hzz2l2v_*.{dat,root} " + svnDir + "/.")
      os.system("cp -rd " + localDir+"/hzz2l2v_*.{dat,root} " + svnDir + "/.")
