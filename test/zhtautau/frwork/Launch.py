#!/usr/bin/env python

import string
import os
import sys
import LaunchOnCondor

FarmDirectory = "FARM"
JobName = "lltautauToysFirstTry"
os.system("make clean; make")
LaunchOnCondor.Jobs_RunHere = 0
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
for i in range(0,11):
   LaunchOnCondor.Jobs_FinalCmds = ["mv justForDebug_"+str(i) + " " + os.getcwd()+"/." ]
   LaunchOnCondor.SendCluster_Push(["BASH", os.getcwd()+"/Analyzer  --File1 "+os.getcwd()+"/plotter_2016_02_09.root --File2 "+os.getcwd()+"/plotter_2016_02_09.root --Index "+str(i)+" --Ntoy 100 --Data --Dir justForDebug_"+str(i)+" --NoXserver" ])
   #LaunchOnCondor.SendCluster_Push(["BASH", os.getcwd()+"/Analyzer  --File1 "+os.getcwd()+"/../plotter_Approval.root --File2 "+os.getcwd()+"/../plotter_Approval_tree.root --Index "+str(i)+" --Ntoy 1 --Dir justForDebug_"+str(i)+" --NoXserver" ])
LaunchOnCondor.SendCluster_Submit()

