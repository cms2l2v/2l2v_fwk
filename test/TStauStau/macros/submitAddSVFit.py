#!/usr/bin/env python

import os
import stat

templateFile = open("addSVFit_template.sh", "r")
baseDir = "./test"
outDir = "/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_fixed"

template = templateFile.read();

if baseDir != "":
  baseDir = baseDir + "/"
  if not os.path.exists(baseDir):
    os.mkdir(baseDir)

if outDir != "":
  outDir = outDir + "/"
  if not os.path.exists(outDir):
    os.mkdir(outDir)

listOfFiles = open("files.txt", "r")
submitCMD = []
for file in listOfFiles.readlines():
  file = file.strip()
  dir, filename = os.path.split(file)
  print "Processing file "+filename

  Info = {"inFile":file, "outFile":outDir+filename}

  outfile = open(baseDir+"submit_"+filename+".sh", "w")
  outfile.write(template.format(**Info))
  outfile.close()

  st = os.stat(baseDir+"submit_"+filename+".sh")
  os.chmod(baseDir+"submit_"+filename+".sh", st.st_mode | stat.S_IEXEC)

  submitCMD.append("qsub "+baseDir+"submit_"+filename+".sh")

templateFile.close()

print "Submitting jobs:"
for cmd in submitCMD:
  os.system(cmd)
