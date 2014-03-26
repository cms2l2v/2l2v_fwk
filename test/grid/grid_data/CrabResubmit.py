#!/usr/bin/env python
import os,sys
import getopt
import commands

workingDir=''
JobsToResubmit=''
JobsToSubmit=''
NJobsToSubmit=0
NJobsToResubmit=0

os.system('rm CrabResubmit.sh')
os.system('rm CrabResubmit.log')
inputFiles = sys.argv[1:]
ReadingData=False
print inputFiles
for input in inputFiles:
	f = open(input,'r')
	for line in f :
                 #identify which is the working directory
                 if(line.find('working directory')>=0):
			linesplits = line.split('/')	
			workingDir = linesplits[len(linesplits)-2]
                        JobsToResubmit=''
                        JobsToSubmit=''
                        NJobsToSubmit=0
			print '###################### %40s ######################' % workingDir

                 #identify the begining of a data block
                 if(line.find('ID    END STATUS            ACTION       ExeExitCode JobExitCode E_HOST')>=0):
                        ReadingData=True
                        continue

                 #if out of a datablock, skip all the lines
                 if(not ReadingData):continue

		 #split the line: ID    END STATUS            ACTION       ExeExitCode JobExitCode E_HOST
		 linesplits = [line[0:5], line[6:9], line[10:27], line[28:41], line[42:52], line[53:64], line[65:len(line)-1] ]

		 #Data block start either by the JobId or by '-----', if it's not the case we reached the end of the datablock
		 if(not (linesplits[0][0].isdigit() or linesplits[0].find('-----')>=0)):
                        print 'END' + workingDir
                        ReadingData=False
                        #Finalize the current project, basically simply printout the resubmit command:
                        NJobsToSubmit=500;
                        NJobsToResubmit=500;
                        #continue

                 #check that we do not need to submit/resubmit jobs
                 if(NJobsToSubmit>=499):
                    if(len(JobsToSubmit)>0):
                       os.system('echo "crab -submit ' + JobsToSubmit + ' -c ' + workingDir+'" >> CrabResubmit.sh')
                    JobsToSubmit=''
                    NJobsToSubmit=0

                 if(NJobsToResubmit>=499):
                    if(len(JobsToResubmit)>0):
                       os.system('echo "crab -resubmit ' + JobsToResubmit + ' -c ' + workingDir+'" >> CrabResubmit.sh')
                    JobsToResubmit=''
                    NJobsToResubmit=0


                 #Job status line should start by JobID
                 if(not linesplits[0][0].isdigit()):continue

                  ############################################### HERE IS THE TREATMENT OF EACH STATUS LINE
                 print linesplits
      
                 #resubmit all aborted jobs or CannotSubmit jobs
                 if(linesplits[1].find('Y')>=0 and (linesplits[2].find('Aborted')>=0 or linesplits[2].find('CannotSubmit')>=0) ):
                        if len(JobsToResubmit)>0: JobsToResubmit+=','
                        JobsToResubmit += str((int(linesplits[0])))
                        NJobsToResubmit+=1;
			continue;

                 #submit jobs in Created State
                 if(linesplits[2].find('Created')>=0):
                        if len(JobsToSubmit)>0: JobsToSubmit+=','
                        JobsToSubmit += str((int(linesplits[0])))
                        NJobsToSubmit+=1
                        continue;

		 #resubmit all jobs finished (END STATUS!=Y) and with exit code != 0
		 if(linesplits[1].find('Y')>=0 and linesplits[4][0].isdigit() and (int(linesplits[4])!=0 or int(linesplits[5])!=0)):
			if len(JobsToResubmit)>0: JobsToResubmit+=','
			JobsToResubmit += str((int(linesplits[0])))
                        NJobsToResubmit+=1;
			continue;

                 #resubmit all jobs finished (END STATUS!=Y) and with exit code != 0
                 if(linesplits[1].find('Y')>=0 and linesplits[5][0].isdigit() and (int(linesplits[5])!=0)):
                        if len(JobsToResubmit)>0: JobsToResubmit+=','
                        JobsToResubmit += str((int(linesplits[0])))
                        NJobsToResubmit+=1;
                        continue;


                 #resubmit all jobs in ready or scheduled state (for too long) 
#                 if(linesplits[1].find('N')>=0 and (linesplits[2].find('Ready')>=0 or linesplits[2].find('Scheduled')>=0)):
#                        if len(JobsToResubmit)>0: JobsToResubmit+=','
#                        JobsToResubmit += str((int(linesplits[0])))
#                        continue;

	f.close()
