This directory holds the scripts and other ancillary information for the TStauStau analysis.


Last Update: 11/06/2014 - By Cristóvão

How-to run  (Work in progress)
-------------
This framework, like many others, runs on its own nTuples, which should be produced in a first step.
In a second step, the analysis itself is run on all the desired nTuples, producing plots and whatever else is desired.
In a next step, the desired plots are combined to produce the final plots/histograms/tables.

#### nTuple production
The nTuples are produced using the GRID infrastructure, submitting jobs with crab via multicrab or directly with crab.

The subdirectories "grid_data", "grid_mc" and "grid_signal" hold the necessary scripts to perform this step, "grid_data" for the data datasets, "grid_mc" for the MC datasets and "grid_signal" for the signal MC datasets.

Inside each subdirectory there is a "crab.cfg" file, this file holds the basic crab configurations to submit the jobs, such as the output storage element, how many events each job should process and so on.
The "multicrab.cfg" file, holds the configurations for multicrab, the first section configures the common properties for all sets of jobs (one set of jobs will be created for each dataset) it starts by referencing the "crab.cfg" file, and then sets whatever properties are needed to be different from the ones defined in the "crab.cfg" file, such as the output directory and for data, the luminosity mask to use (Golden JSON for instance).
Then, several sections are defined in "multicrab.cfg", one for each dataset, specifying which dataset is to be processed and any other dataset specific options desired.
On a lxplus machine, in order to use the crab and multicrab executables, instead of running the command `cmsenv`, one should issue the following commands:
```bash
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh 
cmsenv
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh
```
Or, on an NCG machine, with the following commands:
```bash
cmsenv
source /exper-sw/cmst3/cmssw/crab/crab.sh
```

##### Job creation and submission
In a first step, the cmssw config files should be created and placed in their respective places. In the parent directory to this one, there should be a set of files named runObjectProducer_*_cfg.py. For each of these files, the following command should be run `edmConfigDump runObjectProducer_*_cfg.py >> grid_*/runObjectProducer_*_cfg.py`. Nb: The data object producer should be placed in the "grid_data" folder, the mc object producer should be placed in the "grid_mc" folder and the "grid_signal" folder and the fastsim object producer should be placed in the "grid_signal" folder.
After producing the python cfg files, some small modifications need to be done to them. The edmConfigDump utility prints out some messages at the beginning of the file which should be removed, up to the line:
```python
import FWCore.ParameterSet.Config as cms
```
The utility also inserts some non-existant process to the cmssw path, which should be removed, search for process.none and remove all instances. Finally, before continuing, the config file should be tested using cmsRun to insure everything is ok.

Next, the jobs should be created, using the command `multicrab -create` inside each subdirectory can achieve this.
In order to start running the jobs, they must be submitted to the GRID, this can be done either with the command `multicrab -submit`, which will submit all jobs of all datasets in one go, or it can be done with the command `crab -submit -c [DIR]`, where DIR is the directory created for one of the datasets, this will only submit the jobs for that dataset. (Notice that if a certain dataset has more than 500 jobs, the jobs should be split in groups of 500  and each group submitted independently, either for the `multicrab` or the `crab` command, ie: `crab -submit 1-500 -c [DIR]` for the first 500 jobs)
The jobs should be monitored regularly to insure everything is running fine, use the commands `multicrab -status` or `crab -status -c [DIR]` (depending on if you want to check on all or just one of the datasets) to get the status of the jobs.
Jobs that have terminated successfully will have a status of 0, other statuses, generally indicate an error and can normally be recovered by resubmitting.
Other useful commands for job management are:
- `crab -kill [jobs] -c [DIR]` - This command will kill the listed jobs of the given dataset, the jobs option is optional for this command and all the following ones and when it is not present, the operation is applied to all jobs. Jobs that have been cancelled, must be killed before resubmitting.
- `crab -getoutput [jobs] -c [DIR]` - This command retrieves the output of the listed jobs for the given dataset, the output is the std::cout and std::err as well as any other reports from the jobs, not necessarily the root files, this depends on the options given to crab. Sometimes it is necessary to run this command on some jobs before resubmitting them.
- `crab -resubmit [jobs] -c [DIR]` - This command resubmits the listed jobs for the given dataset.
- `crab -report -c [DIR]` - This command generates a report on the processed jobs. This is important to generate the pileup distribution for data datasets.

*(There should be versions of all these commands for multicrab, but for these fine-grained controls, I find it more useful to control on which dataset I am running it on, so I use these versions much more often.)*

The jobs can fail due to too much "Wall clock time", this means that the jobs have been runnning for too long. If the number of events the jobs are running on has been optimized, you might just have been unlucky and resubmitting the job should solve the problem. However, sometimes the job really does have too many events to be processed in the given time. In this case, no matter how many times you resubmit, the jobs will fail. If the fraction of jobs in this status is significant, it is probably a better idea to change the "crab.cfg" file, reducing the number of events per job (for data this is defined as the lumis per job) and resubmitting the whole dataset. If the fraction of jobs in this situation is small, you could consider accepting the already processed jobs and using the "task_missingLumiSummary.json" generated by the report crab command as the lumifilter for a new set of jobs on that dataset, in order to fill in the missing lumisections.
(**ToDo: More details on how to do this**)

Once all jobs and events have been processed, the resulting output must be merged (so we have a more reasonable number of files). The merged nTuples must also be small enough so that when running our analysis on each file, it does not take too much time (lxbatch has several queues, we want to use the 8nh queue, which has a limit of 8 "normalised" hours).
However, before merging, we will create the pileup distributions and get the integrated luminosity for data, these should be used later in order to perform the pileup reweighting of the MC. For this, make sure you have ran the report crab command on each of the data datasets, which should create three files inside the directory "[DIR]/res/":
- inputLumiSummaryOfTask.json   - This one holds the lumis that were requested to be run on (ie: the lumi mask given to crab).
- task_missingLumiSummary.json  - This one, which was already mentioned, holds the missing lumisections, there shouldn't be any missing lumisections.
- lumiSummary.json              - This one holds all the lumis that are present, it is the one we will be using, either for obtaining the PU distribution or the integrated luminosity.

##### PU distribution calculation
(*Instructions to produce PU distribution:* [TWiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Calculating_Your_Pileup_Distribu))

The following actions have to be repeated for each data dataset, where DIR is the directory for that dataset. The commands given assume the current directory is "grid_data".
  1. Make sure the file "[DIR]/res/lumiSummary.json" was created
  2. Get the file with the pileup per lumisection, it is made available centrally in the directory "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/". I used the command `cp /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-208686_All_2012_pixelcorr.txt ./` - This step only needs to be done once, not once for each dataset... :stuck_out_tongue:
  3. Calculate the pileup distribution with `pileupCalc.py`. In my case I used the command `pileupCalc.py -i [DIR]/res/lumiSummary.json --inputLumiJSON pileup_JSON_DCSONLY_190389-208686_All_2012_pixelcorr.txt --calcMode observed --minBiasXsec 69400 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram.root` - nb: for minBiasXsec you should use 68000 or 69400 depending if you are using 2011 or 2012 data, respectively. Also, the number of bins can and should be customized, from my experience, 100 bins has been adequate, also note that the number of bins corresponds to the max number of bins, we want one bin for each number of vertexes, since this value is an integer, however the script should be able to handle a non integer bin size.
  4. **ToDo: Convert the histogram to a list to be placed in config file for analysis**

##### Integrated Luminosity calculation
(*Instructions to get integrated luminosity:* [TWiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiCalc))

The CMSSW is probably not yet set up to produce the integrated luminosity information, so run the commands:
```bash
git clone https://github.com/cms-sw/RecoLuminosity-LumiDB.git $CMSSW_BASE/src/RecoLuminosity/LumiDB
cd $CMSSW_BASE/src/RecoLuminosity/LumiDB
git checkout V04-02-10
scram b
```
  1. Calculate the integrated luminosity with `pixelLumiCalc.py`. In my case I used the command `pixelLumiCalc.py overview -i [DIR]/res/lumiSummary.json >> [DIR].lumi` (The output was saved in a file with the same name as the directory so there is an association between the two)

##### nTuple merging
The merging of the nTuples should be done on a machine with direct access to the output of the jobs, in my case and with the default options in the crab and multicrab configuration files, this is on a NCG machine.
A CMSSW environment with the framework will have to be replicated there as well, in case you didn't start working from there directly.

In order to merge the nTuples, the json file with the list of processes will be used. There might not be enough space to directly process all the samples in one go, so please check if this is the case. If so, you will need to split the json into parts, processing one at a time, when doing this, pay attention to the json syntax, in particular the commas between elements, never leave a comma dangling at the end of a line if there is no next element.
