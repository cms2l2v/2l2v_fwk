#!/usr/bin/env python
"""
Driver for making figures 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import os
import subprocess


def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print './fig.py label'
        sys.exit()
        
    label = args[0]    
    test = option_exists(args, '-t')
    inDir = 'results/%s/' % label 
    outDir = 'plots/%s' % label 
    outFile = 'plots/%s/plotter.root' % label 
    cmd = 'runPhoJetPlotter --inDir %s --outDir %s --outFile %s' % (
        inDir, outDir, outFile)
    output = proc_cmd(cmd, test=test)
    if output:
        print output

# ----------------------------------------------------------
# Supporting Functions
# ----------------------------------------------------------
def option_exists(args, key): 
    exists = False 
    if key in args:
        exists = True
    return exists


def proc_cmd(cmd, test=False, verbose=1, procdir=None, shell=False):
    if test:
        sys.stdout.write(cmd+'\n')
        return 
    cwd = os.getcwd()
    if procdir is not None:
        os.chdir(procdir)
    args = cmd.split()
    if isinstance(cmd, list) or shell:
        args = cmd 
    process = subprocess.Popen(args, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=shell)
    stdout = process.communicate()[0]
    if 'error' in stdout:
        sys.stdout.write(stdout)
    if procdir is not None:
        os.chdir(cwd)
    return stdout


if __name__ == '__main__':
    main()

    
