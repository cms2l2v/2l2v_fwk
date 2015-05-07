#!/usr/bin/env python
"""
Driver for making figures 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import os
import subprocess

label_json_dict = {
    'pho': '../../data/phojet/phys14_samples.json',
    'dy': '../../data/phojet/phys14_dy_samples.json'
    }

options_dict = {
    '-h': 'print out the help message',
    '-t': 'only print out the command name',
    }

def print_usage():
    sys.stdout.write('NAME\nfig.py - draw figures\n')
    sys.stdout.write('\nSYNOPSIS\n\tfig.py LABEL version [-t] [-h] \n ')
    sys.stdout.write('\nLABEL\n')
    for label, json in label_json_dict.items():
        sys.stdout.write('\t%-5s  %-40s\n' % (label, json))

    sys.stdout.write('\nOPTIONS\n')
    for opt, comment in options_dict.items():
        sys.stdout.write('\t%-5s  %-40s\n' % (opt, comment))

    sys.stdout.write('\nEXAMPLE\n')
    sys.stdout.write('\t./fig.py dy v1.1 \n')
    
    sys.stdout.write('\nAUTHOR\n\tXin Shi <Xin.Shi@cern.ch>\n')

    
def main():
    args = sys.argv[1:]
    if option_exists(args, '-h') or len(args) < 2:
        return print_usage()
        
    label = args[0]
    if label not in label_json_dict:
        return print_usage()
        
    test = option_exists(args, '-t')
    version = args[1] 
    inDir = 'results/%s/' % version 
    outDir = 'plots/%s' % version 
    outFile = 'plots/%s/plotter_%s.root' % (version, label)
    json = label_json_dict[label]
    
    cmd = 'runPhoJetPlotter --inDir %s --outDir %s --outFile %s --json %s' % (
        inDir, outDir, outFile, json)
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

    
