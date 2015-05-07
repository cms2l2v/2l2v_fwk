#!/usr/bin/env python
"""
Save the output files 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import os
import shutil

options_dict = {
    '-h': 'print out the help message',
}


def print_usage():
    sys.stdout.write('NAME\nsav.py - save files from crab\n')
    sys.stdout.write('\nSYNOPSIS\n\tsav.py src dst [-h] \n ')
    sys.stdout.write('\nOPTIONS\n')
    for opt, comment in options_dict.items():
        sys.stdout.write('\t%-5s  %-40s\n' % (opt, comment))

    sys.stdout.write('\nEXAMPLE\n')
    sys.stdout.write('\t./sav.py crab_phojet/crab_GJet_Pt-15to3000/results/ results/v1.1/MC13TeV_GJet_Pt-15to3000 \n')
    
    sys.stdout.write('\nAUTHOR\n\tXin Shi <Xin.Shi@cern.ch>\n')


def main():
    args = sys.argv[1:]
    if option_exists(args, '-h') or len(args) < 1:
        return print_usage()
        
    src = args[0]    
    dst = args[1]
    shutil.move(src, dst) 


# ----------------------------------------------------------
# Supporting Functions
# ----------------------------------------------------------
def option_exists(args, key): 
    exists = False 
    if key in args:
        exists = True
    return exists

if __name__ == '__main__':
    main()

    
