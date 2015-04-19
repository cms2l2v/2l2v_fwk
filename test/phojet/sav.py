#!/usr/bin/env python
"""
Save the output files 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import os
import shutil 

def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print './sav.py label'
        sys.exit()
        
    src = args[0]    
    test = option_exists(args, '-t')
    dst = args[1]
    shutil.move(src, dst) 

if __name__ == '__main__':
    main()

    
