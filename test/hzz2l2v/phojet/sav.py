#!/usr/bin/env python
"""
Save the output files 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"

import sys
import os
import shutil
import filecmp

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
    sys.stdout.write('\t./sav.py crab_GJet_Pt-15To6000/results results/v1.0 \n')
    
    sys.stdout.write('\nAUTHOR\n\tXin Shi <Xin.Shi@cern.ch>\n')


def main():
    args = sys.argv[1:]
    if option_exists(args, '-h') or len(args) < 1:
        return print_usage()

    src = args[0]    
    dst = args[1]
    mctruthmode = 0
    dtag = ''
    if 'GJet_Pt-15To6000' in src:
        mctruthmode = 22
        dtag = 'MC13TeV_GJets_15To6000'
    
    #shutil.copy2(src, dst) 
    for root, dirs, files in os.walk(src):
        for f in files:
            dstname = convert_dstname(f, dtag, mctruthmode)
            srcfile = os.path.join(src, f)
            dstfile = check_and_join(dst, dstname)
            shutil.copy2(srcfile, dstfile)
            check_update_status(dstfile, 1)
              
    
# ----------------------------------------------------------
# Supporting Functions
# ----------------------------------------------------------
def option_exists(args, key): 
    exists = False 
    if key in args:
        exists = True
    return exists

def convert_dstname(f, dtag, mctruthmode):
    jobId = int(f.split('.')[0].split('_')[1])
    dstname = '%s_%s_filt%s.root' %(dtag, jobId+1, mctruthmode)
    return dstname 

def make_tmpfile(f):
    path, name = os.path.split(f)
    tmpname = '.tmp_' + name
    tmpfile = os.path.join(path, tmpname)
    return tmpfile


def check_update_status(f, verbose=0):
    tmpfile = make_tmpfile(f)
    if not os.access(tmpfile, os.F_OK):
        message = 'created %s ...\n' %f
        shutil.copy2(f, tmpfile)
    elif filecmp.cmp(f, tmpfile, shallow=False):
        message = 'up-to-date: %s\n' % f
    else:
        message = 'updated %s ...\n' %f
    if verbose > 0 :
        sys.stdout.write(message)
    return message


def check_and_join(filepath, filename=None):
    if not os.access(filepath, os.F_OK):
        sys.stdout.write('Creating dir %s ...' % filepath)
        os.makedirs(filepath)
        os.chmod(filepath, 0777)
        sys.stdout.write(' OK.\n')

    if filename == None:
        return
    
    file_ = os.path.join(filepath, filename)
    if os.access(file_, os.F_OK) :
        tmpfile = make_tmpfile(file_)
        shutil.copy2(file_, tmpfile)
    return file_



if __name__ == '__main__':
    main()

    
