#!/usr/bin/env python
from subprocess import call, Popen, PIPE
from os import path

EOSCMD = '/afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select'

def get_eoslslist(directory, tag='', prepend='root://eoscms/'):
    '''Takes a directory on eos (starting from /store/...) and returns
    a list of all files with root://eoscms//eos/cms/ prepended'''
    data = Popen([EOSCMD, 'ls', directory], stdout=PIPE)
    out,err = data.communicate()

    if not directory.startswith('/'):
        directory = '/'+directory

    a = []
    for line in out.split('\n'):
        if len(line.split()) == 0: continue
        if not line.strip().endswith('.root'): continue
        a.append(prepend + directory + '/' + line)

    ## strip the list of files
    if tag != '':
        b = [s for s in a if tag in s]
        return b
    ## if no tag given, run over all files
    else:
        return a

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
             for i in range(wanted_parts) ]


def mergeAndMove(packedargs):
    (sample_name, n, filelist, output_dir, options) = packedargs
    # hadd files to temp dir
    haddtarget = path.join(options.tmpDir, sample_name,
                           '%s_%d.root'%(sample_name,n))
    print 'merging to', haddtarget
    if not options.dry:
        call(['hadd', '-f', '-k', haddtarget] + filelist)
    else:
        print ['hadd', '-f', '-k', haddtarget] + filelist


    # Create target dir on eos
    if not options.dry:
        call([EOSCMD, 'mkdir', '-p', output_dir])
    else:
        print [EOSCMD, 'mkdir', '-p', output_dir]

    # Move files to target directory
    print 'moving to', output_dir
    if output_dir.startswith('/eos/cms/'):
        cmsoutdir = output_dir[8:]
    else:
        cmsoutdir = output_dir
    cmdlist = ['cmsStage', '-f', haddtarget, cmsoutdir]
    if options.dry:
        cmdlist.insert(2, '-d')
    call(cmdlist)
    return True

if __name__ == "__main__":
    from optparse import OptionParser
    usage = """
    \t %prog inputdir nsplit outputdir
    \t %prog MC8TeV_TTJets 10 MC8TeV_TTJets_merged --eosBaseDir /eos/cms/store/cmst3/user/stiegerb/Mar26/

    inputdir and outputdir are relative to eosBaseDir

    This will split the list of files into 'nsplit' groups, merge them
    to a temporary directory, then move them to the output directory
    and finally clean up.
    """
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--eosBaseDir",
                      default="/eos/cms/store/cmst3/user/stiegerb/Mar26/",
                      action="store", type="string", dest="eosBaseDir",
                      help=("Base directory on eos [default: %default]"))
    parser.add_option("-t", "--tmpDir",
                      default="/tmp/stiegerb/",
                      action="store", type="string", dest="tmpDir",
                      help=("Temporary local directory for merging"
                            "[default: %default]"))
    parser.add_option("-d", "--dry",
                      action="store_true", dest="dry",
                      help=("Dry running"))
    (options, args) = parser.parse_args()


    if len(args)==3:
        input_dir = path.join(options.eosBaseDir,args[0])
        sample_name = args[0]
        nsplit = int(args[1])
        output_dir = path.join(options.eosBaseDir,args[2])

        inputfiles = get_eoslslist(input_dir)
        step = len(inputfiles)//nsplit
        jobs = split_list(inputfiles, 10)
        print ('Found %d input files, dividing them into %d groups of %d '
               'or %d' % (len(inputfiles), nsplit,
                          len(jobs[0]), len(jobs[0])+1))

        # Create temp directory
        if not options.dry:
            call(['mkdir', '-p', path.join(options.tmpDir,sample_name)])

        tasklist = []
        for n,job in enumerate(jobs):
            print 20*'-'
            print 'Group %d (%d files):'%(n+1, len(job))
            for filename in job:
                print filename

            output_filename = '%s_%d.root'%(sample_name,n)
            haddtarget = path.join(options.tmpDir, sample_name,
                                   output_filename)
            print 'merging to', haddtarget
            print 'moving to', path.join(output_dir, output_filename)
            tasklist.append((sample_name, n, job, output_dir, options))

        raw_input("Press Enter to continue...")

        from multiprocessing import Pool
        pool = Pool(min(nsplit,8))
        pool.map(mergeAndMove, tasklist)


        # Clean up:
        print 20*'-'
        print 'Cleaning up tmp dir...'
        if not options.dry:
            call(['rm', '-r', path.join(options.tmpDir,sample_name)])
        else:
            print ['rm', '-r', path.join(options.tmpDir,sample_name)]

        print 20*'='
        print 'DONE'
        print 20*'='
        exit(0)

    else:
        parser.print_help()
        exit(-1)
