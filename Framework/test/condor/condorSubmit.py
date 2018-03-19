#!/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_3_3/external/slc6_amd64_gcc630/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

import sys, os
from os import system, environ
sys.path = [environ["CMSSW_BASE"] + "/src/Exploration/condor/",] + sys.path

from samples import SampleCollection
import optparse 
import subprocess

# Parse command line arguments
parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('-n',  dest='numfile', type='int', default = 5, help="number of files per job")
parser.add_option ('-d',  dest='datasets', type='string', default = '', help="List of datasets, comma separated")
parser.add_option ('-l',  dest='dataCollections', action='store_true', default = False, help="List all datacollections")
parser.add_option ('-L',  dest='dataCollectionslong', action='store_true', default = False, help="List all datacollections and sub collections")
parser.add_option ('-c',  dest='noSubmit', action='store_true', default = False, help="Do not submit jobs.  Only create condor_submit.txt.")
parser.add_option ('--explore',  dest='explore', default = 'b', help="ExploreTopTagger (t), ExploreBackground (b), ExploreEventSelection (s)")

options, args = parser.parse_args()

# Prepare the list of files to transfer
mvaFileName = ""
with file(environ["CMSSW_BASE"] + "/src/Exploration/TopTagger.cfg") as meowttcfgFile:
    for line in meowttcfgFile:
        if "modelFile" in line:
            mvaFileName = line.split("=")[1].strip().strip("\"")
            break

filestoTransfer = [environ["CMSSW_BASE"] + "/src/Exploration/MyAnalysis", 
                   environ["CMSSW_BASE"] + "/src/Exploration/%s" % mvaFileName,
                   environ["CMSSW_BASE"] + "/src/Exploration/TopTagger.cfg",
                   environ["CMSSW_BASE"] + "/src/TopTagger/TopTagger/test/libTopTagger.so"
                   ]

def makeExeAndFriendsTarball(filestoTransfer, fname):
    system("mkdir -p %s" % fname)
    for fn in filestoTransfer:
        system("cd %s; ln -s %s" % (fname, fn))
        
    tarallinputs = "tar czvf %s.tar.gz %s --dereference"% (fname, fname)
    print tarallinputs
    system(tarallinputs)
    system("rm -r %s" % fname)


if not options.dataCollections and not options.dataCollectionslong:
    makeExeAndFriendsTarball(filestoTransfer, "exestuff")
    system("tar --exclude-caches-all --exclude-vcs -zcf ${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=src --exclude=tmp")

submitFile = """Universe   = vanilla
Executable = run_Exploration_condor.tcsh
Transfer_Input_Files = CMSSW_9_3_3.tar.gz, exestuff.tar.gz
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT

"""

sc = SampleCollection()
if options.dataCollections or options.dataCollectionslong:
    scl = sc.sampleCollectionList()
    for sampleCollection in scl:
        sl = sc.sampleList(sampleCollection)
        print sampleCollection
        if options.dataCollectionslong:
            sys.stdout.write("\t")
            for sample in sl:
                sys.stdout.write("%s  "%sample[1])
            print ""
            print ""
    exit(0)

datasets = []
if options.datasets:
    datasets = options.datasets.split(',')
else:
    print "No dataset specified"
    exit(0)

fileParts = [submitFile]
nFilesPerJob = options.numfile
for ds in datasets:
    ds = ds.strip()
    print ds
    # create the directory
    if not os.path.isdir("output-files/%s" % ds):
        os.makedirs("output-files/%s" % ds)

    for s, n in sc.sampleList(ds):
        print "s:", s, ", n:", n
        print "\t%s"%n
        f = open(s)
        if not f == None:
            count = 0
            for l in f:
                if '.root' in l:
                    count = count + 1
            for startFileNum in xrange(0, count, nFilesPerJob):
                fileParts.append("transfer_output_remaps = \"MyAnalysis_%s_%s.root = output-files/%s/MyAnalysis_%s_%s.root\"\n" % (n, startFileNum, ds, n, startFileNum))
                fileParts.append("Arguments = %s %i %i %s %s\n"%(n, nFilesPerJob, startFileNum, s, options.explore))
                fileParts.append("Output = log-files/MyAnalysis_%s_%i.stdout\n"%(n, startFileNum))
                fileParts.append("Error = log-files/MyAnalysis_%s_%i.stderr\n"%(n, startFileNum))
                fileParts.append("Log = log-files/MyAnalysis_%s_%i.log\n"%(n, startFileNum))
                fileParts.append("Queue\n\n")

            f.close()

fout = open("condor_submit.txt", "w")
fout.write(''.join(fileParts))
fout.close()

if not options.noSubmit: 
    system('mkdir -p log-files')
    system("echo 'condor_submit condor_submit.txt'")
    system('condor_submit condor_submit.txt')

