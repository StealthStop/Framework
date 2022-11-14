#!/bin/bash 

cfgDir=${CMSSW_BASE}/src/Framework/Framework/cfg

#Default options
sampleSCfgFile="sampleSets_UL_v2.cfg"
sampleCCfgFile="sampleCollections_UL_v2.cfg"

sampleSCfgFileSkim='sampleSets_UL_v1_skim.cfg'
sampleCCfgFileSkim='sampleCollections_UL_v1_skim.cfg'

function print_help {
    echo "Usage:"
    echo "getSamplesCfg.sh                                                 | Makes default softlinks"
    echo "getSamplesCfg.sh skims                                           | Makes default skim softlinks"
    echo "getSamplesCfg.sh -s [SampleSetConfig] -c [SampleCollectionConfg] | Makes softlinks to given SampleSet and SampleCollection (path already included)"
    echo ""
}

while getopts "h?s:c:" opt; do
    case "$opt" in
    h|\?)
        print_help
        exit 0
        ;;
    s)  sampleSCfgFile=$OPTARG
        ;;
    c)  sampleCCfgFile=$OPTARG
        ;;
    esac
done

subcommand=$1
case "$subcommand" in
    skims)
        sampleSCfgFile=$sampleSCfgFileSkim
        sampleCCfgFile=$sampleCCfgFileSkim
esac

ln -s -f $cfgDir/$sampleSCfgFile sampleSets.cfg
ln -s -f $cfgDir/$sampleCCfgFile sampleCollections.cfg
echo "Made soft link: "$sampleSCfgFile
echo "Made soft link: "$sampleCCfgFile
