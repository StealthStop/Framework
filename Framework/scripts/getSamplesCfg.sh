#!/bin/bash 

cfgDir=${CMSSW_BASE}/src/Framework/Framework/cfg
sampleSCfgFile="sampleSets_v2.cfg"
sampleCCfgFile="sampleCollections_v2.cfg"

function print_help {
    echo "Usage:"
    echo "getSamplesCfg.sh -s [SampleSetConfig] -c [SampleCollectionConfg]"
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

ln -s $cfgDir/$sampleSCfgFile sampleSets.cfg
ln -s $cfgDir/$sampleCCfgFile sampleCollections.cfg

echo "Made soft link: "$sampleSCfgFile
echo "Made soft link: "$sampleCCfgFile
