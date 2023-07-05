#!/bin/bash

export X509_USER_PROXY=$(voms-proxy-info -path)

usage()
{
    echo "Usage:"
    echo ""
    echo "  sh condor/runScoutingPlots_onCondor.sh input_dir output_dir"
    echo ""
    echo "The output_dir will be created in /ceph/cms/store/user/$USER/ScoutingRun3Output/"
    echo "Control the jobs to be run by editing runScoutingPlots_onCondor.sub"
    echo ""
    exit
}

if [ -z $1 ]; then usage; fi

export SCOUTINGINPUTDIR=$1
export SCOUTINGOUTPUTDIR=$2

mkdir -p condor/plotting_logs
mkdir -p /ceph/cms/store/user/$USER/ScoutingRun3Output/$SCOUTINGOUTPUTDIR

sh condor/create_package.sh

condor_submit condor/runScoutingPlots_onCondor.sub
