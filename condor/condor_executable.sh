#!/bin/bash

SCRAMARCH=slc7_amd64_gcc10
CMSSWVERSION=CMSSW_12_6_0_patch1

INDIR=$1
OUTDIR=$2
INDEX=$3
PACE=$4

function stageout {
    COPY_SRC=$1
    COPY_DEST=$2
    retries=0
    COPY_STATUS=1
    until [ $retries -ge 10 ]
    do
        echo "Stageout attempt $((retries+1)): env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}"
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}
        COPY_STATUS=$?
        if [ $COPY_STATUS -ne 0 ]; then
            echo "Failed stageout attempt $((retries+1))"
        else
            echo "Successful stageout with $retries retries"
            break
        fi
        retries=$[$retries+1]
        echo "Sleeping for 5m"
        sleep 5m
    done
    if [ $COPY_STATUS -ne 0 ]; then
        echo "Removing output file because gfal-copy crashed with code $COPY_STATUS"
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-rm --verbose ${COPY_DEST}
        REMOVE_STATUS=$?
        if [ $REMOVE_STATUS -ne 0 ]; then
            echo "Uhh, gfal-copy crashed and then the gfal-rm also crashed with code $REMOVE_STATUS"
            echo "You probably have a corrupt file sitting on hadoop now."
            exit 1
        fi
    fi
}

ulimit -s unlimited
source /cvmfs/cms.cern.ch/cmsset_default.sh

tar xvf package.tar.gz
cd ScoutingRun3/
export SCRAM_ARCH=${SCRAMARCH} && scramv1 project CMSSW ${CMSSWVERSION}
cd ${CMSSWVERSION}/src
cmsenv

cd ../.. # Get back to ScoutingRun3/
echo $INDIR $INDEX $PACE $OUTDIR
python3 plotScouting.py $INDIR $INDEX $PACE $OUTDIR

for FILE in $(ls $OUTDIR);
do
  echo "File $FILE to be copied..."
  echo ""
  COPY_SRC="file://`pwd`/${OUTDIR}/$FILE"
  COPY_DEST="davs://redirector.t2.ucsd.edu:1095/store/user/$USER/ScoutingRun3Output/${OUTDIR}/$FILE"
  stageout $COPY_SRC $COPY_DEST
done;
