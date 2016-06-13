#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

cd $2/src
eval `scramv1 runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}


xrdcp root://cmseos.fnal.gov/$(echo $5 | sed 's|/eos/uscms||') .

./CS $1 -1 $3 $4 "condor"

rm $(echo $5 | sed 's|.*/||')

