#!/bin/bash

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
python /home/wluszczak/ara/secondaries/scripts/generate_evts.py $1 $2
