#!/bin/sh

export DIRPWD=$(pwd)

export JOB=$(basename $1 .inp)
export RESDIR=$2 

export OUTPUT_d77=${JOB}.log

sh "${DIR_d77}"/run_d77.sh >& ${OUTPUT_d77}
#sed -e s///g $1.log > temp
#mv temp $1.log
#mv $1.log $2 



