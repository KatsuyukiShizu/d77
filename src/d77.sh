#!/bin/sh

# This script is part of d77.

# d77 is free software and can be redistributed and/or modified
# under the terms of the GNU General Public License v3.0
# as published by the Free Software Foundation.
# https://www.gnu.org/licenses/gpl-3.0.html

# For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

export DIRPWD=$(pwd)

export JOB=$(basename $1 .inp)
export RESDIR=$2 

export OUTPUT_d77=${JOB}.log

sh "${DIR_d77}"/run_d77.sh >& ${OUTPUT_d77}
#sed -e s///g $1.log > temp
#mv temp $1.log
#mv $1.log $2 



