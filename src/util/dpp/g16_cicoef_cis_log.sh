#!/bin/sh

FILENAME=$(basename "$0")
echo
echo "----- Entering $FILENAME -----"

# Input files  
G16_CIS_LOG_FILE=$1
THRESHOLD_CI_COEF=$2


#################################################
#  Extracting CI coefficients from G16LOG_FILE  #
#################################################

echo
echo 'Reading CI coefficients from '$G16_CIS_LOG_FILE
echo 'g16_log' > FILE_TYP

grep --quiet 'THRESHOLD_CICOEF' $DPPINP
CHECK_THRESHOLD_CICOEF=$(echo $?)
if [ $CHECK_THRESHOLD_CICOEF = 1 ] ; then
  THRESHOLD_CICOEF='0.0D0'
else
  THRESHOLD_CICOEF=$(grep THRESHOLD_CICOEF $DPPINP | awk '{print $3}')
fi
echo $THRESHOLD_CICOEF > THRESHOLD_VALUES_CI
echo 'THRESHOLD_CICOEF: '$THRESHOLD_CICOEF 


NUM_ELEC_ALPHA=$(grep 'beta electrons' $G16_CIS_LOG_FILE | \
                 awk '{print $1}')
NUM_ELEC_BETA=$(grep 'beta electrons' $G16_CIS_LOG_FILE | \
                awk '{print $4}')
NUM_ELEC=$((NUM_ELEC_ALPHA+NUM_ELEC_BETA))
echo $NUM_ELEC > G16_CIS 

NUM_BASIS=$(grep 'NBasis' $G16_CIS_LOG_FILE | awk 'END{print $2}')
echo $NUM_BASIS >> G16_CIS 

NUM_EX_STATES=$(grep 'Excited State' $G16_CIS_LOG_FILE | \
                awk 'END{print NR}')
echo $NUM_EX_STATES >> G16_CIS 

grep 'Excited State' $G16_CIS_LOG_FILE | \
awk '{print $5}' > EX_ENERGY  

awk /'Excited State   1'/,/'Copying SCF densities'/ \
$G16_CIS_LOG_FILE > TEMP_G16_CIS

if [ $NUM_EX_STATES -ge 9 ] ; then
  echo ' Excited State  '$((NUM_EX_STATES+1))\
  >> TEMP_G16_CIS
else
  echo ' Excited State   '$((NUM_EX_STATES+1))\
  >> TEMP_G16_CIS
fi


for i in `seq 1 $NUM_EX_STATES`
do
  if [ $i -ge 10 ] ; then
    INITIAL=' Excited State  '$i
  else
    INITIAL=' Excited State   '$i
  fi

  if [ $((i+1)) -ge 10 ] ; then
    FINAL=' Excited State  '$((i+1))
  else
    FINAL=' Excited State   '$((i+1))
  fi

  awk /"$INITIAL"/,/"$FINAL"/ TEMP_G16_CIS > TEMP_CICOEF$i
  grep '[0-9].[0-9][0-9][0-9][0-9][0-9]' TEMP_CICOEF$i | \
  grep -v 'State' | grep -v 'Total' | grep -v 'Leave' | \
  sed -e 's/A/ /g' > TEMP_CICOEF$i\_2

  TEMP_EXMULT=$(head -1 TEMP_CICOEF$i | awk '{print $4}')
  if [ `echo $TEMP_EXMULT | grep 'Singlet'` ]; then
    EXMULT='Singlet'
  elif
   [ `echo $TEMP_EXMULT | grep 'Triplet'` ]; then
    EXMULT='Triplet'
  elif
   [ `echo $TEMP_EXMULT | grep '<S**2>=0.750'` ]; then
    EXMULT='DOUBLET'
  else
    EXMULT='Unknown'
  fi  
  echo $EXMULT >> G16_CIS

  NUM_ELEC_CONFIG=$(cat TEMP_CICOEF$i\_2 | wc -l)
  echo $NUM_ELEC_CONFIG >> G16_CIS 
  cat TEMP_CICOEF$i\_2 >> G16_CIS 
done

echo 'Generating input data for Density 77'
"$DIR_d77_dpp_f90"/dpp.exe
echo 'Done'
rm -f FILE_TYP

mv G16_CIS CICOEF EX_ENERGY INP_ELEC

echo
echo "----- Leaving $FILENAME -----"
echo
exit
