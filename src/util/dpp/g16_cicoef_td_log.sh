#!/bin/sh

FILENAME=$(basename "$0")
echo
echo "----- Entering $FILENAME -----"

# Input files  
G16_TD_LOG_FILE=$1

grep --quiet 'THRESHOLD_CICOEF' $DPPINP
CHECK_THRESHOLD_CICOEF=$(echo $?)
if [ $CHECK_THRESHOLD_CICOEF = 1 ] ; then
  THRESHOLD_CICOEF='0.0D0'
else
  THRESHOLD_CICOEF=$(grep THRESHOLD_CICOEF $DPPINP | awk '{print $3}')
fi
echo $THRESHOLD_CICOEF > THRESHOLD_VALUES_CI


##########################################################
#  Extracting CI coefficients from G16_TD_LOG_FILE_FILE  #
##########################################################

echo 
echo 'Reading CI coefficients > '$THRESHOLD_CICOEF' from '$G16_TD_LOG_FILE

NUM_ELEC_ALPHA=$(grep 'beta electrons' $G16_TD_LOG_FILE | \
                 awk '{print $1}')
NUM_ELEC_BETA=$(grep 'beta electrons' $G16_TD_LOG_FILE | \
                awk '{print $4}')
NUM_ELEC=$((NUM_ELEC_ALPHA+NUM_ELEC_BETA))
echo $NUM_ELEC > G16_TD_X 

NUM_BASIS=$(grep 'NBasis' $G16_TD_LOG_FILE | awk 'END{print $2}')
echo $NUM_BASIS >> G16_TD_X 

NUM_EX_STATES=$(grep 'Excited State' $G16_TD_LOG_FILE | \
                awk 'END{print NR}')
echo $NUM_EX_STATES >> G16_TD_X 
echo $NUM_EX_STATES > G16_TD_Y 

grep 'Excited State' $G16_TD_LOG_FILE | \
awk '{print $5}' > EX_ENERGY  

awk /'Excited State   1'/,/'SavETr'/ \
$G16_TD_LOG_FILE > TEMP_G16_TD

if [ $NUM_EX_STATES -ge 9 ] ; then
  echo ' Excited State  '$((NUM_EX_STATES+1))\
  >> TEMP_G16_TD
else
  echo ' Excited State   '$((NUM_EX_STATES+1))\
  >> TEMP_G16_TD
fi

if [ $NUM_ELEC_ALPHA -eq $NUM_ELEC_BETA ]; then
  echo 'g16_td_log' > FILE_TYP
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

    awk /"$INITIAL"/,/"$FINAL"/ TEMP_G16_TD > TEMP_CICOEF$i
    grep '[0-9].[0-9][0-9][0-9][0-9][0-9]' TEMP_CICOEF$i \
    > TEMP_CICOEF$i\_2
    grep '>' TEMP_CICOEF$i\_2 > TEMP_X$i
    grep '<' TEMP_CICOEF$i\_2 > TEMP_Y$i

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
    echo $EXMULT >> G16_TD_X
    echo $EXMULT >> G16_TD_Y
    
    NUM_ELEC_CONFIG_X=$(cat TEMP_X$i | wc -l)
    echo $NUM_ELEC_CONFIG_X >> G16_TD_X
    NUM_ELEC_CONFIG_Y=$(cat TEMP_Y$i | wc -l)
    echo $NUM_ELEC_CONFIG_Y >> G16_TD_Y
    cat TEMP_X$i >> G16_TD_X
    cat TEMP_Y$i >> G16_TD_Y
  done
else
  echo 'g16_td_log_unres' > FILE_TYP
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

    awk /"$INITIAL"/,/"$FINAL"/ TEMP_G16_TD > TEMP_CICOEF$i
    grep '[0-9].[0-9][0-9][0-9][0-9][0-9]' TEMP_CICOEF$i \
    > TEMP_CICOEF$i\_2
    grep '>' TEMP_CICOEF$i\_2 > TEMP_X$i
    grep '<' TEMP_CICOEF$i\_2 > TEMP_Y$i
    
    grep A TEMP_X$i > TEMP_X$i\_ALPHA
    NUM_ALPHA_X=$(awk 'END{print NR}' TEMP_X$i\_ALPHA)
    grep B TEMP_X$i > TEMP_X$i\_BETA
    NUM_BETA_X=$(awk 'END{print NR}' TEMP_X$i\_BETA)

    grep A TEMP_Y$i > TEMP_Y$i\_ALPHA
    NUM_ALPHA_Y=$(awk 'END{print NR}' TEMP_Y$i\_ALPHA)
    grep B TEMP_Y$i > TEMP_Y$i\_BETA
    NUM_BETA_Y=$(awk 'END{print NR}' TEMP_Y$i\_BETA)

    TEMP_EXMULT=$(head -1 TEMP_CICOEF$i | awk '{print $4}')
    if [[ $(`echo "$TEMP_EXMULT" | grep "Singlet"`) -eq 1  ]]; then
      EXMULT='Singlet'
    elif [[ $(`echo "$TEMP_EXMULT" | grep "Triplet"`) ]]; then
      EXMULT='Triplet'
    elif [[ $(`echo "$TEMP_EXMULT" | grep "<S**2>=0.750"`) ]]; then
      EXMULT='DOUBLET'
    else
      EXMULT='Unknown'
    fi  
    echo $EXMULT >> G16_TD_X
    echo $EXMULT >> G16_TD_Y
  
    NUM_ELEC_CONFIG_X=$(cat TEMP_X$i | wc -l)
    echo $NUM_ALPHA_X >> G16_TD_X
    cat TEMP_X$i\_ALPHA >> G16_TD_X 
    echo $NUM_BETA_X >> G16_TD_X
    cat TEMP_X$i\_BETA >> G16_TD_X 

    NUM_ELEC_CONFIG_Y=$(cat TEMP_Y$i | wc -l)
    echo $NUM_ALPHA_Y >> G16_TD_Y
    cat TEMP_Y$i\_ALPHA >> G16_TD_Y 
    echo $NUM_BETA_Y >> G16_TD_Y
    cat TEMP_Y$i\_BETA >> G16_TD_Y 

    sed -e 's/A/ /g' G16_TD_X > G16_TD_X_temp
    sed -e 's/B/ /g' G16_TD_X_temp > G16_TD_X_temp_2
    mv G16_TD_X_temp_2 G16_TD_X

    sed -e 's/A/ /g' G16_TD_Y > G16_TD_Y_temp
    sed -e 's/B/ /g' G16_TD_Y_temp > G16_TD_Y_temp_2
    mv G16_TD_Y_temp_2 G16_TD_Y
  done
fi

sed -e 's// /g' G16_TD_X > G16_TD_X_temp
mv G16_TD_X_temp G16_TD_X
sed -e 's/->/  /g' G16_TD_X > G16_TD_X_temp
mv G16_TD_X_temp G16_TD_X

sed -e 's// /g' G16_TD_Y > G16_TD_Y_temp
mv G16_TD_Y_temp G16_TD_Y
sed -e 's/<-/  /g' G16_TD_Y > G16_TD_Y_temp
mv G16_TD_Y_temp G16_TD_Y

# Generating CICOEF
echo 'Generating input data for Density 77'
"$DIR_d77_dpp_f90"/dpp.exe
echo 'Done'
rm -f FILE_TYP

mv G16_TD_X G16_TD_Y XY XY_COUNT EX_ENERGY INP_ELEC

echo
echo "----- Leaving $FILENAME -----"
echo
exit
