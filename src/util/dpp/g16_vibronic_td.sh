#!/bin/sh

FILENAME=$(basename "$0")
echo
echo "----- Entering $FILENAME -----"

rm -f TEMP_CONTROL_ELEC
echo $QC_PROGRAM_ELEC  >> TEMP_CONTROL_ELEC
echo $QC_PROGRAM_VIB  >> TEMP_CONTROL_ELEC
echo $QC_METHOD  >> TEMP_CONTROL_ELEC
echo $PROPERTY   >> TEMP_CONTROL_ELEC

# Input files  
G16_TD_FCHK_FILE=$1
G16_TD_LOG_FILE=$2
G16_FREQ_FCHK_FILE=$3

sh "$DIR_d77_dpp"/g16_elec_fchk.sh $G16_TD_FCHK_FILE
sh "$DIR_d77_dpp"/g16_cicoef_td_log.sh $G16_TD_LOG_FILE
sh "$DIR_d77_dpp"/g16_freq_fchk.sh $G16_FREQ_FCHK_FILE

echo  
echo '----- Leaving '$FILENAME' -----'
echo  
exit
