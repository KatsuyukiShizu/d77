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
G16_CIS_FCHK_FILE=$1
G16_CIS_LOG_FILE=$2
G16_FREQ_FCHK_FILE=$3

sh "$DIR_d77_dpp"/g16_elec_fchk.sh $G16_CIS_FCHK_FILE
sh "$DIR_d77_dpp"/g16_cicoef_cis_log.sh $G16_CIS_LOG_FILE
sh "$DIR_d77_dpp"/g16_freq_fchk.sh $G16_FREQ_FCHK_FILE

echo
echo '----- Leaving '$FILENAME' -----'
echo
exit
