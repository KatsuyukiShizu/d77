#!/bin/sh

FILENAME=$(basename "$0")
echo
echo "----- Entering $FILENAME -----"

# Input file
G16_FREQ_FCHK_FILE=$1
#echo 'freqfchk'

#################################################################
#  Extracting Cartesian Hessian matrix from G16_FREQ_FCHK_FILE  #
#################################################################
echo
echo 'Reading frequencies and vibrational modes from '$G16_FREQ_FCHK_FILE

grep --quiet 'Number of Normal Modes' $G16_FREQ_FCHK_FILE
VIBMODE=$(echo $?)
echo $VIBMODE >> TEMP_CONTROL_VIB 'VIBMODE: 0 = Yes; 1 = No.'
if [ $VIBMODE = 1 ] ; then
  echo 'No vibrational modes were detected in '$G16_FREQ_FCHK_FILE
  exit
else
  echo 'Vibrational modes were detected in '$G16_FREQ_FCHK_FILE

  NUMATM=$(grep 'Number of atoms' $G16_FREQ_FCHK_FILE | awk '{print $5}') #20200225
  echo $NUMATM >> TEMP_CONTROL_VIB 'Number of atoms' 

  NUMMODE=$(grep 'Number of Normal Modes' $G16_FREQ_FCHK_FILE | \
            awk 'END{print}' | awk '{print $NF}')
  echo $NUMMODE >> TEMP_CONTROL_VIB \
  'Number of normal modes in reference system'

  NUMCOORDMODE=$(grep 'Vib-Modes' $G16_FREQ_FCHK_FILE | \
                 awk 'END{print}' | awk '{print $NF}')
  echo $NUMCOORDMODE >> TEMP_CONTROL_VIB\
  'Number of elements of vibrational vectors in reference system'
fi

awk /'Vib-E2'/,/'Vib-Modes'/ $G16_FREQ_FCHK_FILE |\
grep -v 'N' > TEMP_FREQ

awk /'Vib-Modes'/,/'ClPar MaxAn'/ $G16_FREQ_FCHK_FILE |\
grep -v 'N' | grep -v 'I' > TEMP_VIBMODE

awk /'Real atomic weights'/,/'Atom fragment'/ $G16_FREQ_FCHK_FILE |\
grep -v 'N' > TEMP_ATMWT

echo 'vib_fchk' > FILE_TYP
"$DIR_d77_dpp_f90"/dpp.exe 

if [ -d INP_VIB ] ; then
  rm -rf INP_VIB
  mkdir INP_VIB
else
  mkdir INP_VIB
fi

mv ATMWT FREQ VIBMODE \
   CONTROL* \
   INP_VIB

echo
echo "----- Leaving $FILENAME -----"
echo
exit
