#!/bin/sh

FILENAME=$(basename "$0")
echo  
echo "----- Entering $FILENAME -----"

# Input file  
G16_ELEC_FCHK_FILE=$1


#######################################################################
#  Extracting ground-state electronic strcure from G16_ELEC_FCHK_FILE  #
#######################################################################

echo
echo 'Reading ground-state electronic strcure from '$G16_ELEC_FCHK_FILE
echo 'elec_fchk' > FILE_TYP

NUMATM=$(grep 'Number of atoms' $G16_ELEC_FCHK_FILE | \
         awk 'END{print}' | awk '{print $NF}')
echo $NUMATM >> TEMP_CONTROL_ELEC 'Number of atoms' 

NUMCOORD=$(grep 'Current cartesian coordinates' $G16_ELEC_FCHK_FILE | \
           awk 'END{print}' | awk '{print $NF}')
echo $NUMCOORD >> TEMP_CONTROL_ELEC 'Current cartesian coordinates'

CHARGREF=$(awk /'Number of atoms'/,/'Atomic numbers'/ $G16_ELEC_FCHK_FILE | \
           grep 'Charge' | \
           awk 'END{print}' | awk '{print $NF}')
echo $CHARGREF >> TEMP_CONTROL_ELEC 'Charge in reference system'

MULTREF=$(awk /'Number of atoms'/,/'Atomic numbers'/ $G16_ELEC_FCHK_FILE | \
          grep 'Multiplicity' | \
          awk 'END{print}' | awk '{print $NF}')
echo $MULTREF >> TEMP_CONTROL_ELEC \
'Multiplicity in reference system'

NUMELECREF=$(grep 'Number of electrons' $G16_ELEC_FCHK_FILE | \
               awk 'END{print}' | awk '{print $NF}')
  echo $NUMELECREF >> TEMP_CONTROL_ELEC \
  'Number of electrons in reference system'

NUMELECAREF=$(grep 'Number of alpha electrons' $G16_ELEC_FCHK_FILE | \
              awk 'END{print}' | awk '{print $NF}')
echo $NUMELECAREF >> TEMP_CONTROL_ELEC \
'Number of alpha electrons in reference system'

NUMELECBREF=$(grep 'Number of beta electrons' $G16_ELEC_FCHK_FILE | \
              awk 'END{print}' | awk '{print $NF}')
echo $NUMELECBREF >> TEMP_CONTROL_ELEC \
'Number of beta electrons in reference system'

NUMBASIS0=$(grep 'Number of basis functions' $G16_ELEC_FCHK_FILE | \
           awk 'END{print}' | awk '{print $NF}')
echo $NUMBASIS0 >> TEMP_CONTROL_ELEC \
'Number of basis functions'

NUMBASIS=$(grep 'Number of independent functions' $G16_ELEC_FCHK_FILE | \
           awk 'END{print}' | awk '{print $NF}')
echo $NUMBASIS >> TEMP_CONTROL_ELEC \
'Number of independent basis functions'

NUMSHELL=$(grep 'Number of contracted shells' $G16_ELEC_FCHK_FILE | \
           awk 'END{print}' | awk '{print $NF}')
echo $NUMSHELL >> TEMP_CONTROL_ELEC \
'Number of contracted shells'

grep --quiet P\(S=P\) $G16_ELEC_FCHK_FILE 
SP=$(echo $?)
echo $SP >> TEMP_CONTROL_ELEC \
'SP-type shells: 0 = Yes; 1 = No.'
NUMPRMEXP=$(grep 'Primitive exponents' $G16_ELEC_FCHK_FILE | \
            awk 'END{print}' | awk '{print $NF}')
echo $NUMPRMEXP >> TEMP_CONTROL_ELEC \
'Number of primitive exponents'

grep --quiet 'Beta MO coefficients' $G16_ELEC_FCHK_FILE
BETAMOREF=$(echo $?)
echo $BETAMOREF >> TEMP_CONTROL_ELEC \
'Beta MOs in reference system: 0 = Yes; 1 = No.'
if [ $BETAMOREF = 1 ] ; then
  echo 'No beta MO coefficients were detected'
  NUMCOEFAREF=$(grep 'Alpha MO coefficients' $G16_ELEC_FCHK_FILE | \
                awk 'END{print}' | awk '{print $NF}')
  echo $NUMCOEFAREF >> TEMP_CONTROL_ELEC\
  'Number of alpha MO coefficients in reference system'
else 
  echo 'Beta MO coefficients were detected'
  NUMCOEFAREF=$(grep 'Alpha MO coefficients' $G16_ELEC_FCHK_FILE | \
                awk 'END{print}' | awk '{print $NF}')
  echo $NUMCOEFAREF >> TEMP_CONTROL_ELEC\
  'Number of alpha MO coefficients in reference system'
  NUMCOEFBREF=$(grep 'Beta MO coefficients' $G16_ELEC_FCHK_FILE | \
                awk 'END{print}' | awk '{print $NF}')
  echo $NUMCOEFBREF >> TEMP_CONTROL_ELEC\
  'Number of beta MO coefficients in reference system'
fi

awk /'Atomic numbers'/,/'Nuclear charges'/ $G16_ELEC_FCHK_FILE |\
grep -v 'N' > TEMP_ATMNUM

  awk /'Nuclear charges'/,/'Current cartesian coordinates'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_NUCCHARG

  awk /'Current cartesian coordinates'/,/'Number of symbols'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_COORD

  awk /'Nuclear ZEff'/,/'Nuclear ZNuc'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_ZEFF_SOC

  awk /'Shell types'/,/'Number of primitives per shell'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_SHELLTYP

  awk /'Number of primitives per shell'/,/'Shell to atom map'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_NUMPGFSHELL

  awk /'Shell to atom map'/,/'Primitive exponents'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_SHELL2ATM

  awk /'Primitive exponents'/,/'Contraction coefficients'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_CNTEXP

if [ $SP = 1 ] ; then
  echo 'No SP-type shells were detected'
  awk /'Contraction coefficients                   R'/,\
/'Coordinates of each shell'/ $G16_ELEC_FCHK_FILE | \
  grep -v 'N' > TEMP_CNTCOEF0
else
  echo 'SP-type shells were detected'
  awk /'Contraction coefficients                   R'/,\
/'P\(S=P\) Contraction coefficients'/ $G16_ELEC_FCHK_FILE | \
  grep -v 'N' > TEMP_CNTCOEF0
  awk /'P\(S=P\) Contraction coefficients'/,\
/'Coordinates of each shell'/ $G16_ELEC_FCHK_FILE | \
  grep -v 'N' > TEMP_CNTCOEFP
fi


if [ $BETAMOREF = 1 ] ; then
  awk /'Alpha MO coefficients'/,/'Orthonormal basis'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_COEFAREF
else
  awk /'Alpha MO coefficients'/,/'Beta MO coefficients'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_COEFAREF
  awk /'Beta MO coefficients'/,/'Orthonormal basis'/ $G16_ELEC_FCHK_FILE |\
  grep -v 'N' > TEMP_COEFBREF
fi

echo 'Generating input data for Density 77'
"$DIR_d77_dpp_f90"/dpp.exe
echo 'Done'
rm -f FILE_TYP


if [ -e COEFBREF ] ; then
  :
else
  cp COEFAREF COEFBREF
fi

mv NUCCHARG COORD ZEFF_SOC \
   SHELLTYP NUMPGFSHELL SHELL2ATM \
   CNTEXP CNTCOEF* ATMORB \
   COEFAREF COEFBREF\
   INP_ELEC

mv ATMNUM CONTROL* \
   INP_ELEC

echo  
echo "----- Leaving $FILENAME -----"
echo  
exit
