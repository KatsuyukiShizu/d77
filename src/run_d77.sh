#!/bin/sh

# This script is part of d77.

# d77 is free software and can be redistributed and/or modified
# under the terms of the GNU General Public License v3.0
# as published by the Free Software Foundation.
# https://www.gnu.org/licenses/gpl-3.0.html

# For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp


shopt -s nocasematch

EXE=d77

if [ ! -d ${RESDIR} ] ; then
  mkdir ${RESDIR}
else
  rm -rf ${RESDIR}
  mkdir ${RESDIR}
fi
cd ${RESDIR}

echo 
echo " --------------- ${EXE} execution script ---------------"

###
echo 
echo " Host = `hostname`"
echo " OS   = `uname`"
echo " Date = `date`"
echo 
echo " d77 is running in directory: $(pwd)"
###

export INPUT_d77=${JOB}.inp
if [ -f "${DIRPWD}/${JOB}.inp" ] ; then
  cp "${DIRPWD}/${JOB}.inp" INPUT
else
  echo " Missing ${JOB}.inp"
  echo " Removing ${RESDIR} directory"
  rm -rf "${DIRPWD}/${RESDIR}"
  echo 
  echo " Error termination"
  echo 
  echo " --------------- End of ${EXE} execution script ---------------"
  exit
fi

PROPERTY=$(grep -i Property INPUT | awk '{print $3}')
METHOD=$(grep -i Method INPUT | awk '{print $3}')
RUNTYP=$(grep -i Runtyp INPUT | awk '{print $3}')

INPDIR_ELEC=$(grep -i INPDIR_ELEC INPUT | awk '{print $3}')
INPDIR_VIB=$(grep -i INPDIR_VIB INPUT | awk '{print $3}')
INPDIR_ELFLD=$(grep -i INPDIR_ELFLD INPUT | awk '{print $3}')
INPDIR_SOC_CGF=$(grep -i INPDIR_SOC_CGF INPUT | awk '{print $3}')


case ${RUNTYP} in 

  Calc_int_cgf) 

    if [ -n "${INPDIR_ELEC}" ]; then
      if [ ! -d "${INPDIR_ELEC}" ] ; then
        echo " No ${INPDIR_ELEC} directory"
        echo 
        echo " Error termination"
        echo 
        echo " --------------- End of ${EXE} execution script ---------------"
        exit
      else
        :
      fi
    else
      echo " INPDIR_ELEC is empty"
      echo 
      echo " --------------- End of ${EXE} execution script ---------------"
      exit
    fi
    ln -s $INPDIR_ELEC/* .
  ;; # Calc_int_cgf)

  Density) 

    if [ -n "${INPDIR_ELEC}" ]; then
      if [ ! -d "${INPDIR_ELEC}" ] ; then
        echo " No ${INPDIR_ELEC} directory"
        echo
        echo " Error termination"
        echo
        echo " --------------- End of ${EXE} execution script ---------------"
        exit
      else
        :
      fi
    else
      echo " INPDIR_ELEC is empty"
      echo
      echo " Error termination"
      echo
      echo " --------------- End of ${EXE} execution script ---------------"
      exit
    fi
    ln -s $INPDIR_ELEC/* .

    case $PROPERTY in

      Dipole)
        :
      ;;
      Vc)

        if [ -n "${INPDIR_VIB}" ]; then
          if [ ! -d "${INPDIR_VIB}" ] ; then
            echo " No ${INPDIR_VIB} directory"
            echo
            echo " Error termination"
            echo
            echo " --------------- End of ${EXE} execution script ---------------"
            exit
          else
            :
          fi
        else
          echo " INPDIR_VIB is empty"
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
        ln -s $INPDIR_VIB/* .

      ;;
    esac

  ;;

  Int_pgf) 

    if [ -n "${INPDIR_ELEC}" ]; then
      if [ ! -d "${INPDIR_ELEC}" ] ; then
        echo " No ${INPDIR_ELEC} directory"
        echo 
        echo " Error termination"
        echo 
        echo " --------------- End of ${EXE} execution script ---------------"
        exit
      else
        :
      fi
    else
      echo " INPDIR_ELEC is empty"
      echo 
      echo " Error termination"
      echo 
      echo " --------------- End of ${EXE} execution script ---------------"
      exit
    fi
    ln -s $INPDIR_ELEC/* .

    case $PROPERTY in

      Dipole)
        :
      ;;  
      Vc)

        if [ -n "${INPDIR_VIB}" ]; then
          if [ ! -d "${INPDIR_VIB}" ] ; then
            echo " No ${INPDIR_VIB} directory"
            echo
            echo " Error termination"
            echo
            echo " --------------- End of ${EXE} execution script ---------------"
            exit
          else
            :
          fi
        else
          echo " INPDIR_VIB is empty"
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
        ln -s $INPDIR_VIB/* .

        if [ -n "${INPDIR_ELFLD}" ]; then
          if [ ! -d "${INPDIR_ELFLD}" ] ; then
            echo " No ${INPDIR_ELFLD} directory"
            echo
            echo " Error termination"
            echo
            echo " --------------- End of ${EXE} execution script ---------------"
            exit
          else
            :
          fi
        else
          echo " INPDIR_ELFLD is empty"
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
        ln -s ${INPDIR_ELFLD}/ELFLD_CGF_ATM .

      ;;
      Soc)

        if [ -n "${INPDIR_SOC_CGF}" ]; then
          :
        else
          echo 'No INPDIR_SOC_CGF directry'
          exit
        fi
        ln -s ${INPDIR_SOC_CGF}/SOC_CGF_X .
        ln -s ${INPDIR_SOC_CGF}/SOC_CGF_Y .
        ln -s ${INPDIR_SOC_CGF}/SOC_CGF_Z .
      ;;

      *)
        echo "Invalid PROPERTY: '${PROPERTY}"
        exit
      ;;
    esac 
  ;; # Int_pgf)
  
  *) 
    echo "Invalid Runtyp: ${RUNTYP}" 
    exit
  ;;
esac


###
echo
echo ' INPUT file       = '${DIRPWD}/${JOB}.inp
echo ' RESULTS directry = '${DIRPWD}/${RESDIR}
###

echo
echo
echo
"${DIR_d77_EXE}/d77.exe"

array=(\
ATMNUM CNTEXP SHELLTYP ATMORB COEFAREF NUCCHARG ZEFF_SOC ATMWT COEFBREF \
NUMPGFSHELL CICOEF EX_ENERGY CONTROL_ELEC QCHEMCI G16_CIS \
XY G16_TD_X G16_TD_Y XY_COUNT \
CNTCOEF0 CONTROL_VIB FREQ NAC \
VIBMODE CNTCOEFP COORD HESSIAN SHELL2ATM \
)

for FILENAME in "${array[@]}"
do
  if [ -f "${FILENAME}" ]; then
    unlink "${FILENAME}"
  else
    :
  fi
done

case ${RUNTYP} in
  Int_pgf)

    array=(\
    ELFLD_CGF_ATM SOC_CGF_X SOC_CGF_Y SOC_CGF_Z \
     \
    )
    
    for FILENAME in "${array[@]}"
    do
      if [ -f "${FILENAME}" ]; then
        unlink "${FILENAME}"
      else
        :
      fi
    done
  ;;
  *)
    :
esac

mv INPUT ${JOB}.inp 
mv ../${JOB}.log . 
echo 
echo " --------------- End of ${EXE} execution script ---------------"
exit







