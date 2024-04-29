#!/bin/sh

# This script is part of d77.

# d77 is free software and can be redistributed and/or modified
# under the terms of the GNU General Public License v3.0
# as published by the Free Software Foundation.
# https://www.gnu.org/licenses/gpl-3.0.html

# For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp


shopt -s nocasematch

if [ -n "${DIR_PWD}" ]; then
  :
else
  echo
  echo " DIR_PWD directory is not specified."
  echo
  echo " --------------- End of ${EXE} execution script ---------------"
  exit
fi

EXE=d77
echo " --------------- ${EXE} execution script ---------------"

###
#echo 
#echo " Host = `hostname`"
#echo " OS   = `uname`"
#echo " Date = `date`"
#echo 
#echo " d77 is running in directory: $(pwd)"
###

if [ -f "${DIR_PWD}/${JOB}.inp" ] ; then
  RUNTYP=$(grep -i Runtyp "${DIR_PWD}/${JOB}.inp" | awk '{print $3}')
else
  echo
  echo " Missing ${JOB}.inp"
  echo " Removing ${DIR_RES} directory"
  rm -rf "${DIR_PWD}/${DIR_RES}"
  echo
  echo " Error termination"
  echo
  echo " --------------- End of ${EXE} execution script ---------------"
  exit
fi

case ${RUNTYP} in

  Cube)
    :
  ;;

  *)

    if [ ! -d ${DIR_RES} ] ; then
      mkdir ${DIR_RES}
    else
      rm -rf ${DIR_RES}
      mkdir ${DIR_RES}
    fi
    cd ${DIR_RES}
  ;;
esac

cp "${DIR_PWD}/${JOB}.inp" INPUT

PROPERTY=$(grep -i Property INPUT | awk '{print $3}')
METHOD=$(grep -i Method INPUT | awk '{print $3}')
RUNTYP=$(grep -i Runtyp INPUT | awk '{print $3}')
CUBE_OP=$(grep -i Cube_op INPUT | awk '{print $3}')

DIR_INP_ELEC=$(grep -i DIR_INP_ELEC INPUT | awk '{print $3}')
DIR_INP_VIB=$(grep -i DIR_INP_VIB INPUT | awk '{print $3}')
DIR_INP_ELFLD=$(grep -i DIR_INP_ELFLD INPUT | awk '{print $3}')
DIR_INP_SOC_CGF=$(grep -i DIR_INP_SOC_CGF INPUT | awk '{print $3}')

CUBE_RES=$(grep -i CUBE_RES INPUT | awk '{print $3}')
CUBE_1=$(grep -i CUBE_1 INPUT | awk '{print $3}')
CUBE_2=$(grep -i CUBE_2 INPUT | awk '{print $3}')

case ${RUNTYP} in 

  Calc_int_cgf) 

    if [ -n "${DIR_INP_ELEC}" ]; then
      if [ -d "../${DIR_INP_ELEC}" ] ; then
        :
      else  
        echo 
        echo " No ${DIR_INP_ELEC} directory"
        echo 
        echo " Error termination"
        echo 
        echo " --------------- End of ${EXE} execution script ---------------"
        exit
      fi
    else
      echo 
      echo " DIR_INP_ELEC directory is not specified."
      echo 
      echo " --------------- End of ${EXE} execution script ---------------"
      exit
    fi
    ln -s "../${DIR_INP_ELEC}"/* .
  ;; # Calc_int_cgf)

  Density) 

    if [ -n "${DIR_INP_ELEC}" ]; then
      if [ -d "../${DIR_INP_ELEC}" ] ; then
        :
      else  
        echo 
        echo " No ${DIR_INP_ELEC} directory"
        echo
        echo " Error termination"
        echo
        echo " --------------- End of ${EXE} execution script ---------------"
        exit
      fi
    else
      echo 
      echo " DIR_INP_ELEC directory is not specified."
      echo
      echo " Error termination"
      echo
      echo " --------------- End of ${EXE} execution script ---------------"
      exit
    fi
    ln -s "../${DIR_INP_ELEC}"/* .

    case ${PROPERTY} in

      Dipole)
        :
      ;;
      Vc)

        if [ -n "${DIR_INP_VIB}" ]; then
          if [ -d "../${DIR_INP_VIB}" ] ; then
            :
          else
            echo 
            echo " No ${DIR_INP_VIB} directory"
            echo
            echo " Error termination"
            echo
            echo " --------------- End of ${EXE} execution script ---------------"
            exit
          fi
        else
          echo 
          echo " DIR_INP_VIB directory is not specified."
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
        ln -s "../${DIR_INP_VIB}"/* .

      ;;
    esac

  ;;

  Int_pgf) 

    if [ -n "${DIR_INP_ELEC}" ]; then
      if [ -d "../${DIR_INP_ELEC}" ] ; then
        :
      else
        echo 
        echo " No ${DIR_INP_ELEC} directory"
        echo 
        echo " Error termination"
        echo 
        echo " --------------- End of ${EXE} execution script ---------------"
        exit
      fi
    else
      echo 
      echo " DIR_INP_ELEC directory is not specified."
      echo 
      echo " Error termination"
      echo 
      echo " --------------- End of ${EXE} execution script ---------------"
      exit
    fi
    ln -s "../${DIR_INP_ELEC}"/* .

    case ${PROPERTY} in

      Dipole)
        :
      ;;  
      Vc)

        if [ -n "${DIR_INP_VIB}" ]; then
          if [ -d "../${DIR_INP_VIB}" ] ; then
            :
          else
            echo 
            echo " No ${DIR_INP_VIB} directory"
            echo
            echo " Error termination"
            echo
            echo " --------------- End of ${EXE} execution script ---------------"
            exit
          fi
        else
          echo 
          echo " DIR_INP_VIB directory is not specified."
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
        ln -s "../${DIR_INP_VIB}"/* .

        if [ -n "${DIR_INP_ELFLD}" ]; then
          if [ -d "../${DIR_INP_ELFLD}" ] ; then
            :
          else
            echo 
            echo " No ${DIR_INP_ELFLD} directory"
            echo
            echo " Error termination"
            echo
            echo " --------------- End of ${EXE} execution script ---------------"
            exit
          fi
        else
          echo 
          echo " DIR_INP_ELEC directory is not specified."
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
        ln -s "../${DIR_INP_ELFLD}"/ELFLD_CGF_ATM .

      ;;
      Soc)

        if [ -n "${DIR_INP_SOC_CGF}" ]; then
          if [ -d "../${DIR_INP_SOC_CGF}" ]; then
            :
          else
            echo
            echo " No ${DIR_INP_SOC_CGF} directory"
            echo
            echo " Error termination"
            echo
            echo " --------------- End of ${EXE} execution script ---------------"
            exit
          fi
        else
          echo 
          echo " DIR_INP_SOC_CGF directory is not specified."
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
        ln -s "../${DIR_INP_SOC_CGF}"/SOC_CGF_X .
        ln -s "../${DIR_INP_SOC_CGF}"/SOC_CGF_Y .
        ln -s "../${DIR_INP_SOC_CGF}"/SOC_CGF_Z .
      ;;

      *)
        echo 
        echo "Invalid PROPERTY: '${PROPERTY}"
        exit
      ;;
    esac 
  ;; # Int_pgf)
  
  Cube) 
    
    if [ -n "${CUBE_RES}" ]; then
      rm -f CUBE_RES ${CUBE_RES}
      echo ${CUBE_RES} > CUBE_FILE_NAME
      
      if [ -n "${CUBE_1}" ]; then

        if [ -f "${CUBE_1}" ]; then
          rm -f CUBE_1
          ln -s ${CUBE_1} CUBE_1
          echo ${CUBE_1} >> CUBE_FILE_NAME
        else
          echo
          echo " No ${CUBE_1} file"
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
 
      else
        echo
        echo " CUBE_1 file is not specified."
        echo
        echo " Error termination"
        echo
        echo " --------------- End of ${EXE} execution script ---------------"
        exit
      fi
   
    else
      echo
      echo " CUBE_RES file is not specified."
      echo
      echo " Error termination"
      echo
      echo " --------------- End of ${EXE} execution script ---------------"
      exit
    fi
  
    case ${CUBE_OP} in
    
      Add | Sub | Mul)

        if [ -n "${CUBE_2}" ]; then
  
          if [ -f "${CUBE_2}" ]; then
            rm -f CUBE_2
            ln -s ${CUBE_2} CUBE_2
            echo ${CUBE_2} >> CUBE_FILE_NAME
          else
            echo
            echo " No ${CUBE_2} file"
            echo
            echo " Error termination"
            echo
            echo " --------------- End of ${EXE} execution script ---------------"
            exit
          fi
  
        else
          echo
          echo " CUBE_2 file is not specified."
          echo
          echo " Error termination"
          echo
          echo " --------------- End of ${EXE} execution script ---------------"
          exit
        fi
      ;;

      *)
        :
      ;;
    esac  
  ;;

  *) 
    echo  
    echo "Invalid Runtyp: ${RUNTYP}" 
    exit
  ;;
esac


###
#echo
#echo ' INPUT file       = '${DIR_PWD}/${JOB}.inp
#echo ' RESULTS directry = '${DIR_PWD}/${DIR_RES}
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
  
  Cube)

    mv CUBE_RES ${CUBE_RES}
    
    array=(\
    CUBE_FILE_NAME CUBE_1 CUBE_2 ${JOB}.log
     \
    )

    for FILENAME in "${array[@]}"
    do
      if [ -f "${FILENAME}" ]; then
        rm -f "${FILENAME}"
      else
        :
      fi
    done

  ;;
  *)
    :
esac

mv INPUT ${JOB}.inp 

case ${RUNTYP} in

  Cube)
    :
  ;;

  *)
    mv ../${JOB}.log . 
  ;;
esac

echo 
echo " --------------- End of ${EXE} execution script ---------------"
exit







