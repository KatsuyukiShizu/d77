#!/bin/sh
shopt -s nocasematch

FILENAME=$(basename "$0")
echo
echo "----- Entering ${FILENAME} -----"

PROPERTY=PROPERTY_DEFAULT_VALUE  
QC_METHOD=QC_METHOD_DEFAULT_VALUE

# Searching PROPERTY
grep --quiet 'PROPERTY' ${DPPINP}
CHECK_PROPERTY=$(echo $?)
if [ ${CHECK_PROPERTY} = 1 ] ; then # Missing PROPERTY
  echo ${WARNING}
  echo "Missing PROPERTY in ${DPPINP}"
  echo ' '
  echo ${MESSAGE_ERROR_TERM}
  exit
else
  PROPERTY=$(grep PROPERTY ${DPPINP} | awk '{print $3}')
fi

echo "PROPERTY: ${PROPERTY}"

grep --quiet 'QC_METHOD' ${DPPINP}
CHECK_QC_METHOD=$(echo $?)
if [ ${CHECK_QC_METHOD} = 1 ] ; then # Missing QC_METHOD
  echo ${WARNING}
  echo "Missing QC_METHOD in ${DPPINP}"
  echo ' '
  echo ${MESSAGE_ERROR_TERM}
  exit
else
  QC_METHOD=$(grep QC_METHOD ${DPPINP} | awk '{print $3}')
fi

echo "QC_METHOD: ${QC_METHOD}"
echo ${QC_METHOD}  >> TEMP_CONTROL_ELEC  
echo ${PROPERTY}   >> TEMP_CONTROL_ELEC  

DATA_ELEC=$(grep DATA_ELEC ${DPPINP} | awk '{print $3}')
DATA_CICOEF=$(grep DATA_CICOEF ${DPPINP} | awk '{print $3}')
DATA_VIB=$(grep DATA_VIB ${DPPINP} | awk '{print $3}')

case ${QC_METHOD} in

  td)

    case ${PROPERTY} in
  
      dipole | soc)
        sh "${DIR_d77_dpp}"/g16_elec_only_td.sh \
           ${DATA_ELEC} ${DATA_CICOEF}
      ;;

      vc)
        sh "${DIR_d77_dpp}"/g16_vibronic_td.sh \
           ${DATA_ELEC} ${DATA_CICOEF} ${DATA_VIB}
      ;;


      *)
        echo 'Invalid PROPERTY: '${PROPERTY}
        echo ' '
        echo ${MESSAGE_ERROR_TERM}
        exit
      ;;

    esac # ${PROPERTY}

  ;; # td)

  cis)

    case ${PROPERTY} in

      dipole | soc)
        sh "${DIR_d77_dpp}"/g16_elec_only_cis.sh \
           ${DATA_ELEC} ${DATA_CICOEF}
      ;;

      vc)
        sh "${DIR_d77_dpp}"/g16_vibronic_cis.sh \
           ${DATA_ELEC} ${DATA_CICOEF} ${DATA_VIB}
      ;;


      *)
        echo 'Invalid PROPERTY: '${PROPERTY}
        echo ' '
        echo ${MESSAGE_ERROR_TERM}
        exit
      ;;

    esac # ${PROPERTY}

  ;; # cis)

  *)
    echo 'Invalid QC_METHOD: '${QC_METHOD}
    echo ' '
    echo ${MESSAGE_ERROR_TERM}
    exit
  ;; # *)

esac # ${QC_METHOD} 

echo
echo "----- Leaving ${FILENAME} -----"
echo

exit

