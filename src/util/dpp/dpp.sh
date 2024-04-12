#!/bin/sh
shopt -s nocasematch

export DPPINP='dpp.inp' # Input file
export MESSAGE_ERROR_TERM='----- dpp process was terminated abnormally -----'
export WARNING='--- Warning ---'
export PROPERTY='PROPERTY_DEFAULT_VALUE'
export QC_METHOD='QC_METHOD_DEFAULT_VALUE'
export QC_PROGRAM_ELEC='QC_PROGRAM_ELEC_DEFAULT_VALUE'
export QC_PROGRAM_VIB='QC_PROGRAM_VIB_DEFAULT_VALUE'

echo  
echo '-----------------------------------' 
echo 'Starting data preprocessing process' 
echo '-----------------------------------' 

if [ -f ${DPPINP} ] ; then
  :
else
  echo ${WARNING}
  echo "Missing ${DPPINP}"
  echo  
  echo ${MESSAGE_ERROR_TERM} 
  exit
fi

# Removing existing temporary files
echo  
echo 'Removing the existing temporary files' 
echo ' ' 
rm -f ${DPPOUT} 
sh "${DIR_d77_dpp}/rm_temp_files.sh"


grep --quiet 'PROPERTY' ${DPPINP}
CHECK_PROPERTY=$(echo $?)
if [ ${CHECK_PROPERTY} = 1 ] ; then
  echo ${WARNING}
  echo "Missing PROPERTY in ${DPPINP}"
  echo ' '
  echo ${MESSAGE_ERROR_TERM}
  exit
else
  export PROPERTY=$(grep PROPERTY ${DPPINP} | awk '{print $3}')
fi


echo
echo 'Removing the existing INP_ELEC directory'
echo 'Creating new INP_ELEC directory'
if [ -d INP_ELEC ] ; then
  rm -rf INP_ELEC
  mkdir INP_ELEC
else
  mkdir INP_ELEC
fi


grep --quiet 'QC_METHOD' ${DPPINP}
CHECK_QC_METHOD=$(echo $?)
if [ ${CHECK_QC_METHOD} = 1 ] ; then
  echo ${WARNING}
  echo 'Missing QC_METHOD in '${DPPINP}
  echo ' '
  echo ${MESSAGE_ERROR_TERM}
  exit
else
  export QC_METHOD=$(grep QC_METHOD ${DPPINP} | awk '{print $3}')
fi

grep --quiet 'QC_PROGRAM_ELEC' ${DPPINP}
CHECK_QC_PROGRAM_ELEC=$(echo $?)
if [ ${CHECK_QC_PROGRAM_ELEC} = 1 ] ; then # Missing QC_PROGRAM_ELEC
  echo ${WARNING}
  echo 'Missing QC_PROGRAM_ELEC in '${DPPINP}
  echo ' ' 
  echo ${MESSAGE_ERROR_TERM}
  exit
else 
  export QC_PROGRAM_ELEC=$(grep QC_PROGRAM_ELEC ${DPPINP} | awk '{print $3}')
fi

# Searching QC_PROGRAM_VIB
grep --quiet 'QC_PROGRAM_VIB' ${DPPINP}
CHECK_QC_PROGRAM_VIB=$(echo $?)
if [ ${CHECK_QC_PROGRAM_VIB} = 1 ] ; then # Missing QC_PROGRAM_VIB
  echo ${WARNING}
  echo 'Missing QC_PROGRAM_VIB in '${DPPINP}
  echo ' '
  echo ${MESSAGE_ERROR_TERM}
  exit
else 
  export QC_PROGRAM_VIB=$(grep QC_PROGRAM_VIB ${DPPINP} | awk '{print $3}')
fi

echo "PROPERTY: ${PROPERTY}"
echo "QC_METHOD: ${QC_METHOD}"
echo "QC_PROGRAM_ELEC: ${QC_PROGRAM_ELEC}"
echo "QC_PROGRAM_VIB: ${QC_PROGRAM_VIB}"


case ${QC_PROGRAM_ELEC} in

  g16)
    sh "${DIR_d77_dpp}"/dpp_g16.sh
  ;;

  *)
    echo ' '
    echo "Invalid QC_PROGRAM_ELEC: ${QC_PROGRAM_ELEC}"
    echo ' '
    echo ${MESSAGE_ERROR_TERM}
    exit
  ;;
esac


echo
echo 'Removing temporary files'
"${DIR_d77_dpp}"/rm_temp_files.sh
echo 'Done'

echo  
echo '---------------------------------------' 
echo 'Data preprocessing process was finished' 
echo '---------------------------------------' 

exit
