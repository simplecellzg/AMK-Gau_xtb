#!/bin/bash
if [ "$3" = "mopac" ]; then
   mopac ${2}/assoc${1}.mop 2>/dev/null
elif [ "$3" = "qcore" ]; then
   entos.py ${2}/assoc${1}.qcore > ${2}/assoc${1}.out 2>&1 
   if [ -f ${2}/assoc${1}_opt.xyz ]; then
      cat ${2}/assoc${1}_opt.xyz >> ${2}/assoc${1}.out 
   else
      echo Error >> ${2}/assoc${1}.out 
   fi
elif [ "$3" = "orca" ]; then
   orca_path=$(whereis orca | awk '{print $2}')
   mkdir ${2}/assoc${1}_orca
   # copy orca input file in scratch dir
   cp ${2}/assoc${1}.inp ${2}/assoc${1}_orca/
   cd ${2}/assoc${1}_orca
  # ${orca_path}  ${2}/assoc${1}_orca/assoc${1}.inp >${2}/assoc${1}_orca/assoc${1}.out 2>/dev/null
   cd ${2}
   sed '2d' ${2}/assoc${1}_orca/assoc${1}.xyz |sed '1a \ ' >${2}/assoc${1}_opt.xyz
   #mv ${2}/assoc${1}_orca/assoc${1}.xyz ${2}/assoc${1}_opt.xyz
   mv ${2}/assoc${1}_orca/assoc${1}.out ${2}/assoc${1}.out
   rm -rf ${2}/assoc${1}_orca
fi
