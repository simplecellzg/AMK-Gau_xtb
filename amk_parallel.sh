#!/bin/bash
# default sbatch FT2
#SBATCH --output=amk_parallel-%j.log
#SBATCH --time=04:00:00
# partition selection

#_remove_this_in_ft_SBATCH -p shared --qos=shared
#SBATCH -c 1 --mem-per-cpu=2048
#SBATCH -n 8

# SBATCH --partition=cola-corta,thinnodes
# SBATCH -c 1
# SBATCH -n 24


#exe=$(basename $0)
# under batchs systems the scripts are copied to a generic script (in slurm slurm_script)
#export code_path="/root/docker_amk_g09/workdir_amk_orca/code"
source ${code_path}/utils.sh

#current working dir
cwd=$PWD
sharedir=${AMK}/share
exe="amk_parallel.sh"
# Printing the references of the method
print_ref

#Checking the input of this script and also the presence of some files
re='^[0-9]+$'
if [ $# -eq 2 ]; then
   if [[ ! $2 =~ $re ]]; then usage "Second argument must be a number";fi
   noj1=$(( $(sort -nr <(find . -maxdepth 1 -type d -print | grep 'batch' | sed 's@batch@@;s@./@@') | head -n1) +1 ))
   nojf=$(( $noj1 + $2 -1))
elif [ $# -eq 3 ]; then
   if [[ ! $3 =~ $re ]]; then usage "Third argument must be a number";fi
   noj1=$2
   nojf=$3
else
   usage "At least 2 arguments required"
fi
inputfile=$1
###Define inputfile and create symbolic link
define_inputfile
###reading input file
read_input
###checks. amk_parallel only with md>0
amk_parallel_check
###keywords check: molecule and rate
keywords_check
### ${molecule}.xyz and/or ${frA}.xyz ${frB}.xyz file check
xyzfiles_check
##
sampling_calcs
##
amkscript=0
print_method_screening
##
echo ""
echo "CALCULATIONS START HERE"
#####for sampling 31 get the association complexes
if [ $sampling -eq 31 ]; then exec_assoc ; fi
###amk_parallel specific setups
amk_parallel_setup
###ExtForce sampling
if [ $sampling -eq 3 ]; then
##select the initial minimum
   minsf=${tsdirll}/MINs/SORTED/MINlist_sorted
   minrf=${tsdirll}/mins.inp
   mindb=${tsdirll}/MINs/SORTED/mins.db
   rm -rf ${tsdirll}/ts_bonds.inp
   if [ -f $minsf ];then
      min=$(awk 'NR==FNR{e[NR]=$1;++n};NR>FNR{ok=1;for(i=1;i<=n;i++) {d=$NF-e[i];if($3 ~/min0/ || d*d < .001) ok=0}; if(ok==1) {print $2;exit} } ' $minrf $minsf)
      emin=$(awk 'NR==FNR{e[NR]=$1;++n};NR>FNR{ok=1;for(i=1;i<=n;i++) {d=$NF-e[i];if($3 ~/min0/ || d*d < .001) ok=0}; if(ok==1) {print $NF;exit} } ' $minrf $minsf)
      if [ -z $min ]; then
         echo All minima have been explored
	 sqlite3 ${tsdirll}/track.db "insert into track (nts) values (-1);"
         exit 1 
      fi
      echo $emin >> $minrf
      echo $natom > ${molecule}.xyz
      echo ""    >> ${molecule}.xyz  
      sqlite3 $mindb "select geom from mins where name='MIN$min'" >> ${molecule}.xyz
      minp=$(echo $emin | awk '{printf " E_MIN = %8.2f\n",$1}')
   else
      minp=$(echo "0.00" | awk '{printf " E_MIN = %8.2f\n",$1}')
      echo "0.000000" > $minrf 
   fi
##
   ofn=extforce.log
   echo " ExtForce:$minp  " 
   echo " ExtForce:$minp  " >> $ofn 

   noj1=1
   re='^[0-9]+$'
   nchan=$(ExtForce.py) 
   if ! [[ $nchan =~ $re ]] ; then
      echo " Error: No neighbors for " $nchan
      echo " Error: No neighbors for " $nchan >> $ofn 
      exit 1 
   else 
      noj2=$nchan
   fi
   nojf=$((nforces*noj2))
   rm -rf bfgs.log none.*
###
   ncf=$(awk 'BEGIN{ncf=0};{nci=0;for(i=1;i<=NF;i++) if($i=="f") nci=1;if(nci==1)++ncf};END{print ncf}' ${tsdirll}/ts_bonds.inp )
   nco=$((noj2-ncf))
   ntri=$((nco+nforces*ncf))
###
   if [ $nojf -eq 0 ]; then
      echo "  No channels    "
      echo "  No channels    " >> $ofn 
      echo "" 
      echo ""  >> $ofn
      exit 1
   else
      echo "   $noj2 channels"
      echo "    $ntri trajs  "  
      echo "   $noj2 channels"  >> $ofn 
      echo "    $ntri trajs  "  >> $ofn 
   fi 
   echo "" 
   echo ""  >> $ofn
fi
###
echo ""
echo "A progress bar will pop up very shortly"
#The last job is tors.sh
nojf=$((nojf + 1))
#Then we submit the nojf-noj1+1 trajectory jobs
#ft2 slurm
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
  if (( $nojf-$noj1+1 < $SLURM_NTASKS )); then 
    echo "WARNING: Number of trajectory jobs ($nojf-$noj1+1) lower than allocated tasks ($SLURM_NTASKS)."
  fi
fi
doparallel "${code_path}/runGP.sh {1} $molecule $cwd $nojf" "$(seq $noj1 $nojf)" 
#doparallel "${code_path}/runGP.sh {1} $molecule $cwd $nojf $code_path" "$(seq $noj1 $nojf)" 

