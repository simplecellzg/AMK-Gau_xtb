#!/bin/bash
# default sbatch FT2
#SBATCH --output=llcalcs-%j.log
#SBATCH --time=04:00:00
# partition selection

#_remove_this_in_ft_SBATCH -p shared --qos=shared
#SBATCH -c 1 --mem-per-cpu=2048
#SBATCH -n 32

# SBATCH --partition=cola-corta,thinnodes
# SBATCH -c 1
# SBATCH -n 48

# first  arg is inputfile
# second arg is nbatches (200 is a good number)
# third  arg is niter  
#if no arguments are provided, then a gui pops up
#convergence is achieved when no new TSs are found in three consecutive iterations
convergence=0
exe="llcalcs_TS.sh"
cwd="$PWD"
iter=0
sharedir=${AMK}/share
source ${code_path}/utils.sh
# Printing the references of the method
print_ref
#
if [ $# -eq 0 ]; then
   FILE="$(zenity --file-selection --filename="$PWD/*.dat" --file-filter="*.dat" --title="Select the input file" 2> /dev/null)"
   inputfile="$(basename $FILE)"
   echo "Selected input file: $inputfile"
   answer="$(zenity --forms --title="llcalcs.sh GUI" --text="Add input data" \
      --add-entry="Number of tasks" \
      --add-entry="Number of iterations" \
      --add-entry="Max number of running tasks" 2>/dev/null | awk 'BEGIN{FS="|"};{print $1,$2,$3}' )"
   nbatch=$(echo "$answer" | awk '{print $1}')
   niter=$(echo "$answer" | awk '{print $2}')
   runningtasks=$(echo "$answer" | awk '{print $3}')
elif [ $# -eq 3 ]; then
   inputfile=$1
   nbatch=$2
   niter=$3
   if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
      runningtasks=$SLURM_NTASKS
   else
     echo "With three arguments it must be run under the SLURM batch system:"
     echo "sbatch $exe inputfile ntasks niter"
     exit 1
   fi
elif [ $# -eq 4 ]; then
   inputfile=$1
   nbatch=$2
   niter=$3
   runningtasks=$4
else
   echo You must provide zero or four arguments:
   echo "nohup $exe >llcalcs.log 2>&1 &"
   echo or
   echo "nohup $exe inputfile ntasks niter runningtasks >llcalcs.log 2>&1 &" 
   exit 1
fi
export runningtasks

###Are we in the right folder?
if [ ! -f $inputfile ];then
   echo "$inputfile is not in this folder"
   exit 1
fi
if [ -z $nbatch ] || [ -z $niter ]; then
   echo "Number of batches and/or number of iterations have not been set"
   exit 1
fi
#EMN. If nbatch=0 do not run dynamics.
if [ $nbatch -eq 0 ]; then niter=1 ; fi
#EMN
read_input
###
echo ""
echo "Number of iterations  = $niter"
echo "Tasks per iteration   = $nbatch"
echo ""
###checks and writing stuff
xyzfiles_check 
###
sampling_calcs
###
amkscript=0
print_method_screening
###
echo ""
echo "CALCULATIONS START HERE"
echo ""
iter=1
#set interactive mode to 0
inter=0
export inter
echo $$ > .script.pid
#
if [ $sampling -eq 3 ]; then echo " CK calculations " > ChemKnow.log ; fi
while [ $iter -le $niter ]; do
   export iter
   echo "======================="
   echo "      Iter: ${iter}/${niter}"
   echo "======================="
   echo "$iter/$niter" > iter.txt
   if [ $nbatch -gt 0 ]; then
      echo "   Running TS search"
      start=$(date +%s.%N)
      ${code_path}/amk_parallel_rt.sh $inputfile $nbatch >/dev/null
      end=$(date +%s.%N)    
      tt=$( echo "$end - $start" | bc -l | awk '{printf "%4.0f",$1}')
      echo "   time: $tt s"
   fi
##check that in ExtForce nts is not -1 (-->all minima explored). If that is the case-->stop
   nts=$(sqlite3 ${tsdirll}/track.db "select nts from track" | awk '{a=$1};END{print a+1-1}')
   if [ $nts -lt 0 ]; then
      echo "     All MIN used   "
      echo "    Stop iters here " 
      #rm -rf batch*
      break 
   fi
##check that tslistll file exists
   if [ ! -f $tslistll ]; then
      echo ""
      echo "        ERROR:      "
      echo "    tslist is empty "
      echo "   Check: batch dirs"
      if [ $sampling -eq 3 ]; then echo "   or ChemKnow.log" ; fi 
      echo ""
      exit
   fi
##Doing IRCs
   echo "   IRC screening calcs  "
   irc.sh screening > /dev/null
###check that the nts is greater than the previous value
   if [ $iter -eq 1 ]; then
     ntsp=0
   else
     ((itm1=iter-1))
     ntsp=$(sqlite3 ${tsdirll}/track.db "select nts from track where id=$itm1" | awk '{print $1+1-1}' )
   fi
   ntsa=$(sqlite3 ${tsdirll}/track.db "select nts from track where id=$iter" | awk '{print $1+1-1}' )
   if (( $(echo "$ntsa <= $ntsp" |bc -l) )); then
      echo "      No new TSs    "
      echo "    Skip this iter  "
      echo ""  
      #rm -rf batch*
      ((iter=iter+1))
      ((convergence=convergence+1))
      if [ $convergence -eq 3 ]; then
         echo "    Stop iters here " 
         echo "      Convergence   " 
         break
      fi
      continue 
   fi
   ((iter=iter+1))
done
echo ""
echo "END OF THE TS CALCULATIONS"
echo ""

