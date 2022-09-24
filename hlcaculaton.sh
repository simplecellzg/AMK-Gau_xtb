#/bin/bash
#step 2: heighlevel
exe="hlcalcs.sh"
source ${code_path}/utils.sh
cp Default.Route /opt/g09/
export code_path=/root/docker_amk_g09/workdir/sub_system_hl
tsdirhl=${code_path}/combine_*/hl_g09
cd ${tsdirhl}

# hl caculation
cat /dev/null >hltslist
gjf_num=`ls -l ts*.gjf|wc -c`
runningtasks=1
for gjffile in `ls ts*.gjf`
do
name=`basename $gjffile .gjf`
echo "g09 <${gjffile} >${name}.log" >>hltslist
done

# inputname=`basename $j .xyz`
# cd ./${inputname}/
# ${code_path}/hlcalcs.sh ${inputname}.dat 1 >hlcalcs.log
# echo ${inputname}
# cd ../
# done
# export runningtasks

# system="$(basename $inputfile .dat)"
echo "   Optimizing TSs    "
cp hltslist runTS.dat
#TS.sh $inputfile > /dev/null
#nohup $parallel $COMMAND ::: $TASKS  >/dev/null 2>&1 & 
nohup parallel -j $runningtasks -a runTS.dat >/dev/null 2>&1 & 

#doparallel "${code_path}/runTS.sh {1}" "$(seq $gjf_num)"





