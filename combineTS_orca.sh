#/bin/bash
cwd=$PWD
export code_path=${cwd}
workdir='./combineTS_orca/'
rm -r ${workdir}
mkdir ${workdir}
mkdir ${workdir}/tsdirLL_assoc/
mkdir ${workdir}/dat_xyz/
mkdir ${workdir}/hl_orca/
#mkdir ${workdir}/tsdirLL_assoc/backup
cat /dev/null >${workdir}/tsdirLL_assoc/tslist 
num=0
#for i in `ls -d ${code_path}/combine/assoc*/`
for i in `ls -d ${cwd}/llcaculation/assoc*/`
do 
    rows=0
    basename=`basename $i`
    cp $i/${basename}.dat ${workdir}/dat_xyz/
    cp $i/${basename}.xyz ${workdir}/dat_xyz/
    cp $i/${basename}_ref.xyz ${workdir}/dat_xyz/
    cp $i/${basename}.dat ${workdir}/assoc.dat
    cp $i/${basename}.xyz ${workdir}/assoc.xyz
    cp $i/${basename}_ref.xyz ${workdir}/assoc_ref.xyz
    sed -i '2c molecule  assoc' ${workdir}/assoc.dat
    cp ${workdir}/assoc.dat ${workdir}/amk.dat
    cd $i/tsdirLL*/
    if [ -f tslist ]; then
        for j in `awk '{print $3}' tslist`
        do 
            num=$(echo "${num}+1"| bc )
            rows=$(echo "${rows}+1"| bc )
            line=`awk 'NR=="'${rows}'"{print}' tslist`
            nt=`awk 'NR=="'${rows}'"{print $2}' tslist`
            name=`awk 'NR=="'${rows}'"{print $3}' tslist`
            image=`awk 'NR=="'${rows}'"{print $4}' tslist`
            energy=`awk 'NR=="'${rows}'"{print $5}' tslist`
            w1=`awk 'NR=="'${rows}'"{print $6}' tslist`
            w2=`awk 'NR=="'${rows}'"{print $7}' tslist`
            w3=`awk 'NR=="'${rows}'"{print $8}' tslist`
            w4=`awk 'NR=="'${rows}'"{print $9}' tslist`
            traj=`awk 'NR=="'${rows}'"{print $10}' tslist`
            traj_num=`awk 'NR=="'${rows}'"{print $11}' tslist`
            path=`awk 'NR=="'${rows}'"{print $12}' tslist`
            path_num=`awk 'NR=="'${rows}'"{print $13}' tslist`

            #others=`awk -v OFS='    ' 'NR=="'${rows}'"{$1=$2=$3=""; print $0}' tslist`
            nt=${num}
            name_last=${name#*_}
            name_combine="ts${num}_${name_last}_${basename}"
            path_num_combine="llcaculation/${basename}/${path_num}"
            #echo ${name_combine}
            #echo ${name_last}
            ookk=$(grep -n "ORCA TERMINATED NORMALLY" ./${name}_ts.out |awk -v FS=":" '{print $1}')
            if [ ! -z ${ookk} ]; then
               printf "ts%5s%35s%12s%15s%10s%10s%10s%10s%10s%5s%10s%30s \n" $nt ${name_combine} ${image} ${energy} ${w1}\
                ${w2} ${w3} ${w4} ${traj} ${traj_num} ${path} ${path_num_combine}>>${cwd}/${workdir}/tsdirLL_assoc/tslist 
               #printf "ts%5s%18s%60s \n" $nt ${name_combine} "$others" >>${cwd}/${workdir}/tsdirLL_assoc/tslist 
               # cp ./${name}.out ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}.out
               # cp ./${name}.den ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}.den
               # cp ./${name}.molden ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}.molden
      
               cp ./${name}_ts.out ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_ts.out
               cp ./${name}_irc_trj.xyz ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_irc_trj.xyz
               cp ./${name}_irc.out ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_irc.out
               cp ./${name}_ts.xyz ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_ts.xyz
               cp ./${name}_ts_vib.xyz ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_ts_vib.xyz
              cp ./${name}_ts_NEB-TS.xyz ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_ts_NEB-TS.xyz
               cp ./${name}_ts.hess ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_ts.hess
             # cp ./${name}_ts.interp ${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_ts.interp
               freq_file="${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_ts.out"
   
               cp ./${name}_ts.xyz ${cwd}/${workdir}/hl_orca/${name_combine}_ts.xyz
             # cp ./${name}_ts_NEB-TS.xyz ${cwd}/${workdir}/hl_orca/${name_combine}_ts_NEB-TS.xyz
               charge=$(awk 'BEGIN{ch=0};{if($1=="charge") ch=$2};END{print ch}' ${cwd}/${workdir}/amk.dat)  
               mult=$(awk 'BEGIN{mult=1};{if($1=="mult") mult=$2};END{print mult}' ${cwd}/${workdir}/amk.dat)
               ts_xyz=`basename ${cwd}/${workdir}/hl_orca/${name_combine}_ts.xyz`
             # ts_xyz=`basename ${cwd}/${workdir}/hl_orca/${name_combine}_ts_NEB-TS.xyz`
               sed 's/carga/'$charge'/;s/mult/'$mult'/;s/llts_for_hlts.xyz/'$ts_xyz'/' ${code_path}/optTS_HL_orca >${cwd}/${workdir}/hl_orca/${name_combine}.inp
               #sed 's/example/'${name_combine}'/' ${cwd}/hl_input_opt_g09 >${cwd}/${workdir}/hl_g09/${name_combine}.gjf
               #awk '{print $1,$2,$3,$3}' geom_file >>${cwd}/${workdir}/hl_g09/${name_combine}.gjf
               #get_geom_g09.sh ${freq_file} |awk '{printf "%-10s %12f %12f %12f\n",$1,$2,$3,$4}' >>${cwd}/${workdir}/hl_g09/${name_combine}.gjf
               #get_geom_g09.sh ${freq_file} >>${cwd}/${workdir}/hl_g09/${name_combine}.gjf
               #echo ${freq_file}
               #echo -e "\n\n\n" >>${cwd}/${workdir}/hl_g09/${name_combine}.gjf
               #echo ${num}
            else
               echo "${cwd}/${workdir}/tsdirLL_assoc/${name_combine}_ts.out NOT TERMINATED NORMALLY"
            fi
        done
        echo ${basename}
        cd $cwd 
    else
        echo "${basename} don't have TS"
        cd $cwd
    fi    
done
cd ${cwd}/${workdir}
${cwd}/screening.sh  amk.dat
