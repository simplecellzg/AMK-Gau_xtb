   combine_dir=combine_hl
   ts_numb=`ls ts*.inp | wc -1`
   ts_inp=`ls ts*.inp`
   ls ts*.inp >ts_inp_list
   sed -n '1,5p' ts_inp_list >selected_ts
    for i in `cat selected_ts`
    do
    orca_path=$(whereis orca | awk '{print $2}')
    orca_pltvib_path=$(whereis orca_pltvib | awk '{print $2}')
    named=`basename $i .inp`
    scratch_dir=${named}_scratch
    if [ ! -d ${scratch_dir} ];then
        mkdir ${scratch_dir}  
    else
        rm -rf ${scratch_dir}/*
    fi
    cp ${i} ${scratch_dir}/
    cp ${named}_ts_NEB-TS.xyz ${scratch_dir}/
     current_path=$PWD
     cd ${scratch_dir}
     #echo $PWD
     #echo $batch
     ${orca_path}  ${named}.inp >${named}.out 2>/dev/null
     ${orca_pltvib_path}  ${named}.hess 6 2>/dev/null
     rm -r *.tmp
     cd ../
     cp ${scratch_dir}/${named}.xyz ../${combine_dir}/
     cp ${scratch_dir}/${named}.inp ../${combine_dir}/
     cp ${scratch_dir}/${named}.out ../${combine_dir}/
     cp ${scratch_dir}/${named}.hess.v006.xyz ../${combine_dir}/
done

