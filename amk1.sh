#!/bin/bash
source ${code_path}/utils.sh
#On exit remove tmp files
tmp_files=(ConnMat tmp* ScalMat *.arc *.mop fort.* partial_opt ts_opt *_dyn* *_backup rotate.dat minn black_list* bfgs.log none.out forces.xyz velocities.xyz restraints.xyz energies.txt freq.molden min.xyz ts_opt.xyz ts.xyz min_opt.xyz v0 grad.*)
#trap 'err_report2 $LINENO $gauss_line' ERR
#trap cleanup EXIT INT

##Defining paths and names
cwd=$PWD
sharedir=${AMK}/share
exe=$(basename $0)
if [ $# -eq 0 ]; then usages "One argument is required" ; fi
inputfile=$1

# Printing the references of the method
print_ref
##Define input file and create symbolic link-->amk.dat
define_inputfile
###Reading stuff from inputfile
read_input
###keywords check
keywords_check
### check files and write some stuff 
xyzfiles_check
###Peforming some calcs for the various samplings
sampling_calcs
### print method and screening sections
amkscript=1
print_method_screening
#################################
##  Starting the calculations
#################################
echo ""
echo "CALCULATIONS START HERE"
#####for association and vdw get the association complexes
if [ $sampling -ge 30 ]; then
   exec_assoc
###for association stop here
   if [ $sampling -eq 30 ]; then
      echo ""
      echo "END OF THE CALCULATIONS"
      echo "Check your ${frA}-${frB} structure in file ${molecule}.xyz" 
      echo ""
      exit 
   fi
fi
###ExtForce--> genenerate forces
if [ $sampling -eq 3 ]; then
   ba=$(echo "$nb" | awk 'BEGIN{ba=0};/batch/{ba=1};END{print ba}')
   if [ $ba -eq 0 ]; then
      echo ExtForce sampling must be submitted with llcalcs.sh or amk_parallel.sh
      exit 1    
   else 
      mdc=2
      nba0=$(echo "$nb" | sed 's/batch//')
      nba=$(echo "$nba0" | awk '{a=$1/'$nforces' ;if(a==int(a)) print a; else print int(a+1)}' )
      fstp=$((nba0-nforces*nba))
      fstp0=$((1-nforces))
      force=$((30*nforces+fstp*30))
      nbondsform=$(awk  'NR=='$nba'{n=0;for(i=1;i<=NF;i++) if($i=="f") ++n };END{print n}' ${tsdirll}/ts_bonds.inp )
      if [ $nbondsform -gt 0 ]; then
         awk 'NR=='$nba'{n=0;for(i=1;i<=NF;i++) if($i=="f") print $(i+1),$(i+2),'$force' }' ${tsdirll}/ts_bonds.inp > fort.69
      else
         if (( $(echo "$fstp > $fstp0" |bc -l) )); then
            echo "No bonds formed-->Only the first batch is employed"
	    exit 1 
	 fi
      fi
      nbondsbreak=$(awk 'NR=='$nba'{n=0;for(i=1;i<=NF;i++) if($i=="b") ++n };END{print n}' ${tsdirll}/ts_bonds.inp )
      if [ $nbondsbreak -gt 0 ]; then
         awk 'NR=='$nba'{n=0;for(i=1;i<=NF;i++) if($i=="b") print $(i+1),$(i+2),"110" }' ${tsdirll}/ts_bonds.inp > fort.68
      fi
   fi
   createMat.py ${molecule}.xyz 1 
else
##select starting structure
   ${code_path}/sel_mol.sh $inputfile $multiple_minima
   frag_check
##lift MD-constraint in subsequent iterations (diels_bias for instance)
   if [ -f $kmcfilell ] && [ -f $minfilell ] && [ $mdc -ge 1 ] && [ $ndis -eq 1 ]; then mdc=0 ; fi
fi
##template for the dynamics
generate_dynamics_template
echo "${dytem1}"
#echo "${dytem1}" >${cwd}/${llcaculation}/dynamics_template_aaaaa
##make temporary folders
make_temp_folders
###Opt the starting structure and get e0 and emaxts
# before opt change program_opt =qcore
#program_opt=qcore
#opt_start
cp ${molecule}_ref.xyz opt_start.xyz
####
if [ $mdc -ge 1 ]; then itrajn=5 ; fi
###Loop over the trajectories

for i in $(seq 1 $itrajn) 
do 
  named=${molecule}_dyn${i}
  echo ""
##Empty temporary folders
  rm -rf partial_opt/* ts_opt/* 
####
#This is only for internal dynamics (MOPAC)
####
  if [ $mdc -ge 1 ] && [ $i -gt 1 ]; then
     echo "Searching for reactive events with:"
  else
     echo ""
     echo "+-+-+-+-+-+-+-+-+-          Trajectory ${i}          +-+-+-+-+-+-+-+-+-"
     if [ $md -eq 0 ]; then
       echo "Performing BXDE MD"
       if [ "$program_md" = "qcore" ]; then sed 's/carga/'$charge'/' ${sharedir}/grad > grad.dat ; fi
       bxde.py $inputfile &>  ${named}.log
       if [ ! -f traj.xyz ]; then
          echo "traj.xyz does not exist"
          continue 
       else
          mv traj.xyz coordir/${named}.xyz
          if [ $postp_alg -eq 2 ]; then mv bond_order.txt coordir/${named}.bo ; fi
       fi
     elif [ $md -eq 1 ]; then 
        echo "Performing standard MD"
        echo "$dytem1"     > ${named}.mop
        initialqv_mopac_samp1.sh ${molecule}_freq.out $seed $excite $nlms $lstnm | nm.exe >> ${named}.mop
        mopac ${named}.mop &> ${named}.log 
        if [ ! -f ${named}.xyz ]; then
           echo "${named}.xyz does not exist"
           continue 
        else
           mv ${named}.xyz coordir
        fi
     elif [ $md -eq 2 ]; then
        echo "Performing standard MD"
        if [ "$program_md" = "mopac" ]; then
           echo "$dytem1"     > ${named}.mop
           if [ $sampling -eq 3 ]; then thmass=-1 ; fi
           initialqv_mopac_samp2.sh ${molecule}_freq.out $excite $nlms $lstnm $thmass | termo.exe | sed 's/ 1.d/1.d/g' >> ${named}.mop
           mopac ${named}.mop &> ${named}.log
           if [ ! -f ${named}.xyz ]; then
              echo "${named}.xyz does not exist"
              continue 
           else
              mv ${named}.xyz coordir
           fi
        elif [ "$program_md" = "qcore" ]; then
           rm -rf traj.xyz
           #sed 's/carga/'$charge'/' ${code_path}/MD_qcore_dft > ${named}.qcore
           sed 's/carga/'$charge'/' ${code_path}/MD_qcore > ${named}.qcore
           #sed -i 's/multp/'$mult'/' ${named}.qcore
           sed -i 's/500/'$nfs'/' ${named}.qcore
           awk 'NR!=2{print $0};NR==2{print ""};END{print '"$excite"'"\n0\n0"}' opt_start.xyz | termo.exe | awk 'BEGIN{c=4.5710047e-9;print '$natom'"\n"};NF==3{print $1*c,$2*c,$3*c}' >v0
           ${code_path}/entos.py ${named}.qcore > ${named}.out 2>&1
           if [ ! -f traj.xyz ]; then
              echo "traj.xyz does not exist"
              continue 
           else
              mv traj.xyz coordir/${named}.xyz
           fi
         elif [ "$program_md" = "orca" ]; then
               echo "$dytem1"     > ${named}.inp
               orca_path=$(whereis orca | awk '{print $2}')
               mkdir ${named}_orca
               # copy orca input file in scratch dir
               cp ${named}.inp ${named}_orca/
               cd ${named}_orca
               #echo $PWD
               #echo $batch
               ${orca_path}  ${named}.inp >${named}.out 2>/dev/null
               cd ../
               mv ${named}_orca/pos.xyz coordir/${named}.xyz
               mv ${named}_orca/${named}.out ${named}.out
               #rm -rf ${named}_orca

        fi
     else
       echo "Reading external dynamics results from coordir"
     fi
  fi
  echo "Standard MD end"
  echo "Performing standard NEB"

# 2021.08.19 After MD calculation, run orca neb_ci with first and last frame
# to have minimum energy path which have better ts.
   # get reactant and product frame
    snapshots_mopac.sh coordir/${named}.xyz $irange >coordir/${named}_mat 
    neb_dir=snapshot_neb_orca_${named}
    prd_opt_dir=product_opt_${named}
    if [ ! -d ${neb_dir} ];then
        mkdir ${neb_dir}  
    else
        rm -rf ${neb_dir}/*
    fi    
    natoms_trj=$(awk 'NR==2{print $1}' coordir/${named}_mat)
    steps_num=$(awk 'NR==3{print $1}' coordir/${named}_mat)
    awk 'NR==2{print $1}' coordir/${named}_mat >${neb_dir}/reactant.xyz
    awk 'NR==2{print $1}' coordir/${named}_mat >${neb_dir}/product.xyz
    
    #reactant: first
    echo "#The num 1 step" >>${neb_dir}/reactant.xyz
    start_num=3
    end_num=$(echo "${natoms_trj}+3" | bc )
    awk 'NR>'${start_num}' && NR<='${end_num}'' coordir/${named}_mat >>${neb_dir}/reactant.xyz    

    #product: last
    echo "#The num ${steps_num} step" >>${neb_dir}/product.xyz 
    end_num=$(wc -l coordir/${named}_mat | awk '{print $1}' )
    start_num=$(echo "${end_num}-${natoms_trj}" | bc )
    awk 'NR>'${start_num}' && NR<='${end_num}'' coordir/${named}_mat >>${neb_dir}/product.xyz    
    #
    #optimize product 
    orca_path=$(whereis orca | awk '{print $2}')
    orca_pltvib_path=$(whereis orca_pltvib | awk '{print $2}')
    if [ ! -d ${prd_opt_dir} ];then
        mkdir ${prd_opt_dir}  
    else
        rm -rf ${ prd_opt_dir}/*
    fi
    cp ${neb_dir}/product.xyz ${prd_opt_dir}/
     current_path=$PWD
     # copy orca input file in scratch dir
     sed 's/carga/'$charge'/;s/mult/'$mult'/' ${code_path}/optTS_opt_orca >${prd_opt_dir}/${named}_prd_opt.inp
     #awk '{print $0}' ${code_path}/NEB_orca >ts_opt/${named}_neb.inp
     #cp ${neb_dir}/reactant.xyz ${named}_NEB_orca/ 
     cp ${neb_dir}/product.xyz ${prd_opt_dir}/ 
     #cp ts_opt/${named}_neb.inp ${named}_NEB_orca/   

     cd ${prd_opt_dir}
     #echo $PWD
     #echo $batch
     ${orca_path}  ${named}_prd_opt.inp >${named}_prd_opt.out 2>/dev/null
     cd ../
     cp ${prd_opt_dir}/${named}_prd_opt.xyz ${neb_dir}/product_opt.xyz

    #NEB_clc
    orca_path=$(whereis orca | awk '{print $2}')
    orca_pltvib_path=$(whereis orca_pltvib | awk '{print $2}')
    if [ ! -d ${named}_NEB_orca ];then
       mkdir ${named}_NEB_orca 
    else
       rm -rf ${named}_NEB_orca/*
     fi
     #mkdir ${named}_NEB_orca
     current_path=$PWD
     # copy orca input file in scratch dir
     sed 's/carga/'$charge'/;s/mult/'$mult'/' ${code_path}/optTS_NEB_orca >${named}_NEB_orca/${named}_neb.inp
     #awk '{print $0}' ${code_path}/NEB_orca >ts_opt/${named}_neb.inp
     cp ${neb_dir}/reactant.xyz ${named}_NEB_orca/ 
     cp ${neb_dir}/product_opt.xyz ${named}_NEB_orca/ 
     #cp ts_opt/${named}_neb.inp ${named}_NEB_orca/   

     cd ${named}_NEB_orca
     #echo $PWD
     #echo $batch
     ${orca_path}  ${named}_neb.inp >${named}_neb.out 2>/dev/null
     python3 ${code_path}/neb_snapshots.py ${named}_neb.interp 2>/dev/nul 
     cd ../
     neb_coordir=neb_coordir
     if [ ! -d ${neb_coordir} ];then
        mkdir ${neb_coordir}  
     else
         rm -rf ${neb_coordir}/*
     fi
     mkdir ${neb_coordir}/${named}_ts_neb_frames
     cp ${named}_NEB_orca/${named}_neb.inp ${neb_coordir}/
     cp ${named}_NEB_orca/${named}_neb.out ${neb_coordir}/
     cp ${named}_NEB_orca/${named}_neb.interp ${neb_coordir}/
     cp ${named}_NEB_orca/neb_frames/* ${neb_coordir}/${named}_ts_neb_frames/
    # cp ${named}_NEB_orca/${named}_neb_initial_path_trj.xyz ${neb_coordir}/${named}_neb_MEP_trj.xyz
     cp ${named}_NEB_orca/${named}_neb_MEP_trj.xyz ${neb_coordir}/
     #cp ${named}_NEB_orca/${named}_neb_NEB-CI_converged.xyz ${neb_coordir}/

     # find local max point
     last_itera_line=$(grep -n "Iteration" ${neb_coordir}/${named}_neb.interp |tail -1 |awk -v FS=":" '{print $1}')
     #last_itera_num=$(grep -n "Iteration" NEB_orca.interp |tail -1 |awk -v FS=":" '{print $3}' |sed 's/^[ \t]*//g')
     last_image_line=$(grep -n "Images" ${neb_coordir}/${named}_neb.interp |tail -1 |awk -v FS=":" '{print $1}')
     last_interp_line=$(grep -n "Interp." ${neb_coordir}/${named}_neb.interp |tail -1 |awk -v FS=":" '{print $1}')
     last_interp_line2=$(echo "${last_interp_line}-3" | bc )
     awk 'NR>'${last_image_line}' && NR<='${last_interp_line2}'' ${neb_coordir}/${named}_neb.interp >${neb_coordir}/last_images
     awk 'NR>'${last_interp_line}'' ${neb_coordir}/${named}_neb.interp >${neb_coordir}/last_interp
     snapshots_mopac.sh  ${neb_coordir}/${named}_neb_MEP_trj.xyz $irange > ${neb_coordir}/${named}_neb_MEP_trj_mat 
     awk '{print $3}' ${neb_coordir}/last_images |xargs >${neb_coordir}/energy_row
     energy_row=$(sed 's/ /,/g' ${neb_coordir}/energy_row)
     python3 ${code_path}/local_max.py ${energy_row} >${neb_coordir}/local_max_num
     rm ${neb_coordir}/energy_row
    # rm -rf ${named}_NEB_orca

## change program_opt=gaussian because gaussian should connect with xtb
 # program_opt=gaussian
 # program_md=gaussian
 # echo "postp_alg is ${postp_alg}"  
  if [ $postp_alg -eq 0 ]; then
     echo "End of traj "$i
     echo "Only trajs. No post-processing algorithm applied to search for TSs"
     break
  elif [ $postp_alg -eq 1 ]; then
     postp_file=bbfs     
  elif [ $postp_alg -eq 2 ]; then
     postp_file=bots     
   elif [ $postp_alg -eq 3 ]; then
     postp_file=neb
  fi
###########
#From here everything is common for internal and external dynamics
###########
  if [ $i -eq 1 ]; then
     echo "  *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*"  > ${postp_file}.out
     echo "                $postp_file algorithm results          " >> ${postp_file}.out
     echo "  *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*" >> ${postp_file}.out
     echo ""                                                        >> ${postp_file}.out
  fi
  echo "  Trajectory $i" >> ${postp_file}.out

  if [ $postp_alg -eq 1 ]; then
     if [ $mdc -eq 0 ]; then
        snapshots_mopac.sh coordir/${named}.xyz  $irange | bbfs.exe >> ${postp_file}.out
     elif [ $mdc -ge 1 ]; then
        namedc=${molecule}_dyn1
        irange=$((20-4*(i - 1) ))
        irangeo2=$(echo "scale=0; $irange/2" | bc )
        echo "Time window (fs) = $irange "
        snapshots_mopac.sh coordir/${namedc}.xyz  $irange | bbfs.exe >> ${postp_file}.out
     fi
  elif [ $postp_alg -eq 2 ]; then
     bots.py $natom $cutoff $stdf ${named} >> ${postp_file}.out
  # make neb.out to record path information
  elif [ $postp_alg -eq 3 ]; then
       path_num=`awk '{print NF}' ${neb_coordir}/local_max_num`
       echo "Number of paths      "${path_num} >>${postp_file}.out
       for ii in $(seq ${path_num})
       do
       path_step=`awk -v a="$ii" '{print $a }' ${neb_coordir}/local_max_num`
       echo "Path=    ${ii} Step= ${path_step}" >>${postp_file}.out
      done
   fi
  path_num=$(awk '/Number of paths/{np=$4};END{print np}' ${postp_file}.out )
  if [ ${path_num} -eq 0 ]; then
     echo "This traj has no paths "
     continue
  fi
  
  echo "Npaths=" ${path_num}
  chapath[0]=0
  for ip in $(seq ${path_num})
  do
    if [ $postp_alg -eq 1 ]; then
      coordir=coordir
      gap_num_rt=10
      gap_num_pd=10
  
      #snapshots_mopac.sh ${neb_coordir}/${named}.xyz $irange >${neb_coordir}/${named}_mat 
      ts_num=$(awk '/Path=    '${ip}' Step=/{ttn=$4};END{print ttn}' ${postp_file}.out)
      reactant_num=$(echo "${ts_num}-${gap_num_rt}" | bc )
      product_num=$(echo "${ts_num}+${gap_num_pd}" | bc )
      mtd_dir=snapshot_path_${ip}_${reactant_num}_${product_num}
      if [ ! -d ${mtd_dir} ];then
          mkdir ${mtd_dir}  
      else
          rm -rf ${mtd_dir}/*
      fi
      natoms_trj=$(awk 'NR==2{print $1}' ${coordir}/${named}_mat)
      steps_num=$(awk 'NR==3{print $1}' ${coordir}/${named}_mat)
      awk 'NR==2{print $1}' ${coordir}/${named}_mat >${mtd_dir}/reactant.xyz
      awk 'NR==2{print $1}' ${coordir}/${named}_mat >${mtd_dir}/product.xyz
      awk 'NR==2{print $1}' ${coordir}/${named}_mat >${mtd_dir}//ts_qst3.xyz
      
      #reactant:if reactant number less than 0 ,use first step snapshot
      if [ ${reactant_num} -ge 0 ];then
        echo "#The num ${reactant_num} step" >>${mtd_dir}/reactant.xyz
        start_num=$(echo "${reactant_num}*${natoms_trj}+3" | bc )
        end_num=$(echo "${reactant_num}*${natoms_trj}+${natoms_trj}+3" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${coordir}/${named}_mat >>${mtd_dir}/reactant.xyz
      else
        echo "#The num 1 step" >>${mtd_dir}/reactant.xyz
        start_num=3
        end_num=$(echo "${natoms_trj}+3" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${coordir}/${named}_mat >>${mtd_dir}/reactant.xyz
      fi
  
      #product: if product number beyond max step ,use last step snapshot
      if [ ${product_num} -lt ${steps_num} ];then
        echo "#The num ${product_num} step" >>${mtd_dir}/product.xyz
        start_num=$(echo "${product_num}*${natoms_trj}+3" | bc )
        end_num=$(echo "${product_num}*${natoms_trj}+${natoms_trj}+3" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${coordir}/${named}_mat >>${mtd_dir}/product.xyz
      else
        echo "#The num ${steps_num} step" >>${mtd_dir}/product.xyz 
        end_num=$(wc -l ${neb_coordir}/${named}_neb_MEP_trj_mat | awk '{print $1}' )
        start_num=$(echo "${end_num}-${natoms_trj}" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${coordir}/${named}_mat >>${mtd_dir}/product.xyz
      fi
  
      #ts
        echo "#The num ${ts_num} step" >>${mtd_dir}/ts_qst3.xyz
        start_num=$(echo "${ts_num}*${natoms_trj}+3" | bc )
        end_num=$(echo "${ts_num}*${natoms_trj}+${natoms_trj}+3" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${coordir}/${named}_mat >>${mtd_dir}/ts_qst3.xyz


  elif [ $postp_alg -eq 3 ]; then

      # 2021.0803 use orca NEB-TS method with xtb2, need to make reactant and product xyz file withtout ts xyz file
      # 2021.0819 after orca NEB-CI method with xtb2, need to make reactant and product xyz file withtout ts xyz file 
      # for another NEB-TS calculation
      gap_num_rt=2
      gap_num_pd=2

      #snapshots_mopac.sh ${neb_coordir}/${named}.xyz $irange >${neb_coordir}/${named}_mat 
      ts_num=$(awk '/Path=    '${ip}' Step=/{ttn=$4};END{print ttn}' ${postp_file}.out)
      reactant_num=$(echo "${ts_num}-${gap_num_rt}" | bc )
      product_num=$(echo "${ts_num}+${gap_num_pd}" | bc )
      mtd_dir=snapshot_path_${ip}_${reactant_num}_${product_num}
      if [ ! -d ${mtd_dir} ];then
          mkdir ${mtd_dir}  
      else
          rm -rf ${mtd_dir}/*
      fi
      natoms_trj=$(awk 'NR==2{print $1}' ${neb_coordir}/${named}_neb_MEP_trj_mat)
      steps_num=$(awk 'NR==3{print $1}' ${neb_coordir}/${named}_neb_MEP_trj_mat)
      awk 'NR==2{print $1}' ${neb_coordir}/${named}_neb_MEP_trj_mat >${mtd_dir}/reactant.xyz
      awk 'NR==2{print $1}' ${neb_coordir}/${named}_neb_MEP_trj_mat >${mtd_dir}/product.xyz
      awk 'NR==2{print $1}' ${neb_coordir}/${named}_neb_MEP_trj_mat >${mtd_dir}//ts_qst3.xyz
      
      #reactant:if reactant number less than 0 ,use first step snapshot
      if [ ${reactant_num} -ge 0 ];then
        echo "#The num ${reactant_num} step" >>${mtd_dir}/reactant.xyz
        start_num=$(echo "${reactant_num}*${natoms_trj}+3" | bc )
        end_num=$(echo "${reactant_num}*${natoms_trj}+${natoms_trj}+3" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${neb_coordir}/${named}_neb_MEP_trj_mat >>${mtd_dir}/reactant.xyz
      else
        echo "#The num 1 step" >>${mtd_dir}/reactant.xyz
        start_num=3
        end_num=$(echo "${natoms_trj}+3" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${neb_coordir}/${named}_neb_MEP_trj_mat >>${mtd_dir}/reactant.xyz
      fi

      #product: if product number beyond max step ,use last step snapshot
      if [ ${product_num} -lt ${steps_num} ];then
        echo "#The num ${product_num} step" >>${mtd_dir}/product.xyz
        start_num=$(echo "${product_num}*${natoms_trj}+3" | bc )
        end_num=$(echo "${product_num}*${natoms_trj}+${natoms_trj}+3" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${neb_coordir}/${named}_neb_MEP_trj_mat >>${mtd_dir}/product.xyz
      else
        echo "#The num ${steps_num} step" >>${mtd_dir}/product.xyz 
        end_num=$(wc -l ${neb_coordir}/${named}_neb_MEP_trj_mat | awk '{print $1}' )
        start_num=$(echo "${end_num}-${natoms_trj}" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${neb_coordir}/${named}_neb_MEP_trj_mat >>${mtd_dir}/product.xyz
      fi

      #ts
        echo "#The num ${ts_num} step" >>${mtd_dir}/ts_qst3.xyz
        start_num=$(echo "${ts_num}*${natoms_trj}+3" | bc )
        end_num=$(echo "${ts_num}*${natoms_trj}+${natoms_trj}+3" | bc )
        awk 'NR>'${start_num}' && NR<='${end_num}'' ${neb_coordir}/${named}_neb_MEP_trj_mat >>${mtd_dir}/ts_qst3.xyz
     fi
# Find the highest energy point
    if [ $postp_alg -eq 1 ]; then
       ijc=$(awk '/Joint path=/{if($2=='$ip') ijc=$5};END{print ijc}' ${postp_file}.out)
    else
       ijc=0
    fi
    jp=$((ip - 1))
##If previous path was multiple, continue 
    chapath[$ip]=$ijc
    if [ ${chapath[$jp]} -eq 1 ]; then continue ; fi
##
    if [ $ijc -eq 0 ]; then
       echo "Path" $ip" (Single): $nppp attempt(s)  to locate the ts" 
       ll=$((wrkmode-1))
       dlt=$((wrkmode+1))
       ul=1
    elif [ $ijc -eq 1 ]; then 
       echo "Path" $ip" (Multiple): several attempts to locate the ts" 
       ll=$((1 - irangeo2))
       dlt=$(echo "scale=2; $irange/6" | bc | awk '{print int($1+0.5)}')
       ul=$irangeo2  
    fi
    npo=0
    for itspt in $(seq $ll $dlt $ul)
    do 
       npo=$((npo + 1))
       ctspt=$((100*ip + irangeo2 + itspt))
       echo "$min_template"         > partial_opt/pes$ctspt
       if [ $postp_alg -eq 1 ]; then
          cat partial_opt/fort.$ctspt >> partial_opt/pes$ctspt
       else
          cat partial_opt/fort.$ip >> partial_opt/pes$ctspt
       fi
       if [ "$program_md" = "mopac" ]; then
          mopac partial_opt/pes$ctspt  2> /dev/null
          geo_pes=$(get_geom_mopac.sh partial_opt/pes${ctspt}.out)
          if [ "$geo_pes" = "Error" ]; then continue ; fi
       elif [ "$program_md" = "qcore" ]; then
          echo $natom > min.xyz
          echo "" >> min.xyz
          if [ $postp_alg -eq 1 ]; then
             awk '{print $1,$2,$4,$6}' partial_opt/fort.$ctspt >> min.xyz
             labels=$(awk '{if($3=="0") {printf "%s%s",sep,NR; sep=","}};END{print ""}' partial_opt/fort.$ctspt)
          else
             awk '{print $1,$2,$4,$6}' partial_opt/fort.$ip >> min.xyz
             labels=$(awk '{if($3=="0") {printf "%s%s",sep,NR; sep=","}};END{print ""}' partial_opt/fort.$ip)
          fi 
          sed 's/labels/'"$labels"'/g;s/carga/'$charge'/' ${code_path}/opt_frozen_qcore > partial_opt/pes_qcore
          ${code_path}/entos.py partial_opt/pes_qcore > partial_opt/pes_qcore.out
          if [ ! -f min_opt.xyz ]; then 
             printf "     Pt%2s: failed-->Partial Opt failed\n" $npo
             continue
          fi
        # gaussian don't use partial opt   
       elif [ "$program_md" = "gaussian" ]; then
         echo "gaussian don't use partial opt!"
       elif [ "$program_md" = "orca" ]; then
         echo "orca don't use partial opt!"
       fi
       name=ts${i}_${ip}_${reactant_num}_${product_num}
       fileden=ts_opt/${name}.den
       if [ "$program_opt" = "mopac" ]; then
          echo "$ts_template"                      > ts_opt/$name
          echo "$geo_pes" | awk 'NF==4{print $0}' >> ts_opt/$name
          echo "$freq_template" | sed 's/oldgeo/oldgeo oldens/' >> ts_opt/$name
          mopac ts_opt/$name 2> /dev/null
          file=ts_opt/${name}.out
#If too many variables, run ts int
          if [ $(awk 'BEGIN{f=0};/Too many variables/{f=1};END{print f}' $file) -eq 1 ]; then
             sed -i 's/ts /ts int /g' ts_opt/$name
             mopac ts_opt/$name 2> /dev/null
          fi
###
       elif [ "$program_opt" = "qcore" ]; then
          mv min_opt.xyz ts.xyz
          sed 's/carga/'$charge'/' ${sharedir}/optTS > ts_opt/ts.dat
          ${code_path}/entos.py ts_opt/ts.dat > ts_opt/${name}.out
          if [ ! -f ts_opt.xyz ]; then 
             printf "     Pt%2s: failed-->No XYZ file found for the TS\n" $npo
             echo Error >> ts_opt/${name}.out
             continue 
          else
             cat ts_opt.xyz >> ts_opt/${name}.out
          fi
          file=ts_opt/${name}.out
       elif [ "$program_opt" = "gaussian" ]; then
      #construct g09 input file
          awk '{print $0}' ${code_path}/hl_input_qst3_xtb >ts_opt/${name}_ts.gjf
          #reactant
          echo "" >>ts_opt/${name}_ts.gjf
          echo "reactant" >>ts_opt/${name}_ts.gjf
          echo "" >>ts_opt/${name}_ts.gjf
          echo "${charge} ${mult}" >>ts_opt/${name}_ts.gjf
          awk 'NR>2 {print $0}' ${mtd_dir}/reactant.xyz >>ts_opt/${name}_ts.gjf
          #product
          echo "" >>ts_opt/${name}_ts.gjf
          echo "product" >>ts_opt/${name}_ts.gjf
          echo "" >>ts_opt/${name}_ts.gjf
          echo "${charge} ${mult}" >>ts_opt/${name}_ts.gjf
          awk 'NR>2 {print $0}' ${mtd_dir}/product.xyz >>ts_opt/${name}_ts.gjf
          #ts
          echo "" >>ts_opt/${name}_ts.gjf
          echo "ts" >>ts_opt/${name}_ts.gjf
          echo "" >>ts_opt/${name}_ts.gjf
          echo "${charge} ${mult}" >>ts_opt/${name}_ts.gjf
          awk 'NR>2 {print $0}' ${mtd_dir}/ts_qst3.xyz >>ts_opt/${name}_ts.gjf
          echo -e "\n\n" >>ts_opt/${name}_ts.gjf
          # freq
          awk '{print $0}' ${code_path}/hl_input_freq_xtb >ts_opt/${name}_freq.gjf
          # irc
          awk '{print $0}' ${code_path}/hl_irc_xtb >ts_opt/${name}_irc.gjf
          
         #  g09_input_xtb
         #  echo -e "$tsg09\n\n" > ts_opt/${name}_ts.dat
         #  echo -e "$freqg09\n\n" > ts_opt/${name}_freq.dat
          cp ${code_path}/xtb.sh ./
          cp ${code_path}/genxyz ./
          cp ${code_path}/extderi ./
          # qst3 get ts
          g09 <ts_opt/${name}_ts.gjf >ts_opt/${name}_ts.out 2>/dev/null && gauss_line=$(echo $LINENO)
          # get freq
          g09 <ts_opt/${name}_freq.gjf >ts_opt/${name}_freq.out 2>/dev/null
          
          ts_file=ts_opt/${name}_ts.out
          freq_file=ts_opt/${name}_freq.out
          chk_file=ts_opt/${name}_ts_freq.chk
          irc_file=ts_opt/${name}_irc.out
          cp mol.chk ${chk_file}
          ok=$(awk 'BEGIN{fok=0;ok=0;err=0};/Frequencies/{++nfreq;if($3<0 && $4>0 && nfreq==1) fok=1};/Error termi/{++err};END{if(err==0 && fok==1) ok=1; print ok}' $freq_file)
#         ok=$(awk 'BEGIN{ok=0};/Frequencies/{++nfreq;if($3<0 && $4>0 && nfreq==1) ok=1};END{print ok}' $file)
          if [ $ok -eq 1 ]; then
             xtb_energy_line=$(grep -r 'Recovered energy=' ${freq_file} )
             echo "${xtb_energy_line}" | awk '{print $3}' > tmp_gauss
             #get_energy_g09_${LLcalc}.sh $ts_file 1   > tmp_gauss
             get_freq_g09.sh ${freq_file} >> tmp_gauss
             # get irc
             g09 <ts_opt/${name}_irc.dat >${irc_file} 2>/dev/null
             mv mol.chk ${chk_file}
          else
             printf "     Pt%2s: failed-->EF algorithm was unable to optimize a TS\n" $npo
             mv mol.chk ${chk_file}
	     continue    
          fi

        elif [ "$program_opt" = "orca" ]; then
        #construct orca input file for second NEB-TS calculation
          orca_path=$(whereis orca | awk '{print $2}')
          orca_pltvib_path=$(whereis orca_pltvib | awk '{print $2}')
         if [ ! -d ${name}_NEB_TS_orca ];then
            mkdir ${name}_NEB_TS_orca 
         else
            rm -rf ${name}_NEB_TS_orca/*
          fi
          #mkdir ${name}_NEB_TS_orca
          current_path=$PWD
          # copy orca input file in scratch dir
          sed 's/carga/'$charge'/;s/mult/'$mult'/' ${code_path}/optTS_NEB_TS_orca >${name}_NEB_TS_orca/${name}_ts.inp
          #sed 's/carga/'$charge'/;s/mult/'$mult'/' ${code_path}/optTS_NEB_TS_orca >${name}_NEB_TS_orca/${name}_ts.inp
          #awk '{print $0}' ${code_path}/NEB_TS_orca >ts_opt/${name}_ts.inp
          cp ${mtd_dir}/reactant.xyz ${name}_NEB_TS_orca/ 
          cp ${mtd_dir}/product.xyz ${name}_NEB_TS_orca
          cp ${mtd_dir}/ts_qst3.xyz ${name}_NEB_TS_orca/tsguess.xyz
          #cp ts_opt/${name}_ts.inp ${name}_NEB_TS_orca/
          
          cd ${name}_NEB_TS_orca
          #echo $PWD
          #echo $batch
          ${orca_path}  ${name}_ts.inp >${name}_ts.out 2>/dev/null
          last_image_line=$(awk '/VIBRATIONAL FREQUENCIES/{print NR}' ${name}_ts.out |tail -1 )
          #last_image_line=$(grep -n "VIBRATIONAL FREQUENCIES" ${name}_ts.out |tail -1 |awk -v FS=":" '{print $1}')
          vb=$(awk 'NR>'${last_image_line}'' ${name}_ts.out |awk 'BEGIN{FS=":"}/imaginary mode/{print $1}')
         # vb=$(awk 'NR>'${last_image_line}'BEGIN{FS=":"}/imaginary mode/{print $1}' ${name}_ts.out)
          ${orca_pltvib_path}  ${name}_ts.hess ${vb} 2>/dev/null
          python3 ${code_path}/neb_snapshots.py ${name}_ts.interp 2>/dev/null
          cd ../

          mkdir ts_opt/${name}_ts_neb_frames
          cp ${name}_NEB_TS_orca/${name}_ts.inp ts_opt/
          cp ${name}_NEB_TS_orca/${name}_ts.hess ts_opt/
          cp ${name}_NEB_TS_orca/${name}_ts.hess*.xyz ts_opt/
          cp ${name}_NEB_TS_orca/${name}_ts.out ts_opt/
        #  cp ${name}_NEB_TS_orca/${name}_ts.xyz ts_opt/
        #  cp ${name}_NEB_TS_orca/${name}_ts_trj.xyz ts_opt/
         cp ${name}_NEB_TS_orca/${name}_ts.interp ts_opt/
         cp ${name}_NEB_TS_orca/neb_frames/* ts_opt/${name}_ts_neb_frames/
         cp ${name}_NEB_TS_orca/${name}_ts_MEP_trj.xyz ts_opt/
         cp ${name}_NEB_TS_orca/${name}_ts_NEB-TS_converged.xyz ts_opt/
          #rm -rf ${name}_NEB_TS_orca
         
          freq_file=ts_opt/${name}_ts.out
#          ok=$(awk 'BEGIN{fok=0;ok=0;err=0};/Frequencies/{++nfreq;if($3<0 && $4>0 && nfreq==1) fok=1};/Error termi/{++err};END{if(err==0 && fok==1) ok=1; print ok}' $freq_file)
#         ok=$(awk 'BEGIN{ok=0};/Frequencies/{++nfreq;if($3<0 && $4>0 && nfreq==1) ok=1};END{print ok}' $file)
          if [ ! -z ${vb} ]; then
             awk '/Final Gibbs free energy/{print $6}' ${freq_file} |tail -1 > tmp_gauss
             #get_energy_g09_${LLcalc}.sh $ts_file 1   > tmp_gauss
             vb_start1=$(awk '/VIBRATIONAL FREQUENCIES/{print NR}' ${freq_file} |tail -1 )
             vb_end1=$(awk '/NORMAL MODES/{print NR}' ${freq_file} |tail -1 )
             vb_start2=$(echo "${vb_start1}+4" | bc)
             vb_end2=$(echo "${vb_end1}-4" | bc)
             awk 'NR>'${vb_start2}' && NR<='${vb_end2}'{print $2}' ${freq_file} | sort -nu >> tmp_gauss
             #get_freq_g09.sh ${freq_file} >> tmp_gauss           
             # run irc
             orca_path=$(whereis orca | awk '{print $2}')
             if [ ! -d ${name}_irc_orca ];then
                mkdir ${name}_irc_orca 
             else
                rm -rf ${name}_irc_orca/*
             fi    
             sed 's/carga/'$charge'/;s/mult/'$mult'/' ${code_path}/optTS_irc_orca >${name}_irc_orca/${name}_irc.inp 
             #cp ts_opt/${name}_ts.xyz  ${name}_irc_orca/ts_for_irc.xyz
             cp ts_opt/${name}_ts_NEB-TS_converged.xyz  ${name}_irc_orca/ts_for_irc.xyz
             cd ${name}_irc_orca
             ${orca_path}  ${name}_irc.inp >${name}_irc.out 2>/dev/null
             cd ../
             cp ${name}_irc_orca/${name}_irc.inp ts_opt/
             cp ${name}_irc_orca/${name}_irc.out ts_opt/
             cp ${name}_irc_orca/${name}_irc_IRC_Full_trj.xyz ts_opt/
             
          else
             printf "     Pt%2s: failed-->NEB-TS algorithm was unable to optimize a TS\n" $npo
	     continue    
          fi
       fi

       fe="$(${code_path}/get_ts_properties.sh $freq_file $prog $tight)"
       fi="$(echo "$fe" | awk '{printf "%10.0f",$1}')"
       ei="$(echo "$fe" | awk '{printf "%14.6f",$2}')"
       if [[ ("$fi" -eq -1) ]]; then
          printf "     Pt%2s: failed-->Lowest real freq is negative\n" $npo
          continue
       elif [[ ("$fi" -eq -2) ]]; then
          printf "     Pt%2s: failed-->Sum of 2 lowest real freqs < 10cm-1\n" $npo
          continue
       elif [[ ("$fi" -eq -3) ]]; then
          printf "     Pt%2s: failed-->Stationary point is a minimum\n" $npo
          continue
       elif [[ ("$fi" -eq -4) ]]; then
          printf "     Pt%2s: failed-->EF algorithm was unable to optimize a TS\n" $npo
          continue
      #  elif (( $(echo "$ei > $emaxts" |bc -l) )); then
      #    printf "     Pt%2s: TS optimized but not added-->E=%20s > %20s \n" $npo $ei $emaxts
      #    continue
       fi
       if [[ ("$fi" -ge "$imag") ]]; then
          string="$(echo "$fe" | awk '{printf "%10.0f %10.4f %10.0f %10.0f %10.0f %10.0f",$1,$2,$3,$4,$5,$6}')"
# GLB added lock to tslist so that duplicate numbers are not created
          (
          flock -x 200 || exit 1
          if [ -f "$tslistll" ]; then
             ok=$(${code_path}/diff.sh $string $tslistll $prog)
             if [[ ("$ok" -eq "-1") ]]; then
                nt=$(awk '{nt=$2};END{print nt + 1}' $tslistll )
                namets=ts${nt}_${nb}
                printf "ts%5s%18s%70s traj= %4s Path= %10s\n" $nt $namets "$string" $i $nb  >> $tslistll
                if [ -f ${fileden} ]; then cp ${fileden} ${tsdirll}/${namets}.den ; fi 
                printf "     Pt%2s: TS optimized and added to ts list\n" $npo
                if [ "$program_opt" = "qcore" ]; then mv freq.molden $tsdirll/${namets}.molden ; fi
                if [ "$program_opt" = "mopac" ]; then get_NM_mopac.sh $tsdirll/${namets}.out $tsdirll/${namets} ; fi
                if [ "$program_opt" = "gaussian" ]; then 
                   cp ${ts_file} ${tsdirll}/${namets}_ts.out
                   cp ${freq_file} ${tsdirll}/${namets}_freq.out
                   cp ${chk_file} ${tsdirll}/${namets}_ts_freq.chk
                   cp ${irc_file} ${tsdirll}/${namets}_irc.out
                fi
                if [ "$program_opt" = "orca" ]; then 
                    cp ts_opt/${name}_ts.hess ${tsdirll}/${namets}_ts.hess
                    cp ts_opt/${name}_ts.out ${tsdirll}/${namets}_ts.out
                  #  cp ts_opt/${name}_ts.xyz ${tsdirll}/${namets}_ts.xyz
                    cp ts_opt/${name}_ts.hess.*.xyz ${tsdirll}/${namets}_ts_vib.xyz
                  #  cp ts_opt/${name}_ts.interp ${tsdirll}/${namets}_ts.interp
                    cp ts_opt/${name}_ts_NEB-TS_converged.xyz ${tsdirll}/${namets}_ts_NEB-TS.xyz
                    cp ts_opt/${name}_irc.out ${tsdirll}/${namets}_irc.out
                    cp ts_opt/${name}_irc_IRC_Full_trj.xyz ${tsdirll}/${namets}_irc_trj.xyz
                 fi
             else
                printf "     Pt%2s: TS optimized but not added-->redundant with ts %4s\n" $npo $ok
             fi
          else
             nt=1
             namets=ts${nt}_${nb}
             printf "ts%5s%18s%70s traj= %4s Path= %10s\n" $nt $namets "$string" $i $nb  >> $tslistll
             cp ${ts_file} ${tsdirll}/${namets}_ts.out
             cp ${freq_file} ${tsdirll}/${namets}_freq.out
             cp ${chk_file} ${tsdirll}/${namets}_ts_freq.chk
             cp ${irc_file} ${tsdirll}/${namets}_irc.out
             if [ -f ${fileden} ]; then cp ${fileden} ${tsdirll}/${namets}.den ; fi 
             printf "     Pt%2s: TS optimized and added to ts list\n" $npo
             if [ "$program_opt" = "qcore" ]; then mv freq.molden $tsdirll/${namets}.molden ; fi
             if [ "$program_opt" = "mopac" ]; then get_NM_mopac.sh $tsdirll/${namets}.out $tsdirll/${namets} ; fi
             if [ "$program_opt" = "gaussian" ]; then 
                cp ${ts_file} ${tsdirll}/${namets}_ts.out
                cp ${freq_file} ${tsdirll}/${namets}_freq.out
                cp ${chk_file} ${tsdirll}/${namets}_ts_freq.chk
                cp ${irc_file} ${tsdirll}/${namets}_irc.out
             fi
             if [ "$program_opt" = "orca" ]; then 
                 cp ts_opt/${name}_ts.hess ${tsdirll}/${namets}_ts.hess
                 cp ts_opt/${name}_ts.out ${tsdirll}/${namets}_ts.out
               #  cp ts_opt/${name}_ts.xyz ${tsdirll}/${namets}_ts.xyz
                 cp ts_opt/${name}_ts.hess.*.xyz ${tsdirll}/${namets}_ts_vib.xyz
               #  cp ts_opt/${name}_ts.interp ${tsdirll}/${namets}_ts.interp
                 cp ts_opt/${name}_ts_NEB-TS_converged.xyz ${tsdirll}/${namets}_ts_NEB-TS.xyz
                 cp ts_opt/${name}_irc.out ${tsdirll}/${namets}_irc.out
                 cp ts_opt/${name}_irc_IRC_Full_trj.xyz ${tsdirll}/${namets}_irc_trj.xyz                 
              fi
          fi
          ) 200>>${tslistll}.lock
          #if [ $mdc -eq 1 ]; then exit ; fi
          break
       else
          printf "     Pt%2s: TS optimized but not added-->imag=%4si cm-1 < %4si cm-1\n" $npo $fi $imag
       fi
    done
    done
  done
done


if [ $sampling -ne 30 ]; then 
   echo ""
   echo "END OF THE CALCULATIONS" 
   echo ""
fi

