#/bin/bash
#clean dir
#rm *.log *.dat *.xyz *.pid
for i in `ls -d */`
do
rm -rf $i
done
cwd=$PWD
#define code path
export code_path=${cwd}
#export code_path=/root/docker_amk_g09/workdir_nc_oh/sub_system 
llcaculation='llcaculation'
mkdir ./${llcaculation}

for dddir in `ls assoc*.dat`
do
assocdir=`basename $dddir .dat`
#assocdir=${assoc_dir}
mkdir ./${assocdir}
chmod -R 750 ../code/
cp ../code/* ${cwd}/
cp ../code/* ${cwd}/${assocdir}/
#step 1: assoc TNT with OH-
#enter assocdir
cd ${cwd}/${assocdir}/
#amk.sh vdW.dat > vdW.log
${code_path}/amk1.sh ${assocdir}.dat > ${assocdir}.log
#enter assoc_fragA_fragB dir
cd ${cwd}/${assocdir}/assoc*/
#cp ../out2xyz.sh ./
#step 2: convert out file to xyz file
#if lowlevel is mopac
# for i in `ls DISCNT*`
# do
#     inputname=${i:7}
#     get_geom_mopac.sh $i >${inputname}.xyz
#     echo ${inputname}.xyz
# done
# cd ${cwd}
# if lowlevel is qcore
cat /dev/null >${cwd}/xyz_tmp
#for i in `ls assoc*opt.xyz`
for i in `ls assoc*.xyz`
#for i in `ls assoc*.qcore`
do
    inputname=`basename $i .xyz`
    #inputname=`basename $i .qcore`
    echo ${inputname}.xyz >>${cwd}/xyz_tmp
done
cat ${cwd}/xyz_tmp
cd ${cwd}



#step 3: run amk lowlevel get opt results

cp ${cwd}/amk_template.dat ${cwd}/${llcaculation}/
#cp ${cwd}/${assocdir}/assoc*/*.xyz ./${llcaculation}/
#qcore
for line in $(cat ${cwd}/xyz_tmp)
do
    cp ${cwd}/${assocdir}/assoc*/${line} ./${llcaculation}/${assocdir}_${line}
done
    rm ./${llcaculation}/*_opt.xyz
done
cd ${cwd}/${llcaculation}
#ll caculation for every assocX.xyz file
for j in `ls *.xyz`
do
    inputname=`basename $j .xyz`
    mkdir ${inputname}
    cp $j ./${inputname}/
    cp ./amk_template.dat ./${inputname}/
    cd ./${inputname}/
    #replace molecular name by xyz file name 
    awk '{sub(/amk_mol/,"'$inputname'"); print$0}' amk_template.dat > ${inputname}.dat
    #nohup llcalcs.sh ${inputname}.dat 10 3 8 >llcalcs.log 2>&1 & 
    nohup ${code_path}/llcalcs_TS.sh ${inputname}.dat 1 1 1 >llcalcs_TS.log 2>&1 & 
    cd ${cwd}/${llcaculation}
done

# #cp -r ./ ../backup/
# # hl caculation