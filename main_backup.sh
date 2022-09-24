#/bin/bash
#clean dir
rm *.log *.dat *.xyz *.pid
for i in `ls -d */`
do
rm -rf $i
done
cwd=$PWD
assocdir='assoc'
mkdir ./${assocdir}
chmod -R 750 ../code/
cp ../code/* ${cwd}/
cp ../code/* ${cwd}/${assocdir}/
#step 1: assoc TNT with OH
# enter assocdir
cd ${cwd}/${assocdir}/
amk.sh assoc.dat > assoc.log
# enter assoc_fragA_fragB dir
cd ${cwd}/${assocdir}/assoc*/
#cp ../out2xyz.sh ./
#step 2: convert out file to xyz file
for i in `ls assoc*.out`
do
    inputname=`basename $i .out`
    get_geom_mopac.sh $i >${inputname}.xyz
    echo ${inputname}.xyz
done
#./out2xyz.sh >xyzfile.log
cd ${cwd}

#step 3: run amk lowlevel get opt results
llcaculation='llcaculation'
mkdir ./${llcaculation}
cp ${cwd}/amk_template.dat ${cwd}/${llcaculation}/
cp ${cwd}/${assocdir}/assoc*/*.xyz ./${llcaculation}/
cd ${cwd}/${llcaculation}
# ll caculation for every assocX.xyz file
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
    nohup ${cwd}/llcalcs_TS.sh ${inputname}.dat 2 2 8 >llcalcs_TS.log 2>&1 & 
    cd ${cwd}/${llcaculation}
done

#cp -r ./ ../backup/
# hl caculation
:<<!
for j in `ls *.xyz`
do
inputname=`basename $j .xyz`
cd ./${inputname}/
nohup hlcalcs.sh ${inputname}.dat 1 >hlcalcs.log 2>&1 &
cd ../
done
!


