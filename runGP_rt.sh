#!/bin/bash
##The last job is the tors scan
cwd=$PWD
#code_path=$5
#echo $code_path
if [ $1 -eq $4 ];then
   batch=torsion
else
   batch=batch$1
fi
#Make batch directory
if [ -d ${batch} ]; then  rm -r ${batch} ; fi
mkdir ${batch}
#copy neccesary files.
cp amk.dat $2.xyz ${batch}
echo "tsdirll $3/tsdirLL_$2"  >> ${batch}/amk.dat
##copy frags only if sampling is vdw
copyfrags=$(awk 'BEGIN{cf=0};{if($1=="sampling" && $2=="vdW") cf=1};END{print cf}'  amk.dat)
if [ $copyfrags -eq 1 ]; then
   frA=$(awk '{if($1=="fragmentA") print $2}' amk.dat)
   frB=$(awk '{if($1=="fragmentB") print $2}' amk.dat)
   cp ${frA}.xyz ${frB}.xyz ${batch}
fi
##Launch the cals
if [ $1 -eq $4 ]; then
   (cd ${batch} && tors.sh > scan_tors.log)
else
   (cd ${batch} && ${code_path}/amk1_rt.sh amk.dat >amk1_rt.log)
fi
