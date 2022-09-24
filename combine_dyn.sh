#/bin/bash
cwd=$PWD
export code_path=${cwd}
workdir='./combine_dyn/'
rm -r ${workdir}
mkdir ${workdir}

#mkdir ${workdir}/tsdirLL_assoc/backup

num=0
#for i in `ls -d ${code_path}/combine/assoc*/`
for i in `ls -d ${cwd}/llcaculation/assoc*/`
do 
    rows=0
    basename=`basename $i`
    cp ${i}batch1/coordir/${basename}_dyn1.xyz ${workdir}
done
