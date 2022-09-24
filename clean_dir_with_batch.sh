# find no ts assoc 
# find assoc dir have batch dir
cwd=$PWD
llcaculation='llcaculation'
cat /dev/null > assoclist
cat /dev/null > assoc_batch_list
for i in `ls -d ${cwd}/${llcaculation}/assoc*/`
do
    echo $i >> assoclist
    cd $i 
    if [ -d batch1 ]
    then
        echo $i >>${cwd}/assoc_batch_list
    fi 
    cd ${cwd}
done

#rm dir and xyz file in assoc_batch_list
for line in  `cat assoc_batch_list`
do 
    j=`basename ${line}`
    xyzfile="./${llcaculation}/$j.xyz"
    echo ${line}
    echo ${xyzfile}
    rm -r ${line}
    rm ${xyzfile}
done