# find no ts assoc 
# find assoc dir have batch dir
amkcaculation='amkcaculation'
echo ""> assoclist
echo ""> assoc_batch_list
for i in `ls -d ./${amkcaculation}/assoc*/`
do
echo $i >> assoclist
cd $i 
if [ -d batch1 ]
then
echo $i >>../../assoc_batch_list
fi 
cd ../../
done

#rm dir and xyz file in assoc_batch_list
for line in  `cat assoc_batch_list`
do 
    j=`basename ${line}`
    xyzfile="./${amkcaculation}/$j.xyz"
    echo ${line}
    echo ${xyzfile}
    rm -r ${line}
    rm ${xyzfile}
done