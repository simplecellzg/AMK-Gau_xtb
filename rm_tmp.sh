for i in `ls -d assoc*/`
do
echo $i
cd ${i}batch1
rm *_NEB_orca/*.tmp
rm ts*/*.tmp
rm produnct*/*.tmp
ls *_NEB_orca/*.tmp
ls ts*/*.tmp
ls produnct*/*.tmp
cd ../../
done