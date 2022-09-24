for i in `ls assoc*.out`
do
inputname=`basename $i .out`
get_geom_mopac.sh $i >${inputname}.xyz
echo ${inputname}.xyz
done