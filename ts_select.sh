mkdir ts_select
for i in $(<ts_select.txt)
do
echo $i
cp ts${i}_*.inp ./ts_select
cp ts${i}_*.xyz ./ts_select
done