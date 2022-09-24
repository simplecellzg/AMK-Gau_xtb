rm /root/docker_amk_g09/workdir/sub_system/llcaculation/assoc3/partial_opt/fort*
snapshots_mopac.sh coordir/assoc4_dyn1.xyz 10 |bbfs.exe >bbfs.out 
echo 23 >min.xyz 
echo '' >>min.xyz 
awk '{print $1,$2,$4,$6}' partial_opt/fort.105 >>min.xyz 
labels=$(awk '{if($3=="0") {printf "%s%s",sep,NR; sep=","}};END{print ""}' partial_opt/fort.105)
echo $labels
sed 's/labels/'"$labels"'/g;s/carga/0/' /opt/AutoMeKin/share/opt_frozen > partial_opt/pes_qcore 
qcore -f json partial_opt/pes_qcore > partial_opt/pes_qcore.out 
mv min_opt.xyz ts.xyz 
sed s/carga/0/ optTS > ts_opt/ts.dat 
qcore -f json ts_opt/ts.dat > ts_opt/ttsss.out