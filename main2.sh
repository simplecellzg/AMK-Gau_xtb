#/bin/bash
#clean wrong dir
cwd=$PWD
#./clean_dir_with_batch.sh 
cd ${cwd}
#combine TS results from different assoc dir
./combineTS.sh
cd ${cwd}/combine/


echo "       IRC calcs    "
irc.sh > /dev/null
:<<!
convergence=0
echo "   Optimizing minima"
min.sh  > /dev/null
echo "    rxnetwork calcs "
rxn_network.sh >/dev/null
echo "       KMC calcs    "
kmc.sh > /dev/null
echo "Making final folder: FINAL_LL"
final.sh > /dev/null
echo ""
echo "END OF THE CALCULATIONS"
echo ""
!