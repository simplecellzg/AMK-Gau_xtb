#!/bin/bash
# This script analyzes a geometry optimization file from ORCA.
# Usage: type in the script and give the output file as the first argument
#
# 18.Nov.2013 IS first version

# Check if the output file was passed as the first argument
if [ -z $1 ]; then
    echo "No output file"
    exit 0
fi
file=$1

# Extract the information from a geometry optimization
grep "FINAL SINGLE POINT ENERGY" $file | awk '{print $5}' | nl -v 0 > $file.energy
grep "     Energy change " $file | awk '{print $3 " " $4}' | nl > $file.e.change
grep "MAX step" $file | awk '{print $3 " " $4}' | nl > $file.MAX.step
grep "RMS step" $file | awk '{print $3 " " $4}' | nl > $file.RMS.step
grep -A2 "Energy change" $file | grep "MAX gradient" | awk '{print $3 " " $4}' | nl > $file.MAX.gradient
grep -A2 "Energy change" $file | grep "RMS gradient" | awk '{print $3 " " $4}' | nl > $file.RMS.gradient

rm -f $file.gnuplot

echo "#!/usr/bin/gnuplot -persist" >> $file.gnuplot
echo "set size 1,1" >> $file.gnuplot
echo "set multiplot" >> $file.gnuplot
echo "set xlabel 'optimization step'" >> $file.gnuplot
echo "" >> $file.gnuplot
echo "#first" >> $file.gnuplot
echo "set size 0.5,0.33" >> $file.gnuplot
echo "set origin 0,0" >> $file.gnuplot
echo "set ylabel 'MAX step'" >> $file.gnuplot
echo "plot '$file.MAX.step' u 1:2 w lp t 'MAX step', '$file.MAX.step' u 1:3 w l lw 1 lc 3 t ''" >> $file.gnuplot
echo "" >> $file.gnuplot
echo "#second" >> $file.gnuplot
echo "set origin 0,0.33" >> $file.gnuplot
echo "set ylabel 'RMS step'" >> $file.gnuplot
echo "plot '$file.RMS.step' u 1:2 w lp t 'RMS step', '$file.RMS.step' u 1:3 w l lw 1 lc 3 t ''" >> $file.gnuplot
echo "" >> $file.gnuplot
echo "#third" >> $file.gnuplot
echo "set origin 0.5,0" >> $file.gnuplot
echo "set ylabel 'MAX gradient'" >> $file.gnuplot
echo "plot '$file.MAX.gradient' u 1:2 w lp t 'MAX gradient', '$file.MAX.gradient' u 1:3 w l lw 1 lc 3 t ''" >> $file.gnuplot
echo "" >> $file.gnuplot
echo "#fourth" >> $file.gnuplot
echo "set origin 0.5,0.33" >> $file.gnuplot
echo "set ylabel 'RMS gradient'" >> $file.gnuplot
echo "plot '$file.RMS.gradient' u 1:2 w lp t 'RMS gradient', '$file.RMS.gradient' u 1:3 w l lw 1 lc 3 t ''" >> $file.gnuplot
echo "" >> $file.gnuplot
echo "#fifth" >> $file.gnuplot
echo "set origin 0,0.66" >> $file.gnuplot
echo "set ylabel 'Energy'" >> $file.gnuplot
echo "plot '$file.energy' u 1:2 w lp t 'Energy'" >> $file.gnuplot
echo "" >> $file.gnuplot
echo "#sixth" >> $file.gnuplot
echo "set origin 0.5,0.66" >> $file.gnuplot
echo "set ylabel 'Energy change'" >> $file.gnuplot
echo "plot '$file.e.change' u 1:2 w lp t 'Energy change', '$file.e.change' u 1:3 w l lw 1 lc 3 t ''" >> $file.gnuplot
echo "" >> $file.gnuplot
echo "unset multiplot" >> $file.gnuplot

gnuplot --persist $file.gnuplot
