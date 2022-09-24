set terminal png
set output "energy_metad.png"
set title "meta pic"
set xlabel 'step'
set ylabel 'E'
#plot 'TEST3.relaxscanact.dat' using 1:2
plot 'assoc_assoc1_opt_dyn1-md-ener.csv' using 1:7 with points

