set terminal postscript eps enhanced color font "Helvetica,24"
set output "./output/baseline.eps"
set style data histogram
set style fill solid
set xlabel "Test Function"
set ylabel "Runtime (s)"

plot "./data/baseline_seq.dat" using 2:3:xtic(1) title "AoS" with boxerrorbars, \
     "./data/baseline_seq.dat" using 4:5:xtic(1) title "SoA" with boxerrorbars, \
     "./data/baseline_seq.dat" using 6:7:xtic(1) title "SoA-efficient" with boxerrorbars