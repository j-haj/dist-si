set terminal postscript eps enhanced color
set output "./output/N.eps"
set style data histogram
set style fill solid
set yrange [0:400]
set xlabel "Number of Particles"
set ylabel "Runtime (s)"

plot "./data/n_scale_par.dat" using 2:xtic(1) title "AoS", \
     "./data/n_scale_par.dat" using 4 title "SoA", \
     "./data/n_scale_par.dat" using 6 title "SoA-efficient"