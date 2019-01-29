set terminal postscript eps enhanced color
set output "./output/dim.eps"
set style data histogram
set style fill solid
set xlabel "Dimension"
set ylabel "Runtime (s)"

plot "./data/dim_scale_par.dat" using 2:xtic(1) title "AoS", \
     "./data/dim_scale_par.dat" using 5 title "SoA", \
     "./data/dim_scale_par.dat" using 8 title "SoA-efficient"