set terminal postscript eps enhanced color
set output "./output/all3steps.eps"
set style data histogram
set style fill solid
set xrange [-.25:2.5]
set xtics ("Ackley" 0, "Rastrigin" 1, "Quadratic" 2)
set xlabel "Test Function"
set ylabel "Number of steps"

plot "./data/all3steps.dat" using 2 title "AoS", \
     "./data/all3steps.dat" using 3 title "SoA", \
     "./data/all3steps.dat" using 4 title "SoA-fast"