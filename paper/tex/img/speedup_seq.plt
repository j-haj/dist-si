set terminal postscript eps enhanced color
set output "./output/speedup_seq.eps"
set ylabel "Speed-up"
set xtics ("Ackley" 1, "Rastrigin" 2, "Quadratic" 3)
set yrange [0:2]
set xrange [0:2.05]

plot "./data/speedup_seq.dat" using 2:xtic(1) with linespoints t "SoA", \
     "./data/speedup_seq.dat" using 3:xtic(1) with linespoints t "SoA-fast"