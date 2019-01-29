set terminal eps
set output "par_rastrigin_d100.eps"
set xlabel "Number of particles"
set ylabel "Runtime (s)"

plot "./data/par_rastrigin_100d.dat" using 1:2 with linespoint, \
     "./data/par_rastrigin_100d.dat" using 1:3 with linespoint, \
     "./data/par_rastrigin_100d.dat" using 1:4 with linespoint
