set terminal eps
set output "ackley.dps"

set xrange [-2:2]
set yrange [-2:2]
set ticslevel 0
set isosample 150
set key off
set contour base
splot -20 * exp(-0.2 * sqrt(0.5 * (x**2 + y**2))) - exp(0.5*(cos(2 * 3.14159 * x)) + cos(2 * 3.14159 * y)) + exp(1) + 20