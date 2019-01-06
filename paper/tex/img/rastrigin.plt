set terminal eps
set output "rastrigin.dps"

set xrange [-1.5:1.5]
set yrange [-1.5:1.5]
set ticslevel 0
set isosample 100
set key off
set contour base
splot (20 + x**2 - 10 * cos( 2 * 3.14159 * x) + y**2 - 10 * cos(2 * 3.14159 * y))