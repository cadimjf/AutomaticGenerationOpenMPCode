set terminal png
set output "dat/V_agos.png"
set grid
set xrange[0:1]
plot "dat/euler_1_1.dat" u 1:2 t "euler" w line,"dat/rk2_1_1.dat" u 1:2 t "rk2" w line,"dat/addt_1_1.dat" u 1:2 t "addt" w line, "dat/cvode_1_1.dat" u 1:2 t "cvode" w line
set output
rep
reset


set terminal png
set output "dat/dt_agos.png"
set grid
set xrange[0:1]
plot "dat/euler_1_1.dat" u 1:3 t "euler" w line,"dat/rk2_1_1.dat" u 1:3 t "rk2" w line,"dat/addt_1_1.dat" u 1:3 t "addt" w line,"dat/addt_2_1.dat" u 1:3 t "pycml" w line,"dat/cvode_1_1.dat" u 1:3 t "CVODE" w line
set output
rep
reset
set output
rep
reset

set terminal png
set output "dat/cvode.png"
set grid
set xrange[0:1]
plot "dat/cvode_1_1.dat" u 1:2 t "agos" w line,"dat/cvode_1_1.dat" u 1:2 t "pycml" w line
set output
rep
reset
