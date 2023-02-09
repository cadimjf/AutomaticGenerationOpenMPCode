set terminal png
set output "dat/dt_agos.png"
set xrange[0:3000]
set grid
plot "dat/addt_1_1.dat" u 1:3 t "addt H" w line, "dat/addt_1_2.dat" u 1:3 t "addt H 2t" w line  
set output
rep
reset


set terminal png
set output "dat/V_agos.png"
set grid
set xrange[0:3000]
plot "dat/addt_1_1.dat" u 1:2 t "v addt h" w line,"dat/rk2_1_1.dat" u 1:2 t "v rk f" w line, "dat/addt_1_2.dat" u 1:2 t "v addt h 2t" w line
set output
rep
reset

