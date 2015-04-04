set xrange [-1:501]
set yrange [-1:126]
set zrange [-1:4]
set xlabel"x"
set xlabel"y"
set xlabel"z"
set terminal wxt size 600, 600
set view equal xyz
set view 0, 0, 1.2
set ticslevel 0
do for [ii=0:9] {
set label 1 sprintf('frame = %d', ii) at 3, 1 right front
splot 'default.data.txt' every :::ii::ii using 1:2:3:($4*2.5):($5*2.5):($6*2.5) with vectors
pause 0.1
}
