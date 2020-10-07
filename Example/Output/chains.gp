reset



set term postscript eps color blacktext "Courier, 12"
set output '-change.eps'

set xlabel "Steps"
set ylabel "x1"
#set yr [16500:17500]
#set xr [0:5000]

unset key

set multiplot
plot 'mcmc_fit.dat' using 0:9 with lines lc "black" lw 4


unset multiplot