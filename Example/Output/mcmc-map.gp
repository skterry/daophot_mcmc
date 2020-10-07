reset



set term postscript eps color blacktext "Courier, 12"
set output 'position-change.eps'

set xlabel "x"
set ylabel "y"
set xr [879:908]
set yr [811:836]


unset key

set multiplot

#plot 'uvp2dual_fit.07.chi2change' every 1000 using 0:8 smooth frequency with steps lc "black" lw 6
plot 'mcmc_fit.dat' using 1:2 with lines lc "blue" lw 0.5 #Source Star
plot 'mcmc_fit.dat' using 3:4 with lines lc "red" lw 0.5 #Lens Star

#plot 'source.txt' using 1:2:3:4 with xyerrorbars ps 0.5 lc rgb "black"
#plot 'lens.txt' using 1:2:3:4 with xyerrorbars ps 0.5 lc rgb "black"


unset multiplot