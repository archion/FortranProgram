energy='e:\projects\fortranprogram\data\energy.dat'
fermi='e:\projects\fortranprogram\data\fermi.dat'
# energy='C:\Users\Administrator\Application Data\SSH\temp\energy.dat'
# fermi='C:\Users\Administrator\Application Data\SSH\temp\fermi.dat'
set multiplot layout 1,2
# plot the energy band
unset key
set size square
set xtic 128
# plot for[i=1:4] energy using 0:(column(2*i-1)):(column(i*2)*2+0.1) with points lt 1 pt 7 ps variable
plot for[i=1:4] energy using 0:(column(2*i-1)-(column(i*2)*0.2)):(column(2*i-1)+(column(i*2)*0.2)) with filledcurves lt 1,\
	for[i=1:4] energy using 0:(column(2*i-1)) with line lt 0

# plot the fermi surface
unset xtic
unset key
set size square
set palette rgbformulae 22,13,-31
set cbrange [0:1]
set pm3d map
set pm3d interpolate 0,0
splot [0:3.1416][0:3.1416] fermi