temp='C:\Users\Administrator\Application Data\SSH\temp\temp.dat'
gap='C:\Users\Administrator\Application Data\SSH\temp\gap.dat'
raman='C:\Users\Administrator\Application Data\SSH\temp\raman.dat'
spect='C:\Users\Administrator\Application Data\SSH\temp\spect_0.dat'
energy='C:\Users\Administrator\Application Data\SSH\temp\energy_090.dat'
fermi='C:\Users\Administrator\Application Data\SSH\temp\fermi.dat'
energy='e:\projects\fortranprogram\data\energy_001.dat'
fermi='e:\projects\fortranprogram\data\fermi.dat'
temp='e:\projects\fortranprogram\data\temp.dat'
gap='e:\Projects\FortranProgram\DATA\gap.dat'
raman='e:\Projects\FortranProgram\DATA\raman.dat'
spect='e:\Projects\FortranProgram\DATA\spect_0.dat'
# set multiplot layout 2,1
# plot the energy band
set term wxt 0
unset key
set ytic 0.01
set size square
plot [:][-0.05:0] for[i=1:4] energy using 0:(column(2*i-1)):(column(i*2)) with points lt 1 pt 7 ps variable,\
	 for[i=1:4] energy using 0:(column(2*i-1)) with line lt 0
reset
# # plot for[i=1:4] energy using 0:(column(2*i-1)-(column(i*2)*0.01)):(column(2*i-1)+(column(i*2)*0.01)) with filledcurves lt 1,\
	# # for[i=1:4] energy using 0:(column(2*i-1)) with line lt 0
set term wxt 1
# # plot the fermi surface
# unset xtic
# unset key
# set size square
# set palette rgbformulae 22,13,-31
# set cbrange [0:1]
# set pm3d map
# set pm3d interpolate 0,0
# splot [0:3.1416][0:3.1416] fermi

set term wxt 2
# plot the temprature dependence
set xtic 50
plot temp using 1:2 with linespoints pt 7 lw 5 , temp using 1:3 with linespoints pt 7 lw 5

set term wxt 3
# # plot gap
# set xtic 5
# plot [:] gap using 1:2 with line
set xtic 0.01
set term wxt 4
plot spect with line

set term wxt 5
# # plot Raman
# set xtic 0.02
# plot [:][:] raman using 1:2 with line axis x1y1, raman using 1:3 with line axis x1y2
# plot [:][:] raman using 1:4 with line axis x1y1, raman using 1:5 with line axis x1y2