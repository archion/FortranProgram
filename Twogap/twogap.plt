energy='C:\Users\Administrator\Application Data\SSH\temp\energy.dat'
fermi='C:\Users\Administrator\Application Data\SSH\temp\fermi.dat'
temp='C:\Users\Administrator\Application Data\SSH\temp\temp.dat'
gap='C:\Users\Administrator\Application Data\SSH\temp\gap.dat'
raman='C:\Users\Administrator\Application Data\SSH\temp\raman.dat'
energy='e:\projects\fortranprogram\data\energy.dat'
fermi='e:\projects\fortranprogram\data\fermi.dat'
temp='e:\projects\fortranprogram\data\temp.dat'
gap='e:\Projects\FortranProgram\DATA\gap.dat'
raman='e:\Projects\FortranProgram\DATA\raman.dat'
set multiplot layout 1,3
# plot the energy band
unset key
# set size square
# set yrange [-0.1:0.1]
# set xrange [0:1500]
# plot for[i=1:4] energy using 0:(column(2*i-1)):(column(i*2)) with points lt 1 pt 7 ps variable,\
	 # for[i=1:4] energy using 0:(column(2*i-1)) with points lt 0
# # # plot for[i=1:4] energy using 0:(column(2*i-1)-(column(i*2)*0.01)):(column(2*i-1)+(column(i*2)*0.01)) with filledcurves lt 1,\
	# # # for[i=1:4] energy using 0:(column(2*i-1)) with line lt 0

# # plot the fermi surface
# unset xtic
# unset key
# set size square
# set palette rgbformulae 22,13,-31
# set cbrange [0:1]
# set pm3d map
# set pm3d interpolate 0,0
# splot [0:3.1416][0:3.1416] fermi

# # plot the temprature dependence
# set xtic 50
# plot temp using 1:2 with linespoints pt 7 lw 5 , temp using 1:3 with linespoints pt 7 lw 5

# # plot gap
# set xtic 5
# plot [:] gap using 1:2 with line

# plot Raman
plot [:][0:] raman using 1:2 with line lw 3 
plot [:][0:] raman using 1:3 with line lw 3
plot raman using 1:4 with line lw 3
