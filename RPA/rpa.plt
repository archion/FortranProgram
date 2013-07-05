rpa='C:\Users\Administrator\Application Data\SSH\temp\rpa.dat'
rpa='e:\projects\fortranprogram\data\rpa.dat'

# plot the fermi surface
# unset key
# set size square
# set palette rgbformulae 22,13,-31
# set cbrange [0:1]
# set pm3d map
# set pm3d interpolate 0,0
# splot rpa

unset key
set size square
plot for[i=1:4] rpa using 0:(column(2*i-1)):(column(i*2)) with points lt 1 pt 7 ps variable,\
	 for[i=1:4] rpa using 0:(column(2*i-1)) with points lt 0
# # plot for[i=1:4] energy using 0:(column(2*i-1)-(column(i*2)*0.01)):(column(2*i-1)+(column(i*2)*0.01)) with filledcurves lt 1,\
	# # for[i=1:4] energy using 0:(column(2*i-1)) with line lt 0