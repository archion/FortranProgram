rpa='../data/rpa'
check='../data/check.dat'
i=0
do for [name in "001 002 003 004 005 006 007 008"]{
	i=i+1
	set term x11 i
	set multiplot layout 1,2
	unset key
	set size square
	set palette rgbformulae 22,13,-31
	# set cbrange [0:1]
	set pm3d map
	set pm3d interpolate 0,0
	splot rpa."_".name.".dat" using 1:2:3
	splot rpa."_".name.".dat" using 1:2:4
	unset multiplot
}

#set term wxt 1
#set term x11 1
#unset key
#set size square
#plot for[i=1:4] check using 0:(column(2*i-1)):(column(i*2)) with points lt 1 pt 7 ps variable,\
	 #for[i=1:4] check using 0:(column(2*i-1)) with points lt 0
# plot for[i=1:4] energy using 0:(column(2*i-1)-(column(i*2)*0.01)):(column(2*i-1)+(column(i*2)*0.01)) with filledcurves lt 1,\
	# for[i=1:4] energy using 0:(column(2*i-1)) with line lt 0
