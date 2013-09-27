dsc='../data/dsc.dat'
sdw='../data/sdw.dat'
cdw='../data/cdw.dat'
splot dsc with line
#i=0
#do for [name in "dsc sdw cdw"]{
	#i=i+1
	#set term x11 i
	#set palette rgbformulae 22,13,-31
	#set pm3d
	#set pm3d interpolate 0,0
	#splot name 
#}
