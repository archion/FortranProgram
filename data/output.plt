set pm3d map
set pm3d corners2color c1
set size square
set palette gray
splot "output.dat" matrix
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
