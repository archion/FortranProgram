reset
#custom
#multiplot
ix=1
iy=1
sx=5*ix
sy=3*iy
xlb=" "
ylb=" "
ml=0.1/sx
mb=0.1/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font "Times Bold,18" size sx,sy enhanced
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics out nomirror scale 0.5 offset 0,0.5
set ytics out nomirror scale 0.5 offset 0.5,0
#set border front lw 5
#set grid xtics ytics
#
set size gx,gy
set size ratio -1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set multiplot
set mxtics 5
set mytics 5
unset xlabel
unset ylabel
unset key
do for[i=0:(ix*iy-1)]{
	oxi=i%ix
	oyi=i/ix+1
	set origin ml+oxi*gx,mb+(iy-oyi)*gy
	set pm3d map
	set pm3d corners2color c2
	#set logscale cb
	set cbrange [:]
	#set size square
	set palette rgbformulae 22,13,-31
	if(oxi!=0){
		set format y ""
	}
	if(oyi!=iy){
		set format x ""
	}
	if(i==(ix*iy-1)){
		#set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 samplen 0.5 opaque# autotitle
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	splot [:][:] for[j=0:0] \
	"-" index 6 using ($1-$4/2):($2-$5/2):($3-$6/2):(abs($4)+abs($5)+abs($6)<0.001?0.1:$4):($5):($6):($7) w vectors nohead palette lw 5 title "",\
	#"-" index 2 using 1:2:3:4 w point pt 7 ps 0.5 palette title "",\
	#"-" index 0 using 1:2:3:4 w point pt 7 ps 0.5 palette title "",\
	#"-" index 1 using ($1-$4/2*sgn($7)):($2-$5/2*sgn($7)):($3-$6/2*sgn($7)):(abs($4)+abs($5)+abs($6)<0.001?0.1:$4*sgn($7)):($5*sgn($7)):($6*sgn($7)):(abs($7)) w vectors nohead palette lw 5 title "",\
	#"-" index 0 using 1:2:3 w point pt 7 ps 0.5 title "",\
	#"-" index 4 every :::j::j using ($1-$4/2):($2-$5/2):($3-$6/2):(abs($4)+abs($5)+abs($6)<0.001?0.1:$4):5:6:7 w vectors head palette lw 3 title ""
	#"-" index 5 every :::j::j using ($1-$4/2*sgn($7)):($2-$5/2*sgn($7)):($3-$6/2*sgn($7)):(abs($4)+abs($5)+abs($6)<0.001?0.1:$4*sgn($7)):($5*sgn($7)):($6*sgn($7)):(abs($7)) w vectors head palette lw 2 title ""
	#"-" index 5 every :::j::j using ($1-$4/2*sgn($8)):($2-$5/2*sgn($8)):($3-$6/2*sgn($8)):(abs($4)+abs($5)+abs($6)<0.001?0.1:$4*sgn($8)):($5*sgn($8)):($6*sgn($8)):(abs($8)) w vectors head palette lw 2 title ""
	#"-" index 5 every :::j::j using ($1-$4/2):($2-$5/2):($3-$6/2):(abs($4)+abs($5)+abs($6)<0.001?0.1:$4):5:6:8 w vectors head palette lw 3 title ""
	#"-" index 4 every :::j::j using ($1-$4/2):($2-$5/2):($3-$6/2):4:5:6 w vectors head filled lt 2 title ""
	#"-" index 4 using 1:2:3:4 w point pt 5 ps 0.3 palette title "",\
	#"-" index 0 using 1:2:3 w point pt 7 ps 0.2 lc 0 title "",\
	#"-" index 4 using 1:2:3:($4*$5) w point pt 5 ps 0.3 palette title "",\
	#"-" index 2 every :::j::j using 1:2:3 w point pt 5 ps 0.3 palette title "",\
	#"-" index 3 every :::j::j using 1:2:3 w point pt 5 ps 0.3 palette title "",\
	#"-" index 0 every :::j::j using 1:2:4 w point pt 7 ps 1 palette title "",\
	#"-" index 0 every :::j::j using 1:2:3 w labels font ",7" tc "white" title "",\

	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
