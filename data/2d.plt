reset
#custom
#multiplot
ix=1
iy=1
sx=5*ix
sy=5*iy
xlb=""
ylb=""
ml=0.2/sx
mb=0.2/sy
mr=0.2/sx
mt=-0.5/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font "Times Bold,18" size sx,sy enhanced
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics out nomirror scale 0.5 offset 0,0.5
set ytics out nomirror scale 0.5 offset 0.5,0
set border front lw 5
#set grid xtics ytics
#
set size gx,gy
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
unset xtics
unset ytics
do for[i=0:(ix*iy-1)]{
	oxi=i%ix
	oyi=i/ix+1
	set origin ml+oxi*gx,mb+(iy-oyi)*gy
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
	set pm3d map
	#set pm3d corners2color c1
	set size square 
	set palette rgbformulae 22,13,-3
	#set pm3d interpolate 0,0
	#set view 70,340
	#set hidden3d
	#set log cb
	#set cbrange [word("0.5 -1",i%2+1):word("1 1",i%2+1)]
	#set logscale x
	#set logscale y
	#splot [:][:][:] "-" u 1:3:4:($5*200) w lp pt 7 ps var,"-" u 1:3:4:(sprintf("{/=8 %0.3f}",column(4))) w labels 
	plot "-" u 1:(-$2):3 w p pt 7 palette
	#splot [:][:] "-" index i using 1:2:($3*((-1)**i)**($1+$2)) matrix w image
	#splot [-0.5:19.5][-0.5:19.5] "-" index i using 1:2:3 matrix w image
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
