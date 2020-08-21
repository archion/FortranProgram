reset
#custom
#multiplot
ix=1
iy=1
sx=4*ix
sy=3.7*iy
#sx=4*ix
#sy=3.72*iy
xlb=" "
ylb=" "
ml=0.1/sx
mb=0.1/sy
mr=0.5/sx
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
#set y2tics out nomirror scale 0.5 offset -1,0
set border front lw 5
unset xtics
unset ytics
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
	#set cbtics 0.05
	#set size square
	set palette rgbformulae 22,13,-31
	if(oxi!=0){
		set format y ""
	}
	if(oyi!=iy){
		set format x ""
	}
	#if(i==(ix*iy-1)){
		set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 samplen 0.5 opaque# autotitle
	#}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	#set y2range [0:]
	#set logscale x
	#set logscale y
	set style arrow 1 head noborder size screen 0.025,20,130 lw 7 lt 8
	#plot [:][:] "-" index 0 u 1:2:(20*(($4+$5)/2-0.8)) w p pt 7 ps var lw 5 lc rgb "#aeafac" title "" , "-" index 1 using ($1-$4/2):($2-$5/2):4:5:7 w vectors nohead lw 8 palette title "", "-" index 0 u ($1-(1*($4-$5))):($2-(1*($4-$5))):(2*($4-$5)):(2*($4-$5)) w vectors arrowstyle 1 title ""
	#set cbrange [0.8:0.7]
	#set cbrange [1.2:1.3]
	#plot [:][:] "-" index 0 using "x":"y":(column("n")) w p pt 7 palette title ""
	set cbrange [:]
	plot [:][:] "-" index 0 using "x":"y":(column("S")) w p pt 7 palette title ""
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
