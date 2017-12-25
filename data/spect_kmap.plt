reset
#custom
#multiplot
ix=1
iy=1
sx=3*ix
sy=3*iy
xlb=" "
ylb=" "
ml=0.5/sx
mb=0.2/sy
mr=0.2/sx
mt=0.2/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
#set term wxt
set term eps font "Times Bold,18" size sx,sy enhanced
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
#set xtics out ("" 0, "" pi, "" -pi) nomirror scale 1.5 offset 0,0.5
#set ytics out ("" 0, "" pi, "" -pi) nomirror scale 1.5 offset 0.5,0
set xtics out nomirror scale 1.5 offset 0,0.5
set ytics out nomirror scale 1.5 offset 0.5,0
set border front lw 5
#set grid xtics ytics
#
set size gx,gy
#set size ratio -1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set multiplot
#set mxtics 0
#set mytics 0
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
	#set arrow from pi,0 to 0,pi nohead lc rgb 'red'
	#set arrow from 0,pi to -pi,0 nohead lc rgb 'red'
	#set arrow from -pi,0 to 0,-pi nohead lc rgb 'red'
	#set arrow from 0,-pi to pi,0 nohead lc rgb 'red'
	#set style fill  transparent solid 0.35 noborder
	#splot [-pi-0.001:pi+0.001][-pi-0.001:pi+0.001] for[j=0:2:2] "-" index 0  using ((j-1)*column(1+j)):((j-1)*column(2+j)):((j-1)*column(5+j/2)):(10*(-$8)) w p pt 7 ps var
	#plot [-pi-0.001:pi+0.001][-pi-0.001:pi+0.001] for[j=0:2:2] "-" index 0  using ((j-1)*column(1+j)):((j-1)*column(2+j)):((-$8)) w circle
	#plot [-pi-0.001:pi+0.001][-pi-0.001:pi+0.001] for[j=-1:1] for[k=-1:1] "-" index 0  using ($1+pi*((j*k==0 && (j!=0 || k!=0)) ? 1/0 : j)):($2+pi*((j*k==0 && (j!=0 || k!=0)) ? 1/0 : k)):(0.1*$4) w p pt 7 ps 1 palette
	plot [-pi:pi][-pi:pi] for [j=0:2:2] "-"  index 0 using (-(j-1)*column(1+j)):(-(j-1)*column(2+j)):(-8*$8) w p pt 6 lw 3 ps var
	#plot for [j=0:1] "-" index 0 u (column(9+j)):(-(0*2-1)*column(5+j)):(-10*$8) w p pt 6 lw 1 ps var
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
