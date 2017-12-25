reset
#custom
#multiplot
ix=5
iy=1
sx=4*ix
sy=3*iy
xlb=" "
ylb=" "
ml=0.6/sx
mb=0.4/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font ",18" size sx,sy
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics out nomirror
set ytics out nomirror
#set y2tics out nomirror
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
do for[i=0:(ix*iy-1)]{
	oxi=i%ix
	oyi=i/ix+1
	set origin ml+oxi*gx,mb+(iy-oyi)*gy
	if(oxi!=0){
		#set format y ""
	}
	if(oxi!=ix-1){
		unset y2tics
		set format y2 ""
	}
	if(oyi!=iy){
		set format x ""
	}
	if(i==(ix*iy-1)){
		set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque autotitle# columnhead
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	#set logscale x
	#set yrange [:]
	#set y2range [:]
	#set xrange [:]
	#plot "-" index i u (100-$1)/100.:2 with lp pt 7 ps 0.5 title "" , "-" index i u (100-$1)/100.:4 with lp pt 7 ps 0.5 title ""
	#plot "-" index i u (100-$1)/100.:7:8 with errorbars pt 7 ps 0.3 title ""
	#set pm3d map
	#set pm3d corners2color c2
	set cbrange [:]
	set size square
	set palette rgbformulae 22,13,-31
	#splot [-1:9.5][-1:10.5] "-" u 2:1:3 matrix with image
	plot [:][0:] "-" u 0:2*i+1:2*i+2 w errorbars ps 0.5 pt 7 title ""
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
