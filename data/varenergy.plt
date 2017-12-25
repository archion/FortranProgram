reset
#custom
#multiplot
ix=1
iy=1
sx=4*ix
sy=3.7*iy
#sx=4*ix
#sy=3.72*iy
xlb="x"
ylb=" "
ml=0.5/sx
mb=0.4/sy
mr=0.2/sx
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
	set cbrange [:-0.39]
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
	if(i==0){
		plot [:][:] "-" u ($2-(-($8)/sqrt(($8)**2+($9)**2)*0.015)/2):($3-(-($9)/sqrt(($8)**2+($9)**2)*0.015)/2):(-($8)/sqrt(($8)**2+($9)**2)*0.015):(-($9)/sqrt(($8)**2+($9)**2)*0.015):4 w vectors head filled palette title ""
	}
	#set logscale x
	#plot [:][:] for[j=0:0] "-" index 0 using 5:8:10 with errorbars pt 7 ps 0.1
	if(i==1){
		plot [:][:] "-" u 2:3:4 w point pt 7 ps 0.5 palette title "",\
		#plot [0:0.4][-0.68:-0.56] for[j=0:0] "-" index 0 using 1:($10*3+5.65*$1):11 with errorbars pt 7 ps 0.3
		#plot [0:0.5][-2:0] for[j=0:0] "-" index 0 using 1:($10/(1-$1)):11 with errorbars pt 7 ps 0.3
		#plot [0:0.4][:] for[j=0:0] "-" index 0 using 1:9:10 with errorbars pt 7 ps 0.3
		#plot [0:0.3][0:0.25] for[j=0:0] "-" index 0 using 1:($17/2.) pt 7 ps 0.3
		#plot [:][-0.68:-0.56] for[j=0:0] "-" index 0 using 1:(3*($5+1./(3*4)*(2-4*$1))+5.65*$1):6 with errorbars pt 7 ps 0.2
		#plot [:][:] for[j=0:0] "-" index 0 using 3:4 pt 7 ps 1
	}
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
