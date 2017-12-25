reset
#custom
#multiplot
ix=1
iy=1
sx=4*ix
sy=3*iy
xlb=" "
ylb=" "
ml=0.5/sx
mb=0.5/sy
mr=0.2/sx
mt=0.2/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font "Times Bold,18" size sx,sy enhanced
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics 10 out nomirror scale 0.5 offset 0,0.5
set ytics 10 out nomirror scale 0.5 offset 0.5,0
set border front lw 5
unset border
unset xtics
unset ytics
unset ztics
#set grid xtics ytics
#
set size gx,gy
set size ratio -1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set multiplot
set mxtics 10
set mytics 10
unset xlabel
unset ylabel
#set grid ytics mytics
unset key
do for[i=0:(ix*iy-1)]{
	oxi=i%ix
	oyi=i/ix+1
	set origin ml+oxi*gx,mb+(iy-oyi)*gy
	set view 0,0
	set view equal xyz
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
	splot [:][:]\
	 for[j=0:0] "-" index 1 using ($1):($2+j*sqrt(3)*4):3 w line title "",\
	 for[j=0:0] "-" index 0 using ($1):($2+j*sqrt(3)*4):3 w point pt 7 ps 1 lc "red" title "",\
	 for[j=0:0] "-" index 0 using ($1):($2+j*sqrt(3)*4):3:4 w labels font ",5" tc "white" title "",\
	 for[j=0:0] "-" index 2 using ($1):($2+j*sqrt(3)*4):3 w point pt 5 ps 1 lc "blue" title "",\
	 for[j=0:0] "-" index 2 using ($1):($2+j*sqrt(3)*4):3:4 w labels font ",5" tc "white" title ""
	#plot [:][:] for[j=0:0] "-" index i every :::j::j using 1:2 w point pt 6 ps 1 lc 1 title "", "-" index 1 every :::j::j using 1:2 w line title ""
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
