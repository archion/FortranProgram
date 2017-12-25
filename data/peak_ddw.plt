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
mb=0.3/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font "Times Bold,18" size sx,sy enhanced
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics 0.5 out nomirror offset 0,0.5
set ytics 0.05 out nomirror offset 0.5,0
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
		#set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque autotitle
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	plot [:][:] for[k=0:2] "-" index 0 using ((column(-1)+1)/15.):(column(7+k)!=0&column(2+k)>5e-2?$1:1/0):(column(2+k)) with p pt 7 ps var lc k+1 title "", "-" index 1 using ((column(0)+1)/15.):1 with line lw 5 lc 0 notitle
	unset label
	unset format
}
#data
