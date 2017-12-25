reset
#custom
#multiplot
ix=1
iy=1
sx=5*ix
sy=3*iy
xlb="T"
ylb="gap size"
ml=0.8/sx
mb=0.6/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font ",18" size sx,sy
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics 0.05 out nomirror
set ytics  out nomirror
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
		#set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque autotitle columnhead
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	set style fill transparent solid 0.3 border
	#set style rect fc lt -1 fs solid 0.15 noborder
	#set obj rect from 0.14, graph 0 to 0.175, graph 1 behind
	#plot [0.05:0.25][0:0.05] "-" using (1-$1):($2/2.):(1e7) smooth acsplines lc 1, "-" u (1-$1):($3/2.):(1e8) smooth acsplines lc 3 , "-" u (1-$1):($4/2.):(1e7) smooth acsplines lc 3 
	plot [0.0:][0:0.1] "-" using (1-$1):2 lc 1, "-" u (1-$1):3 lc 3 , "-" u (1-$1):4 lc 3 
    unset label
	unset format
}
#data
