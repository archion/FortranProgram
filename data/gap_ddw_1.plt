reset
#custom
#multiplot
ix=1
iy=1
sx=4*ix
sy=3*iy
xlb=" "
ylb=" "
ml=0.8/sx
mb=0.6/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font ",18" size sx,sy enhance
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics 0.5 out nomirror
set ytics 0.05 out nomirror
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
		set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque autotitle #columnhead
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	plot [:][0:] for[j=0:3] "-" index 0 using ($0/10.):(-column(2+j*2)):(column(3+j*2)<1e-8?1e-1:column(3+j*2)) with p pt 7 ps var lc j+1 notitle, for[j=0:1] "-" index 0 using ($0/10.):(column(14+j)) with line lw 5 lc j+4 notitle
	unset label
	unset format
}
#data
