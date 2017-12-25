reset
#custom
#multiplot
ix=16
iy=1
sx=2*ix
sy=3*iy
xlb=" "
ylb=" "
ml=0.1/sx
mb=0.3/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font "Times Bold,18" size sx,sy enhanced dash
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics 0.1 out nomirror scale 0.5 offset 0,0.5
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
do for[i=0:(ix*iy-1)]{
	oxi=i%ix
	oyi=i/ix+1
	set origin ml+oxi*gx,mb+(iy-oyi)*gy
	#if(oxi!=0){
		set format y ""
	#}
	if(oyi!=iy){
		set format x ""
	}
	if(i==0){
		#set key font ",16" at screen 0.05+ml+0.01/sx,1-mt-0.02/sy horizontal maxcols 1 samplen 0.5# opaque #autotitle
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	plot [:0.39][0:] for[j=0:16] "-" index i/2 every :::j::j using 5:(column(9+4*(i%2))) with line lw 8 title " "#, for[j=0:3] "-" index 0 every :::word("0 6 9",i+1)::word("0 6 9",i+1) using 1:(column(10+j)==1?column(2+j)-word("8e-5 1.5e-4 1.7e-4",i+1):1/0) with points pt 9 ps 0.7 lc 0#,"-" index 0 every :::i::i using 1:(column(4+j)>0?column(2+j):1/0):(sprintf("{/=10 %0.3f}",column(1))) with labels notitle
	unset key
	unset label
	unset format
}
#data
