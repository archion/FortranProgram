reset
#custom
#multiplot
ix=3
iy=1
sx=2*ix
sy=2*iy
xlb=" "
ylb=" "
ml=0.19/sx
mb=0.16/sy
mr=0.05/sx
mt=0.05/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font "Times Bold,10" size sx,sy enhanced
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics ("X" 256) out nomirror scale 0.5 offset 0,0.5
set ytics 0.1 out nomirror scale 0.5  offset 0.5,0
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
	set xzeroaxis
	plot [:][-0.2:0.2] for[j=0:3] "-" index 0 every  :::word("0 6 9",i+1)::word("0 6 9",i+1) using 0:(column(1+j*2)):(column(2+j*2)) with lp pt word("7 6 6 7",j+1) pi 3 ps var lc word("3 1 3 1",j+1) , for[j=0:3] for[k=0:3] "-" index 0 every  :::word("0 6 9",i+1)::word("0 6 9",i+1) using ($0):(column(9+j)!=0?column(1+k*2):1/0):(sprintf("{/=8 %0.4f}",column(1+k*2))) with labels textcolor lt k+1 rotate by 90 notitle
	unset label
	unset format
}
#data
