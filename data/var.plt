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
set xtics out nomirror scale 0.5 offset 0,0.5
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
	if(oxi!=0){
		set format y ""
	}
	if(oyi!=iy){
		#set format x ""
	}
	if(i==(ix*iy-1)){
		#set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 samplen 0.5 opaque# autotitle
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	#splot [:][:] for[j=18:18] "-" index 0 every :::j::j using 0:2:14 with lp title ""
	#i=0
	#if(i%2!=0){
	#plot [:][-0.45564:-0.45] for[j=38:38] "-" index 0 every :::i/2::i/2 using 0:6:7 with errorbars lw 1 title "","-" index 0 every :::i/2::i/2 using 0:6 w line
	#plot [:][:] for[j=38:38] "-" index 0 every :::i::i using 0:6:7 with errorbars lw 1 title "","-" index 0 every :::i::i using 0:6 w line
	#plot [:][:] for[j=38:38] "-" index 0 every :::i::i using 0:4:5 with errorbars lw 1 title "","-" index 0 every :::i::i using 0:4 w line
	#plot [:][:] for[j=38:38] "-" index 0 every :::i::i using 0:10:11 with errorbars lw 1 title "","-" index 0 every :::i::i using 0:10 w line
	#plot [:][:] for[j=38:38] "-" index 0 every :::i::i using 0:10:11 with errorbars lw 1 title "","-" index 0 every :::i::i using 0:10 w line
	#plot [:][:] for[j=38:38] "-" index 0 every :::i::i using 0:(abs($10)) w line
	#plot [-2:][:] for[j=10:10] "-" index 0 every :::i::i using 0:j:j+1 with errorbars lw 1 title "","-" index 0 every :::i::i using 0:j w line
	plot [-2:][:] for[j=1:1] "-" index 0 every :::i::i using 0:j:(column(j+1)) w err lw 1 ps 0.5 pt 7 title ""
	#}else{
	#plot [:][:] for[j=38:38] "-" index 0 every :::i/2::i/2 using 0:5 with lp
	#}
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
