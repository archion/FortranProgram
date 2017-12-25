reset
#custom
#multiplot
ix=3
iy=1
sx=1.7*ix
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
set term eps font "Times Bold,18" size sx,sy enhanced
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics 0.05 out nomirror scale 0.5 offset 0,0.5
set ytics 0.0012 out nomirror scale 0.5 offset 0.5,0
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
	if(i==(ix*iy-1)){
		#set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque autotitle
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	#plot [word("0.13 0.05 0",i+1):word("0.32 0.28 0.24",i+1)][word("0.045 0.02 0",i+1):] for[j=0:11] "-" index i every :::j::j using 1:($5+(11-j)*0.02) with line lw 5 lc j  title "",  for[j=0:11] "-" index i every :::j::j using 1:($13!=0?$5+(11-j)*0.02-word("8e-3 8e-3 8e-3",i+1):1/0) with points pt 9 ps 0.7 lc 0#,"-" index 0 every :::i::i using 1:(column(4+j)>0?column(2+j):1/0):(sprintf("{/=10 %0.3f}",column(1))) with labels notitle
	plot [0:word("0.09 0.18 0.18",i+1)][0:] for[j=0:11] "-" index i every :::j::j using 1:($9+(11-j)*0.0002) with line lw 5 lc j  title "",  for[j=0:11] "-" index i every :::j::j using 1:($17!=0?$9+(11-j)*0.0002-word("8e-5 1.5e-4 1.7e-4",i+1):1/0) with points pt 9 ps 0.7 lc 0#,"-" index 0 every :::i::i using 1:(column(4+j)>0?column(2+j):1/0):(sprintf("{/=10 %0.3f}",column(1))) with labels notitle
	unset label
	unset format
}
#data
