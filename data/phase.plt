reset
#custom
#multiplot
ix=1
iy=1
sx=3.9*ix
sy=3*iy
xlb="x"
ylb="T"
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
set xtics 0.05 out nomirror scale 0.5 offset 0,0.5
set ytics out nomirror scale 0.5 offset 0.5,0
set border front lw 5
#set grid xtics ytics
#
set size gx,gy
#set size ratio -1
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
sp='< awk -v RS="\n\n\n" -v cmd="sort -g -k1 -k2" ''{print | cmd; close(cmd); print ""}'' "-" | awk ''NF < 3 {print;next};$1 != prev {printf "\n"; prev=$1};{print}'''
do for[i=0:(ix*iy-1)]{
	oxi=i%ix
	oyi=i/ix+1
	set origin ml+oxi*gx,mb+(iy-oyi)*gy
	set pm3d map
	set pm3d corners2color c2
	#set logscale cb
	set cbrange [:]
	#set size square
	set palette rgbformulae 22,13,-31
	if(oxi!=0){
		set format y ""
	}
	if(oyi!=iy){
		set format x ""
	}
	if(i==(ix*iy-1)){
		set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 samplen 0.5 opaque# autotitle
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	plot [:0.3][0:0.08] for [i=0:2] for[j=0:1] "-" index i using (column(1+j*2)):((column(2+j*2))>2e-4?(column(2+j*2)):"nan") w lp pt 7 ps 0.5 lw 5 title word("SC rSC DDW rDDW SDW rSDW",i*2+j+1)
	#splot [:][:] sp index 0 using 1:2:(-$5)
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
