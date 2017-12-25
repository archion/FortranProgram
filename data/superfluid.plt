reset
#custom
#multiplot
ix=1
iy=1
sx=5*ix
sy=3.2*iy
xlb=" "
ylb=" "
ml=0.6/sx
mb=0.4/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font "Times Bold,18" size sx,sy enhanced
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics out nomirror scale 0.8 offset 0,0.5
set ytics out nomirror scale 0.8 offset 0.5,0
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
		set key font ",14" at screen 1-mr-0.10/sx,1-mt-0.10/sy horizontal maxcols 1 samplen 0.8 opaque# autotitle
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	if(i==1){
	plot [0:][0:] for[j=0:18] "-" index 0 every :::j::j using 1:2  with lp pt 7 ps 0.3 title "".j
	}
	if(i==0){
	plot [0:][0:] for[j=0:5] "-" index 0 every :::0::0 using 1:(column(2+j))  with lp pt 7 ps 0.3 title "".j
	#plot [:][0:] for[j=4:2] "-" index 0 every :::j::j using ($1/word("0.8125E-03 0.3896E-02 0.8510E-02 0.1064E-01 0.1680E-01 0.1775E-01 0.2032E-01 0.2210E-01 0.2597E-01 0.2824E-01 0.3207E-01 0.2726E-01 ",j+1)):($2/word("0.3294E-01 0.4953E-01 0.6924E-01 0.8087E-01 0.1119E+00 0.1236E+00 0.1395E+00 0.1572E+00 0.1650E+00 0.1689E+00 0.1909E+00 0.2344E+00 ",j+1)) with l lw 4 title word("0.05 0.07 0.09 0.1 0.12 0.125 0.13 0.135 0.145 0.15 0.18 0.25",j+1)
	#plot [:][0:] for[j=0:2] "-" index 0 every :::j::j using ($1/word("0.1680E-01 0.2210E-01 0.3207E-01 ",j+1)):($2/word("0.1119E+00 0.1572E+00 0.1909E+00",j+1)) smooth bezier with l lw 5 lc j+6 title word("0.120 0.135 0.180",j+1) 
	#plot [:][:] for[j=0:11] "-" index 0 every :::j::j using ($1):($2-word("0.3294E-01 0.4953E-01 0.6924E-01 0.8087E-01 0.1119E+00 0.1236E+00 0.1395E+00 0.1572E+00 0.1650E+00 0.1689E+00 0.1909E+00 0.2344E+00 ",j+1)) with lp pt 7 ps 1 lw 4 title word("0.05 0.07 0.09 0.1 0.12 0.125 0.13 0.135 0.145 0.15 0.18 0.25",j+1)
	}
	if(i==2){
	plot [:][:] for[j=0:8] "-" index 0 every :::j::j using ($1):($2-word("0.1129E+00 0.1243E+00 0.1399E+00 0.1573E+00 0.1651E+00 0.1689E+00 0.1909E+00 0.2045E+00",j+1)) with lp pt 7 ps 0.5 title word("0.12 0.125 0.13 0.135 0.145 0.15 0.18 0.2",j+1)
	}
	if(i==3){
	plot [:][0:] for[j=0:8] "-" index 0 every :::j::j using 1:($2)  with lp pt 7 ps 0.3  title word("0.12 0.125 0.13 0.135 0.145 0.15 0.18 0.2",j+1)
	}
	unset label
	unset format
}
unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
