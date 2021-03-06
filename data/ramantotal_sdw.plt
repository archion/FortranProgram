reset
#custom
#multiplot
ix=1
iy=2
sx=4*ix
sy=2*iy
xlb="T/Tc"
ylb="Normalized gap"
ml=0.8/sx
mb=0.6/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font ",18" size sx,sy
set output "ramantotal_ddw."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics out nomirror
set ytics out nomirror
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
		set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque
	}
	set label sprintf("(%c)",97+i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	plot [:1.5][0.2:1.3] for[j=0:14:1] "ramantotal_ddw.dat" index j using 1:i+2:(100.) lc j smooth acsplines notitle, for[j=0:14:1] '' index j every 5 using 1:i+2 with point lc j pt j+1 ps 0.5 title sprintf("%1.3f",0.1+0.005*j)
	unset label
	unset format
}
#data
