reset
#custom
#multiplot
ix=1
iy=2
sx=4*ix
sy=3*iy
xlb=" "
ylb=" "
ml=0.43/sx
mb=0.4/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font ",18" size sx,sy
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics 0.5 out nomirror
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
		set key font ",16" at graph 0.5,0.4 horizontal maxcols 2 opaque autotitle samplen 1 # columnhead
	}
	#set label sprintf("(%c)",97+i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	#plot [:1.5][0.2:1.3] for[j=0:14:1] "-" index j using 1:i+2:(100.) smooth acsplines w line lc j+1 lw 4 title sprintf("%1.3f",0.1+0.005*j)
	#plot [:1.1][0.3:1.3] for[j=0:9] "-" index j using 3:i+4:(1e3) smooth acsp w lines lw 4 title word("0.1 0.12 0.125 0.13 0.135 0.14 0.145 0.15 0.16 0.2",j+1)
	plot [:1.1][0.45:1.15] for[j=1:9] "-" index j using 2:i+3 smooth bezier w l lw 4 title word("0.1 0.12 0.125 0.13 0.135 0.14 0.145 0.15 0.16 0.2",j+1)
	unset label
	unset format
}
#data
