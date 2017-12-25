reset
#custom
#multiplot
ix=1
iy=1
sx=5*ix
sy=4*iy
xlb=""
ylb=""
ml=0.8/sx
mb=0.6/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font ",18" size sx,sy
set output "rpa."."eps"
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
		set key font ",16" at screen 1-mr-0.01/sx,1-mt-0.01/sy horizontal maxcols 1 opaque autotitle columnhead
	}
	#set label sprintf("(%1.3f)",0.11+0.005*i) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	set label word("对角 水平 0.51 0.7",i+1) at graph 0.1/(gx*sx),1.-0.15/(gy*sy)
	set pm3d map corners2color c1 interpolate 1,1
	set size square
	set palette gray
	set hidden3d
	set cbrange [:]
	splot [:][:] "rpa.dat" index i matrix w image
	#plot [:] "rpa.dat" u 1:($2/256.) w line, "rpa.dat" u 1:(-$3) w line
	unset label
	unset format
}
	unset multiplot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
