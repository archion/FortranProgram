reset
#custom
#multiplot
ix=1
iy=1
sx=6*ix
sy=4*iy
xlb=" "
ylb=" "
ml=0.5/sx
mb=0.5/sy
mr=0.1/sx
mt=0.1/sy
gx=(1.-ml-mr)/ix
gy=(1.-mb-mt)/iy
#term
set term eps font ",18" size sx,sy
set output "-."."eps"
set label ylb font ",20" rotate by 90 center at character 0.8,screen 0.5*(1+mb-mt)
set label xlb font ",20" center at screen 0.5*(1+ml-mr),character 0.5
set xtics out nomirror
set ytics out nomirror
#
set size gx,gy
set size ratio -1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set multiplot
#set mxtics 5
#set mytics 5
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
	#set tics scale 0,0.001
	#set xtics 1
	#set ytics 1
	#set mxtics 2
	#set mytics 2
	#set grid front mxtics mytics lw 1.5 lt -1 lc rgb 'white'
	plot [:][:] "-" index 0  using ($1-4*$4*sin($6)/2.):($2-4*$4*cos($6)/2.):(4*$4*sin($6)):(4*$4*cos($6)) w vectors lw 7 title "", "-" index 0 u 1:2 w p pt 7 ps 0.5, "-" index 2 using ($1-$4/2.):($2-$5/2.):($4):($5):7 w vectors nohead lw 5 palette
}
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data
