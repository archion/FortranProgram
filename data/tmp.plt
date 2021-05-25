reset
sx=(exist("sx")?sx:4)
sy=(exist("sy")?sy:4)
flabel="(a)"
lmg=7.
rmg=1.0
umg=2.5
dmg=2.5
fontsize=15
ticfontsize=floor(fontsize*0.7)
#term
c2sx=0.1235/12*fontsize
c2sy=0.2206/12*fontsize
set term eps rounded fontscale 1 font "Helvetica,".fontsize lw 5 size sx+(lmg+rmg)*c2sx,sy+(dmg+umg)*c2sy 
set output "-."."eps"
unset key
unset tics
set ticslevel 0
#set grid xtics ytics
set border 8+4+2+1+16+32+64+128+256+512+1024+2048 front
set margin lmg,rmg,dmg,umg
set size ratio -1

set pm3d corners2color c2
set palette rgbformulae 22,13,-31
#set logscale cb
#set cbrange [0.64:0.68]
set cbtics 0.05 scale 0.5 offset 0,character 2.1 font ",".(ticfontsize)
set colorbox horiz user origin graph 0,1+c2sy/sy*0.5 size graph 1,character 0.6
#set cbtics scale 0.5 offset character -0.7, 0 font ",".(ticfontsize)
#set colorbox user origin graph 1+c2sx/sx*0.5,0 size character 1, graph 1

set xrange [:]
#set xrange [-1:5]
#set yrange [1:1e-5]
set yrange [:]
#set yrange [-1:5]
set zrange [:]
#set logscale x
#set logscale y
set cbrange [-1:1]
a=-0.771
b=0.0001/5
i=10
plot "-" index 0 every :::i::i u 1:2:($3==0?1/0:($3==2?0:$3)) with point pt 7 palette
set key font ",".ticfontsize at graph 1.00,0.95,1 horizontal maxcols 1 spacing 0.8 samplen 0.5 #opaque autotitle
#if(exist(flabel)) set label flabel font ",".fontsize front center textcolor rgb "black" at graph 0+c2sx/sx*2,c2sy/sy,1+2*c2sy/sy
if(exists("zlb")){
	set key at screen 1,1,1
}

set mxtics 5
set mytics 5
set mx2tics 5
set my2tics 5
set mztics 5
set format "%h"
if(exists("zlb")){
	set ztics out scale 0.7, 0.3 offset 0.7,0 font ",".(ticfontsize)
	set zlabel zlb offset character 3,0 
}
if(exists("x1lb")) {
	set xlabel x1lb offset character 0,1.2 
	#set label x1lb center at graph 0.5,screen 0 front offset character 0, 0.7
	set xtics out nomirror scale 0.7, 0.3 offset 0,0.5 font ",".(ticfontsize)
	if(exists("zlb")) {
		set xlabel offset character 0,0
		set xtics offset 0,0 
	}
} 
if(exists("x2lb")) {
	set x2label x2lb offset character 0,-1.2 
	#set label x2lb center at graph 0.5,screen 1 front offset character 0, -0.7
	set x2range [GPVAL_X_MIN:GPVAL_X_MAX]
	set x2tics out nomirror scale 0.7, 0.3 offset 0,-0.5 font ",".(ticfontsize)
	if(exists("zlb")) {
		set x2label offset character 0,0
		set x2tics offset 0,0
	}
}
if(exists("y1lb")) {
	set ylabel y1lb offset character 3,0 
	#set label y1lb center at screen 0,graph 0.5 front offset character 1,0 rotate by 90
	set ytics out mirror scale 0.7, 0.3 offset 0.7,0 font ",".(ticfontsize)
	if(exists("zlb")) {
		set ylabel offset character 0,0
		set ytics offset 0,0
	}
}
if(exists("y2lb")) {
	set y2label y2lb offset character -3,0 rotate by -90
	#set label y2lb center at screen 1,graph 0.5 front offset character -1,0 rotate by -90
	set y2range [GPVAL_Y_MIN:GPVAL_Y_MAX]
	set y2tics out mirror scale 0.7, 0.3 offset -0.7,0 font ",".(ticfontsize)
	if(exists("zlb")) {
		set y2label offset character 0,0
		set y2tics offset 0,0
	}
}


set output "-."."eps"
set multiplot
replot
set colorbox horiz user origin graph -100,-100
set size 0.6,0.6
set origin .3, .3
set size ratio -1
unset xrange
unset yrange
unset xtic
unset ytic
unset border
unset xlabel
unset ylabel
plot for[i=89:89] "-" index 0 every ::(i-1)*id[1]::(i*id[1]-1) u "x":"y":"nc" with p pt 7 palette notitle "N_{c}"#, for[i=1:1] "-" index 0 every ::(i-1)*id[1]::(i*id[1]-1) u "x":"y":(($0)+1) with labels font ",".ticfontsize/2 notitle
#data
