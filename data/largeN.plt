reset
sx=4.3
sy=3.12
x1lb="{/Symbol w}"
y1lb="A_B"
#x2lb="X"
#y2lb="Y"
#zlb="Z"
lmg=5.5
rmg=8.5
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
#set size ratio -1

set pm3d corners2color c2
set palette rgbformulae 22,13,-31
#set logscale cb
#set cbrange [0.8:100]
set cbrange [:]
set cbtics scale 0.5 offset 0,character 2.1 font ",".(ticfontsize)
set colorbox horiz user origin graph 0,1+c2sy/sy*0.5 size graph 1,character 0.6
#set cbtics scale 0.5 offset character -0.7, 0 font ",".(ticfontsize)
#set colorbox user origin graph 1+c2sx/sx*0.5,0 size character 1, graph 1

set xrange [-.1:.1]
set yrange [:]
set zrange [:]
#set logscale x
#set logscale y
#a=0.229-1.
#a=0.106-1.
#a=0.03992865381623536-1.
#a=0.414862855232646-1.
#a=0.394862855232646-1.
#a=0.1064525448624235-1.
#a=0.22940981529676907-1.
#a=0.3286489218401046-1.
#a=0.06309878156026809-1.
#a=0.543604996816506-1.
#a=0.014569749254921006-1.
#a=0.15-1.
a=0.45
b=0.000001/0.11
#fit [x=0.00000001:0.00001] b*x**a  "-" index 70 u 1:(abs($2)/10000.) via a,b
plot for[i=0:0] for[j=1:2] "-" index i u ((column("Tk")-1e-4)<1e-11?column("omega"):column("omega")):(-column(word("iGphi0 rhophi",j))) with l notitle#, 10000.*b*x**a w line dt 2 title sprintf("f(x)= %1.2fx^{%1.3f}",b*10000.,a)#,10000.*b*x**(-0.894) w line dt 2 title sprintf("f(x)= %1.2fx^{%1.3f}",b*10000.,-0.894)
#plot "-" index 0 u 2:5 with l notitle
set key font ",".(ticfontsize-2) at graph 0.84,0.95,1 horizontal maxcols 1 spacing 0.9 samplen 1.5 autotitle #opaque
set label "(b)" font ",".fontsize front center textcolor rgb "black" at graph 0+c2sx/sx*2,1-c2sy/sy,1+2*c2sy/sy
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
replot
if(GPVAL_TERM eq "qt"){
	pause -1
}
#data

