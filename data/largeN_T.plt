reset
sx=4.3
sy=3.12
x1lb="T"
y1lb="Y"
#x2lb="X"
#y2lb="Y"
#zlb="Z"
lmg=8.5
rmg=5.5
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

set xrange [:]
set yrange [:]
set zrange [:]
set logscale x
set logscale y
a_=-0.9
b_=1.
n=4
array a[n]
array b[n]
#do for[i=1:n] {
	#fit [x=1e-8:1e-6] b_*x**a_  "-" index 0 every :::i-1::i-1 u "T":(abs(column("rGB0"))) via a_,b_
	#a[i]=a_
	#b[i]=b_
#}
#plot for[i=1:n] "-" every :::i-1::i-1 u (column("T")):(abs(column("rGB0"))) with l lc i notitle "".i, for[i=1:n] b[i]*x**a[i] w l lc i dt 2 title sprintf("{/Symbol a}=%1.3f",1+a[i])
unset logscale y
set xrange [:]
set yrange [:]
plot for[i=1:9] for[j=1:1] "-" every :::i-1::i-1 u (column("T")):(abs(column(word("S1 S2",j)))) with l lc i dt j notitle "".j#, log(2.) w l dt 2 lc 0
set key font ",".fontsize at graph 1,0.98,1 horizontal maxcols 1 spacing 1.0 samplen 1.5 #opaque autotitle
set label "(a)" font ",".fontsize front center textcolor rgb "black" at graph 0+c2sx/sx*2,1-c2sy/sy,1+2*c2sy/sy
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

