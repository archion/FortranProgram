reset
sx=2.5
sy=3
x1lb="{/Symbol w}"
#y1lb="Im[ln(-G_f^{-1})]"
y1lb="Im(G_B)"
flabel=""
#x2lb="X"
#y2lb="Y"
#zlb="Z"
lmg=5.5
rmg=1.0
umg=0.5
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
#set xrange [-2e-1:2e-1]
#set yrange [-10:10]
#set yrange [1e-3:2000]
set zrange [:]
set logscale x
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
#a=-1
#fit [x=0.00000001:0.00001] b*x**a  "-" index 70 u 1:(abs($2)/10000.) via a,b
array y[6]=[120,174,176,176,176,176]
array am[6]=[1,0.85,1,1,1,1]
array id[5]=[1,2,3,5,6]
#array a[3]=[0.-1.,0.03-1.,0.41-1.]
array a[3]=[0.25-1.,0.15-1.,0.41-1.]
array b[3]=[1.5,0.5,1.0]

plot for[i=117:120:1] "-" index 0 every :::i-1::i-1 u (am[i%2+3]*column("omega")):(-column("iGB")*am[i%2+1]) with l lw 2 dt 2-i%2 lc i%2*2 notitle word("T1 T2 LM C MCK C' LM'",i%4+1)#,x**(0.3+0.15)
#plot for[i=1:80:2] "-" index 0 every :::i-1::i-1 u (column("omega")):(-column("iGf")) with l lw 2 lc i notitle word("LM C MCK C' LM'",i)#,for[i=0:80:2] "-" index 0 every :::i-1::i-1 u (column("omega")/column("T")):(-0.048*column("T")**0.23*column("iGf")) with l lw 2 lc i title word("LM C MCK C' LM'",i)
#plot for[i=1:1] "-" index id[i]-1 every :::y[i]::y[i] u (abs(column("omega"))<2e11?column("omega")*am[i]:1/0):(-column("rSEf")) with l lw 2 lc i title word("LM C MCK C' LM'",i), for[i=1:1] "-" index id[i]-1 every :::y[i]+1::y[i]+1 u (abs(column("omega"))<2e11?column("omega")*am[i]:1/0):(-column("rSEf")) with l lw 2 dt 2 lc i+1 notitle
#plot "-" index 0 u 2:5 with l notitle
set key font ",".(ticfontsize-1) at graph 0.39,0.4,1 horizontal maxcols 1 spacing 1 samplen 1.5 autotitle #opaque
set label flabel font ",".fontsize front center textcolor rgb "black" at graph 0+c2sx/sx*2,1-c2sy/sy,1+2*c2sy/sy
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
	set xtics 1000 out nomirror scale 0.7, 0.3 offset 0,0.5 font ",".(ticfontsize)
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
	set ylabel y1lb offset character 4.5,0 
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

