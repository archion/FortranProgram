reset
sx=4.3
sy=3.12
x1lb="{/Symbol w}"
y1lb="Intensity"
#x2lb="X"
#y2lb="Y"
#zlb="Z"
lmg=5.5
rmg=1.5
umg=1.5
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
#set logscale x
#set logscale y
#a=-0.771
#a=-0.785
a=0.25-1
#a=-1./1.5
b=0.0001/0.11
#fit [x=0.00000001:0.00001] b*x**a  "-" index 70 u 1:(abs($2)/10000.) via a,b
#plot for[i=72*2-2:72*2-1] for[j=1:1] "-" index i u "omega":((column(word("Af",j)))*.046**(i%2)) with l notitle
set style line 1 dashtype 3
set style line 2 dashtype 2
#plot for[i=1:1] for[j=1:2] "-" index i every :::(j-1)::(j-1) u (column("omega")):((column(word("iGB rGB iGf rGf",j)))) with l ls i dt j title word("iGf rGf",j)
#plot for [k=71:71] for[i=k*2-2:k*2-1] for[j=1:1] "-" index i u (column("omega")):(1.05**(i%2)*column(word("iGf iGB",j))*column("iGB")) with l title word("iGf rGf",j)
#plot for [k=1:4:2] for[i=k*1-1:k*1-1] for[j=2:2] "-" index i u (column("omega")):((-1)**(k/2)*(column(word("iGf iGf",j)))) with l title word("iGf rGf",j)
array cl[1]=["iGf"]
#array fact[2]=[0.,0.27]
array fact[2]=[0,0]
#set logscale x
#set xrange [-0.1:0.1]
#set yrange [-10:10]
#set logscale y
#set logscale x
omg=-4.72870288e-01
set yrange [:]
set xrange [:]
plot for[i=1:1] for[j=1:2] "-" index i every :::(j-1)::(j-1) u (column("omega")):((column(word("iGb rGB iGf rGf",4))-0.0*(j-1))) with l ls j dt j title word("iGf rGf",j)
#set xrange [-omg-0.01:-omg+0.01]
#set xrange [0.47:0.475]
#plot for[i=1:4] for[j=1:5] "-" index i-1 every :::(j-1)::(j-1) u ((column("x"))):((column("".word("Gb Gf1 Gf2 Ga rt2",5)))) w l lc j
array dis=[-10.,-4.5191247050730e-08,-10.-4.5191247050730e-08]
#set xrange [10-0.000001:10]
#set xrange [-10:10]
#set xrange [0.5:1]
#set yrange [-1e1:1e1]
#set logscale y
#plot "-" index 0 u ($1):"iGf1" w p pt 7 ps 0.5 lc 0, for[i=1:|dis|] "-" index 0 u ($1-dis[i]):(0.1*i) w p pt 7 ps 0.5 lc 0, for[i=1:1] for[j=2:2] "-" index i every :::(j-1)::(j-1) u ((column(3))):(0.4) w p pt 7 ps 0.5 lc j,
#plot for[j=7:15:1] for[i=1:1] "-" index 1 every :::(j-1)::(j-1) u "y":(column("rt")) w l pt 7 ps 0.5 lc i
#plot for[i=1:1] for[j=1:1:1] for[k=1:2] "-" index i every :::(j-1)::(j-1) u (column(word("x z",k))):(-column(word("rt2",1))) w lp pt 7 ps 0.5 lc k
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

