reset
sx=4.
sy=3.12
x1lb="[{/Symbol p}/sin({/Symbol p}{/Symbol t}T)]^{{/Symbol a}_{f}}"
#x1lb="{/Symbol t}"
#y1lb="{/Symbol c}({/Symbol t})"
#y1lb="T^{-{/Symbol a}_{f}}G_{f}({/Symbol t})"
y1lb="T^{-{/Symbol a}_{f}}G_{f}({/Symbol t})"
flabel=""
#x2lb="X"
#y2lb="Y"
#zlb="Z"
lmg=8.
rmg=0.5
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

set xrange [:10.]
#set xrange [:0.15]
set yrange [:]
set zrange [:]
unset logscale x
unset logscale y
#a=0.769
a=0.41
b=0.0001/1.1*tan(pi*(1.+a)/2.)
#fit [x=0.00000001:0.00001] b*x**a  "-" index 70 u 1:(abs($2)/10000.) via a,b
array fact[4]=[1.,0.006,1.,1.]
#array a[4]=[0.55,1.45,0.29,0.67]
#array a[4]=[0.77,0.15,0.41,0.03]
#array b[4]=[0.38,0.13,0.53,0.15]
n=3
array a[n]
array b[n]
array c[n]
a_=0.5
b_=1.
c_=0.45
j=8
do for[i=n:n] {
	#fit [x=1e-8:1e-6] b_*x**a_  "-" index i-1 every :::j-1::j-1 u (pi*column("T")/sin(pi*column("tau")*column("T"))):(-column("Gft")) via a_,b_
	#fit [x=1.8:5.5] b_*x+a_  "-" index i-1 every :::44-1::44-1 u ((pi/sin(pi*column("tau")*column("T")))**(0.4)):(-column("T")**(-0.4)*column("Gft")) via a_,b_
		a[i]=a_
		b[i]=b_
}
#fit [x=1.5:10.] b_*x+a_  "-" index 0 every :::34-1::34-1 u ((pi/sin(pi*column("tau")*column("T")))**(0.45)):(-column("T")**(-0.45)*column("Gft")) via a_,b_
a[1]=a_
b[1]=b_
#fit [x=25:35] b_*x+a_  "-" index 0 every :::29-1::29-1 u ((pi/sin(pi*column("tau")*column("T")))**(0.45)):(-column("T")**(-0.45)*column("Gft")) via b_,a_
a[2]=a_
b[2]=b_
b_=-1.
#fit [x=1.5:30] (b[2])*x+a_*x**b_  "-" index 0 every :::29-1::29-1 u ((pi/sin(pi*column("tau")*column("T")))**(0.45)):(-column("T")**(-0.45)*column("Gft")) via a_,b_
a[3]=a_
b[3]=b_

set logscale x
set logscale y
#plot for[i=1:4] for[j=9:1:-1] "-" index i-1 every :::j-1::j-1 u (pi*column("T")/sin(pi*column("tau")*column("T"))):(-column("Gft")) with l lc j lw 3 notitle# ,for[i=n:n]  1./(x**(-a[i])/b[i]+2.) w l dt 2 lc 0 title sprintf("{/Symbol a}=%1.4f",a[i])
#plot for[i=1:1] for[j=44:20:-5] "-" index i-1 every :::j-1::j-1 u ((pi/sin(pi*column("tau")*column("T")))**(0.45)):(-column("T")**(-0.45)*column("Gft")) with l lc j lw 3 title "T_{".j."}", for[i=2:1:-1]  b[i]*x+a[i] w l lw 2 dt 2 lc 2+i title sprintf("%1.2fx%+1.2f",b[i],a[i]), for[i=3:3]  (b[2])*x+a[i]*x**b[i] w l lw 2 dt 3 lc 7 title sprintf("%1.2fx%+1.2fx^{%1.2f}",b[2],a[i],b[i])
#plot for[i=n:n] for[j=44:20:-5] "-" index i-1 every :::j-1::j-1 u ((pi*column("T")/sin(pi*column("tau")*column("T")))**(0.338)):(-column("Gft")) with l lc j lw 3 notitle#, x**(1) w l dt 2# ,for[i=n:n]  1./(x**(-a[i])/b[i]+2.) w l dt 2 lc 0 title sprintf("{/Symbol a}=%1.4f",a[i])
unset logscale x
unset logscale y
set xrange [1.2:2.0]
#plot for[i=1:1] for[j=9:1:-1] "-" index i-1 every :::j-1::j-1 u (column("tau")*column("T")):(-column("Gft")*column("T")**(-1./1.5+0.00)) with l lc j lw 3 notitle#, x**(1./1.5) w l dt 2# ,for[i=n:n]  1./(x**(-a[i])/b[i]+2.) w l dt 2 lc 0 title sprintf("{/Symbol a}=%1.4f",a[i])
array fact[2]=[1.,0.691]
plot for[i=1:2] for[j=9:9:-1] "-" index i-1 every :::j-1::j-1 u ((pi/sin(pi*column("tau")*column("T")))**0.2):(-fact[i]*column("Gft")*(column("T")**(-0.2))) with lp lc i lw 1 notitle, 0.385*x w l dt 2 lc 0# ,for[i=n:n]  1./(x**(-a[i])/b[i]+2.) w l dt 2 lc 0 title sprintf("{/Symbol a}=%1.4f",a[i])
#plot for[i=n:n] for[j=5:9] "-" index i-1 every :::j-1::j-1 u (0.5-column("tau")*column("T")):(-column("Gft")/(column("T")**(0.2/2.))) with l lc j lw 1 notitle# ,for[i=n:n]  1./(x**(-a[i])/b[i]+2.) w l dt 2 lc 0 title sprintf("{/Symbol a}=%1.4f",a[i])
set key font ",".ticfontsize at graph 0.95,0.55,1 horizontal maxcols 1 spacing 0.9 samplen 1.5 #opaque autotitle
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
	set ylabel y1lb offset character 1,0 
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

