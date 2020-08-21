reset
sx=3
sy=3
#x1lb="{/Symbol k}"
#x1lb=""
#x1lb="T/T@_K^0"
x1lb="T"
#y1lb="Im[ln(-G_f^{-1})]"
#y1lb="0.6951-s_{imp}"
y1lb="s_{imp}"
#x2lb="X"
#y2lb="entropy"
#zlb="Z"
flabel=" "
lmg=4.5
rmg=1.0
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
set cbrange [:]
set cbtics scale 0.1 offset 0,character 2.1 font ",".(ticfontsize)
set colorbox horiz user origin graph 0,1+c2sy/sy*0.5 size graph 1,character 0.6
#set cbtics scale 0.5 offset character -0.7, 0 font ",".(ticfontsize)
#set colorbox user origin graph 1+c2sx/sx*0.5,0 size character 1, graph 1

#set xrange [:]
set xrange [:]
set yrange [:]
#set yrange [1:1e-5]
#set yrange [:0.156]
set zrange [:]
set logscale x
set logscale y
#a=-0.771
#b=0.0001/5
#set style line 1 linetype dl

set style line 2 dashtype 3
set style line 3 dashtype 2
#set label "MCK" at 2e-8,0.585 font ",".ticfontsize
#set label "C'" at 2e-8,0.705 font ",".ticfontsize
#set label "C" at 3e-8,0.688 font ",".ticfontsize
#set label "LM" at 4e-5,0.685 font ",".ticfontsize
#set label "LM'" at 2e-8,0.732 font ",".ticfontsize
#set arrow from 5.5e-8,0.687 to 1e-6,0.675
#set arrow from 3e-5,0.685 to 3e-6,0.675
array fact[3]=[5,5,1]
array id[5]=[1,2,3,4,6]
n=5
array a[n]
array b[n]
array c[n]
a_=0.066
b_=0.1
c_=0.70345
do for[i=n:n] {
	#fit [x=1e-9:1e-5] b_*x**a_+c_  "-" index 0 every :::i-1::i-1 u "T":"Sn"  via a_,b_,c_
		a[i]=a_
		b[i]=b_
		c[i]=c_
}
unset logscale y
#unset logscale x
array fact[5]=[1.,0.65,0.6,0.6,0.5]
array Tk[5]=[0.146906338116241,0.09094957347610863,0.04090102787182847,0.008271824915486082,0.00006842308743245633]
array Jk[5]=[1,0.8,0.6,0.4,0.2]
array Jk[11]=[2.0,1.8,1.6,1.4,1.2,1,0.8,0.6,0.4,0.2]
fcut=1.
Anorm=0.521387293498214
#Anorm=0.465670688757744
#Anorm=0.243747036329501
plot for[j=1:1] for[i=1:11:1] "-" index j-1 every :::i-1::i-1 u (column("T")/(1.*exp(-0./(Anorm*Jk[i])))):(fact[1]*(column(word("Sn Sa Snp",1))-0*0.7441)) w l lw 2 lc (i+1)/2 dt 1 title word("J=2.0 1.8 1.6 1.4 1.2 1 0.8 0.6 0.4 0.2",i)#, log(2) w l lc 0 dt 2
#plot for[i=1:1] for[j=1:1] "-" index 0  every :::j-1::j-1 u "T":(column(word("Sn Sa Snp",i))-0.62455) with l lw 2 lc i title word("Sn Sa C' LM'",i),"-" u "r":(column("Sn")-column("Sa")) with l lw 2 lc 0 dt 2 title "Sn-Sa"
#plot for[i=n:n] b[i]*x**a[i] w l lc i dt 2 title sprintf("#x^{%1.4f}+%1.4f",a[i],c[i]), for[i=n:n] "-" index 0 every :::i-1::i-1 u "T":((column(word("Sn Sa Snp",1))-c[i])) with l lw 2 lc i title word("k=0.5 1 1.5",i)

#unset logscale y
#unset logscale x
#array e[8]=[1.,1.,1.,1.,1.2,1.5,1.8,2.0]
#plot a u (($0)+1):($2) w lp pt 7 title "exponents", c u (($0)):($2) w lp pt 7 axis x1y2 title "entropy", e u ($0+1):(2./(1.+$2)) w lp pt 7
#plot for[i=1:*] for[j=1:1] "-" index 0 every :::i-1::i-1 u "T":(column(word("Sn Sa Snp",j))) with l lw j lc i dt j title word("LM C MCK MCK C' LM'",i)#, for[i=1:6] for[j=2:2] "-" index 0 every :::i-1::i-1 u "T":(column(word("Q gW Snp",j))) with l lw j lc i dt j notitle word("LM C MCK MCK C' LM'",i),
#plot for[i=2:2] for[j=1:8] "-" index 0 every :::i-1::i-1 u "T":(column(word("Q gW GB g f B c phi",j))) with l lc j lw 2 dt j/5+1 title word("Q gW GB g f B c phi",j)
#plot "-" index 0 u ((abs(column("T")-1e-5)<1e-11)?column("J"):1/0):(column(word("Sn Sa",1))) with p pt 7 ps .5 notitle, "-" index 0 u ((abs(column("T")-1e-7)<1e-11)?column("J"):1/0):(column(word("Sn Sa",1))) with p pt 7 ps .5 notitle, "-" index 1 u ((abs(column("T")-1e-7)<1e-11)?column("J"):1/0):(column(word("Sn Sa",1))) with p pt 7 ps .5 notitle, "-" index 1 u ((abs(column("T")-5e-10)<1e-11)?column("J"):1/0):(column(word("Sn Sa",1))) with p pt 7 ps .5 notitle
set key font ",".ticfontsize at graph 0.5,0.45 vertical maxrows 5 spacing 0.8 samplen 1.0 #opaque autotitle
set label flabel font ",".fontsize front center textcolor rgb "black" at graph 0+c2sx/sx*2,c2sy/sy,1+2*c2sy/sy
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
	set xtics 100 out nomirror scale 0.7, 0.3 offset 0,0.5 font ",".(ticfontsize)
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
	set ylabel y1lb offset character 2.5,0 
	#set label y1lb center at screen 0,graph 0.5 front offset character 1,0 rotate by 90
	set ytics out nomirror scale 0.7, 0.3 offset 0.7,0 font ",".(ticfontsize)
	if(exists("zlb")) {
		set ylabel offset character 0,0
		set ytics offset 0,0
	}
}
if(exists("y2lb")) {
	set y2label y2lb offset character -3,0 rotate by -90
	#set label y2lb center at screen 1,graph 0.5 front offset character -1,0 rotate by -90
	set y2range [GPVAL_Y2_MIN:GPVAL_Y2_MAX]
	set y2tics out nomirror scale 0.7, 0.3 offset -0.7,0 font ",".(ticfontsize)
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
