reset
sx=3.
sy=3
#x1lb="{/Symbol w}/T"
x1lb="T"
#y1lb="Im[ln(-G_f^{-1})]"
#y1lb="{/Symbol c}"
y1lb=" "
#y1lb="T^{0.2}Im(G_B)"
flabel=""
#x2lb="X"
#y2lb="Y"
#zlb="Z"
lmg=7.5
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
set yrange [:]
#set yrange [1e-5:]
set zrange [:]
set logscale x
set logscale y
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
a=0.15-1.
#a=-1
#fit [x=0.00000001:0.00001] b*x**a  "-" index 11 u "omega":"iGf" via a,b
array y[6]=[120,174,176,176,176,176]
array am[10]=[8,0.338,1,1,1,1,1,1,1,1]
array id[5]=[1,2,3,5,6]
#array a[3]=[0.-1.,0.03-1.,0.41-1.]
#array a[3]=[0.25-1.,0.15-1.,0.41-1.]
#array b[3]=[1.5,0.5,1.0]

#plot for[j=1:1] for[i=1:8] for[k=1:-1:-2] "-" index j-1 every :::i-1::i-1 u (k*(column("omega")/column("T"))):(-column("lgGf")/(exp((column("omega")/column("T")))+1.)) with l lw 2 lc k<0?i:i+1 dt k+2<2?k+2:2 notitle word("T1 T2 LM C MCK C' LM'",i%4+1)#,5*x**(-0.2)
#plot for[j=1:1] for[i=20:25] for[k=1:-1:-2] "-" index j-1 every :::i-1::i-1 u (k*(column("omega")/column("T"))):(-column("lgGf")/(exp(column("omega")/column("T"))+1.)) with l lw 2 lc k<0?i:i+1 dt k+2<2?k+2:2 notitle word("T1 T2 LM C MCK C' LM'",i%4+1)#,5*x**(-0.2)
#plot for[i=1:40] "-" index 0 every :::i-1::i-1 u "omega":(abs(column("iGB"))) with lp pt 7 ps 0.5 lc i dt i%2+1 notitle,x**-0.2# word("T1 T2 LM C MCK C' LM'",i%4+1)#, for[i=1:10] "-" index 1 every :::i-1::i-1 u ((column("omega"))):(-column("lgGf")) with l lw 2 lc i dt 0 notitle word("T1 T2 LM C MCK C' LM'",i%4+1)#,0.4*x**(-0.75)
#plot for[i=1:20] "-" index 0 every ::0:i-1:257*2-1:i-1 u "T":((-column("rGB"))) with p pt 7 ps 0.5 lc 1 notitle# word("T1 T2 LM C MCK C' LM'",i%4+1)#, for[i=1:10] "-" index 1 every :::i-1::i-1 u ((column("omega"))):(-column("lgGf")) with l lw 2 lc i dt 0 notitle word("T1 T2 LM C MCK C' LM'",i%4+1)#,0.4*x**(-0.75)
#plot for[i=1:10] for[j=1:3] "-" index 0 every :::i-1::i-1 u ((column("omega"))):(-column(word("iGf iGfa iGft",j))) with l lw 1 lc i dt j notitle word("T1 T2 LM C MCK C' LM'",i%4+1)
#plot for[i=1:10] "-" index 0 every :::i-1::i-1 u ((column("omega"))):(abs(column("iGft")-column("iGfa"))) with l lw 1 lc i notitle word("T1 T2 LM C MCK C' LM'",i%4+1)
#plot for[j=1:4] for[i=31:31] "-" index 40 every :::i-1::i-1 u (column("omega")):(abs(column(word("iSEf iSEB iGf iGB",j)))) with l lw 1 lc j title word("iSEf iSEB Af AB",j),x**(-0.666) with l dt 2 lc 4 title "x^{-0.666}",x**(-0.333) with l dt 2 lc 3 title "x^{-0.333}",x**(0.666) with l dt 2 lc 2 title "x^{0.666}",x**(0.333) with l dt 2 lc 1 title "x^{0.333}"
aa=0.2
#plot for[j=1:4] for[i=1:40] "-" index 0 every :::i-1::i-1 u (column("T")):(abs(column("omega"))<1.3e-11&&column("omega")>0?abs(column(word("iSEf rSEB iGf rGB rSEphi",j))+(j==2?1./0.51:0.)):1/0) with p pt 6 ps 1.0 lw 1. lc j notitle,0.5*x**(-aa) with l dt 2 lc 4 title "rGB",1.5*x**(aa-1.) with l dt 2 lc 3 title "iGf",2*x**(aa) with l dt 2 lc 2 title "rSEB",0.6*x**(1.-aa) with l dt 2 lc 1 title "iSEf"
#plot for[j=1:38] for[i=30:30] "-" index j-1 every :::i-1::i-1 u (1./column("J")):(abs(column("omega"))<1.3e-11&&column("omega")>0?abs(column(word("iSEf iSEB iGf iGB rSEphi",5))):1/0) with p pt 7 ps 0.8 lc 1 notitle, 0.25*exp(1./(0.5)*x)
#set arrow from 5.41e-4, graph 0 to 5.41e-4, graph 1 nohead dt 2
#plot "-" index 0 u "omega":"iGphi0" w l
#plot for[j=1:1] for[i=1:9:1] "-" index j-1 every 2:::i-1::i-1 u (column("omega")/column("T")):(-column("T")**0.6666*column("iGB")) with l ps 1.0 pt 7 lw 1 lc i title word("T/T@_K^{0}=10^{-1} 10^{-2} 10^{-3} 10^{-4} 10^{-5} 10^{-6} 10^{-7} 10^{-8} 10^{-9} T2 LM C MCK C' LM'",i)#,99*x**(-0.666)#,for[j=1:1] for[i=1:1:1] "-" index 0 every :::i-1::i-1 u (column("omega")/column("T")):(column("omega")/column("T")) with p ps 0.2 pt 7 lw 1 lc 3 title word("omega 10^{-2} 10^{-3} 10^{-4} 10^{-5} 10^{-6} 10^{-7} 10^{-8} 10^{-9} T2 LM C MCK C' LM'",i)
plot for[j=1:1] for[i=1:12*10:1] "-" index j-1 every 2:::i-1::i-1 u (column("omega")):(-column("iGf")) with l ps 1.0 pt 7 lw 1 lc i notitle word("T=10^{-1} 10^{-2} 10^{-3} 10^{-4} 10^{-5} 10^{-6} 10^{-7} 10^{-8} 10^{-9} 10^{-10} 10^{-11} 10^{-12}",i),for[j=1:1] for[i=1:12*10:1] 'largeN_omega_my_a_0.dat' index j-1 every 2:::i-1::i-1 u (column("omega")):(abs(column("iGf"))) with l dt 2 ps 1.0 pt 7 lw 1 lc i notitle word("T/T@_K^{0}=10^{-1} 10^{-2} 10^{-3} 10^{-4} 10^{-5} 10^{-6} 10^{-7} 10^{-8} 10^{-9}",i),99*x**(0.666-1) notitle
#plot for[i=1:80:2] "-" index 0 every :::i-1::i-1 u (column("omega")):(-column("iGf")) with l lw 2 lc i notitle word("LM C MCK C' LM'",i)#,for[i=0:80:2] "-" index 0 every :::i-1::i-1 u (column("omega")/column("T")):(-0.048*column("T")**0.23*column("iGf")) with l lw 2 lc i title word("LM C MCK C' LM'",i)
#plot for[i=1:1] "-" index id[i]-1 every :::y[i]::y[i] u (abs(column("omega"))<2e11?column("omega")*am[i]:1/0):(-column("rSEf")) with l lw 2 lc i title word("LM C MCK C' LM'",i), for[i=1:1] "-" index id[i]-1 every :::y[i]+1::y[i]+1 u (abs(column("omega"))<2e11?column("omega")*am[i]:1/0):(-column("rSEf")) with l lw 2 dt 2 lc i+1 notitle
#plot "-" index 0 u 2:5 with l notitle
set key font ",".(ticfontsize-1) at graph 0.39,0.49,1 horizontal maxcols 1 spacing 1 samplen 1.5 autotitle #opaque
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

