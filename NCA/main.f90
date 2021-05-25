include "../largeN/fn.f90"
include "../lib/utility.f90"
include "../largeN/mbroyden.f90"
module param
	use M_fn
	use M_solution
	use M_const
	use M_utility
	use M_matrix
	use mkl_service
	real(wp) :: fcut=1.0_wp,fbeta=15._wp
	real(wp) :: bcut=15._wp,min_omega=1e-10_wp
	integer, parameter :: nlog=13*2,nl=24&
		!,nadd(7)=[50,20,30,100,100,30,10] !fcut, ed, igs(1) igs(2) ed+U bcut
		,nadd(5)=[50,20,30,30,10] !fcut, ed, igs(1) igs(2) ed+U bcut
		!,nadd(7)=[50,20,20,30,30,30,10] !fcut, ed, igs(1) igs(2) ed+U bcut
	integer, parameter :: n=(nl*nlog+sum(nadd))*2,nk=32*8,D=3
	real(wp) :: Tk,beta,Tkunit=1._wp,boson=0.05_wp,etab=0.005_wp,eta=0.01_wp
	real(wp) :: omega(n),tau(n),omega_ek(600)&
	!,igs(3)=[0.1_wp,1.1_wp,1._wp-0.1_wp]
	,igs(1)=[0.1_wp]
	real(wp) :: r=nan,V2(1,1)=nan,ed(1)=nan,gap=nan,scgap=nan,U=nan,GM=nan,g2=nan
	integer :: ch=1
	logical :: ocd=&
		!.true.
		.false.
contains
	subroutine set_omega_grid()
		real(wp) :: gap_,U_
		gap_=merge(gap,scgap,isnan(scgap))
		U_=merge(0._wp,U,isnan(U))
		!call set_grid(omega(n/2+1:),[reset,add_log_linear,(add_linear,i=1,size(nadd))],from=[min_omega,max(fcut-0.1_wp,0.05_wp),max(phicut-0.01_wp,0.005_wp),1.2_wp,bcut/2._wp],to=[bcut-1e-6,max(fcut+0.7_wp/fbeta,0.1_wp),max(phicut+0.7_wp/phibeta,0.01_wp),bcut/2._wp,bcut],n=[nlog,nl, nadd])
		if(isnan(gap_)) then
			call set_grid(omega(n/2+1:),[reset,add_log_linear,(add_linear,i=1,size(nadd))],&
				from=[min_omega,0.5_wp,abs(ed)*0.5_wp,abs(ed+U_)*0.5_wp,igs*0.8_wp,bcut/2._wp],&
				to=[bcut-1e-6_wp,3._wp,abs(ed)*1.5_wp,abs(ed+U_)*1.5_wp,igs*1.2_wp,bcut],n=[nlog,nl, nadd])
		else
			call set_grid(omega(n/2+1:),[reset,add_log_linear,(add_linear,i=1,size(nadd))],&
				from=[min_omega,gap_*0.9_wp,abs(ed)*0.5_wp,abs(ed+U_)*0.5_wp,igs*0.5_wp*gap_,bcut/2._wp],&
				to=[bcut-1e-6_wp,gap_*1.3_wp,abs(ed)*1.5_wp,abs(ed+U_)*1.5_wp,igs*gap_*1.5_wp,bcut],n=[nlog,nl, nadd])
		endif
		!call set_grid(omega(n/2+1:),[reset,add_log_linear],from=[min_omega],to=[bcut],n=[nlog,nl])
		do i=n/2+1,n-1
			if(omega(i+1)<omega(i)) then
				write(*,*)"error, omega is not increase"
				stop
			elseif(omega(i+1)-omega(i)<1e-14_wp) then
				write(*,*)"warning, omega is too close",omega(i),omega(i+1)
				omega(i)=0.5_wp*(omega(max(i-1,n/2+1))+omega(i+1))
			endif
		enddo
		omega(1:n/2)=-omega(n:n/2+1:-1)

		call set_grid(omega_ek,[reset,add_linear],from=[-9._wp],to=[9._wp],n=[size(omega_ek)])
		omega_ek=omega_ek-(omega_ek(1)+omega_ek(size(omega_ek)))*0.5_wp
	end subroutine
end module
include "NCA_global.f90"
program main
	use global
	implicit none
	integer :: i,j,e
	complex(wp) :: SEsave(n,size(ed),size(ed))
	logical :: flag,last
	real(wp) :: a_f,a_b,tmp(1)
	compute => compute_

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	!call omp_set_schedule(omp_sched_static,0)

	!call mkl_set_num_threads(24)
	call omp_set_num_threads(mkl_get_max_threads())
	!call omp_set_num_threads(1)
	
	open(in_init,file="../data/init.dat")
	open(out_Gw,File="../data/NCA_omega.dat")
	open(out_Gt,File="../data/NCA_tau.dat")
	open(out_T,File="../data/NCA_T.dat")
	open(233,File="../data/check.dat")

	r=0._wp
	ch=1
	!ed=-0.15_wp
	ed=-0.6_wp
	!ed=-0.55_wp
	!ed=-0.02_wp*pi
	!ed=0.0_wp
	!ed=-0.02_wp
	!gap=2.07e-6_wp
	!scgap=0.01_wp
	!gap=0.01_wp
	!U=-ed(1)*2._wp
	U=1.2_wp
	GM=0.15_wp
	!GM=0.08_wp
	!g2=0.005_wp
	!ed(2)=ed(1)+0.05_wp
	niter(SC1)=100; iter_err(SC1)=1e-5_wp; iter_rate(SC1)=5e-1_wp
	niter(SC2)=300; iter_err(SC2)=1e-5_wp; iter_rate(SC2)=3e-3_wp


	!ed=-0.01_wp
	!!U=0.102_wp
	!U=0.005_wp
	!GM=0.01_wp
	!Tkunit=sqrt(U*GM*0.5_wp)*exp(pi*ed*(ed+U)/(2._wp*U*GM))
	!!gap=4e-4_wp*Tkunit
	!gap=1e-2_wp*Tkunit


	inode(tau_grid,1:1)=[PMT]
	inode(rhoc,1:1)=[PMT]
	inode(rhophi,1:1)=[PMT]
	inode(SEb,1:1)=[tSEb]
	inode(SEf,1:1)=[tSEf]
	inode(Gf,1:4)=[SEf,tSEf,lambda,PMT]
	inode(Gb,1:4)=[SEb,tSEb,lambda,PMT]
	inode(Gbt,1:3)=[Gb,tau_grid,PMT]
	inode(Gft,1:3)=[Gf,tau_grid,PMT]
	inode(SEc,1:3)=[Gd,rhoc,PMT]
	inode(Gc,1:2)=[SEc,rhoc]

	inode(tSEf,1:1)=[SC1]
	inode(lambda,1:1)=[SC1]
	!if(isnan(U)) then
		!inode(Gd,1:2)=[Gb,Gf]
		!inode(tSEb,1:3)=[Gf,rhoc,PMT]
		!inode(Q,1:2)=[Gf,Gb]
		!inode(SC1,1:5)=[rhoc,Gb,Q,tSEf,lambda]
		!inode(Fimp,1:2)=[lambda,Q]
		!underscore=evaluate(graph=[[Gbt,Gft],Fimp,Gc,SEc,Gd,[SC1,Q,Gb,SEb,tSEb,Gf,SEf,tSEf,lambda],rhoc,tau_grid,PMT])
	!else
		inode(Ga,1:4)=[SEa,tSEa,lambda,PMT]
		inode(SEa,1:1)=[tSEa]
		inode(Gd,1:3)=[Gb,Gf,Ga]
		inode(Q,1:3)=[Gf,Gb,Ga]
		if(ocd) then
			inode(tSEa,1:4)=[Gf,Gb,rhoc,PMT]
			!inode(tSEa,1:1)=[SC1]
			inode(tSEb,1:1)=[SC1]
			inode(SC1,1:8)=[rhoc,Gb,Ga,Gf,Q,tSEb,tSEf,lambda]
			inode(Fimp,1:2)=[lambda,Q]
			underscore=evaluate(graph=[[Gbt,Gft],Fimp,Gc,SEc,Gd,[SC1,Q,Ga,SEa,tSEa,Gb,SEb,Gf,SEf,tSEb,tSEf,lambda],rhoc,tau_grid,PMT])
			!underscore=evaluate(graph=[[Gbt,Gft],Fimp,Gc,SEc,Gd,[SC1,Q,Ga,SEa,Gb,SEb,Gf,SEf,tSEa,tSEb,tSEf,lambda],rhoc,tau_grid,PMT])
		else
			inode(tSEa,1:4)=[Gf,rhoc,rhophi,PMT]
			inode(tSEb,1:3)=[Gf,rhoc,PMT]
			inode(SC1,1:7)=[rhoc,rhophi,Gb,Ga,Q,tSEf,lambda]
			inode(Fimp,1:2)=[lambda,Q]
			inode(optcond,1:2)=[Gd,PMT]
			underscore=evaluate(graph=[[Gbt,Gft],optcond,Fimp,Gc,SEc,Gd,[SC1,Q,Ga,SEa,tSEa,Gb,SEb,tSEb,Gf,SEf,tSEf,lambda],rhophi,rhoc,tau_grid,PMT])
		endif
	!endif

	!inode(tSEf,:)=0
	!inode(Ga,1:4)=[SEa,tSEa,lambda,PMT]
	!inode(SEa,1:1)=[tSEa]
	!inode(tSEa,1:4)=[Gf,rhoc,rhophi,PMT]
	!inode(tSEb,1:3)=[Gf,rhoc,PMT]
	!inode(Q,1:2)=[Gf,Gb]
	!inode(tSEf,1:1)=[SC1]
	!inode(lambda,1:1)=[SC2]
	!inode(SC1,1:3)=[rhoc,Gb,tSEf]
	!inode(SC2,1:2)=[Q,lambda]
	!inode(Fimp,1:2)=[lambda,Q]
	!underscore=evaluate(graph=[[Gbt,Gft],Fimp,Gc,SEc,Gd,[SC1,[SC2,Q,Ga,SEa,tSEa,Gb,SEb,tSEb,Gf,lambda],SEf,tSEf],rhophi,rhoc,tau_grid,PMT])


	!inode(tSEf,:)=0
	!inode(Q,1:2)=[Gf,Gb]
	!inode(tSEf,1:1)=[SC1]
	!inode(tlambda,1:2)=[SEb,SEf]
	!inode(SC1,1:5)=[tlambda,rhoc,Gb,Q,tSEf]
	!inode(Gf,1:4)=[SEf,tSEf,tlambda,PMT]
	!inode(Gb,1:4)=[SEb,tSEb,tlambda,PMT]
	!inode(Fimp,1:2)=[tlambda,Q]
	!underscore=evaluate(graph=[[Gbt,Gft],Fimp,Gc,SEc,Gd,[SC1,Q,Gb,SEb,tSEb,Gf,SEf,tlambda,tSEf],rhoc,tau_grid,PMT])

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!inode(tSEf,:)=0
	!inode(tSEb,:)=0
	!inode(Q,1:2)=[Gf,Gb]
	!inode(tSEf,1:1)=[SC1]
	!inode(tSEb,1:1)=[SC1]
	!inode(lambda,1:1)=[SC1]
	!inode(SC1,1:7)=[rhoc,Gb,Gf,Q,tSEf,tSEb,lambda]
	!inode(Fimp,1:2)=[lambda,Q]
	!underscore=evaluate(graph=[[Gbt,Gft],Fimp,Gc,SEc,Gd,[SC1,Q,Gb,Gf,SEb,SEf,tSEb,tSEf,lambda],rhoc,tau_grid,PMT])

	!inode(tSEf,:)=0
	!inode(tSEb,:)=0
	!inode(Q,1:2)=[Gf,Gb]
	!inode(tSEf,1:1)=[SC1]
	!inode(tSEb,1:1)=[SC1]
	!inode(SC1,1:5)=[rhoc,Gb,Gf,tSEf,tSEb]
	!inode(lambda,1:1)=[SC2]
	!inode(SC2,1:2)=[Q,lambda]
	!inode(Fimp,1:2)=[lambda,Q]
	!underscore=evaluate(graph=[[Gbt,Gft],Fimp,Gc,SEc,Gd,[SC1,[SC2,Q,Gb,Gf,lambda],SEb,SEf,tSEb,tSEf],rhoc,tau_grid,PMT])

	!inode(tSEf,:)=0
	!inode(tSEb,:)=0
	!inode(Q,1:2)=[Gf,Gb]
	!inode(tSEf,1:1)=[SC1]
	!inode(tSEb,1:1)=[SC1]
	!inode(tlambda,1:2)=[SEb,SEf]
	!inode(SC1,1:7)=[tlambda,rhoc,Gb,Gf,Q,tSEf,tSEb]
	!inode(Gf,1:4)=[SEf,tSEf,tlambda,PMT]
	!inode(Gb,1:4)=[SEb,tSEb,tlambda,PMT]
	!inode(Fimp,1:2)=[tlambda,Q]
	!underscore=evaluate(graph=[[Gbt,Gft],Fimp,Gc,SEc,Gd,[SC1,Q,Gb,Gf,tlambda,SEf,SEb,tSEb,tSEf],rhoc,tau_grid,PMT])




	null=compute(omega_grid)
	call self%io([out_info],[out_info])
	
	do i=1,n
		self%tSEf(i,:,:)=cmplx(0._wp,-ff(-omega(i))*ff(omega(i)),kind=wp)
	enddo
	self%tSEb=cmplx(0._wp,-ff(-omega)*ff(omega),kind=wp)
	self%tSEa=cmplx(0._wp,-ff(-omega)*ff(omega),kind=wp)
	self%lambda=0._wp
	null=compute(lattice)
	null=compute(rhoc)
	do 
		!V2=2._wp/pi*for_in(x=[GM],id=2) ! T0=0.000187
		!V2=for_in(x=[0.1_wp],id=2)/pi*2._wp ! T0=0.000187
		!V2=for_in(x=[0.001_wp],id=2)/pi ! T0=0.000187
		!V2=for_in(x=[0.0098_wp],id=2)/pi ! T0=0.000187
		!Tkunit=0.1_wp*exp(pi/2._wp*(-4._wp))
		!gap=for_in(x=[0.045_wp],id=2)*Tkunit ! T0=0.000187
		!V2=2._wp/pi*GM
		!V2=2._wp/pi*GM
		!V2(1,2)=V2(1,2)*0.5_wp
		!V2(2,1)=V2(2,1)*0.5_wp
		!V2(2,2)=V2(2,2)*0.5_wp
		V2(1,1)=0.1_wp
		!V2(1,2)=0.05_wp
		!V2(2,1)=0.05_wp
		!V2(2,2)=0.15_wp
		V2=V2/pi!*2._wp
		write(*,*)"Kondo Temperature and gap is: ",Tkunit,gap/Tkunit
		!g2=for_in(x=[0:80:10]*0.00025_wp,id=2)/pi ! T0=0.000187
		!if(any(isnan([g2,V2]))) then
			!exit
		!endif

		flag=.true.
		Tk=for_in(id=1)
		call self%io([out_param],[out_param])
		SEsave(1,1,1)=nan
		do
			!Tkunit=exp(-1._wp/(self%rhoc(n/2)*Jk))
			Tk=for_in(x=[&
				![100:15:-5]*1e-1_wp&
				!,[100:15:-5]*1e-2_wp&
				!,[100:20:-10]*1e-3_wp&
				!,[100:20:-10]*1e-4_wp&
				!,[100]*1e-5_wp&
				![100:15:-5]*1e-1_wp&
				![15:100:5]*1e-1_wp&
				!,[15:100:5]*1e-0_wp&
				!,[15:100:5]*1e1_wp&
				!,[15:100:5]*1e2_wp&
				![100:20:-50]*1e2_wp&
				!,[100:20:-50]*1e1_wp&
				!,[100:20:-50]*1e-0_wp&
				![100:20:-10]*1e-1_wp&
				![100:20:-10]*1e-2_wp&
				[20:50:+10]*1e-2_wp&
				,[100:20:-10]*1e-3_wp&
				,[100:20:-10]*1e-4_wp&
				,[100:20:-10]*1e-5_wp&
				,[100:20:-20]*1e-6_wp&
				,[100:20:-20]*1e-7_wp&
				,[100]*1e-8_wp&
				!,[100:20:-20]*1e-8_wp&
				!,[100:20:-20]*1e-9_wp&
				!,[100:20:-20]*1e-10_wp&
				!,[100]*1e-11_wp&
				![100:20:-50]*1e-2_wp&
				!,[100:20:-50]*1e-3_wp&
				!,[100:20:-50]*1e-4_wp&
				!,[100:20:-50]*1e-5_wp&
				!,[100:20:-50]*1e-6_wp&
				!,[100:20:-50]*1e-7_wp&
				!,[100:20:-50]*1e-8_wp&
				!,[100:20:-50]*1e-9_wp&
				!,[100:20:-50]*1e-10_wp&
				!,[100]*1e-11_wp&
				!,[100:20:-10]*1e-9_wp&
				!,[100:20:-10]*1e-10_wp&
				!,[100]*1e-11_wp&
				!,[100:20:-10]*1e-10_wp&
				!,[100:20:-10]*1e-11_wp&
				!,[100:20:-10]*1e-13_wp&
				!,[100:20:-10]*1e-14_wp&
				!,[100:20:-10]*1e-15_wp&
				!,[100:20:-10]*1e-16_wp&
				!,[100:20:-10]*1e-15_wp&
				!,[100:100:-50]*1e-16_wp&
				!,[100:100:-50]*1e-17_wp&
				!,[100:100:-50]*1e-18_wp&
				!,[100]*1e-11_wp&
				!,[100:20:-20]*1e-12_wp&
				!,[96:10:-2]*1e-7_wp&
				!,[96:10:-2]*1e-8_wp&
				!,[96:10:-2]*1e-9_wp&
				!,[96:10:-2]*1e-10_wp&
				!,[99:10:-1]*1e-11_wp&
				!,[99:10:-1]*1e-12_wp&
				],id=1)!*Tkunit
			if(isnan(Tk)) then
				write(*,*)"End of T!!!"
				exit
			endif
			!if(Tk<1.101e-3_wp) then
				!Tk=Tk+0.0005_wp
				!!niter(SC1)=-10
				!!iter_err(SC1)=1e-5_wp;
				!iter_rate(SC1)=3e-3_wp
			!endif
			beta=1d0/Tk


			tau(1:n/2)=omega(n/2+1:n)/(3._wp*omega(n)-omega(n-1))/Tk
			tau(n/2+1:)=1._wp/Tk-tau(n/2:1:-1)

			if(any(abs(Tk/Tkunit-[1e4_wp,1e3_wp,1e2_wp,1e1_wp,1e0_wp,1e-1_wp,1e-2_wp,1e-3_wp,1e-4_wp,1e-5_wp,1e-6_wp,1e-7_wp,1e-8_wp,1e-9_wp,1e-10_wp,1e-11_wp,1e-12_wp,1e-13_wp])<1e-15_wp)) then
				self%conv=evaluate(nodes=[&
					optcond,&
					Fimp,Gbt,Gft,Gc,Gd],updated_nodes=[onode(SC1,1:onode(SC1,0)),rhoc,PMT])
			else
				self%conv=evaluate(nodes=[Fimp,Gbt,Gft,Gc,Gd],updated_nodes=[onode(SC1,1:onode(SC1,0)),rhoc,PMT])
			endif

			if(self%conv) then
				call self%io([out_T],[out_T])
				if(isnan(abs(SEsave(1,1,1)))) then
					SEsave=self%SEf
				endif
				if(any(abs(Tk/Tkunit-[1e4_wp,1e3_wp,1e2_wp,1e1_wp,1e0_wp,5e-1_wp,1e-1_wp,1e-2_wp,1e-3_wp,1e-4_wp,1e-5_wp,1e-6_wp,1e-7_wp,1e-8_wp,1e-9_wp,1e-10_wp,1e-11_wp,1e-12_wp,1e-13_wp])<1e-15_wp)) then
				!if(Tk<1e-3_wp) then
					call self%io([out_Gw,out_Gt],[out_Gw,out_Gt])
					write(*,*)integrate(omega,A=-self%Gd%im/pi)
				endif
				!if(abs(Tk-1e-5_wp)<1e-9_wp) then
					!!call get_convolution2(omega,ga=-2._wp*V2**2*self%Gc0([n:1:-1]),gb=self%Gc0([n:1:-1]),A=self%Gf,B=self%Gf,C=self%Gb,rt=self%tSEa)
					!!call get_convolution2(omega,ga=-2._wp*V2**2*self%Gc0,gb=self%Gc0,A=self%Gf,B=self%Gf,C=self%Ga,rt=self%tSEb)
					!call get_convolution2(omega,ga=-2._wp*V2**2*self%Gc0([n:1:-1]),gb=self%Gc0,A=self%Gb,B=self%Ga,C=self%Gf,rt=self%tSEf)
				!endif
			endif
			!iter_rate(SC1)=3e-3_wp
		enddo
		if(flag) then
			write(out_T,"(x)")
			write(out_Gw,"(x)")
			write(out_Gt,"(x)")
		endif
		self%SEf=SEsave
		exit
	enddo
end program
