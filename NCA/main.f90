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
	real(wp) :: fcut=1.0_wp,fbeta=50._wp
	real(wp) :: bcut=10._wp,min_omega=1e-12_wp
	integer, parameter :: nlog=33,nl=14,nadd(2)=[80,10]
	integer, parameter :: n=(nl*nlog+sum(nadd))*2
	real(wp) :: omega(n),tau(n)
	real(wp) :: r=nan,kp=nan,V2=nan,ed=nan,gap=nan
	real(wp) :: Tk,beta,Tkunit=1._wp
contains
	subroutine set_omega_grid()
		!call set_grid(omega(n/2+1:),[reset,add_log_linear,(add_linear,i=1,size(nadd))],from=[min_omega,max(fcut-0.1_wp,0.05_wp),max(phicut-0.01_wp,0.005_wp),1.2_wp,bcut/2._wp],to=[bcut-1e-6,max(fcut+0.7_wp/fbeta,0.1_wp),max(phicut+0.7_wp/phibeta,0.01_wp),bcut/2._wp,bcut],n=[nlog,nl, nadd])
		call set_grid(omega(n/2+1:),[reset,add_log_linear,(add_linear,i=1,size(nadd))],from=[min_omega,0.5_wp,bcut/2._wp],to=[bcut-1e-6_wp,3._wp,bcut],n=[nlog,nl, nadd])
		!call set_grid(omega(n/2+1:),[reset,add_log_linear],from=[min_omega],to=[bcut],n=[nlog,nl])
		omega(1:n/2)=-omega(n:n/2+1:-1)
	end subroutine
end module
include "NCA_global.f90"
program main
	use global
	implicit none
	integer :: i,j,e
	complex(wp) :: SEsave(n)
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
	
	open(in_init,file="../data/init.dat")
	open(out_Gw,File="../data/NCA_omega.dat")
	open(out_Gt,File="../data/NCA_tau.dat")
	open(out_T,File="../data/NCA_T.dat")
	open(233,File="../data/check.dat")

	inode(tau_grid,1:1)=[PMT]
	inode(rhoc,1:1)=[PMT]
	inode(Gb,1:4)=[SEb,rhoc,lambda,PMT]
	inode(tSEb,1:3)=[Gf,rhoc,PMT]
	inode(tSEf,1:3)=[Gb,rhoc,PMT]
	inode(SEb,1:1)=[tSEb]
	inode(SEf,1:1)=[tSEf]
	inode(Gd,1:2)=[Gb,Gf]
	inode(Gf,1:3)=[SEf,lambda,PMT]
	inode(Gbt,1:3)=[Gb,tau_grid,PMT]
	inode(Gft,1:3)=[Gf,tau_grid,PMT]
	inode(SEc,1:3)=[Gd,rhoc,PMT]
	inode(Gc,1:2)=[SEc,rhoc]

	inode(tSEf,:)=0
	inode(Q,1:2)=[Gf,Gb]
	inode(tSEf,1:1)=[SC1]
	inode(lambda,1:1)=[SC1]
	inode(SC1,1:5)=[rhoc,Gb,Q,tSEf,lambda]
	underscore=evaluate(graph=[[Gbt,Gft],Gc,SEc,Gd,[SC1,Q,Gb,SEb,tSEb,Gf,SEf,tSEf,lambda],rhoc,tau_grid,PMT])

	!inode(tSEf,:)=0
	!inode(Q,1:2)=[Gf,Gb]
	!inode(tSEf,1:1)=[SC1]
	!inode(lambda,1:1)=[SEf]
	!inode(SC1,1:4)=[rhoc,Gb,Q,tSEf]
	!underscore=evaluate(graph=[[Gbt,Gft],SEc,Gd,[SC1,Q,Gb,SEb,tSEb,Gf,lambda,SEf,tSEf],rhoc,tau_grid,PMT])


	null=compute(omega_grid)
	call self%io([out_info],[out_info])
	niter(SC1)=200; iter_err(SC1)=1e-6_wp; iter_rate(SC1)=3e-3_wp
	niter(SC2)=200; iter_err(SC2)=1e-6_wp; iter_rate(SC2)=3e-3_wp


	r=0._wp
	kp=0.5_wp
	ed=-0.9_wp
	!ed=-0.04_wp
	!ed=0.0_wp
	!ed=0.02_wp
	!gap=0.1_wp
	
	self%tSEf=cmplx(0._wp,-ff(-omega)*ff(omega),kind=wp)
	self%lambda=0._wp
	do 
		!V2=for_in(x=[0.0063662_wp],id=2) ! T0=0.000187
		V2=for_in(x=[0.1_wp/pi],id=2) ! T0=0.000187
		!Tkunit=0.1_wp*exp(pi/2._wp*(-4._wp))
		if(any(isnan([V2,kp,r]))) then
			exit
		endif

		flag=.true.
		Tk=for_in(id=1)
		call self%io([out_param],[out_param])
		null=compute(rhoc)
		!SEsave(1)=nan
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
				[100:15:-5]*1e-2_wp&
				,[100:15:-5]*1e-3_wp&
				,[100:15:-5]*1e-4_wp&
				,[100:15:-5]*1e-5_wp&
				,[100:15:-5]*1e-6_wp&
				,[100:15:-5]*1e-7_wp&
				,[100:15:-5]*1e-8_wp&
				,[100:15:-5]*1e-9_wp&
				,[100:15:-5]*1e-10_wp&
				,[100]*1e-11_wp&
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
				],id=1)*Tkunit
			if(isnan(Tk)) then
				write(*,*)"End of T!!!"
				exit
			endif
			beta=1d0/Tk

			tau(1:n/2)=omega(n/2+1:n)/(3._wp*omega(n)-omega(n-1))/Tk
			tau(n/2+1:)=1._wp/Tk-tau(n/2:1:-1)

			self%conv=evaluate(nodes=[Gbt,Gft,Gc,Gd],updated_nodes=[lambda,tSEf,PMT])

			if(self%conv) then
				call self%io([out_T],[out_T])
				!if(isnan(abs(SEsave(1)))) then
					!SEsave=self%SEf
				!endif
				if(any(abs(Tk/Tkunit-[1e4_wp,1e3_wp,1e2_wp,1e1_wp,1e0_wp,1e-1_wp,1e-2_wp,1e-3_wp,1e-4_wp,1e-5_wp,1e-6_wp,1e-7_wp,1e-8_wp,1e-9_wp,1e-10_wp,1e-11_wp,1e-12_wp,1e-13_wp])<1e-15_wp)) then
					call self%io([out_Gw,out_Gt],[out_Gw,out_Gt])
				endif
			endif
		enddo
		if(flag) then
			write(out_T,"(x)")
			write(out_Gw,"(x)")
			write(out_Gt,"(x)")
		endif
		!self%SEf=SEsave
	enddo
end program
