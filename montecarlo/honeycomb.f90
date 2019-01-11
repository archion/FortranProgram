module coarray
	integer, parameter :: ica1=1
	integer :: ica2
end module
include "../lib/utility.f90"
include "../lib/hamilton_final_1.f90"
include "vmc_utility.f90"
module global
	use vmc_utility
	use ifport, ifort_qsort => qsort
	use omp_lib
	use mkl_service
	implicit none
	real(8) :: t(3)=[1._wp,0.1_wp/2.8_wp,0.07_wp/2.8_wp]
	!real(8) :: t(1)=[1._wp]
	real(8), parameter :: V=0._wp,U=3.6_wp
	type(t_mc) :: mc[ica1,*]
	integer :: isc
contains
	subroutine initial()
		integer :: i,l,idx,sb,tp,tmp(3)
		real(8) :: q(2,ica1)=reshape([&
			pi,0._wp&
			],[2,ica1])
!***************************parameter setting*************************
		is_project=.false.
		is_ph=.true.
		dx=0.3_wp
		call init_random_seed()

!****************************lattice**********************************
		latt%is_all=.true.
		latt%a1=[1.5_wp,sqrt(3._wp)/2._wp,0._wp]
		latt%a2=[0._wp,sqrt(3._wp),0._wp]
		latt%c1=latt%a1
		latt%c2=latt%a2
		!latt%T1=[3._wp,0._wp,0._wp]*4._wp
		!latt%T2=[0._wp,sqrt(3._wp),0._wp]*7._wp
		latt%T1=latt%a1*12._wp
		latt%T2=latt%a2*12._wp
		latt%bdc=[1._wp,1._wp,0._wp]
		allocate(latt%rsb(2,3))
		latt%rsb(1,:)=[0._wp,0._wp,0._wp]
		latt%rsb(2,:)=[1._wp,0._wp,0._wp]
		call latt%gen_latt(size(t))
		call latt%gen_brizon(brizon)
		if(this_image()==1) then
			call check_lattice(101)
			write(*,*)"Total site number is: ",latt%Ns
		endif
		Ns=latt%Ns

!*********************************************************************
		allocate(Ham%var(-10:10),Hmf%var(-10:10),Hja%var(-10:10),mc%sphy%var(-10:10),mc%dphy%var(-10:10))

!*************************meanfield***********************************
		!! cp
		!idx=Hmf%add(nb=0,ca=[(c("i",i,+1),c("i",i,-1),c("i",i,+2),c("i",i,-2),i=1,2)],n=2,sg=[1._wp,-1._wp,1._wp,-1._wp],label="cp",is_var=.false.)
		!Hmf%var(idx)%bd=-1._wp
		!!Hmf%var(idx)%val=-1._wp
		!Hmf%var(idx)%val=1.0218_wp

		! d+id
		idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",2,-2),c("j",2,+1),c("i",1,-2),c("i",1,+2),c("j",2,-1),c("j",2,+2),c("i",1,-1)],n=2,cg=[.false.,.false.,.true.,.true.],is_var=.false.)
		isc=idx
		do i=1,size(Hmf%var(idx)%bd)
			if(abs(latt%nb(1)%bd(i)%dr(2))<1e-6_wp) then
				Hmf%var(idx)%bd(i)=1._wp
			elseif(latt%nb(1)%bd(i)%dr(1)*latt%nb(1)%bd(i)%dr(2)>0._wp) then
				Hmf%var(idx)%bd(i)=exp(img*4._wp/3._wp*pi)
			else
				Hmf%var(idx)%bd(i)=exp(img*2._wp/3._wp*pi)
			endif
		enddo
		Hmf%var(idx)%val=1.9e-2_wp


		idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",2,-1),c("j",2,+1),c("i",1,-1),c("i",1,+2),c("j",2,-2),c("j",2,+2),c("i",1,-2)],n=2,sg=[+1._wp,+1._wp,-1._wp,-1._wp],is_var=.false.)
		Hmf%var(idx)%bd=-1._wp
		Hmf%var(idx)%val=1._wp

		if(size(t)>1) then
			idx=Hmf%add(nb=2,ca=[(c("i",i,+1),c("j",i,-1),c("j",i,+1),c("i",i,-1),c("i",i,+2),c("j",i,-2),c("j",i,+2),c("i",i,-2),i=1,2)],n=2,sg=[+1._wp,+1._wp,-1._wp,-1._wp,+1._wp,+1._wp,-1._wp,-1._wp],is_var=.false.)
			Hmf%var(idx)%bd=-1._wp
			Hmf%var(idx)%val=t(2)/t(1)
		endif

		if(size(t)>2) then
			idx=Hmf%add(nb=3,ca=[c("i",1,+1),c("j",2,-1),c("j",2,+1),c("i",1,-1),c("i",1,+2),c("j",2,-2),c("j",2,+2),c("i",1,-2)],n=2,sg=[+1._wp,+1._wp,-1._wp,-1._wp],is_var=.false.)
			Hmf%var(idx)%bd=-1._wp
			Hmf%var(idx)%val=t(3)/t(1)
		endif

!************************jastrow**************************************
		idx=Hja%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2),c("i",2,+1),c("i",2,-1),c("i",2,-2),c("i",2,+2)],n=4)
		Hja%var(idx)%bd=-1._wp
		Hja%var(idx)%val=5.84e-1_wp

!*************************Hamiltionian********************************
		idx=Ham%add(nb=1,ca=[c("i",1,+1),c("j",2,-1),c("j",2,+1),c("i",1,-1),c("i",1,-2),c("j",2,+2),c("j",2,-2),c("i",1,+2)],n=2)
		Ham%var(idx)%bd=-1._wp
		Ham%var(idx)%val=t(1)

		if(size(t)>1) then
			idx=Ham%add(nb=2,ca=[(c("i",sb,+1),c("j",sb,-1),c("j",sb,+1),c("i",sb,-1),c("i",sb,-2),c("j",sb,+2),c("j",sb,-2),c("i",sb,+2),sb=1,2)],n=2)
			Ham%var(idx)%bd=-1._wp
			Ham%var(idx)%val=t(2)
		endif

		if(size(t)>2) then
			idx=Ham%add(nb=3,ca=[c("i",1,+1),c("j",2,-1),c("j",2,+1),c("i",1,-1),c("i",1,-2),c("j",2,+2),c("j",2,-2),c("i",1,+2)],n=2)
			Ham%var(idx)%bd=-1._wp
			Ham%var(idx)%val=t(3)
		endif

		! hubbard
		idx=Ham%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2),c("i",2,+1),c("i",2,-1),c("i",2,-2),c("i",2,+2)],n=4)
		Ham%var(idx)%bd=1._wp
		Ham%var(idx)%val=U

		!! t-J
		!is_project=.true.
		!idx=Ham%add(nb=1,ca=[c("i",1,+1),c("i",1,+2),c("j",2,-2),c("j",2,-1),c("j",2,+1),c("j",2,+2),c("i",1,-2),c("i",1,-1)],n=4,V=DJ/2._wp,label="J")
		!Ham%var(idx)%bd=1._wp
		!Ham%var(idx)%val=1._wp

		!idx=Ham%add(nb=1,ca=[c("i",1,+1),c("i",1,-1),c("j",2,+1),c("j",2,-1),c("i",1,-2),c("i",1,+2),c("j",2,-2),c("j",2,+2),c("i",1,+1),c("i",1,-1),c("j",2,-2),c("j",2,+2),c("i",1,-2),c("i",1,+2),c("j",2,+1),c("j",2,-1)],n=4,sg=[V+DJ/4._wp,V+DJ/4._wp,V-DJ/4._wp,V-DJ/4._wp],label="4term")
		!Ham%var(idx)%bd=1._wp
		!Ham%var(idx)%val=1._wp


!*************************static measurement**************************
		!idx=mc%sphy%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2),c("i",2,+1),c("i",2,-1),c("i",2,-2),c("i",2,+2)],n=2,sg=[1._wp,-1._wp,1._wp,-1._wp],V=1._wp/Ns,label="1s")
		!do i=1,size(mc%sphy%var(idx)%bd)
			!if(abs(latt%nb(0)%bd(i)%sb(1))==1) then
				!mc%sphy%var(idx)%bd(i)=1._wp
			!else
				!mc%sphy%var(idx)%bd(i)=-1._wp
			!endif
		!enddo
		!mc%sphy%var(idx)%val=0._wp

		!idx=mc%sphy%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2),c("i",2,+1),c("i",2,-1),c("i",2,-2),c("i",2,+2)],n=2,sg=[1._wp,1._wp,1._wp,1._wp],V=1._wp/Ns,label="1n")
		!mc%sphy%var(idx)%bd=1._wp
		!mc%sphy%var(idx)%val=0._wp

		!idx=mc%sphy%add(nb=1,ca=[c("i",1,+1),c("j",2,-2),c("j",2,+1),c("i",1,-2)],n=2,V=1._wp/(Ns*Ns),label="2sc",extdat=[real(8)::ichar("r"),1e-7_wp])
		!do i=1,size(mc%sphy%var(idx)%bd)
			!if(abs(latt%nb(1)%bd(i)%dr(2))<1e-6_wp) then
				!mc%sphy%var(idx)%bd(i)=1._wp
			!elseif(latt%nb(1)%bd(i)%dr(1)*latt%nb(1)%bd(i)%dr(2)>0._wp) then
				!mc%sphy%var(idx)%bd(i)=exp(img*4._wp/3._wp*pi)
			!else
				!mc%sphy%var(idx)%bd(i)=exp(img*2._wp/3._wp*pi)
			!endif
		!enddo
		!mc%sphy%var(idx)%val=0._wp


!************************dynamic measurement**************************
		!idx=mc%dphy%add(nb=0,ca=[c("i",1,+1),c("i",1,+2),c("i",2,+1),c("i",2,+2)],n=2,label="2s",extdat=[q(:,this_image(mc,1)),0._wp])
		!mc%dphy%var(idx)%bd=1._wp
		!mc%dphy%var(idx)%val=0._wp

		call Hmf%init()
		call Hja%init()
		call Ham%init()
		call mc%sphy%init()
		call mc%dphy%init()
	end subroutine
end module
program main
	use global
	implicit none
	logical :: f
	integer :: i,j,l
	real(8) :: et,E
	character(10) :: hostname
	i=hostnm(hostname)
	ica2=num_images()/ica1
	if(num_images()/=ica2*ica1) stop "plz check the coarray number"
	write(*,*)"runing ",this_image()," in ",hostname
	sync all
	if(this_image()==1) then
		write(*,*)"coarray info: ",ica1,ica2,"out of",num_images()
		open(101,file="../data/lattice.dat")
		open(111,file="../data/tmp.dat")
		open(20,file="../data/var.dat")
		open(30,file="../data/phyri.dat")
		open(40,file="../data/phyvar.dat")
		open(50,file="../data/matrix.dat",access="stream")
		open(70,file="../data/spect.dat")
		open(71,file="../data/spect_kmap.dat")
		open(80,file="../data/band.dat")
	endif

	call omp_set_nested(.false.)
	!call omp_set_max_active_levels(1)
	call omp_set_dynamic(.false.)
	call mkl_set_dynamic(0)
	call omp_set_schedule(omp_sched_static,0)

	call mkl_set_num_threads(1)
	call omp_set_num_threads(mkl_get_max_threads())

	call initial()
	if(this_image()==1) then
		write(50)size(brizon%k,1),brizon%k,brizon%nk
	endif
	otime=0._wp

	mc%ne(1)=Ns/2+nint(Ns/8._wp)
	if(is_ph) then
		mc%ne(2)=Ns-mc%ne(1)
	else
		mc%ne(2)=mc%ne(1)
	endif

	call mc%init(.true.)
	if(this_image()==1) then
		call Hmf%band(80,[0._wp,0._wp,0._wp],brizon%Ta(1,:),100)
		call Hmf%band(80,brizon%Ta(1,:),(brizon%Ta(1,:)+brizon%Ta(2,:))/2._wp,100)
		call Hmf%band(80,(brizon%Ta(1,:)+brizon%Ta(2,:))/2._wp,[0._wp,0._wp,0._wp],100)
	endif
	stop

	mc%hot=1024*16*8
	mc%step=nint(sqrt(real(Ns)))

	Hmf%var(isc)%val=1e-4_wp
	mc%samp=1024*16*8
	call mc%do_var(20)
	mc%samp=mc%samp*8*8*8*2
	mc%sg=1
	do
		!Hmf%var(isc)%val=for_in([1e-3_wp,3e-3_wp,7e-3_wp,2e-2_wp],id=1)
		!Hmf%var(isc)%val=for_in([1e-2_wp],id=1)
		if(isnan(Hmf%var(isc)%val(1))) exit
		call mc%init(.true.)
		call mc%do_vmc()
		if(this_image()==1) then
			write(*,*)Hmf%var(isc)%val,mc%E,mc%err
		endif
		exit
		!Hmf%var(isc)%val=for_in([1e-3_wp,2e-3_wp,3e-3_wp,5e-3_wp,7e-3_wp,1e-2_wp,2e-2_wp,3e-2_wp],id=1)
		!Hmf%var(isc)%val=for_in([2e-3_wp,5e-3_wp,1e-2_wp,3e-2_wp],id=1)
		!if(isnan(Hmf%var(isc)%val(1))) exit
	enddo

end program
