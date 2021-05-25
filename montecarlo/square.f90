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
	real(wp) :: t(1)=[1._wp]
	!real(8) :: t(1)=[1._wp]
	real(wp), parameter :: V=0._wp,U=4._wp,DJ=0.3_wp
	type(t_mc) :: mc[ica1,*]
	integer :: ivar
contains
	subroutine initial()
		integer :: i,l,idx,sb,tp
		real(8) :: q(2,ica1)=reshape([&
			pi,0._wp&
			],[2,ica1])
!***************************parameter setting*************************
		is_project=.false.
		is_ph=.true.
		dx=0.1_wp
		!call init_random_seed()

!****************************lattice**********************************
		latt%is_all=.true.
		latt%a1=[1._wp,0._wp,0._wp]
		latt%a2=[0._wp,1._wp,0._wp]
		latt%c1=latt%a1
		latt%c2=latt%a2
		latt%c1=[1._wp,1._wp,0._wp]
		latt%c2=[-1._wp,1._wp,0._wp]
		!latt%c1=[4._wp,1._wp,0._wp]
		!latt%c2=[0._wp,2._wp,0._wp]
		!latt%c1=[8._wp,0._wp,0._wp]
		!latt%c2=[0._wp,1._wp,0._wp]
		latt%T1=[1._wp,0._wp,0._wp]*10._wp
		latt%T2=[0._wp,1._wp,0._wp]*10._wp
		latt%bdc=[1._wp,-1._wp,0._wp]
		allocate(latt%rsb(1,3))
		latt%rsb(1,:)=[0._wp,0._wp,0._wp]
		call latt%gen_latt(size(t))
		call latt%gen_brizon(brizon)
		if(this_image()==1) then
			call check_lattice(101)
			write(*,*)"Total site number is: ",latt%Ns
		endif
		Ns=latt%Ns

!*********************************************************************
		allocate(Ham%var(-10:10),Hmf%var(-10:10),Hja%var(-10:10),mc%sphy%var(-10:10),mc%dphy%var(-10:10))
		!! cp
		!idx=Hmf%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,+2),c("i",1,-2)],n=2,sg=[1._wp,-1._wp],label="cp",is_var=.false.)
		!Hmf%var(idx)%bd=-1._wp
		!Hmf%var(idx)%val=-0._wp

		! d-wave sc
		idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2),c("j",1,+2),c("i",1,-1),c("i",1,+2),c("j",1,-1)],n=2,label="sc",is_var=.true.)
		do i=1,size(Hmf%var(idx)%bd)
			if(abs(latt%nb(1)%bd(i)%dr(1))>1e-6_wp) then
				Hmf%var(idx)%bd(i)=1._wp
			else
				Hmf%var(idx)%bd(i)=-1._wp
			endif
		enddo
		Hmf%var(idx)%val=1e-1_wp

		!! PDW
		!idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2),c("j",1,+2),c("i",1,-1),c("i",1,+2),c("j",1,-1)],n=2,label="sc")
		!do i=1,size(Hmf%var(idx)%bd)
			!Hmf%var(idx)%bd(i)=merge(1._wp,-1._wp,abs(latt%nb(1)%bd(i)%dr(1))>1e-6_wp)*cos(sum([2._wp*pi/8._wp,0._wp,0._wp]*(latt%nb(1)%bd(i)%r)))
			!!if(abs(latt%nb(1)%bd(i)%r(1)-3.5_wp)<1e-5_wp) then
				!!Hmf%var(idx)%bd(i)=0._wp
			!!elseif(abs(latt%nb(1)%bd(i)%r(1)+8-7.5_wp)<1e-5_wp) then
				!!Hmf%var(idx)%bd(i)=0._wp
			!!elseif(latt%nb(1)%bd(i)%r(1)>3.5_wp) then
				!!Hmf%var(idx)%bd(i)=-merge(1._wp,-1._wp,abs(latt%nb(1)%bd(i)%dr(1))>1e-6_wp)
			!!else
				!!Hmf%var(idx)%bd(i)=merge(1._wp,-1._wp,abs(latt%nb(1)%bd(i)%dr(1))>1e-6_wp)
			!!endif
		!enddo
		!Hmf%var(idx)%val=3._wp

		!! ddw
		!idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("j",1,+2),c("i",1,-2),c("i",1,+2),c("j",1,-2)],n=2,sg=[1._wp,-1._wp,-1._wp,1._wp])
		!do i=1,size(Hmf%var(idx)%bd)
			!if(abs(latt%nb(Hmf%var(idx)%nb)%bd(i)%dr(mod(nint(sum(latt%nb(0)%bd(latt%nb(Hmf%var(idx)%nb)%bd(i)%i(1))%r)),2)+1))>1e-6_wp) then
				!Hmf%var(idx)%bd(i)=img
			!else
				!Hmf%var(idx)%bd(i)=-img
			!endif
		!enddo
		!Hmf%var(idx)%val=1e-1_wp

		!! sdw
		!idx=Hmf%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,+2),c("i",1,-2)],n=2)
		!do i=1,size(Hmf%var(idx)%bd)
			!if(mod(nint(sum(latt%nb(0)%bd(latt%nb(0)%bd(i)%i(1))%r)),2)==0) then
				!Hmf%var(idx)%bd(i)=1._wp
			!else
				!Hmf%var(idx)%bd(i)=-1._wp
			!endif
		!enddo
		!Hmf%var(idx)%val=1.3e0_wp
		!ivar=idx

		!! bond order
		!do l=2,size(t)
			!idx=Hmf%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,-2)],n=2,sg=[+1._wp,+1._wp,-1._wp,-1._wp])
			!Hmf%var(idx)%bd=-1._wp
			!Hmf%var(idx)%val=-t(l)
		!enddo

		! hp
		do l=1,size(t)
			idx=Hmf%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,-2)],n=2,sg=[+1._wp,+1._wp,-1._wp,-1._wp],is_var=.false.)
			Hmf%var(idx)%bd=-1._wp
			Hmf%var(idx)%val=t(l)/t(1)
		enddo

!************************jastrow**************************************
		idx=Hja%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=4)
		Hja%var(idx)%bd=-1._wp
		Hja%var(idx)%val=0.9_wp

		!!!*******************Vzz term(exp(Siz*Sjz))********************* 
		!idx=Hja%add(nb=1,ca=[c("i",1,+1),c("i",1,-1),c("j",1,+1),c("j",1,-1),c("i",1,+1),c("i",1,-1),c("j",1,-2),c("j",1,+2),&
		!&c("i",1,-2),c("i",1,+2),c("j",1,+1),c("j",1,-1),c("i",1,-2),c("i",1,+2),c("j",1,-2),c("j",1,+2)],n=4,sg=[+1._wp,-1._wp,-1._wp,+1._wp],is_var=.true.)
		!Hja%var(idx)%bd=-1._wp
		!Hja%var(idx)%val=.2_wp
		!ivar=idx

!*************************Hamiltionian********************************
		do l=1,size(t)
			idx=Ham%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,-2),c("j",1,+2),c("j",1,-2),c("i",1,+2)],n=2,label="tp")
			Ham%var(idx)%bd=-1._wp
			Ham%var(idx)%val=t(l)
		enddo

		!! t-J interaction
		!idx=Ham%add(nb=1,ca=[c("i",1,+1),c("i",1,+2),c("j",1,-2),c("j",1,-1),c("j",1,+1),c("j",1,+2),c("i",1,-2),c("i",1,-1)],n=4,V=DJ/2._wp,label="J")
		!Ham%var(idx)%bd=1._wp
		!Ham%var(idx)%val=1._wp

		!idx=Ham%add(nb=1,ca=[c("i",1,+1),c("i",1,-1),c("j",1,+1),c("j",1,-1),c("i",1,-2),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,+1),c("i",1,-1),c("j",1,-2),c("j",1,+2),c("i",1,-2),c("i",1,+2),c("j",1,+1),c("j",1,-1)],n=4,sg=[V+DJ/4._wp,V+DJ/4._wp,V-DJ/4._wp,V-DJ/4._wp],label="4term")
		!Ham%var(idx)%bd=1._wp
		!Ham%var(idx)%val=1._wp

		! hubbard interaction
		idx=Ham%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=4)
		Ham%var(idx)%bd=1._wp
		Ham%var(idx)%val=U


!*************************static measurement**************************
		!idx=mc%sphy%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=2,sg=[1._wp,-1._wp],V=1._wp/Ns,label="1s")
		!do i=1,size(mc%sphy%var(idx)%bd)
			!if(mod(nint(sum(latt%nb(0)%bd(latt%nb(0)%bd(i)%i(1))%r)),2)==0) then
				!mc%sphy%var(idx)%bd(i)=1._wp
			!else
				!mc%sphy%var(idx)%bd(i)=-1._wp
			!endif
		!enddo

		!idx=mc%sphy%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=2,V=1._wp/Ns,label="1n")

		!idx=mc%sphy%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2)],n=2,V=1._wp/(Ns*Ns),label="2sc",extdat=[real(8)::ichar("r"),1e-7_wp])
		!do i=1,size(mc%sphy%var(idx)%bd)
			!if(abs(latt%nb(1)%bd(i)%dr(1))>1e-6_wp) then
				!mc%sphy%var(idx)%bd(i)=1._wp
			!else
				!mc%sphy%var(idx)%bd(i)=-1._wp
			!endif
		!enddo

		!idx=mc%sphy%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2)],n=2,V=1._wp/(sqrt(real(Ns))),label="2rsc",val=[real(8)::(i,i=1,size(latt%nb(1)%bd))],extdat=[real(8)::ichar("r"),1e-7_wp,ichar("f"),merge(1._wp,0._wp,abs(latt%nb(1)%bd(:)%r(1))<1e-6_wp)])

!************************dynamic measurement**************************
		!idx=mc%dphy%add(nb=0,ca=[c("i",1,+1),c("i",1,+2)],n=2,label="2s",extdat=[real(8)::ichar("q"),q(:,this_image(mc,1)),0._wp])

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



	!if(this_image()==1) then
		!call Hmf%band(80,[0._wp,0._wp,0._wp],brizon%Ta(1,:),100)
		!call Hmf%band(80,brizon%Ta(1,:),(brizon%Ta(1,:)+brizon%Ta(2,:))/2._wp,100)
		!call Hmf%band(80,(brizon%Ta(1,:)+brizon%Ta(2,:))/2._wp,[0._wp,0._wp,0._wp],100)
	!endif
	!stop


	mc%ne(1)=Ns/2
	if(is_ph) then
		mc%ne(2)=Ns-mc%ne(1)
	else
		mc%ne(2)=mc%ne(1)
	endif

	mc%step=nint(sqrt(real(Ns)))

	mc%samp=1024*16*8!*4*2

	call mc%init(.true.)
	call mc%do_var(100)

	do 
		!Hmf%var(ivar)%val=for_in([-10:10]*0.1_wp,id=1)
		Hja%var(ivar)%val=for_in([-10:10]*0.1_wp,id=1)
		!if(isnan(Hmf%var(ivar)%val(1))) then
		if(isnan(Hja%var(ivar)%val(1))) then
			exit
		endif
		mc%sg=2
		call mc%init(.true.)
		call mc%do_vmc()
		if(this_image()==1) then
			!write(*,*)Hmf%var(ivar)%val,mc%E,mc%err
			!write(81,*)Hmf%var(ivar)%val,mc%E,mc%err
			!write(*,*)Hja%var(ivar)%val,mc%E,mc%err
			!write(81,*)Hja%var(ivar)%val,mc%E,mc%err
			write(*,"(es12.4$)")Hmf%var(1:)%val(1),Hja%var(1:)%val(1)
			write(*,"(es14.6$)")mc%E
			write(*,"(es9.2$)")mc%err
			write(*,"(*(i3))")int(sign(1._wp,mc%g))
		endif
	enddo
	stop

	call mc%do_var(20)

end program
