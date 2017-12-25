module coarray
	integer, parameter :: ica1=1
	integer :: ica2
end module
include "../lib/hamilton_final_1.f90"
include "vmc_utility.f90"
module global
	use vmc_utility
	use ifport, ifort_qsort => qsort
	use omp_lib
	use mkl_service
	implicit none
	real(8) :: t(2)=[1d0,-0.3d0]
	real(8), parameter :: DJ=0.3d0,V=0d0,U=2.4d0
	type(t_mc) :: mc[ica1,*]
contains
	subroutine initial()
		integer :: i,l,idx,sb,tp
		real(8) :: q(2,ica1)=reshape([&
			pi,0d0&
			],[2,ica1])
!***************************parameter setting*************************
		is_project=.true.
		is_ph=.true.
		dx=0.3d0
		!call init_random_seed()

!****************************lattice**********************************
		latt%is_all=.true.
		latt%a1=[1d0,0d0,0d0]
		latt%a2=[0d0,1d0,0d0]
		latt%c1=latt%a1
		latt%c2=latt%a2
		latt%c1=[1d0,1d0,0d0]
		latt%c2=[-1d0,1d0,0d0]
		!latt%c1=[4d0,1d0,0d0]
		!latt%c2=[0d0,2d0,0d0]
		latt%c1=[8d0,0d0,0d0]
		latt%c2=[0d0,1d0,0d0]
		latt%T1=[1d0,0d0,0d0]*8d0
		latt%T2=[0d0,1d0,0d0]*8d0
		latt%bdc=[1d0,1d0,0d0]
		allocate(latt%rsb(1,3))
		latt%rsb(1,:)=[0d0,0d0,0d0]
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
		! cp
		idx=Hmf%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,+2),c("i",1,-2)],n=2,sg=[1d0,-1d0],label="cp")
		Hmf%var(idx)%bd=-1d0
		Hmf%var(idx)%val=-0.361803398874990d0

		! d-wave sc
		idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2),c("j",1,+2),c("i",1,-1),c("i",1,+2),c("j",1,-1)],n=2,label="sc")
		do i=1,size(Hmf%var(idx)%bd)
			if(abs(latt%nb(1)%bd(i)%dr(1))>1d-6) then
				Hmf%var(idx)%bd(i)=1d0
			else
				Hmf%var(idx)%bd(i)=-1d0
			endif
		enddo
		Hmf%var(idx)%val=2.92d0

		!! PDW
		!idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2),c("j",1,+2),c("i",1,-1),c("i",1,+2),c("j",1,-1)],n=2,label="sc")
		!do i=1,size(Hmf%var(idx)%bd)
			!Hmf%var(idx)%bd(i)=merge(1d0,-1d0,abs(latt%nb(1)%bd(i)%dr(1))>1d-6)*cos(sum([2d0*pi/8d0,0d0,0d0]*(latt%nb(1)%bd(i)%r)))
			!!if(abs(latt%nb(1)%bd(i)%r(1)-3.5d0)<1d-5) then
				!!Hmf%var(idx)%bd(i)=0d0
			!!elseif(abs(latt%nb(1)%bd(i)%r(1)+8-7.5d0)<1d-5) then
				!!Hmf%var(idx)%bd(i)=0d0
			!!elseif(latt%nb(1)%bd(i)%r(1)>3.5d0) then
				!!Hmf%var(idx)%bd(i)=-merge(1d0,-1d0,abs(latt%nb(1)%bd(i)%dr(1))>1d-6)
			!!else
				!!Hmf%var(idx)%bd(i)=merge(1d0,-1d0,abs(latt%nb(1)%bd(i)%dr(1))>1d-6)
			!!endif
		!enddo
		!Hmf%var(idx)%val=3d0

		!! ddw
		!idx=Hmf%add(nb=1,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("j",1,+2),c("i",1,-2),c("i",1,+2),c("j",1,-2)],n=2,sg=[1d0,-1d0,-1d0,1d0])
		!do i=1,size(Hmf%var(idx)%bd)
			!if(abs(latt%nb(Hmf%var(idx)%nb)%bd(i)%dr(mod(nint(sum(latt%nb(0)%bd(latt%nb(Hmf%var(idx)%nb)%bd(i)%i(1))%r)),2)+1))>1d-6) then
				!Hmf%var(idx)%bd(i)=img
			!else
				!Hmf%var(idx)%bd(i)=-img
			!endif
		!enddo
		!Hmf%var(idx)%val=1d-1

		!! sdw
		!idx=Hmf%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,+2),c("i",1,-2)],n=2)
		!do i=1,size(Hmf%var(idx)%bd)
			!if(mod(nint(sum(latt%nb(0)%bd(latt%nb(0)%bd(i)%i(1))%r)),2)==0) then
				!Hmf%var(idx)%bd(i)=1d0
			!else
				!Hmf%var(idx)%bd(i)=-1d0
			!endif
		!enddo
		!Hmf%var(idx)%val=1d-1

		! bond order
		do l=2,size(t)
			idx=Hmf%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,-2)],n=2,sg=[+1d0,+1d0,-1d0,-1d0])
			Hmf%var(idx)%bd=-1d0
			Hmf%var(idx)%val=-t(l)
		enddo

		! hp
		do l=1,size(t)
			idx=Hmf%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,-2)],n=2,sg=[+1d0,+1d0,-1d0,-1d0],is_var=.false.)
			Hmf%var(idx)%bd=-1d0
			Hmf%var(idx)%val=t(l)
		enddo

!************************jastrow**************************************
		!idx=Hja%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=4)
		!Hja%var(idx)%bd=-1d0
		!Hja%var(idx)%val=0d0

!*************************Hamiltionian********************************
		do l=1,size(t)
			idx=Ham%add(nb=l,ca=[c("i",1,+1),c("j",1,-1),c("j",1,+1),c("i",1,-1),c("i",1,-2),c("j",1,+2),c("j",1,-2),c("i",1,+2)],n=2,label="tp")
			Ham%var(idx)%bd=-1d0
			Ham%var(idx)%val=t(l)
		enddo

		! t-J interaction
		idx=Ham%add(nb=1,ca=[c("i",1,+1),c("i",1,+2),c("j",1,-2),c("j",1,-1),c("j",1,+1),c("j",1,+2),c("i",1,-2),c("i",1,-1)],n=4,V=DJ/2d0,label="J")
		Ham%var(idx)%bd=1d0
		Ham%var(idx)%val=1d0

		idx=Ham%add(nb=1,ca=[c("i",1,+1),c("i",1,-1),c("j",1,+1),c("j",1,-1),c("i",1,-2),c("i",1,+2),c("j",1,-2),c("j",1,+2),c("i",1,+1),c("i",1,-1),c("j",1,-2),c("j",1,+2),c("i",1,-2),c("i",1,+2),c("j",1,+1),c("j",1,-1)],n=4,sg=[V+DJ/4d0,V+DJ/4d0,V-DJ/4d0,V-DJ/4d0],label="4term")
		Ham%var(idx)%bd=1d0
		Ham%var(idx)%val=1d0

		!! hubbard interaction
		!idx=Ham%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=4)
		!Ham%var(idx)%bd=1d0
		!Ham%var(idx)%val=U


!*************************static measurement**************************
		!idx=mc%sphy%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=2,sg=[1d0,-1d0],V=1d0/Ns,label="1s")
		!do i=1,size(mc%sphy%var(idx)%bd)
			!if(mod(nint(sum(latt%nb(0)%bd(latt%nb(0)%bd(i)%i(1))%r)),2)==0) then
				!mc%sphy%var(idx)%bd(i)=1d0
			!else
				!mc%sphy%var(idx)%bd(i)=-1d0
			!endif
		!enddo

		!idx=mc%sphy%add(nb=0,ca=[c("i",1,+1),c("i",1,-1),c("i",1,-2),c("i",1,+2)],n=2,V=1d0/Ns,label="1n")

		!idx=mc%sphy%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2)],n=2,V=1d0/(Ns*Ns),label="2sc",extdat=[real(8)::ichar("r"),1d-7])
		!do i=1,size(mc%sphy%var(idx)%bd)
			!if(abs(latt%nb(1)%bd(i)%dr(1))>1d-6) then
				!mc%sphy%var(idx)%bd(i)=1d0
			!else
				!mc%sphy%var(idx)%bd(i)=-1d0
			!endif
		!enddo

		!idx=mc%sphy%add(nb=1,ca=[c("i",1,+1),c("j",1,-2),c("j",1,+1),c("i",1,-2)],n=2,V=1d0/(sqrt(real(Ns))),label="2rsc",val=[real(8)::(i,i=1,size(latt%nb(1)%bd))],extdat=[real(8)::ichar("r"),1d-7,ichar("f"),merge(1d0,0d0,abs(latt%nb(1)%bd(:)%r(1))<1d-6)])

!************************dynamic measurement**************************
		idx=mc%dphy%add(nb=0,ca=[c("i",1,+1),c("i",1,+2)],n=2,label="2s",extdat=[real(8)::ichar("q"),q(:,this_image(mc,1)),0d0])

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
	integer :: i,j,l,idx
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
	otime=0d0



	!if(this_image()==1) then
		!call Hmf%band(80,[0d0,0d0,0d0],brizon%Ta(1,:),100)
		!call Hmf%band(80,brizon%Ta(1,:),(brizon%Ta(1,:)+brizon%Ta(2,:))/2d0,100)
		!call Hmf%band(80,(brizon%Ta(1,:)+brizon%Ta(2,:))/2d0,[0d0,0d0,0d0],100)
	!endif
	!stop


	mc%step=Ns
	mc%delay=2
	mc%ne(1)=Ns/2-6
	if(is_ph) then
		mc%ne(2)=Ns-mc%ne(1)
	else
		mc%ne(2)=mc%ne(1)
	endif
	!Hmf%var(1:vn)%val(1)=[-6.2487E-01,2.5764E-01,1.0392E-01,-3.6208E-02] ! dsc+mu+t'+SDW E=-3.2721E-01
	mc%samp=1024*8*16*2
	mc%hot=1024*8*2
	mc%step=nint(sqrt(real(Ns)))
	!call mc%do_var(10)
	!!stop
	mc%num=this_image(mc,1)

	mc%sg=1
	!mc%samp=1024*8*8*16
	mc%samp=1024*8*8*4
	mc%hot=1024*8*8
	mc%step=nint(sqrt(real(Ns)))*2

	critical
	if(this_image(mc,2)==1) then
		do idx=1,Hmf%rg(2)
			if(Hmf%var(idx)%label=="sc") then
				do i=1,size(latt%nb(1)%bd)
					write(30,"(*(es12.4))")latt%nb(1)%bd(i)%r,latt%nb(1)%bd(i)%dr,Hmf%var(idx)%val(Hmf%var(idx)%bd2v(i))*Hmf%var(idx)%bd(i)
				enddo
				write(30,"(x/)")
			endif
		enddo
	endif
	endcritical

	call mc%init(.true.)
	call mc%do_vmc()
	critical
	if(this_image(mc,2)==1) then
		write(*,"(*(es12.4))")mc%dphy%var(1:mc%dphy%rg(2))%extdat(1)/pi,mc%dphy%var(1:mc%dphy%rg(2))%extdat(2)/pi,Hmf%var(1:)%val(1),mc%E,mc%err,mc%dphy%var(1:mc%dphy%rg(2))%val(1)
		do idx=1,mc%sphy%rg(2)
			if(mc%sphy%var(idx)%label=="2rsc") then
				do i=1,size(latt%nb(1)%bd)
					write(30,"(*(es12.4))")latt%nb(1)%bd(i)%r,latt%nb(1)%bd(i)%dr,mc%sphy%var(idx)%val(mc%sphy%var(idx)%bd2v(i))*mc%sphy%var(idx)%bd(i)
				enddo
			else
				write(*,"(es12.4$)")mc%sphy%var(idx)%val(:)
			endif
		enddo
		error stop
	endif
	endcritical

	if(this_image()==1) then
		do l=2,ica1
			mc[1,1]%E=mc[1,1]%E+mc[l,1]%E
		enddo
		mc[1,1]%E=mc[1,1]%E*1d0/ica1
		sync images(*)
	else
		sync images(1)
		mc%E=mc[1,1]%E
	endif
	sync all

	mc%sg=3
	mc%ne=mc%ne+1
	mc%hot=1024*8*8
	mc%samp=1024*8*8*8 !dsc 16x16
	mc%step=nint(sqrt(real(Ns)))*2
	call mc%init(.true.)
	call mc%do_vmc()

	if(this_image()==1) then
		write(*,*)"finished, exporting data....",mc%samp
		do l=1,ica1
			write(50)-1,size(mc[l,1]%psi0),Ns,mc[l,1]%dphy%var(1)%extdat,mc[l,1]%E,mc[l,1]%dphy%var(1)%val(1),mc[l,1]%Emf,mc[l,1]%psi0,mc[l,1]%Ok2,mc[l,1]%Ek2,mc[l,1]%nn(1:,:)
			!call mc%spect(70,0.02d0,[-10d0,10d0],5000)
			!mc[l,1]%Ek2=transpose(conjg(mc[l,1]%Ek2))
			!call mc%spect(70,0.02d0,[-10d0,10d0],5000)
		enddo
		rewind(50)
		read(50)i
		deallocate(brizon%k)
		allocate(brizon%k(i,3))
		read(50)brizon%k
		read(50)brizon%nk
		do l=1,ica1
			read(50)i,j,Ns
			if(allocated(mc%psi0)) deallocate(mc%psi0,mc%Ok2,mc%Ek2,mc%Emf,mc%nn)
			allocate(mc%psi0(j),mc%Ok2(j,j),mc%Ek2(j,j),mc%nn(j,2),mc%Emf(Ns*2))
			read(50)mc%dphy%var(1)%extdat,mc%E,mc%dphy%var(1)%val(1),mc%Emf,mc%psi0,mc%Ok2,mc%Ek2,mc%nn
			write(*,*)mc%dphy%var(1)%extdat,mc%E,Ns,i,j
			call mc%spect(70,0.02d0,[-10d0,10d0],5000)
			!call mc%spect(70,0.1d0,[-10d0,10d0],5000)
			mc%Ek2=transpose(conjg(mc%Ek2))
			call mc%spect(70,0.02d0,[-10d0,10d0],5000)
			!call mc%spect(70,0.1d0,[-10d0,10d0],5000)
			mc%Ek2=0.5d0*(mc%Ek2+transpose(conjg(mc%Ek2)))
			call mc%spect(70,0.02d0,[-10d0,10d0],5000)
			!call mc%spect(70,0.1d0,[-10d0,10d0],5000)
			deallocate(mc%psi0,mc%Ok2,mc%Ek2,mc%Emf,mc%nn)
		enddo
	endif
end program
