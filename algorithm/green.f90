module green
	implicit none
	contains
	subroutine pade(u,z,omg,it,c)
		complex(8), parameter :: img=(0d0,1d0)
		real(8), parameter :: pi=3.14159265358979d0
		real(8) :: Tk,it,rtmp
		complex(8) :: u(:),c(:),p(size(u),size(u)),a(size(u)),DA(0:1),DB(0:1),z(:),ctmp,ctmp1,omg(:)
		integer :: n,m,i,j,k
		k=0
		n=size(u)
		m=size(omg)
		!z=(/(cmplx(0d0,pi*Tk*(2d0*(i-1)+1d0)),i=1,n)/)
		p(1,:)=u
		a(1)=u(1)
		do i=2,n
			do j=i,n
				p(i,j)=(p(i-1,i-1)-p(i-1,j))/((z(j)-z(i-1))*p(i-1,j))
			enddo
			a(i)=p(i,i)
		enddo
		do j=1,m
			DA(0)=0
			DA(1)=a(1)
			DB(0:1)=1d0
			do i=2,n
				ctmp=(omg(j)-z(i-1))*a(i)
				ctmp1=DA(1)+ctmp*DA(0)
				DA(0)=DA(1)
				DA(1)=ctmp1
				ctmp1=DB(1)+ctmp*DB(0)
				DB(0)=DB(1)
				DB(1)=ctmp1
				if ((abs(dreal(DA(1)))>1.d100).or.(abs(dimag(DA(1)))>1.d100)) then
					rtmp=1d0/abs(DA(1))
					DA(:)=DA(:)*rtmp
					DB(:)=DB(:)*rtmp
				endif
				if ((abs(dreal(DB(1)))>1.d100).or.(abs(dimag(DB(1)))>1.d100)) then
					rtmp=1d0/abs(DB(1))
					DA(:)=DA(:)*rtmp
					DB(:)=DB(:)*rtmp
				endif
			enddo
			c(j)=DA(1)/DB(1)
		enddo
	end subroutine
end
