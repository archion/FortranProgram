!Gram-Schmidt Orthonormalization Process
!A INPUT/OUTPUT, complex*16,dimension(n,m), target matrix
!n,m INPUT, integers, dimensions of target matrix

subroutine orth(A,n,m)
	implicit none
	integer,intent(IN) ::n,m
	complex*16,dimension(n,m),intent(INOUT) ::A
	
	integer ::i,j,k
	complex*16,dimension(n) ::b,c
	complex*16 ::m1,m2
	
	if(n<1) return
	if(m<1) return
	if(m>n) return
	
	do i=2,m,1
		b=A(1:n,i)
		do j=1,i-1,1
			m1=0.0D0
			m2=0.0D0
			do k=1,n
				m1=m1+conjg(A(k,j))*A(k,i)
				m2=m2+conjg(A(k,j))*A(k,j)
			enddo
			b=b-(m1/m2)*A(1:n,j)
		enddo
		A(1:n,i)=b
	enddo
	
	
	do i=1,m,1
		m1=0.0D0
		do j=1,n,1
			m1=m1+conjg(A(j,i))*A(j,i)
		enddo
		A(1:n,i)=A(1:n,i)/sqrt(m1)
	enddo
	return
end subroutine orth

subroutine orth1(A,n,m)
	implicit none
	integer,intent(IN) ::n,m
	complex*16,dimension(n,m),intent(inout) ::A
	
	integer ::i,j,k
	complex*16,dimension(n) ::b,c
	complex*16 ::m1
	
	if(n<1) return
	if(m<1) return
	if(m>n) return
	
	do i=1,m,1
		b=A(1:n,i)
		do j=1,i-1,1
			m1=0.0D0
			do k=1,n,1
				m1=m1+conjg(A(k,j))*A(k,i)
			enddo
			b=b-m1*A(1:n,j)
		enddo
		m1=0.0D0
		do k=1,n,1
			m1=m1+abs(b(k))**2
		enddo
		m1=sqrt(m1)
		do k=1,n,1
			b(k)=b(k)/m1
		enddo
		A(1:n,i)=b
	enddo
	return
end subroutine orth1

subroutine orth2(A,n,m,d)
	implicit none
	integer,intent(IN) ::n
	integer,intent(IN) ::m
	complex*16,dimension(n,m),intent(inout) ::A
	
	integer ::i,j,k
	complex*16 ::m1
	double precision,dimension(m),intent(OUT) ::d
	
	if(n<1) return
	if(m<1) return
	if(m>n) return
	
	d=0.0D0
	do i=1,m,1
	
		if(i>1)then
			do j=1,i-1,1
				m1=0.0D0
				do k=1,n,1
					m1=m1+conjg(A(k,j))*A(k,i)
				enddo
				A(1:n,i)=A(1:n,i)-m1*A(1:n,j)
			enddo
		endif
		
		m1=0.0D0
		
		do k=1,n,1
			m1=m1+conjg(A(k,i))*A(k,i)
		enddo
		
		m1=sqrt(m1)
		d(i)=abs(m1)
		A(1:n,i)=A(1:n,i)/m1
	enddo
	
	return
end subroutine orth2

subroutine orth3(A,n,m,d)
	implicit none
	integer,intent(IN) ::n,m
	complex*16,dimension(n,m),intent(INOUT) ::A
	
	integer ::i,j,k
	complex*16,dimension(n) ::b,c
	double precision,dimension(m) ::d
	complex*16 ::m1,m2
	
	if(n<1) return
	if(m<1) return
	if(m>n) return
	
	do i=2,m,1
		b=A(1:n,i)
		do j=1,i-1,1
			m1=0.0D0
			m2=0.0D0
			do k=1,n
				m1=m1+conjg(A(k,j))*A(k,i)
				m2=m2+conjg(A(k,j))*A(k,j)
			enddo
			b=b-(m1/m2)*A(1:n,j)
		enddo
		A(1:n,i)=b
	enddo
	
	d=0.0D0
	do i=1,m,1
		m1=0.0D0
		do j=1,n,1
			m1=m1+conjg(A(j,i))*A(j,i)
		enddo
		d(i)=abs(sqrt(m1))
		A(1:n,i)=A(1:n,i)/sqrt(m1)
	enddo
	return
end subroutine orth3