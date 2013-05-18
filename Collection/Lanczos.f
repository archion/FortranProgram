      implicit none
	
      double precision t,u,h,u1,u2,h3,w1,x1,y1,temp,z,e0,vt,v0,gs,vs
      integer i,j,l,m,p,q,k,o,v,v1,c,x,c1,c2,n,g,n1
      dimension h(36,36),v(36,4),v1(6,2),c(4,3),h3(36,36),u1(36,36)
      dimension w1(36),x1(36),y1(36),u2(36)
      dimension e0(36),vt(36),v0(36,36),vs(36)


      open(1,file="D:\FortranProgram\DATA\LanczosOutPut.txt")
      t=1
      u=6
      
!     basic wavefunction generation
      p=0
      l=0
11    l=l+1
      m=l
12    m=m+1
      p=p+1
      v1(p,1)=l
      v1(p,2)=m
      if(p.gt.6) goto 20
      if(m.lt.4) goto 12
      if(l.lt.3) goto 11
      
20    continue
      
      q=-1
      
21    q=q+1
      p=0
22    p=p+1
      
      k=p+6*q
      v(k,1)=v1(q+1,1)
      v(k,2)=v1(q+1,2)
      v(k,3)=v1(p,1)
      v(k,4)=v1(p,2)
      
      if(p.lt.6) goto 22
      if(q.lt.5) goto 21
      write(1,*) "Initial Vector Space:"
      do 30 q=1,36,1
        print*,v(q,1),v(q,2),v(q,3),v(q,4)
        write(1,*) v(q,1),v(q,2),v(q,3),v(q,4)
30    continue

!     evaluate the hamiltonian
      do 51 p=1,36
        do 51 q=1,36
          h(p,q)=0.0
        continue
51    continue

     
      do 61 m=0,2,1
        do 61 p=1,36,1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
          c(1,1)=v(p,m+1)+1
          c(1,2)=v(p,m+2)
          c(2,1)=v(p,m+1)-1
          c(2,2)=v(p,m+2)
          c(3,1)=v(p,m+1)
          c(3,2)=v(p,m+2)+1
          c(4,1)=v(p,m+1)
          c(4,2)=v(p,m+2)-1
          do 52 i=1,4,1
            do 53 j=1,2,1
              if(c(i,j).lt.1) c(i,j)=4
              if(c(i,j).gt.4) c(i,j)=1
53          continue
          
            if(c(i,1).lt.c(i,2))then
              c(i,3)=1
            else if(c(i,1).eq.c(i,2))then
              c(i,3)=0
            else if(c(i,1).gt.c(i,2))then
              x=c(i,1)
              c(i,1)=c(i,2)
              c(i,2)=x
              c(i,3)=-1
            endif         
52        continue
          do 61 l=1,4
            if(c(l,3).ne.0)then
              do 66 q=1,36,1
                if(v(q,m+1).eq.c(l,1))then
                  if(v(q,m+2).eq.c(l,2))then
                    h(p,q)=h(p,q)+c(l,3)*t
                  else
                    h(p,q)=h(p,q)
                  endif
                else
                  h(p,q)=h(p,q)
                endif
66            continue
            endif
          continue
        continue
61    continue

      do 70 p=1,36,1
        do 70 l=1,2,1
          do 70 m=3,4,1
            if(v(p,l).eq.v(p,m))then
              h(p,p)=h(p,p)+u
            else
              h(p,p)=h(p,p)
            endif
          continue
        continue
70    continue
      write(1,*) "Initail Hamiltonian:"
      write(1,67) ((h(i,j),j=1,36),i=1,36)
67    format(1x,6(6(f4.1," "),"|"),";")
!!!!!!!!!!!!!!start lanczos!!!!!!!!!!!!!!!!!!!!!!!

      do 101 i=1,36
        do 101 j=1,36
          h3(i,j)=0
          u1(i,j)=0
        continue
        y1(i)=0.0
101   continue
      n1=1
      o=0
      g=1
      n=0
      temp=0.0
!!!! initiate random vector !!!!!!!!!!!!!!!!
      
110   z=(o+n1)
      call nrnds(36,u2,z)
      if(n1.gt.1)then
        do 111 j=1,n1-1
          temp=0.0
          do 112 i=1,36
            temp=temp+u1(i,j)*u2(i)
112       continue
          do 113 i=1,36
            u2(i)=u2(i)-temp*u1(i,j)
113       continue
111     continue
      endif
      
      l=0
      do 114 i=1,36
        if(u2(i).ge.0.00001) l=l+1
114   continue
      if(l.eq.0)then
        o=o+1
        goto 110
      endif
      o=0
      temp=0.0
      do 115 i=1,36
        temp=temp+u2(i)*u2(i)
115   continue
      do 116 i=1,36
        u1(i,n1)=u2(i)/sqrt(temp)
116   continue
      
120   n=n+1
!!!!! get beta n-1 !!!!!!!!!      
      if(n.gt.n1) h3(n,n-1)=h3(n-1,n)
!!!! get alpha n !!!!!!!!!!!!
      temp=0.0
      do 121 i=1,36
        y1(i)=0.0
        do 122 j=1,36
          y1(i)=y1(i)+h(i,j)*u1(j,n)
122     continue
        temp=temp+u1(i,n)*y1(i)
121   continue
      h3(n,n)=temp
      e0(n)=h3(n,n)

      temp=0.0
!!!!!!!! get beta n !!!!!!!!!!!
      do 123 i=1,36
        w1(i)=y1(i)-h3(n,n)*u1(i,n)
        if(n.gt.n1) w1(i)=w1(i)-h3(n,n-1)*u1(i,n-1)
        x1(i)=w1(i)
123   continue

      do 127 l=1,n
        do 127 i=1,36
          temp=0.0
          do 129 j=1,36
            temp=temp+u1(j,l)*w1(j)
129       continue
          x1(i)=x1(i)-temp*u1(i,l)
        continue
127   continue

      temp=0.0
      do 130 i=1,36
        temp=temp+x1(i)*x1(i)
130   continue
      write(*,132) n,h3(n,n),n,sqrt(temp)
132   format(1x,"a(",I2,")=",f10.6,"b(",I2,")=",f10.6)
      if(temp.le.0.00001)then
        print*,"!"
        goto 300
      endif
      if(n.eq.36)then
        print*,"B36=",sqrt(temp)
        goto 300
      endif
      h3(n,n+1)=sqrt(temp)
      vt(n)=h3(n,n+1)
!!!!!!get v n+1 !!!!!!!!!
      do 133 i=1,36
        u1(i,n+1)=x1(i)/h3(n,n+1)
133   continue
      goto 120
300   if(n.lt.36)then
        vt(n)=0.0
        print*,"degeneration at n=",n
        n1=n+1
        g=g+1
        goto 110
      endif
      print*,"G=",g
      write(1,*) "Lanczos Matrix:"
      write(1,311) ((h3(i,j),j=1,n),i=1,n)
      write(1,*) "Lanczos Vectors:"
      write(1,311) ((u1(i,j),i=1,36),j=1,n)
      write(1,*) "Degeneration=",g
311   format(1x,36(f5.2," "))

!!!!!!!!! eigenvalues !!!!!!!!!!!!
      do 401 i=1,36
        do 401 j=1,36
          if(i.eq.j)then
            v0(i,j)=1.0
          else
            v0(i,j)=0.0
          endif
        continue
401   continue
      call csstq(36,e0,vt,v0,1e-5,l)
      if(l.ne.0)then
        gs=e0(1)
        do 410 i=1,36
          if(e0(i).lt.gs)then
            gs=e0(i)
            m=i
          endif
410     continue
        do 412 i=1,36
          temp=0.0
          do 413 j=1,36
            temp=temp+u1(i,j)*v0(j,m)
413       continue
          vs(i)=temp
412     continue
        
        write(1,*) "Eigenvalues:"
        write(1,411) (e0(i),i=1,36)
411     format(1x,f10.6)
        write(1,*) "Transformation:"
        write(1,311) ((v0(i,j),j=1,36),i=1,36)
        write(1,*) "Grandstate energy=",gs
        print*,"Grandstate energy=",gs
        write(1,*) "Eigenstate:"
        write(1,411) (vs(i),i=1,36)
      else
        print*,"Error!"
      endif
      close(1)
      end

      
!!!!!!!!!!!!! random number generator!!!!!!      
      subroutine nrnds(n,ak,rk)
      integer m,n,i
      dimension ak(n)
      double precision ak,rk,sk,uk,vk
      sk=65563.0
      uk=2053.0
      vk=13849.0
      do 200 i=1,n
        rk=uk*rk+vk
        m=rk/sk
        rk=rk-m*sk
        ak(i)=rk/sk
200   continue
      return
      end
      
!!!!!! QR !!!!!!!!!!!!!!!!!!

	SUBROUTINE CSSTQ(N,B,C,Q,EPS,L)
	DIMENSION B(N),C(N),Q(N,N)
	DOUBLE PRECISION B,C,Q,D,H,P,R,F,E,S,G
	C(N)=0.0
	D=0.0
	F=0.0
        DO 550 J=1,N
	  IT=0
	  H=EPS*(ABS(B(J))+ABS(C(J)))
	  IF (H.GT.D) D=H
	  M=J-1
510       M=M+1
	  IF (M.LE.N) THEN
            IF (ABS(C(M)).GT.D) GOTO 510
	  END IF
	  IF (M.NE.J) THEN
515         IF (IT.EQ.60) THEN
	      L=0
              WRITE(*,518)
518           FORMAT(1X,'  FAIL')
	      RETURN
	    END IF
	    IT=IT+1
	    G=B(J)
	    P=(B(J+1)-G)/(2.0*C(J))
	    R=SQRT(P*P+1.0)
	    IF (P.GE.0.0) THEN
	      B(J)=C(J)/(P+R)
	    ELSE
	      B(J)=C(J)/(P-R)
	    END IF
	    H=G-B(J)
            DO 520 I=J+1,N
520           B(I)=B(I)-H
	    F=F+H
	    P=B(M)
	    E=1.0
	    S=0.0
            DO 540 I=M-1,J,-1
	      G=E*C(I)
	      H=E*P
	      IF (ABS(P).GE.ABS(C(I))) THEN
	        E=C(I)/P
	        R=SQRT(E*E+1.0)
	        C(I+1)=S*P*R
	        S=E/R
	        E=1.0/R
	      ELSE
	        E=P/C(I)
	        R=SQRT(E*E+1.0)
	        C(I+1)=S*C(I)*R
	        S=1.0/R
	        E=E/R
	      END IF
	      P=E*B(I)-S*G
	      B(I+1)=H+S*(E*G+S*B(I))
              DO 530 K=1,N
	        H=Q(K,I+1)
	        Q(K,I+1)=S*Q(K,I)+E*H
	        Q(K,I)=E*Q(K,I)-S*H
530           CONTINUE
540         CONTINUE
	    C(J)=S*P
	    B(J)=E*P
            IF (ABS(C(J)).GT.D) GOTO 515
	  END IF
	  B(J)=B(J)+F
550     CONTINUE
        DO 580 I=1,N
	  K=I
	  P=B(I)
	  IF (I+1.LE.N) THEN
	    J=I
560         J=J+1
	    IF (J.LE.N) THEN
	      IF (B(J).LE.P) THEN
	        K=J
	        P=B(J)
                GOTO 560
	      END IF
	    END IF
	  END IF
	  IF (K.NE.I) THEN
	    B(K)=B(I)
	    B(I)=P
            DO 570 J=1,N
	      P=Q(J,I)
	      Q(J,I)=Q(J,K)
	      Q(J,K)=P
570         CONTINUE
	  END IF
580     CONTINUE
	L=1
	RETURN
	END

