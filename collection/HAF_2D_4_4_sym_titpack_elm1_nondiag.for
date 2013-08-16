c************* main program: sym_elm1_nondiag **************
c*      2D Heisenberg antiferromagnet, 6*6 spins,          *
c*   H_{bond }=J*\sum_{i,j}[ ( S^{+}_{i}S^{-}_{j}          *
c*                           + S^{-}_{i}S^{+}_{j} ) / 2    *
c*                   + zrtio * S^{z}_{i}S^{z}_{j} ] .      *
c*          Translational symmetry is used .               *
c*               Eigenvalues by lnc1                       *
c*               Eigenvector by inv1                       *
c*      Precision  check and correlation functions         *
c*      Degeneracy check by various initial vectors        *
c***********************************************************
      parameter (n=16,nx=4,ny=4,ibond=2*n)
c     parameter (n=36,nx=6,ny=6,ibond=2*n)
      parameter (idim0=20000)      ! idim0 = estimated idim
c     parameter (idim0=15900000)
      parameter (nbond=1)
      implicit real*8 (a-h,o-z)
      dimension E(4)
c     dimension diag(idim0),nondiag(idim0,ibond)
      real*8,  allocatable, dimension(:  ) ::    diag
      integer, allocatable, dimension(:,:) :: nondiag
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)     ! number of member states
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
c
      ! Lanczos/CG vectors and wave functions
      complex*16 wk(idim0,4)
      complex*16  x(idim0)
c
      ! correlation functions
      dimension sxx(-nx:nx,-ny:ny),szz(-nx:nx,-ny:ny)
c
c     1-----2-----3-----4-----5-----6-----1
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     7-----8-----9----10----11----12-----7
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     13---14----15----16----17----18----13
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     19---20----21----22----23----24----19
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     25---26----27----28----29----30----25
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     31---32----33----34----35----36----31
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     1-----2-----3-----4-----5-----6-----1
c
      data ipair/1, 2,    2, 3,    3, 4,    4, 1,
     &           5, 6,    6, 7,    7, 8,    8, 5,
     &           9,10,   10,11,   11,12,   12, 9,
     &          13,14,   14,15,   15,16,   16,13,
     &           1, 5,    5, 9,    9,13,   13, 1,
     &           2, 6,    6,10,   10,14,   14, 2,
     &           3, 7,    7,11,   11,15,   15, 3,
     &           4, 8,    8,12,   12,16,   16, 4/

c      data ipair/
c     &  1, 2,    2, 3,    3, 4,    4, 5,    5, 6,    6, 1,
c     &  7, 8,    8, 9,    9,10,   10,11,   11,12,   12, 7,
c     & 13,14,   14,15,   15,16,   16,17,   17,18,   18,13,
c     & 19,20,   20,21,   21,22,   22,23,   23,24,   24,19,
c     & 25,26,   26,27,   27,28,   28,29,   29,30,   30,25,
c     & 31,32,   32,33,   33,34,   34,35,   35,36,   36,31,
c     &
c     &  1, 7,    7,13,   13,19,   19,25,   25,31,   31, 1,
c     &  2, 8,    8,14,   14,20,   20,26,   26,32,   32, 2,
c     &  3, 9,    9,15,   15,21,   21,27,   27,33,   33, 3,
c     &  4,10,   10,16,   16,22,   22,28,   28,34,   34, 4,
c     &  5,11,   11,17,   17,23,   23,29,   29,35,   35, 5,
c     &  6,12,   12,18,   18,24,   24,30,   30,36,   36, 6/
      data bondwt/ibond*1.d0/ ! -0.5d0 to 1.0d0 , see 'mltply' !
      data  zrtio/ibond*1.d0/
c
      print *
      print *, '%=========================================%'
      print *, '| 2D HAF Square Lattice, 4*4 Spins, J=1.0 |'
      print *, '%=========================================%'
      print *
c
c*** Allocate dynamic arrays as large as possible
      allocate(   diag(idim0)      )
      allocate(nondiag(idim0,ibond))
c
c     +---------------------------+
c     | Total Momentum = (Kx, Ky) | for 2D Square Lattice
c     +---------------------------+
c
c
      ! spin inversion quantum number: 0, 1
      do 999 nsi = 0, 1
c
      ! x-axis   mirror reflection quantum number: 0, 1
      do 999 nmx = 0, 1
c
      ! y-axis   mirror reflection quantum number: 0, 1
c     do 999 nmy = 0, 0
      do 999 nmy = nmx, nmx
c
      ! diagonal mirror reflection quantum number: 0, 1
      ! compatible with nmx=nmy only
      do 999 nmd = 0, 1
c
      ! translation quantum numbers,
      ! only (Kx,Ky)=(0,0) (pi,pi) (0,pi) (pi,0)
      ! are compatible with SI,Mx,My.
      ! only (Kx,Ky)=(0,0) (pi,pi)
      ! are compatible with SI,Mx,My,Md.
      do 999 nkx = 0, nx/2, 2
c     do 999 nky = 0, ny/2, 1
      do 999 nky = nkx, nkx

         print *
         print *, '+---------------------------------------------+'
         print 200, nsi,nmx,nmy,nmd,2*nkx,nx,2*nky,ny
 200     format( ' | Quantum numbers ='i2','i2 ','i2 ','i2 
     &              ', (',i1,'/',i1,', 'i1,'/',i1,')*pi |')
         print *, '+---------------------------------------------+'
         print *
c
c*** Number of eigenvectors to calculate by 'lncv1'
         nvec=1
c
c*** Configurations with the specified sz 
         call szdysym(n,nx,ny,idim,idim0,
     &                nsi,nmx,nmy,nmd,nkx,nky,
     &                0.0d0,
     &                list1,list2,nmem)
         print *,'[Subspace dimension]'
         print *, idim
         print *
c
c         call cpu_time(timetmp)
c         cputime1=timetmp
c         print *,'Basis coding done !'
c         print *,'Time for coding by szdysym : ',
c     &            cputime1
c         print *,'Time elapsed : ',cputime1
c         print *
c
c*** Initial vector of the Lanczos and inverse iteration
         iv=idim/3
c
c*** Elements of Hamiltonian
         call elm1(n,nx,ny,idim,idim0,
     &             nsi,nmx,nmy,nmd,nkx,nky,
     &             ipair,bondwt,zrtio,ibond,
     &             diag,nondiag,
     &             list1,list2,nmem)
c
c         call cpu_time(timetmp)
c         cputime2=timetmp
c         print *,'Hamiltonian elements done !'
c         print *,'Time for H elements by elm1 : ',
c     &            cputime2-cputime1
c         print *,'Time elapsed : ',cputime2
c         print *
c 
c*** Eigenvalues
         call lnc1(n,nx,ny,idim,idim0,
     &             nsi,nmx,nmy,nmd,nkx,nky,
     &             ipair,bondwt,zrtio,ibond,
     &             nvec,iv,E,itr,wk,idim0,
     &             diag,nondiag,
     &             list1,list2,nmem)
         print 100,E,itr
 100     format(/' [Eigenvalues]  '/2x,4f14.8
     &          //' [Iteration number]'/i8)
         print *
c
c         call cpu_time(timetmp)
c         cputime3=timetmp
c         print *,'Lanczos diagonalization done !'
c         print *,'Time for diagonalization by lnc1 : ',
c     &            cputime3-cputime2
c         print *,'Time elapsed : ',cputime3
c         print *
c
c*** Ground-state eigenvector
c- You may use lncv1 for faster speed and smaller working area -
c----- Note: wk(idim,2) for lncv1, and wk(idim,4) for lnv1 -----
        call lncv1(n,nx,ny,idim,idim0,
     &             nsi,nmx,nmy,nmd,nkx,nky,
     &             ipair,bondwt,zrtio,ibond,
     &             nvec,iv,x,itr,wk,idim0,
     &             diag,nondiag,
     &             list1,list2,nmem)
c- You may alternatively use inv1 to obtain higher precision -
c         call inv1(n,nx,ny,idim,idim0,
c     &             nsi,nmx,nmy,nmd,nkx,nky,
c     &             ipair,bondwt,zrtio,ibond,
c     &             E(1),iv,x,wk,idim0,
c     &             diag,nondiag,
c     &             list1,list2,nmem)
c         print *,'[Eigenvector components (selected)]'
c         print 120,(x(j),j=13,idim,idim/20)
c 120     format('    (',1d18.9,',',1d18.9,')')
c
c         call cpu_time(timetmp)
c         cputime4=timetmp
c         print *,'GS wave functions done !'
c         print *,'Time for GS vector by lncv1 : ',
c     &            cputime4-cputime3
c         print *,'Time elapsed : ',cputime4
c         print *
c
c*** Precision check and correlation functions
c         call check1(n,nx,ny,idim,idim0,
c     &               nsi,nmx,nmy,nmd,nkx,nky,
c     &               ipair,bondwt,zrtio,ibond,
c     &               x,wk,Hexpec,
c     &               diag,nondiag,
c     &               list1,list2,nmem)
c
c         call cpu_time(timetmp)
c         cputime5=timetmp
c         print *,'Precision check done !'
c         print *,'Time for precision check by check1 : ',
c     &            cputime5-cputime4
c         print *,'Time elapsed : ',cputime5
c         print *
c
c*** Spin-spin correlation functions
c         print *,'[Spin-spin correlation functions]'
c          call xcorr(n,nx,ny,idim,idim0,
c     &               nsi,nmx,nmy,nmd,nkx,nky,
c     &               x,sxx,
c     &               list1,list2,nmem)
c          call zcorr(n,nx,ny,idim,idim0,
c     &               nsi,nmx,nmy,nmd,nkx,nky,
c     &               x,szz,list1)
c          do ix=0,nx
c          do iy=0,ny
c          print 130,ix,iy,sxx(ix,iy),szz(ix,iy)
c 130      format('  ix,iy=',2i3,'    sxx=',d18.10,'    szz=',d18.10)
c          end do
c          end do
c
c*** Total spin of ground state
          call stot(n,nx,ny,idim,idim0,
     &              nsi,nmx,nmy,nmd,nkx,nky,
     &              x,stot2,
     &              list1,list2,nmem)
          print 160,stot2
 160      format(/' [Total spin of ground state]'/
     &            '    S(S+1)=',f14.8/)
c
c         call cpu_time(timetmp)
c         cputime6=timetmp
c         print *,'Correlation functions done !'
c         print *,'Time for correlation functions : ',
c     &            cputime6-cputime5
c         print *,'Time elapsed : ',cputime6
c         print *
c
 999  continue

c*** Deallocate large dynamic arrays
      deallocate(   diag)
      deallocate(nondiag)

      end
c
c***********************************************************
c*                                                         *
c*               Basis Coding and Decoding                 *
c*         for Spin-1/2 Systems with Symmetries            *
c*                                                         *
c*               by Yi-Fei Wang, Dec. 2009                 *
c*                                                         *
c*       Ref: H. Nishimori,    TITPACK Ver 2. (1991)       *
c*       Ref: Stephan W. Haas,     PhD Thesis (1995)       *
c*       Ref: Andreas M. Lauchli,  PhD Thesis (2002)       *
c*                                                         *
c***********************************************************
c
c    Bit-Wise Operations on Integers (Fortran 90):
c
c    i = BTEST(i, ipos)       Returns .TRUE. if bit at ipos is 1
c    i = IAND(i, j)           Logical AND
c    i = IBCLR(i, ipos)       Clears one bit at ipos to 0
c    i = IBITS(i, ipos, len)  Extracts a sequence of bits
c    i = IBSET(i, ipos)       Set one bit of i at ipos to 1
c    i = IEOR(i, j)           Exclusive OR
c    i = IOR(i, j)            Inclusive OR
c    i = ISHFT(i, j)          Logical  shift left (right if j<0)
c    i = ISHFTC(i, j, size)   Circular shift left (right if j<0)
c    i = NOT(i)               Logical complement
c
c    The bits are numbered from 0 to 31/63, from right to left.
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c
c============= Geometry-Dependent Subroutines ==============
c
c****** generates members of a momentum symmetry ***********
c****** class for an input state  'istate', if ifind=0 *****
c-------------- 2D square lattice geometry -----------------
c-------------- with mirror reflection and -----------------
c-------------- spin inversion symmetries ------------------
c****** find the representative state and symmetry *********
c****** operations relating it to 'istate', if ifind=1 *****
c
      subroutine symsqlat(istate,ifind,iflag,
     &                    n,nx,ny,
     &                    nsi,nmx,nmy,nmd,nkx,nky,
     &                    nmem0,
     &                    irprsnt,indsym)
c
c    n,nx,ny   @ lattice size, n=nx*ny
c    nkx,nky   @ momentum of x(y)-direction translational
c                operator Tx(y), in units of 2*pi/nx(y),
c                nx(y) is number of sites in x(y)-direction
c    istate    @ a given input state
c    iflag     # =1 only when a new class is successfully
c                generated, else iflag=0
c    nmem0     # number of member states in a class, =< n ,
c                number of classes >= subspace dimension
c    listmem   # stores the configuration belonging to
c                a class generated by the symmetry operation
c    charmem   # store the character of each member state
c
      implicit real*8(a-h,o-z)
      dimension  listmem(576)
      complex*16 charmem(576)
      
      dimension indsym(6)
      complex*16 chrctr
      complex*16 charac,charsum
c
      pi=2.d0*acos(0.d0)
c
      listmem(1) = istate
      charmem(1) = (1.d0,0.d0)
      iflag  = 0
      nmem0  = 1
c
      iseed  = istate
      ismall = istate
c
      do k=1,6
         indsym(k)=0
      end do
c
c*** runs over two spin inversion operations
      do 10 nspinv=0,1
c
c*** runs over two x-axis   mirror reflections
      do 10 nmrx=0,1
c           
c*** runs over two y-axis   mirror reflections
      do 10 nmry=0,1
c           
c*** runs over two diagonal mirror reflections
      do 10 nmrd=0,1
c
c*** runs over all x-axis translational operations
      do 10 ntrx=0,nx-1
c           
c*** runs over all y-axis translational operations
      do 10 ntry=0,ny-1

         if(     nspinv.eq.0
     &       .and. nmrx.eq.0 
     &       .and. nmry.eq.0 
     &       .and. nmrd.eq.0
     &       .and. ntrx.eq.0 
     &       .and. ntry.eq.0   ) go to 10

            itemp=iseed


            if(nspinv.eq.1) then
c*** perform spin inversion operation
            ispinv=ieor(itemp,2**n-1)
            itemp=ispinv
            end if


            if(nmrx.eq.1) then
            ireflx=0
c*** perform x-axis   mirror reflection operation (one loop)
            do iy=1,ny
            call mvbits( itemp, nx*(iy-1 ), nx,
     &                  ireflx, nx*(ny-iy)    )
            end do
            itemp=ireflx
            end if


            if(nmry.eq.1) then
            irefly=0
c*** perform y-axis   mirror reflection operation (two loops)
            do ix=1,nx
            do iy=1,ny
            call mvbits( itemp,     iy*nx-ix  ,  1, 
     &                  irefly, (iy-1)*nx+ix-1    )
            end do
            end do
            itemp=irefly
            end if


            if(nmrd.eq.1) then
            irefld=0
c*** perform diagonal mirror reflection operation (two loops)
            do ix=1,nx
            do iy=1,ny
            call mvbits( itemp, (nx-ix)*ny+ny-iy,  1,
     &                  irefld, (iy-1)*nx+ix-1      )
            end do
            end do
            itemp=irefld
            end if


            do itrx=1,ntrx   ! 'ntrx' x-axis translations
c*** translate 'itemp' to right in x-axis by 1 spacing (better)
            itranx=0
            do iy=1,ny
            itranx=itranx
     &     +iand(ishft(itemp,   1), 2**(iy*nx)-2**((iy-1)*nx+1))
     &     +iand(ishft(itemp,1-nx), 2**((iy-1)*nx)             )
            end do
            itemp=itranx
            end do


c*** translate 'itemp' to up in y-axis by ntry spacing (best)
            itrany=ishftc(itemp, -nx*ntry, n)
            itemp=itrany


            itrans=itemp

c*** if ifind=1, and the input state is greater, then
c    find the corresponding representative state
            if(ifind.eq.1 .and.
     &         itrans.lt.ismall) then
               ismall=itrans
c*** store the symmetry index temporarily
               indsym(1)=nspinv
               indsym(2)=nmrx
               indsym(3)=nmry
               indsym(4)=nmrd
               indsym(5)=ntrx
               indsym(6)=ntry
            else if(ifind.eq.1 .and.
     &         itrans.ge.ismall) then
c*** try out remaining symmetry operations
               go to 10
            else
               continue
            end if

c*** return if ifind=0 and the input state is greater,
c    i.e. it is not a representative state
            if(ifind.eq.0 .and.
     &         itrans.lt.istate) go to 20

c*** store the character temporarily
            charac=chrctr(n,nx,ny,
     &                       nsi, nmx, nmy, nmd, nkx, nky,
     &                    nspinv,nmrx,nmry,nmrd,ntrx,ntry)

c*** store a new generated configuration
            nmem0=nmem0+1
            listmem(nmem0)=itrans
            charmem(nmem0)=charac

 10   continue
c
c*** store the found representative and return (ifind=1)
      if(ifind.eq.1) then
         irprsnt=ismall
         go to 20
      end if
c
c*** continue to eliminate degenerate states   (ifind=0)

c*** check the degeneracy of symmetry-related states
      ndegen=0
      charsum=(0.d0,0.d0)
      do iu=1,nmem0
         if(listmem(iu).eq.istate) then
            ndegen=ndegen+1
            charsum=charsum+charmem(iu)
         end if
      end do
c*** check if the sum of charactors of degenerate states
c    vanishes or not
      if(abs(charsum) .gt. 1.0d-4) then
         nmem0=nmem0/ndegen
c        print *, 'nmem0=',nmem0
      else
         go to 20
      end if
c
c*** a new class is successfully generated
      iflag=1
 20   return
      end
c
c********** function 'chrctr' output character *************
c************* of an input symmetry operation **************
c
      function chrctr(n,nx,ny,
     &                   nsi, nmx, nmy, nmd, nkx, nky,
     &                nspinv,nmrx,nmry,nmrd,ntrx,ntry)
      implicit real*8(a-h,o-z)
      complex*16 chrctr

      pi=2.d0*acos(0.d0)

      phsreal=cos( 2.d0*pi*(nsi*nspinv)/2
     &           + 2.d0*pi*(nmx*nmrx)/2
     &           + 2.d0*pi*(nmy*nmry)/2
     &           + 2.d0*pi*(nmd*nmrd)/2
     &           + 2.d0*pi*(nkx*ntrx)/nx
     &           + 2.d0*pi*(nky*ntry)/ny )
      phsimag=sin( 2.d0*pi*(nsi*nspinv)/2
     &           + 2.d0*pi*(nmx*nmrx)/2
     &           + 2.d0*pi*(nmy*nmry)/2
     &           + 2.d0*pi*(nmd*nmrd)/2
     &           + 2.d0*pi*(nkx*ntrx)/nx
     &           + 2.d0*pi*(nky*ntry)/ny )

      chrctr=dcmplx(phsreal,phsimag)

      end
c
c============ Geometry-Independent Subroutines =============
c
c******** configurations with the specified sz **********
c****************  cf. szdy/TITPACK  ********************
c********** input 'nboson' instead of 'szval' ***********
c************** by Yi-Fei Wang, Jan. 2010 ***************
c
      subroutine szdysym(n,nx,ny,idim,idim0,
     &                   nsi,nmx,nmy,nmd,nkx,nky,
     &                   szval,
     &                   list1,list2,nmem)
c
c    n         @ lattice size
c    idim      # dimension of the Hilbert subspace with
c                given quantum numbers: Sz, kx, ...,
c                number of classes >= subspace dimension
c    idim0     @ estimated dimension, should be >= idim
c    nkx       @ momentum of x-direction translational
c                operator Tx, in units of 2*pi/nx, 
c                nx is number of sites in x-direction
c    nboson    @ total number of hard-core bosons,
c                i.e. total number of up-spins
c    list1(j)  # spin configuration with seiral number j, e.g.
c                j=1, list1(1)=|dduu>=(0011)_binary= 3_decimal
c                j=2, list1(2)=|dudu>=(0101)_binary= 5_decimal
c                j=6, list1(6)=|uudd>=(1100)_binary=12_decimal
c    list2     # inverse table of list1 expressed by the
c                2D (two-table) search method of H. Q. Lin :
c                list2(1,ia)+list2(2,ib)=ja+jb=j ,
c                list1(j)=i=(ia)(ib)_binary .
c
c========================================================
c     This routine is equivalent to szsym,
c     but is faster than szsym for larger n's.
c========================================================
      implicit real*8(a-h,o-z)
      real*8    szval
      dimension list1(idim0),list2(2,0:2**20),j1(65)
      dimension  nmem(idim0)

      dimension indsym(6)
c
c* initialization
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
      iupspn=n/2+mod(n,2)+int(szval+0.001d0)
c     iupspn=nboson
      icnt=0
c     ja=0
c     jb=0
c     ibpatn=0
c
      idim=0
      iflag=0
      icnt=0
c
      do 1 n1=1,iupspn
 1    j1(n1)=n1           ! record sites of up-spins from left
      j1(iupspn+1)=n+1
c
c*** main loop to generate decimal representations:
c    e.g. |s16 s15 ... s02 s01>
c        =|dddduuuudddduuuu>
c        =|0000111100001111>
c        =(0000111100001111)_binary
c        =3855_decimal
 2    continue
      icnt=icnt+1
c
      i=0
      do 4 k=1,iupspn
 4    i=ibset(i,j1(k)-1)  ! set bits at sites with up-spins
c
c     ia=iand(i,irght)
c     ib=iand(i,ilft)/ihfbit
c     if(ib.eq.ibpatn)then
c        ja=ja+1
c     else
c        ibpatn=ib
c        ja=1
c        jb=icnt-1
c     end if
c
c*** generate members of a symmetry class for a state 'i'
c    e.g. |0011> -- |0110> -- |1100> -- |1001> ,
c    or   |0101> -- |1010> .
         call symsqlat(i,0,iflag,
     &                 n,nx,ny,
     &                 nsi,nmx,nmy,nmd,nkx,nky,
     &                 nmem0,
     &                 irprsnt,indsym)
c
c*** store only the representative state with the
c    smallest decimal representation , 
c    e.g. |0011>, or |0101> .
      if(iflag .eq. 1) then
         idim=idim+1
         list1(idim)=i
c        list2(1,ia)=ja
c        list2(2,ib)=jb
         nmem(idim)=nmem0
         if(idim .gt. idim0) then
         print *,' #(E02)# Incorrect idim0 or n given to szdysym'
         stop
         end if
      end if
c
      do 5 k=1,iupspn
 5    if(j1(k+1) .gt. j1(k)+1) go to 6
c
      return              ! all basis states are found
c
 6    j1(k)=j1(k)+1       ! shift one up-spin to right
      if(k .eq. 1) go to 2
c
      do 7 k1=1,k-1
 7    j1(k1)=k1           ! record other up-spins from left
      go to 2
c 
      end
c
c******** finds whether the input state 'istate' is ********
c************** a representative of a class ****************
c
      subroutine search(istate,idim,idim0,iflag,jserl,list1)
c
c    idim      @ dimension of the matrix
c    list1     @ spin configurations generated in 'sz'
c    iflag     # =1 if the input state 'istate' is a
c                representative of a class
c    jserl     # the serial number of 'istate' if iflag=1
c
      implicit real*8(a-h,o-z)
      dimension list1(idim0)
c
      iflag=0
      n2=idim                  ! bisection search method,
      n1=1                     ! slower than 2-table method.
 34   if(n1 .gt. n2) go to 50
         jserl=(n1+n2)/2
         if(istate .gt. list1(jserl)) go to 35
         if(istate .eq. list1(jserl)) go to 36
         n2=jserl-1
         go to 34
 35   n1=jserl+1
      go to 34
 36   iflag=1
 50   continue
      return
      end
c
c******* display of an input state 'istate' with its *******
c*********** decimal and binary representation *************
c
      subroutine statdisp(n,istate,fnorm)
c
c    n         @ lattice size
c    istate    @ a given input state
c    fnorm     @ normalization factor of 'istate"
c
      implicit real*8 (a-h,o-z)
      dimension ibnry(40)
      real*8    fnorm
c
      itmp=istate
      do 10 i=n,1,-1
         ibnry(i)=mod(itmp,2)
         itmp=itmp/2
 10   continue
      write(*,100) istate, fnorm, (ibnry(i),i=1,n)
 100  format(1h ,i14,' ---', f16.10,'    ', 40i1)
      return
      end
c
c***********************************************************
c*                                                         *
c*             Static Correlation Functions                *
c*         for Spin-1/2 Systems with Symmetries            *
c*                                                         *
c*               by Yi-Fei Wang, Dec. 2009                 *
c*                                                         *
c*       Ref: H. Nishimori, TITPACK Ver 2.        (1991)   *
c*       Ref: G.Misguich et al., PRB 60,1064      (1999)   *
c*       Ref: A.Lauchli  et al., PRB 67,100409(R) (2003)   *
c*       Ref: D.N.Sheng  et al., PRB 79,205112    (2009)   *
c*                                                         *
c***********************************************************
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c
c*************** xx correlation function ****************
c****************  cf. TITPACK Ver 2. *******************
c------------- 2D square lattice geometry ---------------
c
      subroutine xcorr(n,nx,ny,idim,idim0,
     &                 nsi,nmx,nmy,nmd,nkx,nky,
     &                 x,sxx,
     &                 list1,list2,nmem)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    npair       @ pair of sites (k,l) <Sx(k)Sx(l)>
c    nbond       @ number of bonds to be investigated
c    x           @ eigenvetor
c    sxx         # xx correlation function
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      dimension   sxx(-nx:nx,-ny:ny)
      dimension   nxx(-nx:nx,-ny:ny)
      complex*16 corr(-nx:nx,-ny:ny)
      complex*16 x(idim0)
c
      dimension indsym(6)
      complex*16 chrctr,charac
      complex*16 ztmp
c
c*** initialization
      do ix=-nx,nx
      do iy=-ny,ny
          nxx(ix,iy)=0
         corr(ix,iy)=(0.d0,0.d0)
      end do
      end do
c
c*** loop over lattice sites
      do 21 isite1=0,n-1
         ix1=mod(isite1,nx)         ! 0,1,...,nx-1
         iy1=isite1/nx              ! 0,1,...,ny-1
      do 22 isite2=0,n-1
         ix2=mod(isite2,nx)         ! 0,1,...,nx-1
         iy2=isite2/nx              ! 0,1,...,ny-1
         ! distance between sites
         ixd=ix2-ix1                ! -(nx-1),...,nx-1
         iyd=iy2-iy1                ! -(ny-1),...,ny-1
         nxx(ixd,iyd)=nxx(ixd,iyd)+1
c
c*** loop over basis states
         do 30 j=1,idim
            istate0=list1(j)
            mask=2**isite1+2**isite2
            ibit=iand(istate0,mask)
            if(ibit.ne.0.and.ibit.ne.mask) then
            istate1=iexch(mask,istate0)
            call symsqlat(istate1,1,iflag0,
     &                    n,nx,ny,
     &                    nsi,nmx,nmy,nmd,nkx,nky,
     &                    nmem0,
     &                    irprsnt,indsym)
c*** bisection search method
            call search(irprsnt,idim,idim0,iflag,
     &                  jnew,list1)
            if(jnew.ne.0 .and. iflag.eq.1) then
               nspinv = indsym(1)
               nmrx   = indsym(2)
               nmry   = indsym(3)
               nmrd   = indsym(4)
               ntrx   = indsym(5)
               ntry   = indsym(6)
               charac=chrctr(n,nx,ny,
     &                          nsi, nmx, nmy, nmd, nkx, nky,
     &                       nspinv,nmrx,nmry,nmrd,ntrx,ntry)
               ztmp=charac*sqrt(dble( nmem(j   ) ))
     &                    /sqrt(dble( nmem(jnew) ))
               corr(ixd,iyd)=corr(ixd,iyd)
     &                      +conjg(x(j))*ztmp*x(jnew)
            end if
            end if
 30     continue

 22   continue
 21   continue 
c
c*** normalization
      do ix=-nx,nx
      do iy=-ny,ny
         sxx(ix,iy)=dreal(corr(ix,iy))
     &             /4.d0
c    &             /dble(nxx(ix,iy))
      end do
      end do
c
      return
      end
c
c*************** zz correlation function ****************
c****************  cf. TITPACK Ver 2. *******************
c------------- 2D square lattice geometry ---------------
c
      subroutine zcorr(n,nx,ny,idim,idim0,
     &                 nsi,nmx,nmy,nmd,nkx,nky,
     &                 x,szz,list1)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    npair       @ pair of sites (k,l) <Sz(k)Sz(l)>
c    nbond       @ number of bonds to be investigated
c    x           @ eigenvetor
c    szz         # zz correlation function
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension list1(idim0)
      dimension szz(-nx:nx,-ny:ny)
      dimension nzz(-nx:nx,-ny:ny)
      complex*16 x(idim0)
c
c*** initialization
      do ix=-nx,nx
      do iy=-ny,ny
         szz(ix,iy)=0.d0
         nzz(ix,iy)=0
      end do
      end do
c
c*** loop over lattice sites
      do 21 isite1=0,n-1
         ix1=mod(isite1,nx)         ! 0,1,...,nx-1
         iy1=isite1/nx              ! 0,1,...,ny-1
      do 22 isite2=0,n-1
         ix2=mod(isite2,nx)         ! 0,1,...,nx-1
         iy2=isite2/nx              ! 0,1,...,ny-1
         ! distance between sites
         ixd=ix2-ix1                ! -(nx-1),...,nx-1
         iyd=iy2-iy1                ! -(ny-1),...,ny-1
         nzz(ixd,iyd)=nzz(ixd,iyd)+1
c
c*** loop over basis states
         do 30 j=1,idim
            istate0=list1(j)
            if(isite1.ne.isite2) then
               mask=2**isite1+2**isite2
               ibit=iand(istate0,mask)
               if(ibit.eq.0.or.ibit.eq.mask) then
               szz(ixd,iyd)=szz(ixd,iyd)
     &                     +(abs(x(j))**2)*0.25d0
               else
               szz(ixd,iyd)=szz(ixd,iyd)
     &                     -(abs(x(j))**2)*0.25d0
               end if
            else   ! at the same site
               szz(ixd,iyd)=szz(ixd,iyd)
     &                     +(abs(x(j))**2)*0.75d0
            end if
 30      continue

 22   continue
 21   continue 
c
c*** normalization
      do ix=-nx,nx
      do iy=-ny,ny
         szz(ix,iy)=szz(ix,iy)
c    &             /dble(nzz(ix,iy))
      end do
      end do
c
      return
      end
c
c******** function "iexch" to exchange two spins ********
c**************** Yi-Fei Wang, Jul. 2010 ****************
c
      function iexch(mask,istate)
c
c    mask        @ mask of 2 spins to be exchanged
c    istate      @ initial state before 2-spin exchange
c    iexch       # final   state after  2-spin exchange
c
      implicit real*8 (a-h,o-z)
c
      ibit=iand(istate,mask)
      if(ibit.eq.0.or.ibit.eq.mask) then
         iexch=istate
      else
         iexch=ieor(istate,mask)
      end if
c
      return
      end
c
c************* total spin of an eigenstate **************
c**************** Yi-Fei Wang, Jan. 2009 ****************
c
      subroutine stot(n,nx,ny,idim,idim0,
     &                nsi,nmx,nmy,nmd,nkx,nky,
     &                x,stot2,
     &                list1,list2,nmem)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    npair       @ pair of sites (k,l) <Sz(k)Sz(l)>
c    nbond       @ number of bonds to be investigated
c    x           @ eigenvetor
c    stot2       # eigenvalue S(S+1) of S^2_{tot}
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      complex*16 x(idim0)
      dimension sxx(-nx:nx,-ny:ny),szz(-nx:nx,-ny:ny)
c
      call xcorr(n,nx,ny,idim,idim0,
     &           nsi,nmx,nmy,nmd,nkx,nky,
     &           x,sxx,
     &           list1,list2,nmem)
      call zcorr(n,nx,ny,idim,idim0,
     &           nsi,nmx,nmy,nmd,nkx,nky,
     &           x,szz,list1)
c
      stot2=0.d0
      do ix=-nx,nx
      do iy=-ny,ny
c        print *, 'sxx(ix,iy),szz(ix,iy)',sxx(ix,iy),szz(ix,iy)
         stot2=stot2+sxx(ix,iy)*2.d0
     &              +szz(ix,iy)
      end do
      end do
c
      return
      end
c
c***********************************************************
c*                                                         *
c*                   TITPACK Ver. 2                        *
c*                                                         *
c*                   February, 1991                        *
c*                                                         *
c*         Copyright (C) Hidetoshi Nishimori               *
c*                                                         *
c***********************************************************
c
c***** Complex H and |Psi>, by Yi-Fei Wang, Dec. 2008 ******
c******** Add symmetries,   by Yi-Fei Wang, Dec. 2009 ******
c
c============  SUBROUTINES / LARGE MATRICES =================
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c*** The following variables are common to all routines
c    n         @ lattice size
c    idim      @ matrix dimension
c    idim0     @ estimated subspace dimension
c    ipair     @ pairs of sites connected by bonds
c    bondwt    @ exchange interaction of each bond Jxy
c    zrtio     @ ratio of Jz to Jxy
c    ibond     @ number of bonds
c
c********* matrix elements for general bond weights ***********
c****************  cf. hamilt/D.Poilblanc  ********************
c
      subroutine elm1(n,nx,ny,idim,idim0,
     &                nsi,nmx,nmy,nmd,nkx,nky,
     &                ipair,bondwt,zrtio,ibond,
     &                diag,nondiag,
     &                list1,list2,nmem)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    ipair       @ pairs of sites connected by bonds
c    bondwt      @ exchange interaction of each bond Jxy
c    zrtio       @ ratio of Jz to Jxy
c    ibond       @ number of bonds
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
c
      complex*16 zi,ztmp
c
      dimension  mask(ibond),wght(ibond)
c
      dimension indsym(6)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
      zi=(0.d0,1.d0)
      pi=2.d0*acos(0.d0)
c
      call datack(ipair,ibond,n)
c
c*** initialization
      do j=1,idim
          diag(j)=0.d0
      end do
      do j=1,idim
      do k=1,ibond
          nondiag(j,k)=0
      end do
      end do
c
c*** store the masks which indicate the position of
c    the two spins to be exchanged
      do k=1,ibond
        isite1=ipair(k*2-1)-1
        isite2=ipair(k*2  )-1
        is1=2**isite1
        is2=2**isite2
        mask(k)=is1+is2
        wght(k)=bondwt(k)*zrtio(k)
      end do
c
      do 10 j=1,idim
c
c*** diagonal contributions,
        do 20 k=1,ibond
          ibit=iand(list1(j),mask(k))
          if(ibit.eq.0 .or. ibit.eq.mask(k)) then ! both down or up
            temp= wght(k)*0.25d0
            diag(j)=diag(j)+temp
          else  ! ibit=is1 or ibit=is2            ! one up, one down
            temp=-wght(k)*0.25d0
            diag(j)=diag(j)+temp
          end if
 20     continue
c
c*** off-diagonal contributions,
c
c    c.f. Weisse-Fehske, Lect.Notes.Phys.739,529(2008).
c    Let P_k be the projector to a fixed-k subspace :
c      P_k=(1/L)\sum^{L-1}_{m=0}[exp{i*k*m}T^{m}] .
c    The basis of fixed-k subspaces are given by
c     |r_k>=P_k|r>/sqrt{<r|P_k|r>} ,
c    and we discard those |r> with <r|P_k|r>=0 .
c    We note that [H, T]=0 , and therefore [H, P_k]=0 ,
c     <r'_k|H|r_k>
c    =<r'|P_k H P_k|r>/sqrt{<r'|P_k|r'>}/sqrt{<r|P_k|r>}
c    =<r'|P_k H|r>/sqrt{<r'|P_k|r'><r|P_k|r>} (D.Poilblanc)
c    =<r'|H P_k|r>/sqrt{<r'|P_k|r'><r|P_k|r>} (S.W.Haas),
c    i.e. we need to apply P_k only once !
c
         kcol=0
         do 25 k=1,ibond
            ibit=iand(list1(j),mask(k))
            
            if(ibit.ne.0 .and. ibit.ne.mask(k)) then
               kcol=kcol+1
               iexchg=ieor(list1(j),mask(k))
c*** find the representative of 'iexchg'
               call symsqlat(iexchg,1,iflag0,
     &                       n,nx,ny,
     &                       nsi,nmx,nmy,nmd,nkx,nky,
     &                       nmem0,
     &                       irprsnt,indsym)
c*** bisection search method
               call search(irprsnt,idim,idim0,iflag,
     &                     jnew,list1)
c              print *, 'jnew=',jnew
c              print *, 'iflag=',iflag
c              if(jnew .eq.0) print *, 'jnew=0'
c              if(iflag.eq.0) print *, 'iflag=0'
               if(jnew.ne.0 .and. iflag.eq.1) then
               ! store the found  'jnew'
               call mvbits(jnew     ,0,30,nondiag(j,kcol), 0)
               ! store symmetry operations 
               call mvbits(indsym(1),0, 1,nondiag(j,kcol),31)
               call mvbits(indsym(2),0, 1,nondiag(j,kcol),32)
               call mvbits(indsym(3),0, 1,nondiag(j,kcol),33)
               call mvbits(indsym(4),0, 1,nondiag(j,kcol),34)
               call mvbits(indsym(5),0, 3,nondiag(j,kcol),35)
               call mvbits(indsym(6),0, 3,nondiag(j,kcol),38)
               ! 2 for bond, 3 for ring3, 4 for ring4, etc
c              call mvbits(        2,0, 4,nondiag(j,kcol),41)                
               end if

             end if

 25       continue

 10   continue
c
      return
      end
c
c************ eigenvalues by the Lanczos method **************
c         --- dummy routine for simple working area
c
      subroutine lnc1(n,nx,ny,idim,idim0,
     &                nsi,nmx,nmy,nmd,nkx,nky,
     &                ipair,bondwt,zrtio,ibond,
     &                nvec,iv,E,itr,wk,ideclr,
     &                diag,nondiag,
     &                list1,list2,nmem)
c
c    nvec        @ number of eigenvectors to calculate in lncvec
c    iv          @ location of the nonzero element of the initial vector
c    E           # eigenvalues
c    itr         # number of iterations required for convergence
c    wk            working areas
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension E(4)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      complex*16 wk(ideclr,2)
c
      if(iv.le.0.or.iv.gt.idim)then
          print *,' #(E06)# Incorrect iv given to lnc1'
          return
      end if
      if(nvec.lt.0.or.nvec.gt.4)then
          print *,' #(W06)# Wrong value given to nvec in lnc1'
          print *,'         Only the eigenvalues are calculated'
          nvec=0
      end if
c
      call lnc1z(n,nx,ny,idim,idim0,
     &           nsi,nmx,nmy,nmd,nkx,nky,
     &           ipair,bondwt,zrtio,ibond,
     &           nvec,iv,E,itr,wk(1,1),wk(1,2),
     &           diag,nondiag,
     &           list1,list2,nmem)
      end
c
c************ eigenvalues by the Lanczos method
c
      subroutine lnc1z(n,nx,ny,idim,idim0,
     &                 nsi,nmx,nmy,nmd,nkx,nky,
     &                 ipair,bondwt,zrtio,ibond,
     &                 nvec,iv,E,itr,v1,v0,
     &                 diag,nondiag,
     &                 list1,list2,nmem)
c
      implicit real*8(a-h,o-z)
      dimension E(4)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      dimension wk(200,5),iwk(200)
      complex*16 v0(idim0),v1(idim0)
      complex*16 temp1,temp2
      common /vecdat/alpha(200),beta(200),coef(200,5)
c
c*** initialization
      do 10 i=1,idim
          v0(i)=0.0d0
          v1(i)=0.0d0
 10   continue
c
      v1(iv)=1.0d0
c
      call datack(ipair,ibond,n)
c
c*** alpha(1) and beta(1)
      call mltply(n,nx,ny,idim,idim0,
     &            nsi,nmx,nmy,nmd,nkx,nky,
     &            ipair,bondwt,zrtio,ibond,
     &            v1,v0,prdct,
     &            diag,nondiag,
     &            list1,list2,nmem)
      alpha1=prdct
      alpha(1)=alpha1
      beta1=0.d0
      do 50 i=1,idim
 50   beta1=beta1+(abs(v0(i)-alpha1*v1(i)))**2
      beta1=sqrt(beta1)
      beta(1)=beta1
c
C*** iteration
      do 100 i=2,200
          do 110 j=1,idim
            temp1=v1(j)
            temp2=(v0(j)-alpha1*v1(j))/beta1
            v0(j)=-beta1*temp1
            v1(j)=temp2
 110      continue
c
          call mltply(n,nx,ny,idim,idim0,
     &                nsi,nmx,nmy,nmd,nkx,nky,
     &                ipair,bondwt,zrtio,ibond,
     &                v1,v0,prdct,
     &                diag,nondiag,
     &                list1,list2,nmem)
          alpha1=prdct
          alpha(i)=alpha1
          beta1=0.d0
          do 120 j=1,idim
 120      beta1=beta1+(abs(v0(j)-alpha1*v1(j)))**2
          beta1=sqrt(beta1)
          beta(i)=beta1
          if(abs(beta(i)).lt.0.5d-30)then
            print *,' #(E07)# Tridiagonalization unsuccessful in lnc1'
            print *,'         Beta(i) is too small at i=  ',i
            stop
          end if
c
c*** convergence check
          if(i.gt.20.and.mod(i,5).eq.0)then
             call bisec(alpha,beta,i,E,4,eps)
             if(abs((ebefor-E(2))/E(2)).lt.1.0d-13)then
                if(nvec.gt.0)call vec12(E,i,nvec,wk(1,1),wk(1,2),
     &                            wk(1,3),wk(1,4),wk(1,5),iwk)
                itr=i
                return
             end if
             ebefor=E(2)
          end if
          if(i.eq.20)then
             eps=1.d-10
             call bisec(alpha,beta,20,E,4,eps)
             ebefor=E(2)
          end if
 100  continue
c
      print *,' #(W07)# lnc1 did not converge within 200 steps'
      itr=200
      return
      end
c
c************ eigenvector by the Lanczos method *************
c
      subroutine lncv1(n,nx,ny,idim,idim0,
     &                 nsi,nmx,nmy,nmd,nkx,nky,
     &                 ipair,bondwt,zrtio,ibond,
     &                 nvec,iv,x,itr,wk,ideclr,
     &                 diag,nondiag,
     &                 list1,list2,nmem)
c
c    nvec        @ number of eigenvectors to be calculated
c    iv          @ location of the nonzero element of the initial vector
c    x           # eigenvector
c    itr         @ number of interations for convergence
c    wk            working area
c    ideclr      @ declared array size in the main program
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      complex*16 x(ideclr,nvec)
      complex*16 wk(ideclr,2)
c
      if(nvec.le.0.or.nvec.gt.4)then
          print *,'#(W08)# nvec given to lncv1 out of range'
          return
      end if
      call lncv1z(n,nx,ny,idim,idim0,
     &            nsi,nmx,nmy,nmd,nkx,nky,
     &            ipair,bondwt,zrtio,ibond,
     &            nvec,iv,x,ideclr,itr,wk(1,1),wk(1,2),
     &            diag,nondiag,
     &            list1,list2,nmem)
      end
c
c************ eigenvector by the Lanczos method
c
      subroutine lncv1z(n,nx,ny,idim,idim0,
     &                  nsi,nmx,nmy,nmd,nkx,nky,
     &                  ipair,bondwt,zrtio,ibond,
     &                  nvec,iv,x,ideclr,itr,v1,v0,
     &                  diag,nondiag,
     &                  list1,list2,nmem)
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      complex*16 x(ideclr,nvec)
      complex*16 v0(idim0),v1(idim0)
      complex*16 temp1,temp2
      common /vecdat/alpha(200),beta(200),coef(200,5)
c
c*** initialization
      do 10 i=1,idim
          v0(i)=0.0d0
          v1(i)=0.0d0
 10   continue
      do 15 k=1,nvec
      do 15 i=1,idim
 15   x(i,k)=0.d0
c
      v1(iv)=1.0d0
      do 20 k=1,nvec
 20   x(iv,k)=coef(1,k)
c
c*** alpha(1) and beta(1)
      call mltply(n,nx,ny,idim,idim0,
     &            nsi,nmx,nmy,nmd,nkx,nky,
     &            ipair,bondwt,zrtio,ibond,
     &            v1,v0,prdct,
     &            diag,nondiag,
     &            list1,list2,nmem)

      alpha1=alpha(1)
      beta1=beta(1)
      do 40 k=1,nvec
      do 40 j=1,idim
 40   x(j,k)=x(j,k)+coef(2,k)*(v0(j)-alpha1*v1(j))/beta1
c
c*** iteration
      do 100 i=2,itr-1
          do 110 j=1,idim
            temp1=v1(j)
            temp2=(v0(j)-alpha1*v1(j))/beta1
            v0(j)=-beta1*temp1
            v1(j)=temp2
 110      continue
c
          call mltply(n,nx,ny,idim,idim0,
     &                nsi,nmx,nmy,nmd,nkx,nky,
     &                ipair,bondwt,zrtio,ibond,
     &                v1,v0,prdct,
     &                diag,nondiag,
     &                list1,list2,nmem)
          alpha1=alpha(i)
          beta1=beta(i)
          do 130 k=1,nvec
          do 130 j=1,idim
 130      x(j,k)=x(j,k)+coef(i+1,k)*(v0(j)-alpha1*v1(j))/beta1
 100  continue
c
c*** normalization
      do 200 k=1,nvec
          dnorm=0.d0
          do 210 j=1,idim
 210      dnorm=dnorm+(abs(x(j,k)))**2
          dnorm=sqrt(dnorm)
          do 220 j=1,idim
 220      x(j,k)=x(j,k)/dnorm
 200  continue
      return
      end
c
c***************** matrix multiplication ******************
c***** Complex H and |Psi>, by Yi-Fei Wang, Dec. 2008 *****
c******** Add symmetries,   by Yi-Fei Wang, Dec. 2009 *****
c****************  cf. hpsi/D.Poilblanc  ******************
c
      subroutine mltply(n,nx,ny,idim,idim0,
     &                  nsi,nmx,nmy,nmd,nkx,nky,
     &                  ipair,bondwt,zrtio,ibond,
     &                  v1,v0,prdct,
     &                  diag,nondiag,
     &                  list1,list2,nmem)
c
c    v1          @  input vector
c    v0          @# output vector : H*v1+v0(input)=v0+H1*v1+H2*v1
c    prdct       #  <v1*H*v1>=<v1*H1*v1>+<v1*H2*v1>
c    list1,list2 @  spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      complex*16 v0(idim0),v1(idim0)
      complex*16 prd,ztmp
      complex*16 offdg
c
      dimension indsym(6)
      complex*16 chrctr
      complex*16 charac
c
      prd=(0.d0,0.d0)
c
c*** H_{bond}*v1 and <v1*H_{bond}*v1>
      do 10 j=1,idim
c
c*** diagonal contributions, 
c
         v0(j)=v0(j)+diag(j)*v1(j)
c
c*** off-diagonal contributions,
         offdg=(0.d0,0.d0)

         do 25 k=1,ibond
            
            jnew  =0
            nspinv=0
            nmrx  =0
            nmry  =0
            nmrd  =0
            ntrx  =0
            ntry  =0
            ! read the stored  'jnew'
            call mvbits(nondiag(j,k), 0,30,  jnew, 0)
            ! read the stored symmetry operations 
            call mvbits(nondiag(j,k),31, 1,nspinv, 0)
            call mvbits(nondiag(j,k),32, 1,  nmrx, 0)
            call mvbits(nondiag(j,k),33, 1,  nmry, 0)
            call mvbits(nondiag(j,k),34, 1,  nmrd, 0)
            call mvbits(nondiag(j,k),35, 3,  ntrx, 0)
            call mvbits(nondiag(j,k),38, 3,  ntry, 0)

c           print *, 'jnew=',jnew
c           if(jnew.eq.0) print *, 'jnew=0'
            if(jnew.ne.0) then
            charac=chrctr(n,nx,ny,
     &                       nsi, nmx, nmy, nmd, nkx, nky,
     &                    nspinv,nmrx,nmry,nmrd,ntrx,ntry)
            ztmp=bondwt(k)*0.5d0
     &                    *charac  !exchange j, jnew below
     &                    *sqrt(dble( nmem(j   ) ))
     &                    /sqrt(dble( nmem(jnew) ))
            v0(j)=v0(j)+ztmp*v1(jnew)
            offdg=offdg+conjg(v1(j))*ztmp*v1(jnew)
            end if

 25      continue 
 
         prd=prd+diag(j)*(abs(v1(j)))**2+offdg

 10   continue
c
      prdct=dreal(prd)

      end
c
c*************** check of the eigenvector and eigenvalue ************
c********** Complex H and |Psi>, by Yi-Fei Wang, Dec. 2008 **********
c************* Add symmetries,   by Yi-Fei Wang, Dec. 2009 **********
c
      subroutine check1(n,nx,ny,idim,idim0,
     &                  nsi,nmx,nmy,nmd,nkx,nky,
     &                  ipair,bondwt,zrtio,ibond,
     &                  x,v,Hexpec,
     &                  diag,nondiag,
     &                  list1,list2,nmem)
c
c    x           @ eigenvector to be checked
c    v           # H*x=H1*x+H2*x
c    Hexpec      # <x*H*x>=<x*H1*x>+<x*H2*x>
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      complex*16 x(idim0),v(idim0)
      complex*16 zi,prd,ztmp
c
      dimension indsym(6)
      complex*16 chrctr
      complex*16 charac
c
      dnorm=0.d0
      do 5 j=1,idim
 5    dnorm=dnorm+(abs(x(j)))**2
      if(dnorm.lt.1.d-30)then
         print *,' #(W09)# Null vector given to check1'
         return
       end if
      do 8 j=1,idim
 8    v(j)=0.0d0
c 
      do 10 j=1,idim

         v(j)=v(j)+diag(j)*x(j)

         do 25 k=1,ibond
            jnew  =0
            nspinv=0
            nmrx  =0
            nmry  =0
            nmrd  =0
            ntrx  =0
            ntry  =0
            call mvbits(nondiag(j,k), 0,30,  jnew, 0)
            call mvbits(nondiag(j,k),31, 1,nspinv, 0)
            call mvbits(nondiag(j,k),32, 1,  nmrx, 0)
            call mvbits(nondiag(j,k),33, 1,  nmry, 0)
            call mvbits(nondiag(j,k),34, 1,  nmrd, 0)
            call mvbits(nondiag(j,k),35, 3,  ntrx, 0)
            call mvbits(nondiag(j,k),38, 3,  ntry, 0)
            if(jnew.ne.0) then
            charac=chrctr(n,nx,ny,
     &                       nsi, nmx, nmy, nmd, nkx, nky,
     &                    nspinv,nmrx,nmry,nmrd,ntrx,ntry)
            ztmp=bondwt(k)*0.5d0
     &                    *charac
     &                    *sqrt(dble( nmem(j   ) ))
     &                    /sqrt(dble( nmem(jnew) ))
            v(j)=v(j)+ztmp*x(jnew)
            end if
 25      continue 

 10   continue
c
      prd=0.0d0
c     print *,prd
      do 40 j=1,idim
 40   prd=prd+v(j)*conjg(x(j))
      Hexpec=dreal(prd)
      print *
      print 200
 200  format(' ------------ Information from check1 -------------- ')
      print 210,prd/n
 210  format(' <x*H*x>/n =', '(',1d16.8,',',1d16.8,')')
      print 220
 220  format(' H*x(j)/x(j) (j=min(idim/3,13),idim,max(1,idim/20))')
      print 230,(v(j)/x(j),j=min(idim/3,13),idim,max(1,idim/20))
 230  format('    (',1d18.9,',',1d18.9,')')
      print 240
 240  format(' ---------------------------------------------------'/)
      return
      end
c
c****************** inverse iteration ************************
c      --- dummy routine for simple working area
c
      subroutine inv1(n,nx,ny,idim,idim0,
     &                nsi,nmx,nmy,nmd,nkx,nky,
     &                ipair,bondwt,zrtio,ibond,
     &                Eig,iv,x,wk,ideclr,
     &                diag,nondiag,
     &                list1,list2,nmem)
c
c    Eig         @ eigenvalue
c    iv          @ location of the nonzero element of the initial vector
c    x           # eigenvector
c    wk            working area
c    ideclr      @ declared array size in the main program
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      complex*16 x(idim0)
      complex*16 wk(ideclr,4)
c
      call inv1z(n,nx,ny,idim,idim0,
     &           nsi,nmx,nmy,nmd,nkx,nky,
     &           ipair,bondwt,zrtio,ibond,
     &           Eig,iv,x,wk(1,1),wk(1,2),wk(1,3),wk(1,4),
     &           diag,nondiag,
     &           list1,list2,nmem)
      end
c
c****************** inverse iteration
c
      subroutine inv1z(n,nx,ny,idim,idim0,
     &                 nsi,nmx,nmy,nmd,nkx,nky,
     &                 ipair,bondwt,zrtio,ibond,
     &                 Eig,iv,x,b,p,r,y,
     &                 diag,nondiag,
     &                 list1,list2,nmem)
c
c    b            working area for the rhs of (H-E(approx))*x=b
c    p,r,y        working area used in the routine cg1
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      complex*16 b(idim0),x(idim0),r(idim0),y(idim0),p(idim0)
      complex*16 xb
c
      do 10 i=1,idim
 10   b(i)=0.0d0
      b(iv)=1.0d0
      do 20 itr=1,20
        call cg1(n,nx,ny,idim,idim0,
     &           nsi,nmx,nmy,nmd,nkx,nky,
     &           ipair,bondwt,zrtio,ibond,
     &           Eig,x,b,p,r,y,iterat,
     &           diag,nondiag,
     &           list1,list2,nmem)
        if(iterat.gt.idim)then
          xnorm=0.0d0
          do 22 i=1,idim
 22       xnorm=xnorm+(abs(x(i)))**2
          xnorm=sqrt(xnorm)
          do 24 i=1,idim
 24       x(i)=x(i)/xnorm
          print *,' #(W10)# Iterat in cg1 exceeds idim or 500'
          print *,'         Approximate eigenvector returned'
          print *,'         Itration number in inv1 is',itr
          return
        end if
        xnorm=0.0d0
        do 30 i=1,idim
 30     xnorm=xnorm+(abs(x(i)))**2
        xnorm=sqrt(xnorm)
        do 40 i=1,idim
 40     x(i)=x(i)/xnorm
        xb=0.0d0
        do 50 i=1,idim
 50     xb=xb+x(i)*conjg(b(i))
        if(abs(abs(xb)-1.0d0).lt.1.0d-12)then
c          print 100,itr
 100      format('       number of iterations in inv1 :',i5)
          return
        end if
        do 60 i=1,idim
 60     b(i)=x(i)
 20   continue
      print *,' #(W11)# inv1 did not converge'
      return
      end
c
c************** solution of linear equations -- cg method ***********
c********** Complex H and |Psi>, by Yi-Fei Wang, Dec. 2008 **********
c************* Add symmetries,   by Yi-Fei Wang, Dec. 2009 **********
c
      subroutine cg1(n,nx,ny,idim,idim0,
     &               nsi,nmx,nmy,nmd,nkx,nky,
     &               ipair,bondwt,zrtio,ibond,
     &               Eig,x,b,p,r,y,itr,
     &               diag,nondiag,
     &               list1,list2,nmem)
c
c    Eig         @ eigenvalue
c    x           # eigenvector
c    b             working area for the rhs of (H-E(approx))*x=b
c    p,r,y         working area used in the routine cg
c    itr         # number of iterations required for convergence
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension diag(idim0),nondiag(idim0,ibond)
      dimension list1(idim0),list2(2,0:2**20)
      dimension  nmem(idim0)
      complex*16 b(idim0),x(idim0),r(idim0),y(idim0),p(idim0)
      complex*16 zi,rp,yp,alpha,beta
      complex*16 ztmp
c
      dimension indsym(6)
      complex*16 chrctr
      complex*16 charac
c
c*** initialization
      bnorm=0.0d0
      do 8 j=1,idim
          bnorm=bnorm+(abs(b(j)))**2
          r(j)=b(j)
          p(j)=b(j)
          x(j)=0.0d0
 8    continue
c
c*** iteration
      do 200 itr=1,min(600,idim)

        do 30 j=1,idim
 30     y(j)=0.0d0
c
        do 10 j=1,idim

          y(j)=y(j)+(diag(j)-Eig)*p(j)

          do 25 k=1,ibond
            jnew  =0
            nspinv=0
            nmrx  =0
            nmry  =0
            nmrd  =0
            ntrx  =0
            ntry  =0
            call mvbits(nondiag(j,k), 0,30,  jnew, 0)
            call mvbits(nondiag(j,k),31, 1,nspinv, 0)
            call mvbits(nondiag(j,k),32, 1,  nmrx, 0)
            call mvbits(nondiag(j,k),33, 1,  nmry, 0)
            call mvbits(nondiag(j,k),34, 1,  nmrd, 0)
            call mvbits(nondiag(j,k),35, 3,  ntrx, 0)
            call mvbits(nondiag(j,k),38, 3,  ntry, 0)
            if(jnew.ne.0) then
            charac=chrctr(n,nx,ny,
     &                       nsi, nmx, nmy, nmd, nkx, nky,
     &                    nspinv,nmrx,nmry,nmrd,ntrx,ntry)
            ztmp=bondwt(k)*0.5d0
     &                    *charac
     &                    *sqrt(dble( nmem(j   ) ))
     &                    /sqrt(dble( nmem(jnew) ))
            y(j)=y(j)+ztmp*p(jnew)
            end if
 25       continue 

 10     continue
c
        rp=0.0d0
        yp=0.0d0
        do 60 j=1,idim
          rp=rp+r(j)*conjg(p(j))
          yp=yp+y(j)*conjg(p(j))
 60     continue
        alpha=rp/yp
        rnorm=0.0d0
        do 70 j=1,idim
          x(j)=x(j)+alpha*p(j)
          rnorm=rnorm+(abs(r(j)))**2
 70     continue
        rnorm2=0.0d0
        do 80 j=1,idim
          r(j)=r(j)-alpha*y(j)
          rnorm2=rnorm2+(abs(r(j)))**2
 80     continue
        beta=rnorm2/rnorm
        do 90 j=1,idim
 90     p(j)=r(j)+beta*p(j)
        if(mod(itr,5).ne.0)go to 200
        if(sqrt(rnorm2).lt.1.0d-9*sqrt(bnorm))then
c          print 150,itr
 150      format('       number of iterations in cg1     :',i5)
          return
        end if
 200  continue
c     print *,' #(Wxx)# cg1 did not converge'
      return
      end
c
c***********************************************************
c*                                                         *
c*                    TITPACK Ver. 2                       *
c*                                                         *
c*                    February, 1991                       *
c*                                                         *
c*          Copyright (C) Hidetoshi Nishimori              *
c*                                                         *
c***********************************************************
c
c***** Complex H and |Psi>, by Yi-Fei Wang, Dec. 2008 ******
c
c================== COMMON SUBROUTINES ====================
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c
c************* data check of pairs of sites ************
c
      subroutine datack(ipair,ibond,n)
      dimension ipair(ibond*2)
c
      do 10 k=1,ibond
         isite1=ipair(k*2-1)
         isite2=ipair(k*2  )
         if(isite1.le.0.or.isite2.le.0.or.
     &      isite1.gt.n.or.isite2.gt.n)then
            print *,' #(E03)# Incorrect data in ipair'
            print *,'         Location :  ',k*2-1,k*2
            stop
         end if
 10    continue
       end
c
c*********** eigenvalues of a real tridiagonal matrix ***********
c******************** by the bisection method *******************
c
      subroutine bisec(alpha,beta,ndim,E,ne,eps)
c
c    alpha  @ diagonal element
c    beta   @ subdiagonal element
c    ndim   @ matrix dimension
c    E      # eigenvalues
c    ne     @ number of eigenvalues to calculate
c    eps    @ limit of error
c
      implicit real*8 (a-h,o-z)
      dimension alpha(ndim),beta(ndim),E(ne),b2(2000)
c
      if(ndim.gt.2000)then
          print *,' #(E04)# ndim given to bisec exceeds 2000'
          stop
      end if
      if(ne.gt.ndim.or.ne.le.0)then
          print *,' #(E05)# ne given to bisec out of range'
          stop
      end if
c
c*** initial bound
      range=abs(alpha(1))+abs(beta(1))
      do 10 k=2,ndim-1
 10   range=max(range,abs(beta(k-1))+abs(alpha(k))+abs(beta(k)))
      range=max(range,abs(beta(ndim-1))+abs(alpha(ndim)))
      range=-range
c
      b2(1)=0.d0
      do 20 i=2,ndim
 20   b2(i)=(abs(beta(i-1)))**2
c
      epsabs=abs(range)*eps
      do 30 i=1,ne
 30   E(i)=-range
      b=range
c
c*** bisection method
      do 100 k=1,ne
        a=E(k)
        do 110 j=1,100
          c=(a+b)/2.d0
          if(abs(a-b).lt.epsabs)goto 100
          numneg=0
          g=1.d0
          ipass=0
          do 120 i=1,ndim
            if(ipass.eq.0)then
              g=c-alpha(i)-b2(i)/g
              else if(ipass.eq.1)then
                ipass=2
              else
                g=c-alpha(i)
                ipass=0
            end if
c
            if(ipass.eq.0)then
              if(g.le.0.d0)numneg=numneg+1
              if(abs(g).le.abs(b2(i)*epsabs*eps))ipass=1
            end if
 120      continue
          numneg=ndim-numneg
          if(numneg.lt.k)then
            b=c
          else
            a=c
            do 130 i=k,min(numneg,ne)
 130        E(i)=c
          end if
 110    continue
 100   continue
       end
c
c************ eigenvector of a real tridiagonal matrix **********
c******************** by inverse iteration **********************
c            ---  for the large/medium routines
c
      subroutine vec12(E,ndim,nvec,di,bl,bu,bv,cm,lex)
c
c    E(4)       @  4 lowest eigenvalues
c    ndim       @  matrix dimension
c    nvec       @  number of vectors to calculate
c    di - lex      working areas
c
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension di(ndim),bl(ndim),bu(ndim),bv(ndim),cm(ndim),lex(ndim)
      common /vecdat/alpha(200),beta(200),v(200,5)
c
      do 10 k=1,nvec
c
        do 100 j=1,ndim
          di(j)=E(k)-alpha(j)
          bl(j)=-beta(j)
          bu(j)=-beta(j)
 100    continue
c
c*** LU decomposition
        do 110 j=1,ndim-1
         if(abs(di(j)).gt.abs(bl(j)))then
c--- non pivoting
           lex(j)=0
           if(abs(di(j)).lt.1.d-13)di(j)=1.d-13
           cm(j+1)=bl(j)/di(j)
           di(j+1)=di(j+1)-cm(j+1)*bu(j)
           bv(j)=0.d0
         else
c--- pivoting
           lex(j)=1
           cm(j+1)=di(j)/bl(j)
           di(j)=bl(j)
           s=bu(j)
           bu(j)=di(j+1)
           bv(j)=bu(j+1)
           di(j+1)=s-cm(j+1)*bu(j)
           bu(j+1)= -cm(j+1)*bv(j)
         end if
 110    continue
        if(abs(di(ndim)).lt.1.d-13)di(ndim)=1.d-13
c
c--- initial vector
        do 120 j=1,ndim
 120    v(j,k)=1.d0/(float(j)*5.d0)
c
c*** degeneracy check up
        if(k.eq.1)then
           km=k
        else if(abs(E(k)-E(km)).gt.1.d-13)then
           km=k
        else
           do 130 i=km,k-1
             prd=0.d0
             do 140 j=1,ndim
 140         prd=prd+v(j,i)*v(j,k)
             do 150 j=1,ndim
 150         v(j,k)=v(j,k)-prd*v(j,i)
 130       continue
        end if
c
c*** inverse iteration
        do 160 l=1,k-km+3
          if((l.ne.1).or.(k.ne.km))then
c--- forward substitution
            do 170 j=1,ndim-1
              if(lex(j).eq.0)then
                v(j+1,k)=v(j+1,k)-cm(j+1)*v(j,k)
              else
                s=v(j,k)
                v(j,k)=v(j+1,k)
                v(j+1,k)=s-cm(j+1)*v(j,k)
              end if
 170        continue
          end if
c--- backward substitution
          do 180 j=ndim,1,-1
            s=v(j,k)
            if(j.le.ndim-1)s=s-bu(j)*v(j+1,k)
            if(j.le.ndim-2)s=s-bv(j)*v(j+2,k)
            v(j,k)=s/di(j)
 180      continue
c
c*** normalization
          dnorm=0.d0
          do 190 j=1,ndim
 190      dnorm=dnorm+(abs(v(j,k)))**2
          if(dnorm.gt.1.d-13)dnorm=1./sqrt(dnorm)
          do 200 j=1,ndim
 200      v(j,k)=v(j,k)*dnorm
 160    continue
c
 10   continue
      end
c
c******* Orthogonalization of the eigenvectors ************
c
      subroutine orthg(idim,ideclr,ev,norm,idgn,numvec)
c
c   idim    @  matrix dimension
c   ideclr  @  declared array size in the main program
c   ev      @# vectors to be orthogonalized / orthogonalized vectors
c   norm(j) #  norm of the j-th vector returned
c   idgn    #  degree of degenearcy
c   numvec  @  number of vectors to be checked
c
      implicit real*8(a-h,o-z)
      dimension norm(numvec)
      complex*16 ev(ideclr,numvec)
      complex*16 prjct,prd
c
      if(numvec.le.1)then
         print *,' #(W03)# Number of vectors is less than 2 in orthg'
         return
      end if
      do 10 i=1,numvec
        dnorm=0.0d0
        do 20 j=1,idim
 20     dnorm=dnorm+(abs(ev(j,i)))**2
        if(dnorm.lt.1.0d-20)then
           print *,' #(W04)# Null vector given to orthg. Location is',i
           return
        end if
        dnorm=1.0d0/sqrt(dnorm)
        do 25 j=1,idim
 25     ev(j,i)=ev(j,i)*dnorm
 10   continue
      idgn=numvec
      norm(1)=1
c*** orthogonalization
      do 30 i=2,numvec
       norm(i)=1
       do 40 j=1,i-1
         prjct=0.0d0
         do 50 l=1,idim
 50      prjct=prjct+ev(l,i)*conjg(ev(l,j))
         do 60 l=1,idim
 60      ev(l,i)=ev(l,i)-prjct*ev(l,j)
 40    continue
       vnorm=0.0d0
       do 70 l=1,idim
 70    vnorm=vnorm+(abs(ev(l,i)))**2
       if(vnorm.gt.1.0d-15)then
         vnorm=1.0d0/sqrt(vnorm)
         do 80 l=1,idim
 80      ev(l,i)=ev(l,i)*vnorm
        else
         do 90 l=1,idim
 90      ev(l,i)=0.0d0
         idgn=idgn-1
         norm(i)=0
       end if
 30   continue
c*** check orthogonality
      do 100 i=2,numvec
       do 100 j=1,i-1
       prd=0.0d0
       do 110 l=1,idim
 110   prd=prd+ev(l,i)*conjg(ev(l,j))
       if(abs(prd).lt.1.0d-10)go to 100
         print 200,i,j
 200     format(' #(W05)# Non-orthogonal vectors at',2i4)
         print 210,prd
 210     format('         Overlap : ',d14.7)
         print *,'         Unsuccessful orthogonalization'
         return
 100  continue
      return
      end
