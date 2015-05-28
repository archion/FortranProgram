# lib 

##Dependence
- lapack95

##usage

### M_matrix
  - det_ratio_row(c,ci,A,iA,pb)
  - det_ratio_col(c,ci,A,iA,pb)
  - det_ratio_rowcol(c,ci,A,iA,iY,pb)
  - inv_update_row(c,ci,pb,A,iA)
  - inv_update_col(c,ci,pb,A,iA)
  - inv_update_rowcol(c,ci,iY,A,iA)

### M_lattice

input:  
real(8) a1(2): principal axis 1 along x direct  
real(8) a2(2): principal axis 2  
real(8) T1(2): lattice size vector 1, along x direct  
real(8) T2(2): lattice size vector 2  
complex(8) bdc(2): boundary condition for T1 and T2  
real(8), allocate sub(:,2): sublattice position  

call:  
subroutine gen_latt(): give i2r  
subroutine gen_neb(): give neb  

get: 
integer Ns: the site number  
real(8) i2r(Ns,2): map site index to real space coordinate  
neb(Ns):  
neb(i)%nb(j)%bond(:): the (j-1)th neighbor of site i  
neb(i)%nb(j)%bdc(:): the boundary condition of neighbor site
neb(i)%nb(j)%r(:,2): the direction of neighbor site

subroutine test(): export lattice
