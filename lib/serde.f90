module M_serde
	implicit none
	type t_serde(n)
		integer, len :: n
		integer :: shap(n)
		!integer, allocatable :: i2idx(:),idx2i(:)
	contains
		procedure :: get_idx,get_i!,set_i2idx,set_idx2i
	end type
contains
	function get_i(self,idx) result(rt)
		class(t_serde(*)) :: self
		integer :: idx(:)
		integer :: rt
		integer :: prod,i
		rt=1
		prod=1
		do i=1,self%n
			rt=rt+(idx(i)-1)*prod
			prod=prod*self%shap(i)
		enddo
	end function
	function get_idx(self,i) result(idx)
		class(t_serde(*)) :: self
		integer :: i
		integer :: idx(self%n)
		integer :: j,prod
		prod=1
		do j=1,self%n
			idx(j)=mod((i-1)/prod,self%shap(j))+1
			prod=prod*self%shap(j)
		enddo
	end function
end module
