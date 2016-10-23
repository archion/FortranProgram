module M_list
	implicit none
	type t_list
		integer :: num=0
		type(t_node), pointer :: head
		type(t_node), pointer :: tail
	contains
		procedure :: insert
		final :: delete_list
	end type
	type t_node
		class(*), pointer :: dat
		type(t_node), pointer :: prev
		type(t_node), pointer :: next
	end type
contains
	subroutine insert(self,dat)
		class(t_list) :: self
		class(*), intent(in), pointer :: dat
		type(t_node), pointer :: node
		allocate(node)
		node%dat => dat

		if(self%num==0) then
			self%head => node
			self%tail => node
		else
			self%tail%next => node
			node%prev => self%tail
			self%tail => node
		endif

		self%num=self%num+1
	end subroutine
	subroutine delete_list(self)
		type(t_list) :: self
		call delete_node(self%head)
		self%num=0
		nullify(self%head,self%tail)
	end subroutine
	recursive subroutine delete_node(node)
		type(t_node),pointer :: node
		if(associated(node)) then
			if(associated(node%dat)) deallocate(node%dat)
			nullify(node%dat)
			call delete_node(node%next)
			nullify(node%prev)
			deallocate(node)
			nullify(node)
		endif
	end subroutine
end module
