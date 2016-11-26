module M_list
	implicit none
	type t_list
		integer :: num=0
		type(t_node), pointer :: head => null()
		type(t_node), pointer :: tail => null()
	contains
		procedure :: insert
		final :: delete_list,delete_list_all
	end type
	type t_node
		class(*), pointer :: dat => null()
		type(t_node), pointer :: prev => null()
		type(t_node), pointer :: next => null()
	end type
contains
	subroutine insert(self,dat)
		class(t_list) :: self
		class(*), pointer :: dat
		type(t_node), pointer :: node
		allocate(node)
		node%dat => dat
		nullify(dat)

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
		if(associated(self%head)) then
			call delete_node(self%head)
		endif
		self%num=0
		nullify(self%head,self%tail)
	end subroutine
	subroutine delete_list_all(self)
		type(t_list) :: self(:)
		integer :: i
		do i=1,size(self)
			if(associated(self(i)%head)) then
				call delete_node(self(i)%head)
			endif
			self(i)%num=0
			nullify(self(i)%head,self(i)%tail)
		enddo
	end subroutine
	recursive subroutine delete_node(node)
		type(t_node), pointer :: node
		if(associated(node)) then
			if(associated(node%dat)) then
				deallocate(node%dat)
				nullify(node%dat)
			endif
			call delete_node(node%next)
			!nullify(node%prev)
			deallocate(node)
			!nullify(node)
		endif
	end subroutine
end module
