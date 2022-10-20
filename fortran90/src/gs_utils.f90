module gs_utils
use gs_types

contains
subroutine init_all(gs, n, pattern, with_noise, method)
    use gs_numeric_utils, only : throw_exception
    implicit none
    type(grayscott_params), intent(inout) :: gs
    integer, intent(in) :: n, pattern, method
    logical, intent(in) :: with_noise

    allocate(gs%u(n,n), gs%v(n,n))   
    if (method .eq. 2) then 
        allocate(gs%au(n), gs%bu(n), gs%cu(n) )
        allocate(gs%av(n), gs%bv(n), gs%cv(n) )
        gs%uf = 0.5 * gs%dt * gs%mu        
        gs%vf = 0.5 * gs%dt * gs%mv 
        
        gs%au(:) = -gs%uf 
        gs%bu(:) = 1+2*gs%uf 
        gs%cu(:) = -gs%uf                 

        gs%av(:) = -gs%vf 
        gs%bv(:) = 1+2*gs%vf 
        gs%cv(:) = -gs%vf      

    else if (method .eq. 1) then
        allocate(gs%tu(n,n), gs%tv(n,n))   
        if (gs%dt .ne. 1.0) then 
            write(*,*) "EULER method only support Dt = 1, so we will use dt = 1"
            gs%dt = 1.0
        endif                 
    else 
        call throw_exception("Wrong method, valid are [1 and 2]")
    endif 

    if (pattern .eq. 1) then 
        call initialize_p1(gs%u, gs%v, with_noise)
    else 
        write(*,*) "Wrong initial pattern, valid patterns are [1]"
        write(*,*) "Selecting pattern 1"
        call initialize_p1(gs%u, gs%v, with_noise)
    endif 

    gs%method = method
    gs%sol_in_uv = .true.

end subroutine init_all

subroutine add_noise(u, v)
    implicit none
    real(dp), dimension(:,:), intent(inout) :: u, v
    real(dp), dimension(size(u,1), size(u,2)) :: rnds

    call random_init(repeatable=.true., image_distinct=.true.)
    call random_number(rnds)
    u = u + rnds*0.05

    call random_number(rnds)
    v = v + rnds*0.05
end subroutine add_noise

subroutine initialize_p1(u, v, with_noise)
    implicit none
    real(dp), dimension(:,:), intent(inout) :: u, v
    logical, intent(in) :: with_noise

    integer :: r, l 
    l = size(u,1)
    r = l/(l/16)

    u(:,:) = 1.0
    v(:,:) = 0.0

    u(l/2-16:l/2+16, l/2-16:l/2+16) = 0.5
    v(l/2-16:l/2+16, l/2-16:l/2+16) = 0.25

    if (with_noise .eqv. .true.) then 
        call add_noise(u,v)
    end if 

end subroutine initialize_p1

end module gs_utils 

