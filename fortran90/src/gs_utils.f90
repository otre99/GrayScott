module gs_utils
use gs_types

contains
subroutine init_all(gs, n, pattern, with_noise, method)
    implicit none
    type(grayscott_params), intent(inout) :: gs
    integer, intent(in) :: n, pattern, method
    logical, intent(in) :: with_noise

    allocate(gs%u(n,n), gs%v(n,n))   
    if (method .eq. 2) then 
        allocate(gs%au(n), gs%bu(n), gs%cu(n) )
        allocate(gs%av(n), gs%bv(n), gs%cv(n) )
    else 
        allocate(gs%tempu(n,n), gs%tempv(n,n))   
    endif 

    if (pattern .eq. 1) then 
        call initialize_p1(gs%u, gs%v, with_noise)
    else 
        write(*,*) "Wrong initial pattern, valid patterns are [1]"
        write(*,*) "Selecting pattern 1"
        call initialize_p1(gs%u, gs%v, with_noise)
    endif 

    gs%method = method

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

