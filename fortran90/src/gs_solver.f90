module gs_solver
use gs_types
contains

! subroutine laplace2d(a, b)
!     implicit none
!     real(dp), dimension(:,:), intent(in) :: a
!     real(dp), dimension(:,:), intent(out) :: b
!     integer :: n, m  

!     n = size(a, 1)
!     m = size(a, 2)
!     b(2:n-1, 2:m-1) = a(2:n-1,1:m-2) + a(2:n-1,3:m)&
!                     + a(1:n-2,2:m-1) + a(3:n,2:m-1) - 4*a(2:n-1, 2:m-1)
                
!     b(1,1) = a(1,2)+a(1,n)+a(n,1)+a(2,1)-4*a(1,1)
!     b(n,m) = a(1,m)+a(n-1,m)+a(n,1)+a(n,n-1)-4*a(n,m)
!     b(1,m) = a(1,m-1)+a(1,1)+a(n,m)+a(2,m)-4*a(1,m)
!     b(n,1) = a(n,2)+a(n,m)+a(1,1)+a(n-1,1)-4*a(n,1)

!     b(1,2:m-1) = a(1,1:m-2) + a(1,3:m) + a(2,2:m-1) + a(n,  2:m-1) - 4*a(1,2:m-1)
!     b(n,2:m-1) = a(n,1:m-2) + a(n,3:m) + a(1,2:m-1) + a(n-1,2:m-1) - 4*a(n,2:m-1)

!     b(2:n-1,1) = a(1:n-2,1) + a(3:n,1) + a(2:n-1,2) + a(2:n-1,  m) - 4*a(2:n-1,1)
!     b(2:n-1,m) = a(1:n-2,m) + a(3:n,m) + a(2:n-1,1) + a(2:n-1,m-1) - 4*a(2:n-1,m)

! end subroutine laplace2d

!subroutine print_matrix(m)
!     implicit none 
!     real(dp), dimension(:,:), intent(in) :: m
!     integer :: i, k 
!     do i=1, size(m,1)
!         write(*,*) (m(i,k), k=1, size(m,2))
!     enddo 
!end subroutine print_matrix

subroutine private_euler(u, v, tu, tv, mu, mv, f, k)
    implicit none 
    real(dp), dimension(:,:), intent(in) :: u, v 
    real(dp), dimension(:,:), intent(out) :: tu, tv 
    real(dp), intent(in) :: mu, mv, f, k

    real(dp), dimension(size(u,1)) :: uv2 
    integer :: i, n 

    n=size(u, 1)
    uv2=u(:,1)*v(:,1)**2
    ! left line 
    tu(2:n-1,1)=(u(1:n-2,1)+u(3:n,1)+u(2:n-1,2)+u(2:n-1,n)-4*u(2:n-1,1))*mu &
                -uv2(2:n-1)+f*(1-u(2:n-1,1))+ u(2:n-1,1) 
    tv(2:n-1,1)=(v(1:n-2,1)+v(3:n,1)+v(2:n-1,2)+v(2:n-1,n)-4*v(2:n-1,1))*mv &
                + uv2(2:n-1)-(f+k)*v(2:n-1,1)+v(2:n-1,1)
    ! tl
    tu(1,1)=(u(1,n)+u(1,2)+u(n,1)+u(2,1)-4*u(1,1))*mu &
            -uv2(1)+f*(1-u(1,1))+u(1,1)
    tv(1,1)=(v(1,n)+v(1,2)+v(n,1)+v(2,1)-4*v(1,1))*mv &
            +uv2(1)-(f+k)*v(1,1)+v(1,1)
    ! bl 
    tu(n,1)=(u(1,1)+u(n-1,1)+u(n,2)+u(n,n)-4*u(n,1))*mu &
            -uv2(n)+f*(1-u(n,1))+u(n,1)
    tv(n,1)=(v(1,1)+v(n-1,1)+v(n,2)+v(n,n)-4*v(n,1))*mv &
            +uv2(n)-(f+k)*v(n,1)+v(n,1)   

    ! right line                 
    uv2=u(:,n)*v(:,n)**2    
    tu(2:n-1,n)=(u(1:n-2,n)+u(3:n,n)+u(2:n-1,1)+u(2:n-1,n-1)-4*u(2:n-1,n))*mu &
               -uv2(2:n-1)+f*(1-u(2:n-1,n))+u(2:n-1,n)
    tv(2:n-1,n)=(v(1:n-2,n)+v(3:n,n)+v(2:n-1,1)+v(2:n-1,n-1)-4*v(2:n-1,n))*mv &
               +uv2(2:n-1)-(f+k)*v(2:n-1,n)+v(2:n-1,n)
    ! tr 
    tu(1,n)=(u(n,n)+u(2,n)+u(1,n-1)+u(1,1)-4*u(1,n))*mu &
            -uv2(1)+f*(1-u(1,n))+u(1,n)
    tv(1,n)=(v(n,n)+v(2,n)+v(1,n-1)+v(1,1)-4*v(1,n))*mv &
            +uv2(1)-(f+k)*v(1,n)+v(1,n)
    ! br 
    tu(n,n)=(u(1,n)+u(n-1,n)+u(n,n-1)+u(n,1)-4*u(n,n))*mu &
            -uv2(n)+f*(1-u(n,n))+u(n,n)
    tv(n,n)=(v(1,n)+v(n-1,n)+v(n,n-1)+v(n,1)-4*v(n,n))*mv &
            +uv2(n)-(f+k)*v(n,n)+v(n,n)
       
    !$omp parallel do private(uv2)  
    do i=2,n-1         
        uv2=u(:,i)*v(:,i)**2
        ! top
        tu(1,i)=(u(n,i)+u(2,i)+u(1,i-1)+u(1,i+1)-4*u(1,i))*mu &
               -uv2(1)+f*(1-u(1,i))+u(1,i)
        tv(1,i)=(v(n,i)+v(2,i)+v(1,i-1)+v(1,i+1)-4*v(1,i))*mv &
               +uv2(1)-(f+k)*v(1,i)+v(1,i)         
        ! center    
        tu(2:n-1,i)=(u(1:n-2,i)+u(3:n,i)+u(2:n-1,i-1)+u(2:n-1,i+1)-4*u(2:n-1,i))*mu &
                -uv2(2:n-1)+f*(1-u(2:n-1,i))+u(2:n-1,i)
        tv(2:n-1,i)=(v(1:n-2,i)+v(3:n,i)+v(2:n-1,i-1)+v(2:n-1,i+1)-4*v(2:n-1,i))*mv &
               +uv2(2:n-1)-(f+k)*v(2:n-1,i)+v(2:n-1,i)        
        ! bottom 
        tu(n,i)=(u(n-1,i)+u(1,i)+u(n,i-1)+u(n,i+1)-4*u(n,i))*mu &
               -uv2(n)+f*(1-u(n,i))+u(n,i)
        tv(n,i)=(v(n-1,i)+v(1,i)+v(n,i-1)+v(n,i+1)-4*v(n,i))*mv &
               +uv2(n)-(f+k)*v(n,i)+v(n,i)         
    enddo
    !$omp end parallel do    

end subroutine private_euler

subroutine euler_time_step(gs)
    implicit none 
    type(grayscott_params) :: gs 

    if (gs%sol_in_uv .eqv. .true.) then 
        call private_euler(gs%u,gs%v,gs%tu,gs%tv,gs%mu,gs%mv,gs%f,gs%k)
        gs%sol_in_uv = .false.
    else 
        call private_euler(gs%tu,gs%tv,gs%u,gs%v,gs%mu,gs%mv,gs%f,gs%k)
        gs%sol_in_uv = .true.
    endif 
end subroutine euler_time_step

subroutine cn_cyclic(gs, u, v)
    use gs_linalg, only : cyclic
    implicit none 
    real(dp), dimension(:), intent(inout) :: u, v  
    type(grayscott_params), intent(in) :: gs 
    
    real(dp), dimension(size(u,1)) :: ru, rv 
    integer :: n 
    
    n = size(u,1)
    ru(1) = (u(2) + u(n)) * gs%uf + (1 - 2 * gs%uf) * u(1)
    rv(1) = (v(2) + v(n)) * gs%vf + (1 - 2 * gs%vf) * v(1)
    
    ru(2:n-1) = (u(1:n-2) + u(3:n))*gs%uf + (1 - 2 * gs%uf)*u(2:n-1)
    rv(2:n-1) = (v(1:n-2) + v(3:n))*gs%vf + (1 - 2 * gs%vf)*v(2:n-1)
    
    ru(n) = (u(1) + u(n-1)) * gs%uf + (1 - 2 * gs%uf) * u(n)
    rv(n) = (v(1) + v(n-1)) * gs%vf + (1 - 2 * gs%vf) * v(n)
    
    call cyclic(gs%au, gs%bu, gs%cu, gs%cu(n), gs%au(1), ru, u)
    call cyclic(gs%av, gs%bv, gs%cv, gs%cv(n), gs%av(1), rv, v)
   
 end subroutine cn_cyclic


subroutine symrk2_time_step(gs)
    implicit none 
    type(grayscott_params) :: gs 
    real(dp), dimension(size(gs%u,1)) :: uv2, uk1, vk1, uk2, vk2
    integer :: i, n
    n = size(gs%u,1)

    ! xdir 
    !$omp parallel do  
    do i=1, n
        call cn_cyclic(gs, gs%u(i,:), gs%v(i,:))
    enddo
    !$omp end parallel do    
    ! ydir 
    !$omp parallel do  
    do i=1, n
        call cn_cyclic(gs, gs%u(:,i), gs%v(:,i))
    enddo
    !$omp end parallel do    
    ! reaction 
    !$omp parallel do private(uv2, uk1, vk1, uk2, vk2)  
    do i=1, n 
        uv2 = gs%u(i,:) * gs%v(i,:)**2
        uk1 = gs%u(i,:) + 0.5 * gs%dt * (-uv2 + gs%f * (1 - gs%u(i,:)))
        vk1 = gs%v(i,:) + 0.5 * gs%dt * ( uv2 - (gs%f + gs%k) * gs%v(i,:))
      
        uv2 = uk1 * vk1 ** 2
        uk2 = -uv2 + gs%f * (1 - uk1)
        vk2 =  uv2 - (gs%f + gs%k) * vk1
        gs%u(i,:) = gs%u(i,:) + gs%dt * uk2
        gs%v(i,:) = gs%v(i,:) + gs%dt * vk2
        !----------------------------------------
        uv2 = gs%u(i,:) * gs%v(i,:)**2
        uk1 = gs%u(i,:) + 0.5 * gs%dt * (-uv2 + gs%f * (1 - gs%u(i,:)))
        vk1 = gs%v(i,:) + 0.5 * gs%dt * ( uv2 -(gs%f + gs%k) * gs%v(i,:))
      
        uv2 = uk1 * vk1 **2
        uk2 = -uv2 + gs%f * (1 - uk1)
        vk2 =  uv2 - (gs%f + gs%k) * vk1
        gs%u(i,:) = gs%u(i,:) + gs%dt * uk2
        gs%v(i,:) = gs%v(i,:) + gs%dt * vk2
    enddo 
    !$omp end parallel do    
    ! ydir 
    !$omp parallel do  
    do i=1, n
        call cn_cyclic(gs, gs%u(:,i), gs%v(:,i))
    enddo
    !$omp end parallel do    
    ! xdir 
    !$omp parallel do  
    do i=1, n
        call cn_cyclic(gs, gs%u(i,:), gs%v(i,:))
    enddo
    !$omp end parallel do    

end subroutine symrk2_time_step


end module gs_solver 
