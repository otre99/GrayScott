module gs_solver
use gs_types
contains


subroutine laplace2d(a, b)
    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:,:), intent(out) :: b
    integer :: n, m  

    n = size(a, 1)
    m = size(a, 2)
    b(2:n-1, 2:m-1) = a(2:n-1,1:m-2) + a(2:n-1,3:m)&
                    + a(1:n-2,2:m-1) + a(3:n,2:m-1) - 4*a(2:n-1, 2:m-1)
                
    b(1,1) = a(1,2)+a(1,n)+a(n,1)+a(2,1)-4*a(1,1)
    b(n,m) = a(1,n)+a(n-1,m)+a(n,1)+a(n,n-1)-4*a(n,m)
    b(1,m) = a(1,n-1)+a(1,1)+a(n,m)+a(2,n)-4*a(1,m)
    b(n,1) = a(n,2)+a(n,m)+a(1,1)+a(n-1,1)-4*a(n,1)

    b(1,2:m-1) = a(1,1:m-2) + a(1,3:m) + a(2,2:m-1) + a(n,  2:m-1) - 4*a(1,2:m-1)
    b(n,2:m-1) = a(n,1:m-2) + a(n,3:m) + a(1,2:m-1) + a(n-1,2:m-1) - 4*a(n,2:m-1)

    b(2:n-1,1) = a(1:n-2,1) + a(3:n,1) + a(2:n-1,2) + a(2:n-1,  m) - 4*a(2:n-1,1)
    b(2:n-1,m) = a(1:n-2,m) + a(3:n,m) + a(2:n-1,1) + a(2:n-1,m-1) - 4*a(2:n-1,m)

end subroutine laplace2d



subroutine euler_time_step1dt(gs)
    type(grayscott_params) :: gs 
    
    call laplace2d(gs%u, gs%tempu)
    call laplace2d(gs%v, gs%tempv) 
   
    gs%tempu = gs%u + gs%mu*gs%tempu - gs%u*gs%v**2 + gs%f*(1.0-gs%u) 
    gs%tempv = gs%v + gs%mv*gs%tempv + gs%u*gs%v**2 - (gs%f + gs%k)*gs%v 

    gs%u = gs%tempu 
    gs%v = gs%tempv 

end subroutine euler_time_step1dt


!        do k=1, J_
!            if (kfactor_by_lats_(k) .ne. 1) then 
!                call reduce_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))
!                call solve_lon_problem_cn( sp(i)%phi(k,:I_/kfactor_by_lats_(k)),&
!                sp(i)%mulon(k,::kfactor_by_lats_(k)),R_, sinlatsC_(k), tau)
!                call expand_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))                
!            else                 
!                call solve_lon_problem_cn( sp(i)%phi(k,:), sp(i)%mulon(k,:), R_, sinlatsC_(k), tau)
!            endif 
!        end do        
!        !$omp end parallel do    
!    end do    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! SIGMA and FORCING [tn-1, tn+1]   !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do i=1, n_species        
!        sp(i)%phi_np = (2*tau*sp(i)%ff_np+(1.0-tau*sp(i)%sigma_np)*sp(i)%phi_np)/(1+tau*sp(i)%sigma_np)        
!        !$omp parallel do  
!        do k=1, I_
!            sp(i)%phi(:,k) = (2*tau*sp(i)%ff(:,k)+(1.0-tau*sp(i)%sigma(:,k))*sp(i)%phi(:,k))&
!            / (1.0+tau*sp(i)%sigma(:,k))        
!        end do 
!        !$omp end parallel do    
!        sp(i)%phi_sp = (2*tau*sp(i)%ff_sp+(1.0-tau*sp(i)%sigma_sp)*sp(i)%phi_sp)/(1+tau*sp(i)%sigma_sp)                
!    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! AdvDiff Lon direction [tn, tn+1] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do i=1, n_species        
!        !$omp parallel do  
!        do k=1, J_
!            if (kfactor_by_lats_(k) .ne. 1) then 
!                call reduce_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))
!                call solve_lon_problem_cn( sp(i)%phi(k,:I_/kfactor_by_lats_(k)),&
!                sp(i)%mulon(k,::kfactor_by_lats_(k)),R_, sinlatsC_(k), tau)
!                call expand_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))                
!            else                 
!                call solve_lon_problem_cn( sp(i)%phi(k,:), sp(i)%mulon(k,:), R_, sinlatsC_(k), tau)
!            endif 
!        end do        
!        !$omp end parallel do    
!    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! AdvDiff Lat direction [tn, tn+1] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do i=1, n_species        
!        call solve_lat_problem_cn(sp(i)%phi_np, sp(i)%phi, sp(i)%phi_sp, sp(i)%mulat,&
!        sinlatsC_, sinlatsB_, R_, tau)
!    end do

!    do i=1, n_species        
!        !$omp parallel do  
!        do k=1, J_
!            if (kfactor_by_lats_(k) .ne. 1) then 
!                call reduce_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))
!                call expand_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))                
!            endif 
!        end do        
!        !$omp end parallel do    
!     end do    
!     
!end subroutine cn_dim_split_linear



end module gs_solver 
