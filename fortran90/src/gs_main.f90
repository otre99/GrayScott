!&input
!f=0.055
!k=0.062
!mu=0.16
!mv=0.08
!nn=512
!dt=1.0
!sv=-1
!nsteps=20000
!method=1
!pattern=1
!with_noise=.true.
!out_folder="./"
!npar=6
!/
PROGRAM gray_scott_problem
    USE gs_utils
    USE gs_solver
    USE gs_io 

    implicit none 
    integer(i4b) :: nn, sv, nsteps, method, pattern
    logical :: with_noise
    real(dp) :: dt, f, k, mu, mv
    character*256 :: out_folder
    integer :: iter, npar  

    type(grayscott_params) :: gs_params
    !!!!!!!!!!!!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    namelist /input/ f, k, mu, mv, nn, dt, sv, nsteps, method, pattern, with_noise,&
    out_folder, npar  
    open(10,file='namelist')
    read(10,nml=input)
    close(10)

    call omp_set_num_threads(npar)
    gs_params%k = k 
    gs_params%f = f
    gs_params%mu = mu 
    gs_params%mv = mv
    gs_params%dt = dt      
    call init_all(gs_params, nn, pattern, with_noise, method)
    write(*,*) "**********************************************************************************"
    if (method .eq. 1) then 
        write(*,*) "Method: EULER"        
    else 
        write(*,*) "Method: SymRK2"        
    endif 
    write(*,"(a,F8.5,a,F8.5,a,F8.5,a,F8.5)") "  f = ", gs_params%f, "  k = ", gs_params%k, "     mu = ", gs_params%mu, &
                " mv = ", gs_params%mv               
    write(*,"(a,I8,a,F8.5,a,I8,a,I8, a,I8)") "Res = ", nn, " Dt = ", dt, " Nsteps = ", nsteps, " Sv = ", sv, " Npar = ", npar 
    write(*,*) "**********************************************************************************"
    call save_sol(gs_params, 0, out_folder)

    !dt = 0.0
    !do iter=1, nn
    !    do npar=1, nn
    !        gs_params%u(iter, npar) = dt 
    !        gs_params%v(iter, npar) = dt
    !        dt = dt + 1.0/16 
    !    enddo    
    !enddo 

    if (method .eq. 1) then ! EULER 
        do iter=0, nsteps-1
            call euler_time_step(gs_params)

            if (sv .gt. 0 .and. mod(iter, sv) == 0) then 
                write(*,*) "Save solution for iter = ", iter
                call save_sol(gs_params, iter, out_folder)
            endif 
        end do 
    else                    ! SymRK2 
        do iter=0, nsteps-1
            call symrk2_time_step(gs_params)

            if (sv .gt. 0 .and. mod(iter, sv) == 0) then 
                write(*,*) "Save solution for iter = ", iter
                call save_sol(gs_params, iter, out_folder)
            endif 
        end do         
    endif 

    !call print_matrix(gs_params%u)
    call save_matrix_to_file(trim(out_folder)//"UFINAL", gs_params%u)
    call save_matrix_to_file(trim(out_folder)//"VFINAL", gs_params%v) 
    
END PROGRAM gray_scott_problem


    
