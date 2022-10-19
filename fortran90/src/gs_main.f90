PROGRAM gray_scott_problem
    USE gs_utils
    USE gs_solver
    USE gs_io 

    implicit none 
    integer(i4b) :: nn, sv, nsteps, method, pattern
    logical :: with_noise
    real(dp) :: dt, f, k, mu, mv
    character*256 :: out_folder
    integer :: iter 

    type(grayscott_params) :: gs_params
    !!!!!!!!!!!!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    namelist /input/ nn, dt, sv, nsteps, method, pattern, with_noise,&
    f, k, mu, mv, out_folder 
    open(10,file='namelist')
    read(10,nml=input)
    close(10)

    gs_params%k = k 
    gs_params%f = f
    gs_params%mu = mu 
    gs_params%mv = mv
    gs_params%dt = dt      
    call init_all(gs_params, nn, pattern, with_noise, method)
    write(*,*) "f=", gs_params%f, "k=", gs_params%k, "mu=", gs_params%mu, &
                "mv=", gs_params%mv

    call save_sol(gs_params, 0, out_folder)
    do iter=0, nsteps-1
        call euler_time_step1dt(gs_params)

        if (sv .gt. 0 .and. mod(iter, sv) == 0) then 
            write(*,*) "Save solution for iter = ", iter
            call save_sol(gs_params, iter, out_folder)
        endif 

    end do 
    call save_matrix_to_file(out_folder//"UFINAL", gs_params%u)
    call save_matrix_to_file(out_folder//"VFINAL", gs_params%v) 
    
END PROGRAM gray_scott_problem


    
