module gs_io
use gs_types

contains
subroutine save_matrix_to_file(fname, data)
    implicit none    
    character(len=*), intent(in) :: fname
    real(dp), dimension(:,:), intent(in) :: data     

    real(dp) :: rows, cols, n 
    n = size(data)
    rows = size(data,1)
    cols = size(data,2)
    open(unit=1, file=fname, action="write", status='replace', form='unformatted', access='stream')
    write(1) TRANSPOSE(data)     
    close(1)
end subroutine save_matrix_to_file

function gen_outname_u(iter)
    implicit none 
    integer(i4b), intent(in) :: iter    
    character(len=25) :: gen_outname_u
    write(gen_outname_u, '(a,I0.8)') 'UITER',iter       
end function gen_outname_u

function gen_outname_v(iter)
    implicit none 
    integer(i4b), intent(in) :: iter    
    character(len=25) :: gen_outname_v
    write(gen_outname_v, '(a,I0.8)') 'VITER',iter       
end function gen_outname_v

subroutine save_sol(gs, iter, outf)
    implicit none    
    type(grayscott_params), intent(in) :: gs   
    character(len=128), intent(in) :: outf
    integer, intent(in) :: iter 

    call save_matrix_to_file(trim(outf)//gen_outname_u(iter),gs%u)
    call save_matrix_to_file(trim(outf)//gen_outname_v(iter),gs%v)
end subroutine save_sol

end module gs_io
