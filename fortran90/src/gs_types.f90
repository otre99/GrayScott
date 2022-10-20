module gs_types
    integer, parameter :: i4b = selected_int_kind(9)
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
    real(dp), parameter :: pio2=1.57079632679489661923132169163975144209858_dp
    real(dp), parameter :: twopi=6.283185307179586476925286766559005768394_dp
    real(dp), parameter :: sqrt2=1.41421356237309504880168872420969807856967_dp
    real(dp), parameter :: euler=0.5772156649015328606065120900824024310422_dp

    type grayscott_params
        real(dp) :: f, k, mu, mv, dt 
        real(dp) :: uf, vf 
        integer :: method 
        real(dp), dimension(:,:), allocatable :: u, v
        real(dp), dimension(:,:), allocatable :: tu, tv 
        real(dp), dimension(:), allocatable :: au, bu, cu
        real(dp), dimension(:), allocatable :: av, bv, cv
        logical :: sol_in_uv

    end type grayscott_params
end module gs_types 

