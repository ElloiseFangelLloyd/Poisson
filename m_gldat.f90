MODULE m_gldat
    use precision
    implicit none
    real(dp), allocatable, dimension(:,:,:) :: u, u_old, f
    real(dp), allocatable, dimension(:) :: x, y, z
    integer :: iostat, iter
    real(dp) :: h,L, residual
END MODULE m_gldat