program poisson3d

    ! Import different methods
    use precision
    use poisson_methods_m
    use m_gldat
    use m_alloc
    use m_init
  
    integer, save :: itermax = 1000
    real(dp), save :: tolerance = 1.e-4_dp
    integer, save :: algorithm = 2 ! 1 for Jacobi, 2 for Gauss-Seidel
  
   
    call read_namelist()   

    call alloc(N)   !allocates u, u_old, f fields and x,y,z vectors

    call init_bounds(N)   !initialise domain and boundary values

    call init_radiator(N)   !initialise f array

    iter = 0
    residual = 10 !initialise with a nonzero value

    do while (residual > tolerance .and. iter < itermax)    !quit if one of these is untrue

      call iterator(algorithm, N)   
      call calculate_residual(N)
      call update_old_field()

      iter = iter + 1

      !printing block
      IF (MOD(iter,10)==0) THEN
        print*,'Residual=',residual
        print*,"Iteration number: ",iter+1
        END IF

    end do
  
  end program poisson3d
  