module poisson_methods_m

    ! Please do not put any computational variable in this module.
    ! Pass all variables as arguments to the subroutine.
    use precision, only: dp ! use dp as the data-type (and for defining constants!)
    use m_gldat
    implicit none
  
    private
    public :: jacobi
    public :: gauss_seidel
    public :: iterator 
    public :: calculate_residual
    public :: update_old_field
  
  contains
  
    subroutine jacobi(N)

      ! -----------------------------------------------------------------------------------------------------------
      ! Description: Performs one iteration of the Jacobi method for solving Poisson's equation.
      ! Inputs:
      !   - N: Integer representing the size of the domain (number of grid points along each dimension)
      ! Outputs:
      !   - Updates the values of the global array u based on neighboring values and the source term f.
      ! Notes:
      !   - The Jacobi method is an iterative numerical technique for solving linear systems of equations.
      !   - In this context, it is used to update the solution field u at each grid point.
      !   - The update formula involves contributions from neighboring grid points and the source term f.
      !   - The factor 1/6.0 is used to weight the contributions appropriately.
      !   - The loop indices i, j, and k iterate over the interior grid points (excluding boundaries).
      !   - The updated values are stored in the array u.
      ! -----------------------------------------------------------------------------------------------------------
        integer :: N
        integer :: i,j,k

        DO i = 2,N-1
            DO j = 2,N-1
                DO k = 2,N-1
                u(i,j,k) = 1.0/6.0 * ((h**2 * f(i,j,k)) + &
                u_old(i-1,j,k) + u_old(i+1,j,k) + &
                u_old(i,j-1,k) + u_old(i,j+1,k) + &
                u_old(i,j,k-1) + u_old(i,j,k+1))
                END DO 
            END DO 
        END DO  


    end subroutine jacobi

    subroutine gauss_seidel(N)

    ! -----------------------------------------------------------------------------------------------------------
    ! Description: Performs one iteration of the Gauss-Seidel method for solving a partial differential equation.
    ! Inputs:
    !   - N: Integer representing the size of the domain (number of grid points along each dimension)
    ! Outputs:
    !   - Updates the values of the global array u using the Gauss-Seidel update formula.
    ! Notes:
    !   - The Gauss-Seidel method is an iterative numerical technique for solving linear systems of equations.
    !   - In this context, it is used to update the solution field u at each grid point.
    !   - The update formula involves contributions from neighboring grid points and the source term f.
    !   - The factor 1/6.0 is used to weight the contributions appropriately.
    !   - The loop indices i, j, and k iterate over the interior grid points (excluding boundaries).
    !   - The updated values are stored in the array u.
    ! -----------------------------------------------------------------------------------------------------------

        integer :: N,i,j,k
          DO i = 2, N-1
            DO j = 2, N-1
              DO k = 2, N-1
                u(i, j, k) = 1.0d0/6.0d0 * (u(i-1, j, k) + u(i+1, j, k) + &
                                u(i, j-1, k) + u(i, j+1, k) + &
                                u(i, j, k-1) + u(i, j, k+1) +  f(i, j, k)*(h**2))
              END DO
            END DO
          END DO
          
      end subroutine gauss_seidel

    subroutine iterator(algorithm, N)
    ! -----------------------------------------------------------------------------------------------------------
    ! Description: Calls either the Jacobi or Gauss-Seidel method based on the specified algorithm.
    ! Inputs:
    !   - algorithm: Integer representing the chosen algorithm (1 for Jacobi, 2 for Gauss-Seidel)
    !   - N: Integer representing the size of the domain (number of grid points along each dimension)
    ! Works by:
    !   - Calls either the Jacobi or Gauss-Seidel subroutine based on the algorithm choice.
    ! Notes:
    !   - The Jacobi and Gauss-Seidel methods are iterative techniques for solving partial differential equations.
    !   - The subroutine updates the solution field u using the chosen method.
    ! -----------------------------------------------------------------------------------------------------------
        integer :: N, algorithm
        
        if (algorithm == 1) then 
            call jacobi(N)
          else if (algorithm == 2) then 
            call gauss_seidel(N)
          end if

    end subroutine iterator

    subroutine calculate_residual(N) 
    ! -----------------------------------------------------------------------------------------------------------
    ! Description: Computes the residual (error) between the updated solution u and the previous solution u_old.
    ! Inputs:
    !   - N: Integer representing the size of the domain (number of grid points along each dimension)
    ! Outputs:
    !   - Updates the global variable 'residual' with the computed residual value.
    ! Notes:
    !   - The residual is calculated as the square root of the sum of squared differences between u and u_old.
    !   - The loop indices i, j, and k iterate over the interior grid points (excluding boundaries).
    ! -----------------------------------------------------------------------------------------------------------
        integer :: N, i, j, k

        residual = 0.0
        DO i = 2,N-1
            DO j = 2,N-1
                DO k = 2,N-1
                residual = residual + SQRT((u(i,j,k) - u_old(i,j,k))**2)
                END DO 
            END DO 
        END DO 
        residual = residual / N
    end subroutine calculate_residual


    subroutine update_old_field()
    ! -----------------------------------------------------------------------------------------------------------
    ! Description: Updates the previous solution field u_old with the current solution field u.
    ! Inputs:
    !   - None (since it operates on global variables u and u_old)
    ! Outputs:
    !   - Updates the values of u_old to match the current solution u.
    ! Notes:
    !   - This subroutine is called after each iteration to maintain the previous solution.
    ! -----------------------------------------------------------------------------------------------------------
        u_old = u
    end subroutine update_old_field
  
  end module poisson_methods_m
  