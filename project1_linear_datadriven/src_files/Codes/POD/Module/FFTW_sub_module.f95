module FFTW_sub_module
    contains
        subroutine FFTW_sub(SVD_V,n)
            use, intrinsic :: iso_c_binding
            implicit none
            include 'fftw3.f03'
            integer*8 :: plan=0                                         ! Plan for FFTW
            integer, parameter :: dp = selected_real_kind(15, 307)    	! Double precision kind.                                                               //// User inputs

            integer :: m						                        ! Number of rows=n and columns=m of input matrix                                       //// User inputs
            integer :: j, i, io_status                                  ! Dynamic integer variables
            real(dp) :: input(8)
            real(kind=dp), allocatable :: row_in_a(:,:)                 ! FFTW function input real
            complex(kind=dp), allocatable :: row_in(:)                  ! FFTW function input complex
            complex(kind=dp), allocatable :: row_out(:)                 ! FFTW function output
            real(kind=dp),  allocatable :: SVD_V_mean(:,:)              ! VD_VT_mean (Mean of each column)
            complex(kind=dp), allocatable :: Psi(:,:)                   ! Storing FFT output


            integer, intent(in) :: n
            real(dp), dimension(:,:), intent(inout) :: SVD_V
            
            m=n

            allocate(SVD_V_mean(m,1))                              

            call Mean_calculation_row(n, m, SVD_V, SVD_V_mean)


            allocate(row_in(n),row_out(n),Psi(n,m),row_in_a(n,m))


            do i = 1, n                 
                do j = 1, m
                    
                    row_in_a(i, j) = SVD_V(i, j) - SVD_V_mean(j, 1)

                end do
            end do


            print *, ' FFT Start '

            do j = 1, m
                ! Copy column i into row_in for FFT computation
                row_in = (/(dcmplx(row_in_a(i, j), 0.0), i = 1, n)/)  ! i/ row j/ column

                ! Create FFTW plan for 1D FFT
                call dfftw_plan_dft_1d(plan, n, row_in, row_out, FFTW_BACKWARD, FFTW_ESTIMATE)

                ! Execute FFT
                call dfftw_execute(plan)

                ! Store the result in Psi (commented out as Psi is not defined)
                Psi(:,j) = row_out

                ! Destroy FFTW plan
                call dfftw_destroy_plan(plan)

            end do

            SVD_V = Psi


            print *, 'FFT operation done'
            print *, '------------------------------------------'

        end subroutine FFTW_sub


        ! Subroutine Mean
        subroutine Mean_calculation_row(n, m, u, u_mean)
            implicit none
            integer, parameter :: dp = selected_real_kind(15, 307)         ! Double precision kind. (Caution! higher numbers might cause memory issues)
            real(kind=dp) :: sum                                        ! Dynamic real variables
            integer :: t, k                                             ! Dynamic integer variables
            integer, intent(in) :: n, m
            real(kind=dp), intent(in) :: u(n,m)                                  ! Input arguments
            real(kind=dp), intent(out) :: u_mean(1,m)                            ! Output arguments
            
            do t=1, m
                sum = 0.00000000
                do k=1, n
                    sum = sum + u(k, t)
                end do
                u_mean(1, t) = sum / n
            end do
        end subroutine Mean_calculation_row
            
end module FFTW_sub_module