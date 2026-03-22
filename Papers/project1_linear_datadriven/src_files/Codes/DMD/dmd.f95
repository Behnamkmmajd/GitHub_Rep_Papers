program DMD
    !use omp_lib                                                                    ! Include OpenMP library
    use stdlib_sorting, only: sort_index
    use HDF5
    use Sparity
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)                          ! Double precision kind.

    ! User inputs (Cehck input.txt for each case)
    integer ::      &
    nt,             &                                                               ! Number of snapshots
    time_phase,     &                                                               ! t0
    ns,             &                                                               ! Number of grid points	
    nv,             &                                                               ! 2D or 3D problem
    np,             &                                                               ! Number of column-wise partition
    np2,            &                                                               ! Number of row-wise partition
    n,              &                                                               ! Number of rows in Full stacked matrix
    m,              &                                                               ! Number of columns in Full stacked matrix
    nc,             &                                                               ! Number of columns in partitioned matrix m/np
    nc2,            &                                                               ! Number of columns in partitioned matrix n/np2
    rank,           &                                                               ! Rank 
    r_mode                                                                          ! User_rank
    integer :: nx,ny,nz, gama_n                                                     ! Number of columns in partitioned matrix n/np2

    integer, allocatable :: Index_s(:)
    integer, allocatable :: nzval(:)
    logical :: found_nan, found_inf

    real(dp), parameter :: pi = 4.0 * atan(1.0_dp), scall_p=1.0_dp                                  ! Double precision kind.                        

    real(dp) ::         &
    dt, rho, scale_factor                                                                              ! Sampling frequency

    integer(HID_T) :: file_id, group_id, group_id2, group_id3, dataset_id, dataspace_id, group_id4, group_id5, &
    group_id6, group_id7, group_id8
    integer ::  hdferr

    integer(HSIZE_T) :: dims2(2), dims1(1)


    integer :: info, lwork, io_status, t, j, i, ii, counter, tt2, &                 ! Dynamic integer variables (loop counters, etc)
    start_c, end_c, start_r, end_r, s_c, e_c, s_c2, e_c2, cc, r_answer, &
    mem_man, subtract_mean_flag, detrend_flag, project_eig_flag, &
    fd_preprocessing_flag                                                           ! FD preprocessing: position -> velocity
    
    complex(dp), allocatable :: Phi(:,:)                                            ! Spatial modes
    complex(dp), allocatable :: Phi_in(:,:)                                            ! Spatial modes
    complex(dp), allocatable :: Mu(:)                                               ! Eigenvalue of A
    complex(dp), allocatable :: Mu2(:)                                               ! Eigenvalue of A
    complex(dp), allocatable :: Mu_scaled(:)
    complex(dp), allocatable :: A_comp(:,:)                                         ! A, DMD matrix
    complex(dp), allocatable :: Vl(:,:)                                             ! Left eigen vector
    complex(dp), allocatable :: Vr(:,:)                                             ! Right eigen vector
    complex(dp), allocatable :: temp_three(:,:)                                     ! temporary dmd
    complex(dp), allocatable :: temp_four(:,:)                                      ! temporary dmd
    complex(dp), allocatable :: cwork(:)
    complex(dp), allocatable :: SVD_U_Comp(:,:)                                     ! Mode Energy
    complex(dp), allocatable :: SVD_VT_Comp(:,:)                                    ! SVD_VT
    complex(dp), allocatable :: work_Comp(:)                                        ! Workspace for SVD
    complex(dp), allocatable :: Amp_b_complex(:)                                        ! Workspace for SVD
    complex(dp), allocatable :: Phi_pinv(:,:)
    complex(dp), allocatable :: P(:,:)
    complex(dp), allocatable :: Vand(:,:)
    complex(dp), allocatable :: q(:)
    complex(dp), allocatable :: q_diag(:,:)
    complex(dp), allocatable :: Amp(:,:), Amp_p(:,:)

    complex(dp) :: cquery(1)

    real(dp), allocatable :: W(:)                                                   ! Output of eiv function (eigenvalues)
    real(dp), allocatable :: E(:)                                                   ! Output of eiv function (eigenvalues)
    real(dp), allocatable :: work(:)                                                ! Workspace for SVD
    real(dp), allocatable :: Stacked(:)                                             ! Matrix of stacked velocity
    real(dp), allocatable :: u_stacked(:,:) 
    real(dp), allocatable :: u_stacked1(:,:)                                        ! Matrix of stacked velocity
    real(dp), allocatable :: u_stacked2(:,:)                                        ! Matrix of stacked velocity
    real(dp), allocatable :: u_stacked3(:,:)                                        ! Matrix of stacked velocity
    real(dp), allocatable :: u_stacked_full(:,:)                                    ! Matrix of stacked velocity
    real(dp), allocatable :: R(:,:)                                                 ! Eigenvectors
    real(dp), allocatable :: A(:,:)                                                 ! A, DMD matrix
    real(dp), allocatable :: Phi_pod(:,:)                                           ! Spatial modes
    real(dp), allocatable :: Psi_lambda(:,:)                                        ! Stacked mode shapes
    real(dp), allocatable :: LAMBDA(:,:)
    real(dp), allocatable :: SVD_VT(:,:)                                            ! SVD_VT
    real(dp), allocatable :: SVD_VT_full(:,:)                                            ! SVD_VT
    real(dp), allocatable :: SVD_Sigma(:)                                           ! SVD_Sigma
    real(dp), allocatable :: SVD_S(:,:)
    real(dp), allocatable :: SVD_Sigma_full(:)                                           ! SVD_Sigma
    real(dp), allocatable :: SVD_Sigma_diag(:,:)                                    ! SVD_Sigma
    real(dp), allocatable :: SVD_U(:,:)                                             ! Mode Energy
    real(dp), allocatable :: SVD_U_full(:,:)                                             ! Mode Energy
    real(dp), allocatable :: Inv_Sigma(:,:)                                         ! SVD_Sigma
    real(dp), allocatable :: temp_one(:,:)                                          ! temporary dmd
    real(dp), allocatable :: temp_two(:,:)                                          ! temporary dmd
    real(dp), allocatable :: rwork(:)
    real(dp), allocatable :: Amp_b(:)                                               ! Workspace for SVD
    real(dp), allocatable :: Amp_b_ordered(:)
    real(dp), allocatable :: Gain(:)                                                ! Workspace for SVD
    real(dp), allocatable :: Freq(:)                                                ! Workspace for SVD
    real(dp), allocatable :: x(:),y(:),z(:), absMu(:)
    real(dp), allocatable :: G(:,:)
    real(dp), allocatable :: gamval(:), Jval(:), Jval_p(:), SVT(:,:)
    real(dp), allocatable :: x_mean(:)                                              ! Temporal mean of snapshots
    real(dp), allocatable :: x_slope(:)                                             ! Linear trend slope for each DOF
    real(dp), allocatable :: x_initial(:)                                           ! Initial position snapshot (for FD preprocessing reconstruction)
    real(dp), allocatable :: snap_prev(:), snap_next(:)                             ! Workspace for FD subroutine
    real(dp) :: t_bar, t_denom                                                      ! Mean time index and denominator for regression


    real(dp) :: query(1), Ej, total_energy

    real(dp) :: start_time, end_time, start_time2, end_time2, elapsed_time, L_start, L_stop, delta         ! CPU time variables

    character(1) :: jobvl, jobvr

    character(512) :: filename, filename2, filename_results                        ! Name of the files 
    character(1030) :: data_folder_dir                                             ! Where the velocity fields are stored
    character(512) :: data_folder                                                  ! The folder of the case
    character(512) :: results_folder_dir                                           ! Results folder
    character(512) :: data_folder_dir2, jobz, uplo, data_folder_dir3


    found_nan = .false.
    found_inf = .false.
    !******************************************************* Directories
    data_folder_dir2 = '../../database/'                                  ! directory of velocity field
    results_folder_dir = './results/'                                              ! Directory to save the files


    ! Prompt the user for the data folder input
    print *, 'Enter the data folder (e.g., data_1, data_2, ...):'
    read(*, '(A)') data_folder  ! Read user input into data_folder

    ! Construct the full path by concatenating base_dir, user input, and the rest of the path
    data_folder_dir = trim(data_folder_dir2)//trim(data_folder)//"/data/"
    data_folder_dir3 = trim(data_folder_dir2)//trim(data_folder)//"/"

    !_______________________________________________________     Program initialization
    !_______________________________________________________
    print *, 'Program initialization'
    print *, '------------------------------------------'
    print *, 'Data directory : ', data_folder_dir
    print *, 'Result directory : ', results_folder_dir 

    print *, 'Reading the case info : '

    write(filename, '(A)') trim(data_folder_dir3)//"input_DMD.txt"
    open(2000, file=TRIM(filename), status='OLD', iostat=io_status)
    if (io_status /= 0) then
        print *, 'Error opening file: ', filename
        stop
    end if

    read(2000, *)
    read(2000, *)
    read(2000, *)

    
    read(2000, *) filename

    filename_results = trim(results_folder_dir)//"Results_"//trim(filename)//".h5"

    write(filename2, '(A)') trim(results_folder_dir)//"input_"//trim(filename)//".txt"
    open(unit=510, file=filename2)

    print *, ' Case :                                       ', filename
    read(2000, *) nt
    write(510, *) nt
    print *, ' Number of snapshots :                        ', nt
    print *, '______________________'
    read(2000, *) dt
    write(510, *) dt
    print *, ' Sampling Frequency :                         ', 1.0/dt, 'Hz'
    print *, '______________________'
    read(2000, *) time_phase
    write(filename, '(A,I0,A)') 'Stacked_', 1+time_phase, '.csv'
    print *, ' The name of the first snapshot is :          ', trim(filename)
    print *, '______________________'
    read(2000, *) ns
    print *, ' Number of grid points :                      ', ns
    write(510, *) ns
    print *, '______________________'
    read(2000, *) nv
    print *, ' 2D or 3D case :                              ', nv,"D"
    write(510, *) nv
    print *, '______________________'
    read(2000, *) np
    print *, ' Number of column-wise partition :            ', np
    write(510, *) np
    print *, '______________________'
    read(2000, *) np2
    print *, ' Number of row-wise partition :               ', np2
    write(510, *) np2
    print *, '______________________'
    read(2000, *) r_mode
    write(510, *) r_mode
    read(2000, *) r_answer
    if (r_answer==1) then
        print *, ' Truncated mode disabled, User Rank :          ', r_mode
    else
        print *, ' Truncated mode enabled '
    end if
    print *, '______________________'

    read(2000, *) nx
    write(510, *) nx

    read(2000, *) ny
    write(510, *) ny

    read(2000, *) nz
    write(510, *) nz

    read(2000, *) gama_n

    print *, ' Number of Gamma value ', gama_n

    read(2000, *) L_start 

    print *, ' Lower bound ', L_start

    read(2000, *) L_stop 

    print *, ' Upper bound ', L_stop

    read(2000, *) rho

    print *, ' rho ', rho

    ! Skip optional rho_min / rho_max lines (used by amp.f95 only)
    read(2000, *)  ! rho_min
    read(2000, *)  ! rho_max

    read(2000, *) mem_man

    print *, ' Memory management (1=on 0=off) ', mem_man

    ! Preprocessing options (0=off, 1=on)
    read(2000, *, iostat=io_status) subtract_mean_flag
    if (io_status /= 0) subtract_mean_flag = 1   ! Default: mean subtraction ON
    print *, ' Mean subtraction (1=on 0=off) ', subtract_mean_flag

    read(2000, *, iostat=io_status) detrend_flag
    if (io_status /= 0) detrend_flag = 0          ! Default: linear detrending OFF
    print *, ' Linear detrending (1=on 0=off) ', detrend_flag

    read(2000, *, iostat=io_status) project_eig_flag
    if (io_status /= 0) project_eig_flag = 0       ! Default: eigenvalue projection OFF
    print *, ' Eigenvalue projection (1=on 0=off) ', project_eig_flag

    read(2000, *, iostat=io_status) fd_preprocessing_flag
    if (io_status /= 0) fd_preprocessing_flag = 0  ! Default: FD preprocessing OFF
    print *, ' FD preprocessing pos->vel (1=on 0=off) ', fd_preprocessing_flag

    close(2000)

    print *, ' Is everything fine with the above case (Yes=1 / No=0 )?'
    read(*,*) info  ! User answer


    allocate(x(nx),y(ny),z(nz))

    if (info==1) then


        write(filename, '(A)') trim(data_folder_dir)//"x.csv"
        open(2000, file=TRIM(filename), status='OLD', iostat=io_status)
        do i=1, nx
            read(2000, *) x(i)
        end do
        close(2000)

        write(filename, '(A)') trim(data_folder_dir)//"y.csv"
        open(2000, file=TRIM(filename), status='OLD', iostat=io_status)

        do i=1, ny
            read(2000, *) y(i)
        end do
        close(2000)

        write(filename, '(A)') trim(data_folder_dir)//"z.csv"
        open(2000, file=TRIM(filename), status='OLD', iostat=io_status)

        do i=1, nz
            read(2000, *) z(i)
        end do
        close(2000)

        n=nv*ns
        m=nt

        ! Allocate FD workspace and save initial position snapshot
        allocate(snap_prev(n), snap_next(n))   ! Workspace for read_snapshot_fd subroutine
        if (fd_preprocessing_flag == 1) then
            allocate(x_initial(n))
            write(filename, '(A, I0, A)') trim(data_folder_dir)//"Stacked_", 1+time_phase, ".csv"
            open(unit=2000, file=TRIM(filename), status='OLD', iostat=io_status)
            if (io_status /= 0) then
                print *, 'Error: cannot read first snapshot for x_initial: ', trim(filename)
                stop
            end if
            read(2000, *, iostat=io_status) x_initial
            close(2000)
            x_initial = x_initial * scall_p
            print *, 'FD preprocessing: initial position saved for reconstruction'
        end if

        !_______________________________________________________Reading the velocity field and preparation of the stacked matrix
        !_______________________________________________________
        !_______________________________________________________
        !_______________________________________________________Allocation of the variables

        !_______________________________________________________Reading the velocity field and stacking procedure
        !_______________________________________________________
        !_______________________________________________________

        print *, 'Reading Velocity field'
        print *, '------------------------------------------'

        ! Initialize HDF5 library
        call h5open_f(hdferr)

        ! Create a new HDF5 file
        call h5fcreate_f(filename_results, H5F_ACC_TRUNC_F, file_id, hdferr)
        ! Create a group
        call h5gcreate_f(file_id, "Raw", group_id, hdferr)

        call h5gcreate_f(file_id, "Standard_Decomposition", group_id6, hdferr)

        call h5gcreate_f(group_id6, "A", group_id7, hdferr)

        call h5gcreate_f(group_id, "x", group_id4, hdferr)

        call h5gcreate_f(group_id, "xDMD_E", group_id3, hdferr)

        call h5gcreate_f(group_id, "SVD_U", group_id5, hdferr)

        call h5gcreate_f(file_id, "Phi", group_id8, hdferr)


        ! allocate(u_stacked_full(n,m))

        ! do t=1, m

        !     write(filename, '(A, I0, A)') trim(data_folder_dir)//"Stacked_",t+time_phase,".csv"
        !     open(unit=2000, file=TRIM(filename), status='OLD', iostat=io_status)
        !     if (io_status /= 0) then
        !         print *, 'Error opening file: ', filename
        !         stop
        !     end if

        !     read(2000, *, iostat=io_status) u_stacked_full(:,t)

        !     close(2000)

        ! end do


        ! dims2(1) = m 
        ! dims2(2) = n
        ! call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        

        ! ! Create the dataset within the "/Phi" group with double precision
        ! call h5dcreate_f(group_id, "X", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(u_stacked_full), dims2, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)

        ! deallocate(u_stacked_full)

        ! Measure the start time
        call CPU_TIME(start_time2)


        if (mem_man == 1) then

            nc=(m-1)/np
            nc2=n/np2

            allocate(u_stacked1(n, nc),u_stacked2(n,nc),Stacked(n),R(m-1,m-1))

            ! Compute temporal mean and linear trend of all m snapshots
            allocate(x_mean(n), x_slope(n))
            x_mean = 0.0_dp
            x_slope = 0.0_dp
            t_bar = real(m + 1, dp) / 2.0_dp
            t_denom = real(m, dp) * (real(m, dp)**2 - 1.0_dp) / 12.0_dp

            if (subtract_mean_flag == 1) then
                do t = 1, m
                    call read_snapshot_fd(data_folder_dir, t, time_phase, m, n, dt, &
                                         fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                    x_mean = x_mean + Stacked * scall_p
                end do
                x_mean = x_mean / real(m, dp)
                print *, 'Temporal mean computed'
            end if

            if (detrend_flag == 1 .and. subtract_mean_flag == 1) then
                ! Second pass for slope: b_i = sum_k (k - t_bar) * (x_i(k) - x_mean_i) / t_denom
                do t = 1, m
                    call read_snapshot_fd(data_folder_dir, t, time_phase, m, n, dt, &
                                         fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                    x_slope = x_slope + (real(t, dp) - t_bar) * (Stacked * scall_p - x_mean)
                end do
                x_slope = x_slope / t_denom
                print *, 'Linear trend computed'
            end if

            R=0.0_dp

            counter=0

            do counter=1, np

                call CPU_TIME(start_time)

                s_c=(nc*(counter-1))+1
                e_c=counter*nc

                cc=1
                ! Measure the start time

                do i=s_c, e_c
                    call read_snapshot_fd(data_folder_dir, i, time_phase, m, n, dt, &
                                         fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                    u_stacked1(:,cc)=Stacked*scall_p - x_mean - x_slope*(real(i,dp) - t_bar)
                    cc=cc+1
                end do

                ! Measure the end time
                ! call CPU_TIME(end_time)
                ! ! Calculate the elapsed time
                ! elapsed_time = end_time - start_time
                !print '(A, I0, F0.3, A)', 'Elapsed time for reading column-wise partition:' ,elapsed_time, ' seconds'
                ! Measure the start time
                ! call CPU_TIME(start_time)

                do tt2=counter, np

                    if (counter==tt2) then

                        start_c=((tt2-1)*nc)+1
                        end_c=start_c-1+nc
            
                        start_r=((counter-1)*nc)+1
                        end_r=start_r-1+nc
            
                        R(start_r:end_r,start_c:end_c)=matmul(transpose(u_stacked1),u_stacked1)

                    else 

                        start_c=((tt2-1)*nc)+1
                        end_c=((tt2-1)*nc)+nc
            
                        start_r=((counter-1)*nc)+1
                        end_r=((counter-1)*nc)+nc

                        s_c2=(nc*(tt2-1))+1
                        e_c2=tt2*nc

                        cc=1
                        do ii=s_c2, e_c2
                            call read_snapshot_fd(data_folder_dir, ii, time_phase, m, n, dt, &
                                                 fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                            u_stacked2(:,cc)=Stacked*scall_p - x_mean - x_slope*(real(ii,dp) - t_bar)
                            cc=cc+1
                        end do
                        R(start_r:end_r,start_c:end_c)=matmul(transpose(u_stacked1),u_stacked2)
                    end if
                end do
                
                ! Measure the end time
                call CPU_TIME(end_time)
                ! Calculate the elapsed time
                elapsed_time = end_time - start_time
                print   '(A, F10.3, A, I0, A, I0)', 'Hrs reamining for R (Time correlation MATRIX)', &
                ((np-(counter-1))*elapsed_time)/3600.0_dp, &
                '__column-wise Partition__' , counter, '/' ,np

            end do

            deallocate(u_stacked1,u_stacked2)

            ! Completing the other half of the matrix since R is symmetry 
            do i=1, m-1
                do j=i+1, m-1
                    R(j,i)=R(i,j)
                end do
            end do

                ! Measure the end time
            call CPU_TIME(end_time2)
                ! Calculate the elapsed time
            elapsed_time = end_time2 - start_time2
            print *, 'Elapsed time for calculating R: ', elapsed_time/3600.0_dp, 'hours'
            !_________________________________________________________________________________________________________________________
            !_________________________________________________________________________________________________________________________
            !___________________________________________________________________________________________________Covariance matrix R end

            print *, ' R is calculated '
            print *, ' ------------------------------------------ '

            allocate(W(m-1))
            
            print *, ' '
            print *, ' Eigenvalue calculation in progress ............... '

            ! _________________________________________________________________________________________________________________________


            ! _________________________________________________________________________________________________________________________
            ! _________________________________________________________________________________________________________________________
            ! ____________________________________________________________________________________________Deallocation of variables end



            !_________________________________________________________________________________________________________________________
            !_________________________________________________________________________________________________________________________
            !____________________________________________________________________Eigen_value_problem________________Method of snapshot
            jobz = 'V' ! Compute eigenvalues and eigenvectors
            uplo = 'U' ! Upper triangle of A is stored

            ! Query for optimal workspace size
            call dsyev(jobz, uplo, m-1, R, m-1, W, query, -1, info)
            lwork = int(query(1))
            allocate(work(lwork))

            !Call LAPACK routine DSYEV to compute eigenvalues
            call dsyev(jobz, uplo, m-1, R, m-1, W, work, lwork, info)

            deallocate(work)

            !Now R contains the eigenvectors and W contains the eigenvalues

            ! Check for success
            if (info == 0) then
                print *, ' Eigenvalue calculation is complete!  '
                print *, ' ------------------------------------------ '
                print *, ' '
            else
                print *, 'Error: DSYEV did not converge, INFO =', info
                print *, ' ------------------------------------------ '
                print *, ' '
            end if

            ! conversion of  eigenvalue to singular value for DMD of the next part (Sigma^2=lambda )
            allocate(SVD_Sigma(m-1))
            do t=m-1, 1, -1
                SVD_Sigma((m-1)-t+1)=sqrt(abs(W(t)))
            end do


            !Truncated mode
            if (r_answer==1) then
                rank=r_mode
            else
                call truncated_series(m-1,SVD_Sigma,rank)
            end if

            !Truncated VT and order correction 
            allocate(SVD_VT(m-1,rank))
            do j=1, rank
                SVD_VT(:,j)=R(:,(m-1)-j+1)
            end do

            allocate(LAMBDA(m-1,m-1),u_stacked3(nc2,m-1))
            allocate(Phi_pod(nc2,m-1))
            allocate(Psi_lambda(m-1,m-1))


            ! diagonalization of the eigenvalues 1/square root (to calculate Phi_pod)
            LAMBDA(:,:)=0.0_dp
            do t=1, m-1

                if (W(t)==0) then

                    LAMBDA(t,t)= 0 

                else

                    LAMBDA(t,t)=1/sqrt(abs(W(t)))

                end if

            end do

            allocate(Inv_Sigma(rank,rank),SVD_S(rank,rank))

            SVD_S = 0.0_dp

            Inv_Sigma=0.0_dp
            do i=1, rank


                if (SVD_Sigma(i)==0) then

                    Inv_Sigma(i,i) = 0

                else

                    Inv_Sigma(i,i)=1/SVD_Sigma(i)

                end if

                ! Inv_Sigma(i,i)=1/SVD_Sigma(i)


                SVD_S(i,i) = SVD_Sigma(i)

            end do
    
            deallocate(SVD_Sigma)


            allocate(A(rank,rank),temp_one(rank,m-1),temp_two(rank,rank))

            temp_one=0.0_dp

            temp_two=0.0_dp

            ! Calculating Phi
            Psi_lambda=matmul(R,LAMBDA)
            deallocate(LAMBDA,W,R) 


            allocate(SVD_U(nc2,rank)) 
            print *, 'Calculating and writing Phi ... '

            dims2(1) = m-1
            dims2(2) = nc2
            call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

            call CPU_TIME(start_time)

            do i=1, np2
                ! Measure the start time
                call CPU_TIME(start_time2)
                start_r=((i-1)*nc2)+1
                end_r=((i-1)*nc2)+nc2

                do tt2=1, m-1
                    call read_snapshot_fd(data_folder_dir, tt2, time_phase, m, n, dt, &
                                         fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)

                    u_stacked3(:,tt2)=Stacked(start_r:end_r)*scall_p - x_mean(start_r:end_r) &
                        - x_slope(start_r:end_r)*(real(tt2,dp) - t_bar)

                end do

                write(filename, '(A, I0, A)') "X_",i

                call h5dcreate_f(group_id4, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
                call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(u_stacked3), dims2, hdferr)
                call h5dclose_f(dataset_id, hdferr)

                Phi_pod=MATMUL(u_stacked3, Psi_lambda)

                !Truncated SVD_U and order correction 
                do j=1, rank
                    SVD_U(:,j)=Phi_pod(:,(m-1)-j+1)
                end do

                ! Measure the end time

                do tt2=2, m
                    call read_snapshot_fd(data_folder_dir, tt2, time_phase, m, n, dt, &
                                         fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)

                    u_stacked3(:,tt2-1)=Stacked(start_r:end_r)*scall_p - x_mean(start_r:end_r) &
                        - x_slope(start_r:end_r)*(real(tt2,dp) - t_bar)

                end do
        
                temp_one = temp_one + matmul(transpose(SVD_U), u_stacked3)


                dims2(1) = rank
                dims2(2) = nc2
                call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

                write(filename, '(A, I0, A)') "SVD_U_",i

                call h5dcreate_f(group_id5, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
                call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, SVD_U, dims2, hdferr)
                call h5dclose_f(dataset_id, hdferr)

                call CPU_TIME(end_time2)

                ! Calculate the elapsed time
                elapsed_time = end_time2 - start_time2
                print   '(A, F10.3, A, I0, A, I0)', 'Hrs reamining (Phi_POD)', ((np2-(i-1))*elapsed_time)/3600.0_dp, &
                '__Row-wise Partition___' , i, '/' ,np2

            end do

            call CPU_TIME(end_time)

            elapsed_time = end_time - start_time
            print   '(A, F10.3)', 'A amtrix is calculated in', elapsed_time/3600.0_dp

            deallocate(Phi_pod,Psi_lambda)

            temp_two = matmul(temp_one,SVD_VT)
            A        = matmul(temp_two,Inv_Sigma)

            print *, 'SVD complete'
            print *, '------------------------------------------'

        else


            !______________________________________________________________________________________________NO-memory management mode
            !______________________________________________________________________________________________
            !______________________________________________________________________________________________
            !______________________________________________________________________________________________
            !______________________________________________________________________________________________

            nc2=n/np2

            allocate(u_stacked_full(n,m-1))

            do t=1, m-1

                write(filename, '(A, I0, A)') trim(data_folder_dir)//"Stacked_",t+time_phase,".csv"
                open(unit=2000, file=TRIM(filename), status='OLD', iostat=io_status)
                if (io_status /= 0) then
                    print *, 'Error opening file: ', filename
                    stop
                end if

                read(2000, *, iostat=io_status) u_stacked_full(:,t)

                close(2000)

            end do

            print *, "done"

            ! In-place FD conversion: position -> velocity (if enabled)
            if (fd_preprocessing_flag == 1) then
                ! Rolling buffer: snap_prev = x_{t-1}, snap_next = x_t
                snap_prev = u_stacked_full(:,1)   ! x_1
                snap_next = u_stacked_full(:,2)   ! x_2

                ! Forward diff for t=1: v_1 = (x_2 - x_1)/dt
                u_stacked_full(:,1) = (snap_next - snap_prev) / dt

                ! Central diff for t=2..m-2
                do t = 2, m-2
                    ! snap_prev = x_{t-1}, u_stacked_full(:,t+1) = x_{t+1} (original)
                    u_stacked_full(:,t) = (u_stacked_full(:,t+1) - snap_prev) / (2.0_dp * dt)
                    ! Update buffers: snap_prev <- x_t, snap_next <- x_{t+1}
                    snap_prev = snap_next
                    snap_next = u_stacked_full(:,t+1)
                end do

                ! Central diff for t=m-1: need x_m from file
                ! After loop: snap_prev = x_{m-2} (correct)
                write(filename, '(A, I0, A)') trim(data_folder_dir)//"Stacked_", m+time_phase, ".csv"
                open(unit=2000, file=TRIM(filename), status='OLD', iostat=io_status)
                if (io_status /= 0) then
                    print *, 'Error: cannot read m-th snapshot for FD: ', trim(filename)
                    stop
                end if
                read(2000, *, iostat=io_status) snap_next   ! x_m
                close(2000)
                u_stacked_full(:,m-1) = (snap_next - snap_prev) / (2.0_dp * dt)

                print *, 'FD preprocessing applied: position -> velocity (in-place)'
            end if

            ! Compute temporal mean (include all m snapshots)
            allocate(x_mean(n), Stacked(n), x_slope(n))
            x_mean = 0.0_dp
            x_slope = 0.0_dp
            t_bar = real(m + 1, dp) / 2.0_dp
            t_denom = real(m, dp) * (real(m, dp)**2 - 1.0_dp) / 12.0_dp

            if (subtract_mean_flag == 1) then
                do t = 1, m-1
                    x_mean = x_mean + u_stacked_full(:,t)
                end do
                ! Read the m-th snapshot (velocity if FD, position otherwise) for the mean
                call read_snapshot_fd(data_folder_dir, m, time_phase, m, n, dt, &
                                     fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                x_mean = x_mean + Stacked
                x_mean = x_mean / real(m, dp)
                print *, 'Temporal mean computed'
            end if

            if (detrend_flag == 1 .and. subtract_mean_flag == 1) then
                ! Need to read m-th snapshot if not already read
                if (subtract_mean_flag == 1) then
                    ! Stacked already has m-th snapshot from above
                else
                    call read_snapshot_fd(data_folder_dir, m, time_phase, m, n, dt, &
                                         fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                end if
                ! Compute linear trend slope for each DOF
                do t = 1, m-1
                    x_slope = x_slope + (real(t, dp) - t_bar) * (u_stacked_full(:,t) - x_mean)
                end do
                x_slope = x_slope + (real(m, dp) - t_bar) * (Stacked - x_mean)
                x_slope = x_slope / t_denom
                print *, 'Linear trend computed'
            end if

            deallocate(Stacked)

            ! Subtract mean and linear trend from all snapshots
            do t = 1, m-1
                u_stacked_full(:,t) = u_stacked_full(:,t) - x_mean - x_slope*(real(t,dp) - t_bar)
            end do

            !_______________________________________________________Economy SVD
            !_______________________________________________________
            !_______________________________________________________
            allocate(SVD_VT_full(m-1,m-1))                
            allocate(SVD_Sigma_full(m-1))             
            allocate(SVD_U_full(n,m-1))                 

            lwork = max(3*min(n, m-1) + max(n, m-1), 5*min(n, m-1))
            PRINT *, 'Calculated lwork = ', lwork
            allocate(work(max(1, lwork)))

            PRINT *, 'Calling dgesvd with dimensions:'
            PRINT *, 'n = ', n
            PRINT *, 'm = ', m-1
            PRINT *, 'lwork = ', lwork

            info = 0

            call dgesvd('S', 'A', n, m-1, u_stacked_full, n, SVD_Sigma_full, SVD_U_full, n, SVD_VT_full, m-1, work, lwork, info)

            if (info /= 0) then
                print *, 'Error: dgesvd returned error code ', info
                stop
            end if

            deallocate(work)

                        !Truncated mode
            if (r_answer==1) then
                rank=r_mode
            else
                call truncated_series(m-1,SVD_Sigma_full,rank)
            end if


            allocate(SVD_VT(m-1,rank))
            SVD_VT = transpose(SVD_VT_full(1:rank,:))
            deallocate(SVD_VT_full)


            allocate(SVD_Sigma(rank))   
            SVD_Sigma = SVD_Sigma_full(1:rank)
            deallocate(SVD_Sigma_full)

            ! Build diagonal SVD_S matrix (needed for G computation later)
            allocate(SVD_S(rank,rank))
            SVD_S = 0.0_dp
            do j = 1, rank
                SVD_S(j,j) = SVD_Sigma(j)
            end do

            allocate(SVD_U(n,rank))
            SVD_U = SVD_U_full(:,1:rank)
            deallocate(SVD_U_full)

            ! Compute reduced dynamics matrix A_tilde = U'*X2*V*Sigma_inv
            ! Need to re-read snapshots 2..m since dgesvd destroyed u_stacked_full
            allocate(A(rank,rank), temp_one(rank,m-1), temp_two(rank,rank))
            allocate(Inv_Sigma(rank,rank), Stacked(n))
            temp_one = 0.0_dp
            Inv_Sigma = 0.0_dp
            do j = 1, rank
                Inv_Sigma(j,j) = 1.0_dp / SVD_Sigma(j)
            end do

            ! Build temp_one = U' * X2, where X2 = snapshots 2..m (mean/trend subtracted)
            do t = 2, m
                call read_snapshot_fd(data_folder_dir, t, time_phase, m, n, dt, &
                                     fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                Stacked = Stacked - x_mean - x_slope*(real(t,dp) - t_bar)
                ! Accumulate column t-1 of temp_one = U' * x2_t
                temp_one(:,t-1) = matmul(transpose(SVD_U), Stacked)
            end do
            deallocate(Stacked)

            ! A = temp_one * V * Sigma_inv = U'*X2 * V * Sigma_inv
            temp_two = matmul(temp_one, SVD_VT)
            A = matmul(temp_two, Inv_Sigma)
            deallocate(Inv_Sigma)

            ! Re-read X1 snapshots for HDF5 storage (dgesvd destroyed u_stacked_full)
            deallocate(u_stacked_full)
            allocate(u_stacked_full(n,m-1), Stacked(n))
            do t = 1, m-1
                call read_snapshot_fd(data_folder_dir, t, time_phase, m, n, dt, &
                                     fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                u_stacked_full(:,t) = Stacked - x_mean - x_slope*(real(t,dp) - t_bar)
            end do
            deallocate(Stacked)

            ! Save X to HDF5 (partitioned by np2)
            dims2(1) = m-1
            dims2(2) = nc2
            call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
            do i = 1, np2
                start_r = ((i-1)*nc2)+1
                end_r = i*nc2
                write(filename, '(A, I0)') "X_", i
                call h5dcreate_f(group_id4, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
                call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(u_stacked_full(start_r:end_r,:)), dims2, hdferr)
                call h5dclose_f(dataset_id, hdferr)
            end do
            call h5sclose_f(dataspace_id, hdferr)

            ! Save SVD_U partitioned by np2
            dims2(1) = rank
            dims2(2) = nc2
            call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
            do i = 1, np2
                start_r = ((i-1)*nc2)+1
                end_r = i*nc2
                write(filename, '(A, I0)') "SVD_U_", i
                call h5dcreate_f(group_id5, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
                call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, SVD_U(start_r:end_r,:), dims2, hdferr)
                call h5dclose_f(dataset_id, hdferr)
            end do
            call h5sclose_f(dataspace_id, hdferr)

            deallocate(u_stacked_full)

            print *, 'SVD complete'
            print *, '------------------------------------------'


        end if

        ! ! !_________________________________________________________________________________________________________________________
        ! ! !_________________________________________________________________________________________________________________________
        ! ! !________________________________________________________________Eigen_value_problem________________Method of snapshot_end

        !     deallocate(u_stacked)

        print *, 'Truncated series: r= ', rank
        print *, '------------------------------------------'

        write(510, *) rank


        !*******************************************************Economy SVD end
        !*******************************************************

        print *, ' DMD initiated . . . . . '
        print *, '------------------------------------------'

        allocate(A_comp(rank,rank))
        allocate(Mu(rank),Mu2(rank),absMu(rank),Index_s(rank),Vr(rank,rank),Vl(rank,rank))

        allocate(Phi(nc2,rank))

        ! Convert real matrix A to complex matrix A_comp
        DO i = 1, rank
            DO j = 1, rank
                A_comp(i, j) = CMPLX(A(i, j), 0.0_dp, kind=dp)
            END DO
        END DO


        jobvr = 'V' ! Right eigenvectors are computed
        jobvl = 'N' ! Left eigenvectors are not computed

        ALLOCATE(rwork(2*rank))

        ! Query for optimal workspace size
        ! NOTE: LAPACK zgeev signature: (JOBVL,JOBVR,N,A,LDA,W, VL,LDVL, VR,LDVR, ...)
        !       arg7=VL (left eigvecs), arg9=VR (right eigvecs)
        call zgeev(jobvl, jobvr, rank, A_comp, rank, Mu2, Vl, rank, Vr, rank, cquery, -1, rwork, info)
        lwork = int(real(cquery(1), dp))
        allocate(cwork(lwork))

        ! Call LAPACK routine zgeev to compute eigenvalues and right eigenvectors
        call zgeev(jobvl, jobvr, rank, A_comp, rank, Mu2, Vl, rank, Vr, rank, cwork, lwork, rwork, info)

        absMu = real(Mu2)


        call sort_index(absMu, Index_s, reverse=.true.)


        do i=1,rank

            Mu(i) = Mu2(i)

        end do


        !_______________________________________Allocating_Vand matrix
        allocate(Vand(rank,m-1), P(rank,rank), G(rank,m-1) ,q(rank), q_diag(rank,rank), Mu_scaled(rank))

        scale_factor = 1.0_dp  ! Choose an appropriate factor based on the range of Mu

        ! Project eigenvalues onto the unit circle for SPDMD Vandermonde construction.
        ! Eigenvalues with |mu| > 1 cause the Vandermonde matrix to overflow,
        ! corrupting the P matrix (Inf/NaN) and making the ADMM optimization fail.
        ! This preserves the frequency (phase) of each mode while ensuring |mu| = 1.
        ! The original eigenvalues (Mu2) are preserved for diagnostics.
        if (project_eig_flag == 1) then
            do j = 1, rank
                if (abs(Mu(j)) > 1.0_dp) then
                    Mu(j) = Mu(j) / abs(Mu(j))
                end if
            end do
            print *, 'Eigenvalues projected onto unit circle for Vandermonde stability'
        else
            print *, 'Eigenvalue projection DISABLED (project_eig=0)'
        end if

        Vand = 0.0_dp

        Vand(:, 1) = 1.0_dp
        do i = 2, m-1
            do j = 1, rank
                Vand(j, i) = ((Mu(j))**(i-1))
            end do
        end do
        
        ! Perform your matrix operations with Vand as usual
        P = matmul(conjg(transpose(Vr)), vr) * conjg(matmul(Vand, conjg(transpose(Vand))))
        
        ! Rescale the final result
        ! P = P * (scale_factor ** 2)

        G = matmul(transpose(SVD_S),transpose(SVD_VT))

        q_diag = conjg(matmul(matmul(Vand,transpose(G)),Vr))

        do i=1,rank

            q(i) = q_diag(i,i)

        end do

        write(510, *) gama_n

        !________________________________________________________Saving_G
        dims2(1) = rank
        dims2(2) = m-1
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        

        ! Create the dataset within the "/Phi" group with double precision
        call h5dcreate_f(group_id6, "G", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, G, dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)



        !________________________________________________________Saving_q

        dims1(1) = rank
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        

        ! Create the dataset within the "/Phi" group with double precision
        call h5dcreate_f(group_id6, "q_Real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, real(q), dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)


        call h5dcreate_f(group_id6, "q_Imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, aimag(q), dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)



        !________________________________________________________Saving_p
        dims2(1) = rank
        dims2(2) = rank
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        

        ! Create the dataset within the "/Phi" group with double precision
        call h5dcreate_f(group_id6, "P_Real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, real(P), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)


        call h5dcreate_f(group_id6, "P_Imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, aimag(P), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)




        !allocate(Amp(m-1,gama_n), Amp_p(m-1,gama_n), Jval(gama_n), Jval_p(gama_n), nzval(gama_n))
                !__________________________________________________________________________------ Obtaining the matrix of amplitudes
        !call Amplitudes_sp(G, P, q, gamval, rho, Amp, Jval, Amp_p, Jval_p, nzval)


        do i=1, np2

            call CPU_TIME(start_time)

            dims2(1) = nc2 
            dims2(2) = rank
            call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
            write(filename, '(A, I0)') "SVD_U_", i
            call h5dopen_f(group_id5, filename, dataset_id, hdferr)
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, SVD_U, dims2, hdferr)            
            call h5sclose_f(dataspace_id, hdferr)

            Phi = matmul(SVD_U,Vr)

            dims2(1) = rank 
            dims2(2) = nc2
            call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
            
            write(filename, '(A, I0)') "A_Real_", i
            ! Create the dataset within the "/Phi" group with double precision
            call h5dcreate_f(group_id7, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(real(Phi)), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)

            write(filename, '(A, I0)') "A_Imag_", i
    
            call h5dcreate_f(group_id7, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(aimag(Phi)), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)
            call h5sclose_f(dataspace_id, hdferr)

            print *, ' "A" matrix partition ', i, 'calculated and saved'

            call CPU_TIME(end_time)

            elapsed_time = end_time - start_time
            print   '(A, F10.3)', 'Hrs Remaining', (np2-(i-1))*elapsed_time/3600.0_dp

        end do



        ! dims1(1) = gama_n 
        ! call h5screate_simple_f(1, dims1, dataspace_id, hdferr)

        !         ! Create the dataset within the "/Phi" group with double precision
        ! call h5dcreate_f(group_id6, "nzval", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, nzval, dims1, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)


        dims2(1) = m-1
        dims2(2) = rank
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        

        ! Create the dataset within the "/Phi" group with double precision
        call h5dcreate_f(group_id6, "C_Real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(real(Vand)), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)


        call h5dcreate_f(group_id6, "C_Imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(aimag(Vand)), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)


        dims2(1) = rank
        dims2(2) = rank
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        
        ! Create the dataset within the "/Phi" group with double precision
        call h5dcreate_f(group_id6, "Vr_Real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, real(Vr), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)


        call h5dcreate_f(group_id6, "Vr_Imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, aimag(Vr), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        

        ! Create the dataset within the "/Phi" group with double precision
        call h5dcreate_f(group_id6, "SVD_S", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, SVD_S, dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save SVD_VT with correct dimensions (m-1, rank)
        dims2(1) = m-1
        dims2(2) = rank
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        call h5dcreate_f(group_id6, "SVD_VT", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, SVD_VT, dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)


        ! Compute total energy directly from singular values (||Sigma||_F^2 = ||Sigma * VT||_F^2 since V is orthogonal)
        total_energy = sum(SVD_S**2)

        write(510, *) total_energy
        close(510)

        ! dims2(1) = gama_n
        ! dims2(2) = m-1
        ! call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        

        ! Create the dataset within the "/Phi" group with double precision
        ! call h5dcreate_f(group_id6, "B_p_Real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(real(Amp_p)), dims2, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)


        ! call h5dcreate_f(group_id6, "B_p_Imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(aimag(Amp_p)), dims2, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)
        ! call h5sclose_f(dataspace_id, hdferr)


            ! Check for success
        if (info == 0) then
            print *, ' Eigenvalue calculation is complete!  '
            print *, ' Optimized mode amplitudes are calculated!  '
            print *, ' ------------------------------------------ '
            print *, ' '
        else
            print *, 'Error: DSYEV did not converge, INFO =', info
            print *, ' ------------------------------------------ '
            print *, ' '
        end if

        deallocate(A_comp,A,temp_one,temp_two,SVD_U,rwork)


        print *, ' Saving ! '
        print *, '------------------------------------------'
        print *, ' '



        allocate(temp_three(nc2,rank),temp_four(m,rank))

        ! Ensure variables needed for DMD mode computation are allocated
        ! (They may have been deallocated in mem_man==0 path)
        if (.not. allocated(Stacked)) allocate(Stacked(n))
        if (.not. allocated(u_stacked3)) allocate(u_stacked3(nc2,m-1))
        if (.not. allocated(Inv_Sigma)) then
            allocate(Inv_Sigma(rank,rank))
            Inv_Sigma = 0.0_dp
            do j = 1, rank
                Inv_Sigma(j,j) = 1.0_dp / SVD_Sigma(j)
            end do
        end if

        dims2(1) = rank 
        dims2(2) = nc2
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)


        print *, 'Saving DMD modes'


        do i=1, np2
            ! Measure the start time
            call CPU_TIME(start_time)
            start_r=((i-1)*nc2)+1
            end_r=((i-1)*nc2)+nc2

            do tt2=2, m
                call read_snapshot_fd(data_folder_dir, tt2, time_phase, m, n, dt, &
                                     fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)

                u_stacked3(:,tt2-1)=Stacked(start_r:end_r)*scall_p - x_mean(start_r:end_r) &
                    - x_slope(start_r:end_r)*(real(tt2,dp) - t_bar)

            end do

            Phi = matmul(u_stacked3,SVD_VT)
            temp_three  = matmul(Phi, Inv_Sigma)
            Phi         = matmul(temp_three,Vr)


            write(filename, '(A, I0)') "Phi_real_", i
    
            ! Create the dataset within the "/Phi" group with double precision
            call h5dcreate_f(group_id8, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(real(Phi)), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)

            write(filename, '(A, I0)') "Phi_imag_", i
    
            call h5dcreate_f(group_id8, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(aimag(Phi)), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)

            print *, ' Mode Partition ', i, 'Saved'

            call CPU_TIME(end_time)

            elapsed_time = end_time - start_time
            print   '(A, F10.3)', 'Hrs Remaining', (np2-(i-1))*elapsed_time/3600.0_dp


        end do
        call h5sclose_f(dataspace_id, hdferr)
        deallocate(SVD_VT,temp_three,temp_four)


        dims1(1) = rank 
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)

                ! Save projected eigenvalues (used in Vandermonde)
        call h5dcreate_f(group_id, "EigenV_Real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, real(Mu), dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)


        call h5dcreate_f(group_id, "EigenV_Imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, aimag(Mu), dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)

                ! Save original (unprojected) eigenvalues for diagnostics
        call h5dcreate_f(group_id, "EigenV_Orig_Real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, real(Mu2), dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)

        call h5dcreate_f(group_id, "EigenV_Orig_Imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, aimag(Mu2), dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)

        call h5sclose_f(dataspace_id, hdferr)


        print *, 'Data saving is complete!'
        print *, '------------------------------------------'

        !*******************************************************reshaping_end
        !*******************************************************

        ! print *, 'Calculating Intial amplitudes ..... !'
        ! print *, '------------------------------------------'

        !allocate(SVD_Sigma(rank),SVD_VT_Comp(rank,rank),SVD_U_Comp(n,rank), SVD_Sigma_diag(rank,rank), Phi_pinv(rank,n))

        ! !______________________________________________________________ Calculating the initial conditions
        ! lwork = max(1,2*MIN(rank,n)+MAX(rank,n))
        
        ! allocate(work_Comp(max(1, lwork)),rwork(5*min(rank,n)),E(rank),temp_four(rank,n),Phi_in(n,rank),)

        ! PRINT *, 'Calling dgesvd with dimensions:'
        ! PRINT *, 'n = ', n
        ! PRINT *, 'm = ', rank
        ! PRINT *, 'lwork = ', lwork

        ! info = 0

        ! Phi_in = Phi

        ! call zgesvd('S', 'S', n, rank, Phi_in, n, SVD_Sigma, SVD_U_Comp, n, SVD_VT_Comp, rank, work_Comp, lwork, rwork, info)

        ! deallocate(Phi_in)

        ! ! Compute the pseudoinverse: A_pinv = VT^T * diag(1/S) * U^T

        ! SVD_Sigma_diag = 0.0_dp

        ! do i = 1, rank

        !     SVD_Sigma_diag(i,i) = 1.0_dp / SVD_Sigma(i)

        ! end do

        ! temp_four = matmul(conjg(transpose(SVD_VT_Comp)), SVD_Sigma_diag)
        ! ! A_pinv = VT^T * diag(1/S) * U^T
        ! Phi_pinv = matmul(temp_four, conjg(transpose(SVD_U_Comp)))

        ! allocate(Amp_b_complex(rank),Amp_b(rank),Amp_b_ordered(rank))

        ! write(filename, '(A, I0, A)') trim(data_folder_dir)//"Stacked_",1+time_phase,".csv"
        ! open(unit=2000, file=TRIM(filename), status='OLD')
        ! read(2000, *, iostat=io_status) Stacked
        ! close(2000)

        ! Amp_b_complex=matmul(Phi_pinv,Stacked)
        ! Amp_b = abs(Amp_b_complex)
        ! Amp_b = Amp_b/maxval(Amp_b)

        ! Ej = 0.0_dp

        ! do j=1, rank
        !     do i=1,m

        !        Ej = Ej + sum(abs(Phi(:,j)*Amp_b_complex(j)*(Mu(j)**(i-1)))**2)

        !     end do

        !     E(j) = Ej

        !     Ej = 0.0_dp

        ! end do

        ! call sort_index(E, Index_s, reverse=.true.)
        
        ! ! do i = 1, rank

        ! !     Gain(i) = log(abs(Mu(Index_s(i))))/dt
        ! !     Freq(i) = atan2(aimag(Mu(Index_s(i))), real(Mu(Index_s(i))))/(2*pi*dt)
            
        ! !     Amp_b_ordered(i) = Amp_b(Index_s(i))

        ! ! end do

        allocate(Gain(rank),Freq(rank))

        do i = 1, rank

            Gain(i) = log(abs(Mu2(i)))/dt
            Freq(i) = atan2(aimag(Mu2(i)), real(Mu2(i)))/(2*pi*dt)

            ! Amp_b_ordered(i) = Amp_b(i)
            
        end do


        call h5gcreate_f(file_id, "Processed", group_id2, hdferr)

        dims1(1) = rank 

        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        

        ! call h5dcreate_f(group_id, "Amp_not_ordered_real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, real(Amp_b_complex), dims1, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)

        ! call h5dcreate_f(group_id, "Amp_not_ordered_imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, aimag(Amp_b_complex), dims1, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)

        ! call h5dcreate_f(group_id2, "Amp_ordered", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, Amp_b_ordered, dims1, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)

        ! call h5dcreate_f(group_id2, "Amp", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, Amp_b, dims1, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)

        ! call h5dcreate_f(group_id2, "E", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, E, dims1, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)

        call h5dcreate_f(group_id2, "Gain", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, Gain, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)

        call h5dcreate_f(group_id2, "Freq", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, Freq, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)

        ! call h5dcreate_f(group_id2, "Index", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, Index_s, dims1, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        write(filename, '(A)') "x"
        dims1(1) = nx
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(file_id, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, x, dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)
        deallocate(x)

        write(filename, '(A)') "y"
        dims1(1) = ny
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(file_id, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, y, dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)
        deallocate(y)

        write(filename, '(A)') "z"
        dims1(1) = nz
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(file_id, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, z, dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)
        deallocate(z)

        ! Save temporal mean for reconstruction
        dims1(1) = n
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "x_mean", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, x_mean, dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save linear trend slope for reconstruction
        dims1(1) = n
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "x_slope", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, x_slope, dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save t_bar for reconstruction
        dims1(1) = 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "t_bar", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, [t_bar], dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save preprocessing flags for reconstruction
        dims1(1) = 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "subtract_mean_flag", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, [subtract_mean_flag], dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        dims1(1) = 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "detrend_flag", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, [detrend_flag], dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        dims1(1) = 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "project_eig_flag", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, [project_eig_flag], dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save FD preprocessing flag
        dims1(1) = 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "fd_preprocessing_flag", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, [fd_preprocessing_flag], dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save dt for reconstruction (needed for trapezoidal integration)
        dims1(1) = 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "dt", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, [dt], dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save initial position snapshot (if FD preprocessing enabled)
        if (fd_preprocessing_flag == 1) then
            dims1(1) = n
            call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
            call h5dcreate_f(group_id, "x_initial", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, x_initial, dims1, hdferr)
            CALL h5dclose_f(dataset_id, hdferr)
            call h5sclose_f(dataspace_id, hdferr)
            deallocate(x_initial)
        end if

        ! Deallocate FD workspace
        if (allocated(snap_prev)) deallocate(snap_prev)
        if (allocated(snap_next)) deallocate(snap_next)

        deallocate(x_mean, x_slope)


        CALL h5fclose_f(file_id, hdferr)
        CALL h5close_f(hdferr)


        ! !*******************************************************Saving_end
        ! !*******************************************************

        print *, 'Writing files complete !'
        print *, ' Saved in ', results_folder_dir

        print *, ' End of the program '
        print *, '------------------------------------------'
        print *, '------------------------------------------'

    else

        print *, ' End of the program '
        print *, ' Cehck the input folder and case info then try again :)'
        print *, '------------------------------------------'

    end if

end program DMD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++zgeev+++++++++++++++________________SUBROUTINES___________________________++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Subroutine Mean
subroutine truncated_series(n, sigma, r)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)          ! Double precision kind. (Caution! higher numbers might cause memory issues)
    real(kind=dp) :: sum
    real(kind=dp) :: total, tolerance                               ! Dynamic real variables
    integer :: t
    integer, intent(in) :: n
    real(kind=dp), intent(in) :: sigma(n)                           ! Input arguments
    integer, intent(out) :: r                                       ! output arguments

    r=0

    tolerance = 1.0e-4

    sum = 0.000
    do t=1, n
        sum = sum + sigma(t)
    end do
    total = sum



    sum = 0.000
    do t=1, n
        sum = sum + sigma(t)
        if (abs(sum/total - 1.0) < tolerance) then

            r=t

            exit

        end if
    end do

end subroutine truncated_series


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine to read a snapshot as position (fd_flag=0) or compute velocity via finite differences (fd_flag=1)
! snap     : output vector of size n_size (position or velocity depending on fd_flag)
! snap_w1, snap_w2 : workspace vectors of size n_size (only used when fd_flag=1)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_snapshot_fd(data_dir, k, tphase, m_total, n_size, delta_t, &
                             fd_flag, snap, snap_w1, snap_w2, ios)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    character(*), intent(in) :: data_dir
    integer, intent(in) :: k, tphase, m_total, n_size, fd_flag
    real(dp), intent(in) :: delta_t
    real(dp), intent(inout) :: snap(n_size), snap_w1(n_size), snap_w2(n_size)
    integer, intent(out) :: ios
    character(512) :: fn

    if (fd_flag == 0) then
        ! Read position snapshot directly
        write(fn, '(A, I0, A)') trim(data_dir)//"Stacked_", k+tphase, ".csv"
        open(unit=3001, file=trim(fn), status='OLD', iostat=ios)
        if (ios /= 0) return
        read(3001, *, iostat=ios) snap
        close(3001)
    else
        ! Compute velocity via finite differences
        if (k == 1) then
            ! Forward difference: v = (x_2 - x_1) / dt
            write(fn, '(A, I0, A)') trim(data_dir)//"Stacked_", 1+tphase, ".csv"
            open(unit=3001, file=trim(fn), status='OLD', iostat=ios)
            if (ios /= 0) return
            read(3001, *, iostat=ios) snap_w1
            close(3001)
            write(fn, '(A, I0, A)') trim(data_dir)//"Stacked_", 2+tphase, ".csv"
            open(unit=3001, file=trim(fn), status='OLD', iostat=ios)
            if (ios /= 0) return
            read(3001, *, iostat=ios) snap_w2
            close(3001)
            snap = (snap_w2 - snap_w1) / delta_t
        else if (k == m_total) then
            ! Backward difference: v = (x_m - x_{m-1}) / dt
            write(fn, '(A, I0, A)') trim(data_dir)//"Stacked_", m_total-1+tphase, ".csv"
            open(unit=3001, file=trim(fn), status='OLD', iostat=ios)
            if (ios /= 0) return
            read(3001, *, iostat=ios) snap_w1
            close(3001)
            write(fn, '(A, I0, A)') trim(data_dir)//"Stacked_", m_total+tphase, ".csv"
            open(unit=3001, file=trim(fn), status='OLD', iostat=ios)
            if (ios /= 0) return
            read(3001, *, iostat=ios) snap_w2
            close(3001)
            snap = (snap_w2 - snap_w1) / delta_t
        else
            ! Central difference: v = (x_{k+1} - x_{k-1}) / (2*dt)
            write(fn, '(A, I0, A)') trim(data_dir)//"Stacked_", k-1+tphase, ".csv"
            open(unit=3001, file=trim(fn), status='OLD', iostat=ios)
            if (ios /= 0) return
            read(3001, *, iostat=ios) snap_w1
            close(3001)
            write(fn, '(A, I0, A)') trim(data_dir)//"Stacked_", k+1+tphase, ".csv"
            open(unit=3001, file=trim(fn), status='OLD', iostat=ios)
            if (ios /= 0) return
            read(3001, *, iostat=ios) snap_w2
            close(3001)
            snap = (snap_w2 - snap_w1) / (2.0_dp * delta_t)
        end if
    end if
end subroutine read_snapshot_fd