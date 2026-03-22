program pod
    use FFTW_sub_module                                                                                 ! Loading the Module for calculating Pis_m
    use omp_lib                                                                                         ! Include OpenMP library
    use stdlib_sorting, only: sort_index
    use HDF5

    implicit none

    integer, parameter ::       dp = selected_real_kind(15, 307)                                        ! Double precision


    integer(HSIZE_T) ::         dims2(2),                                                           &
                                dims1(1)

    integer(HID_T) ::           file_id,                                                            &
                                group_id,                                                           &
                                group_id2,                                                          &
                                group_id3,                                                          &
                                group_id4,                                                          &
                                dataset_id,                                                         &
                                dataspace_id

    integer ::                  hdferr
!_________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________
    integer ::                  nt,                                                                 &   ! Number of snapshots
                                time_phase,                                                         &   ! t0
                                ns,                                                                 &   ! Number of grid points
                                nv,                                                                 &   ! 2D or 3D problem
                                np,                                                                 &   ! Number of column-wise partition
                                np2,                                                                &   ! Number of row-wise partition
                                nms,                                                                &   ! Rank
                                n,                                                                  &   ! Number of rows in Full stacked matrix
                                m,                                                                  &   ! Number of columns in Full stacked matrix
                                nc,                                                                 &   ! Number of columns in partitioned matrix m/np
                                nc2,                                                                &   ! Number of columns in partitioned matrix n/np2
                                nx,ny,nz,                                                           &   ! Number of columns in partitioned matrix n/np2
                                subtract_mean_flag,                                                 &   ! Mean subtraction flag (1=on 0=off)
                                detrend_flag,                                                       &   ! Linear detrending flag (1=on 0=off)
                                fd_preprocessing_flag,                                              &   ! FD preprocessing: position -> velocity (1=on 0=off)
                                ! Dynamic integer variables (loop counters, etc)
                                info,                                                               &
                                lwork,                                                              &
                                io_status,                                                          &
                                t, j, i, ii,                                                        &
                                counter, tt2,                                                       & 
                                start_c, end_c, start_r,                                            &
                                end_r, s_c, e_c, s_c2,                                              &
                                e_c2, cc, status
    
    integer, allocatable ::     index_m(:)
                                
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________
    
    real(dp), allocatable ::    W(:),                                                               &   ! Output of eiv function (eigenvalues)
                                work(:),                                                            &   ! Workspace for SVD
                                Stacked(:),                                                         &   ! Matrix of stacked velocity
                                u_stacked1(:,:),                                                    &   ! Matrix of stacked velocity
                                u_stacked2(:,:),                                                    &   ! Matrix of stacked velocity
                                u_stacked3(:,:),                                                    &   ! Matrix of stacked velocity
                                R(:,:),                                                             &   ! Eigenvectors
                                Phi(:,:),                                                           &   ! Spatial modes
                                Psi_lambda(:,:),                                                    &   ! Stacked mode shapes
                                LAMBDA(:,:),                                                        &   ! Mode Energy                                                                                                  
                                Freqs(:),                                                           &   ! Time vector (It is actually the FFTW of time vector)
                                x(:),y(:),z(:),                                                     &
                                x_mean(:),                                                          &   ! Temporal mean of snapshots
                                x_slope(:),                                                         &   ! Linear trend slope for each DOF
                                x_initial(:),                                                       &   ! Initial position snapshot (for FD reconstruction)
                                snap_prev(:), snap_next(:),                                         &   ! Workspace for FD subroutine
                                Re_error(:)
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________

    real(dp) ::                 normal,                                                             &   ! Nondimensionalzing the Freq
                                fs,                                                                 &
                                query(1),                                                           &
                                start_time, end_time,                                               &
                                start_time2, end_time2,                                             &
                                elapsed_time,                                                       &   ! CPU time variables
                                summm,                                                              &
                                orgi,                                                               &
                                scalling_P,                                                         &
                                t_bar, t_denom                                                          ! Mean time index and denominator for regression
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________

    character(512) ::           filename, filename_results,                                                           &
                                filename2,                                                          &   ! Name of the files 
                                data_folder_dir,                                                    &   ! Where the velocity fields are stored
                                data_folder,                                                        &   ! The folder of the case
                                results_folder_dir,                                                 &   ! Results folder
                                data_folder_dir2, jobz, uplo, data_folder_dir3
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________
!___________________________________________________________________________________________________________________________________________________________________________


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

    write(filename, '(A)') trim(data_folder_dir3)//"input_POD.txt"
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
    open(unit=503, file=filename2)

    print *, ' Case :                                       ', filename
    read(2000, *) nt
    write(503, *) nt
    print *, ' Number of snapshots :                        ', nt
    print *, '______________________'
    read(2000, *) time_phase
    write(filename, '(A,I0,A)') 'Stacked_', 1+time_phase, '.csv'
    print *, ' The name of the first snapshot is :          ', trim(filename)
    print *, '______________________'
    read(2000, *) ns
    print *, ' Number of grid points :                      ', ns
    write(503, *) ns
    print *, '______________________'
    read(2000, *) nv
    print *, ' 2D or 3D case :                              ', nv,"D"
    write(503, *) nv
    print *, '______________________'
    read(2000, *) np
    print *, ' Number of column-wise partition :            ', np
    write(503, *) np
    print *, '______________________'
    read(2000, *) np2
    print *, ' Number of row-wise partition :               ', np2
    write(503, *) np2
    print *, '______________________'
    read(2000, *) nms
    print *, ' Rank(Number of modes to be saved) :          ', nms
    write(503, *) nms
    print *, '______________________'
    read(2000, *) normal
    print *, ' Nondimensionalized factor fn :               ', normal
    write(503, *) normal

    print *, '______________________'
    read(2000, *) fs
    print *, ' Sampling frequency :                         ', fs
    write(503, *) fs
    print *, '______________________'
    

    read(2000, *) nx
    write(503, *) nx

    read(2000, *) ny
    write(503, *) ny

    read(2000, *) nz
    write(503, *) nz

    read(2000, *) scalling_P
    print *, ' Nondimensionalized factor based on bubble Diameter :               ', scalling_P
    write(503, *) scalling_P

    read(2000, *, iostat=io_status) subtract_mean_flag
    if (io_status /= 0) subtract_mean_flag = 0   ! Default: mean subtraction OFF
    print *, ' Mean subtraction (1=on 0=off) ', subtract_mean_flag
    write(503, *) subtract_mean_flag

    read(2000, *, iostat=io_status) detrend_flag
    if (io_status /= 0) detrend_flag = 0          ! Default: linear detrending OFF
    print *, ' Linear detrending (1=on 0=off) ', detrend_flag
    write(503, *) detrend_flag

    read(2000, *, iostat=io_status) fd_preprocessing_flag
    if (io_status /= 0) fd_preprocessing_flag = 0  ! Default: FD preprocessing OFF
    print *, ' FD preprocessing pos->vel (1=on 0=off) ', fd_preprocessing_flag
    write(503, *) fd_preprocessing_flag

    close(2000)
    close(503)

    

    print *, ' Is everything fine with the above case (Yes=1 / No=0 )?'
    read(*,*) info                                                                                      ! User answer

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
        nc=m/np
        nc2=n/np2

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
            x_initial = x_initial / scalling_P
            print *, 'FD preprocessing: initial position saved for reconstruction'
        end if

        
        !_______________________________________________________
        !_______________________________________________________Allocation of the variables
        allocate(u_stacked1(n, nc),u_stacked2(n,nc),Stacked(n))

        ! Compute temporal mean and linear trend of all m snapshots
        allocate(x_mean(n), x_slope(n))
        x_mean = 0.0_dp
        x_slope = 0.0_dp
        t_bar = real(m + 1, dp) / 2.0_dp
        t_denom = real(m, dp) * (real(m, dp)**2 - 1.0_dp) / 12.0_dp

        if (subtract_mean_flag == 1) then
            do t = 1, m
                call read_snapshot_fd(data_folder_dir, t, time_phase, m, n, &
                    1.0_dp/fs, fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                if (io_status /= 0) then
                    print *, 'Error reading snapshot for mean computation, t=', t
                    stop
                end if
                x_mean = x_mean + Stacked / scalling_P
            end do
            x_mean = x_mean / real(m, dp)
            print *, 'Temporal mean computed'
        end if

        if (detrend_flag == 1 .and. subtract_mean_flag == 1) then
            ! Second pass for slope: b_i = sum_k (k - t_bar) * (x_i(k) - x_mean_i) / t_denom
            do t = 1, m
                call read_snapshot_fd(data_folder_dir, t, time_phase, m, n, &
                    1.0_dp/fs, fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                if (io_status /= 0) then
                    print *, 'Error reading snapshot for slope computation, t=', t
                    stop
                end if
                x_slope = x_slope + (real(t, dp) - t_bar) * (Stacked / scalling_P - x_mean)
            end do
            x_slope = x_slope / t_denom
            print *, 'Linear trend computed'
        end if

        print *, 'Reading Velocity field'
        print *, '------------------------------------------'


        ! Measure the start time
        call CPU_TIME(start_time2)


        allocate(R(m,m))

        R = 0.0_dp

        counter=0

        do counter=1, np

            s_c=(nc*(counter-1))+1
            e_c=counter*nc

            cc=1
            ! Measure the start time
            call CPU_TIME(start_time)

            do i=s_c, e_c

                call read_snapshot_fd(data_folder_dir, i, time_phase, m, n, &
                    1.0_dp/fs, fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                if (io_status /= 0) then
                    print *, 'Error reading snapshot for R-matrix, i=', i
                    stop
                end if

                u_stacked1(:,cc)=Stacked/scalling_P - x_mean - x_slope*(real(i,dp) - t_bar)
                cc=cc+1
            end do

            ! Measure the end time
            call CPU_TIME(end_time)
            ! Calculate the elapsed time
            elapsed_time = end_time - start_time
            print *, 'Elapsed time for reading one stacked: ', elapsed_time, ' seconds'
            ! Measure the start time
            call CPU_TIME(start_time)

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
                        call read_snapshot_fd(data_folder_dir, ii, time_phase, m, n, &
                            1.0_dp/fs, fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                        if (io_status /= 0) then
                            print *, 'Error reading snapshot for R-matrix, ii=', ii
                            stop
                        end if
                        u_stacked2(:,cc)=Stacked/scalling_P - x_mean - x_slope*(real(ii,dp) - t_bar)
                        cc=cc+1
                    end do
                    R(start_r:end_r,start_c:end_c)=matmul(transpose(u_stacked1),u_stacked2)
                end if
            end do
            
            ! Measure the end time
            call CPU_TIME(end_time)
            ! Calculate the elapsed time
            elapsed_time = end_time - start_time
            print *, 'Elapsed time', elapsed_time, ' seconds', 'For partition', counter
        end do

        deallocate(u_stacked1,u_stacked2)

        ! Completing the other half of the matrix since R is symmetry 
        do i=1, m
            do j=i+1, m
                R(j,i)=R(i,j)
            end do
        end do

            ! Measure the end time
        call CPU_TIME(end_time2)
            ! Calculate the elapsed time
        elapsed_time = end_time2 - start_time2
        print *, 'Elapsed time for calculating R: ', elapsed_time, ' seconds'
        !_________________________________________________________________________________________________________________________
        !_________________________________________________________________________________________________________________________
        !___________________________________________________________________________________________________Covariance matrix R end

        print *, ' R is calculated '
        print *, ' ------------------------------------------ '

        allocate(W(m))
        
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
        call dsyev(jobz, uplo, m, R, m, W, query, -1, info)
        lwork = int(query(1))
        allocate(work(lwork))

        !Call LAPACK routine DSYEV to compute eigenvalues
        call dsyev(jobz, uplo, m, R, m, W, work, lwork, info)

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

        allocate(LAMBDA(m,m),u_stacked3(nc2,m))
        allocate(Phi(nc2,m))
        allocate(Psi_lambda(m,m),Re_error(m))


        ! diagonalization of the eigenvalues 1/square root
        LAMBDA(:,:)=0.0_dp
        do t=1, m
            LAMBDA(t,t)=1/sqrt(abs(W(t)))
        end do

        summm = 0.0_dp

        do t=1, m

            summm = summm + sqrt(abs(W(t)))**2

        end do

        orgi = summm

        summm = 0.0_dp

        do t=1, m

            do i=t+1, m

                summm = summm + sqrt(abs(W(m-i+1)))**2

            end do

            Re_error(t) = (summm / orgi)*100.0_dp

            summm = 0.0_dp

        end do

        ! Calculating Phi
        Psi_lambda=matmul(R,LAMBDA)

        deallocate(LAMBDA)

    
        print *, 'Calculating and writing Phi ... '


        ! Initialize HDF5 library
        call h5open_f(hdferr)
        ! Create a new HDF5 file
        call h5fcreate_f(filename_results, H5F_ACC_TRUNC_F, file_id, hdferr)
        call h5gcreate_f(file_id, "Raw", group_id, hdferr)
        call h5gcreate_f(group_id, "Phi_partitioned", group_id2, hdferr)
        call h5gcreate_f(group_id, "X_partitioned", group_id3, hdferr)
        call h5gcreate_f(group_id, "xPOD_partitioned", group_id4, hdferr)

        dims2(1) = m 
        dims2(2) = nc2
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

        
        do i=1, np2
            ! Measure the start time
            call CPU_TIME(start_time2)
            start_r=((i-1)*nc2)
            end_r=((i-1)*nc2)+nc2

            do tt2=1, m

                call read_snapshot_fd(data_folder_dir, tt2, time_phase, m, n, &
                    1.0_dp/fs, fd_preprocessing_flag, Stacked, snap_prev, snap_next, io_status)
                if (io_status /= 0) then
                    print *, 'Error reading snapshot for Phi, tt2=', tt2
                    stop
                end if

                u_stacked3(:,tt2)=Stacked(start_r:end_r)/scalling_P &
                    - x_mean(start_r:end_r) - x_slope(start_r:end_r)*(real(tt2,dp) - t_bar)

                close(2000)

            end do

            Phi=MATMUL(u_stacked3, Psi_lambda)

            write(filename, '(A, I0, A)') "Phi_",i

            call h5dcreate_f(group_id2, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(Phi), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)

            ! Restore original (un-preprocessed) data for saving to X_partitioned
            ! Mean subtraction / detrending is internal to POD; stored snapshots must be original
            do tt2 = 1, m
                u_stacked3(:,tt2) = u_stacked3(:,tt2) + x_mean(start_r:end_r) &
                    + x_slope(start_r:end_r)*(real(tt2,dp) - t_bar)
            end do

            write(filename, '(A, I0, A)') "X_",i

            call h5dcreate_f(group_id3, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(u_stacked3), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)

            call CPU_TIME(end_time2)

            ! Calculate the elapsed time
            elapsed_time = end_time2 - start_time2
            print *, 'Elapsed time for calculating Phi',i, elapsed_time, ' seconds'

        end do
        call h5sclose_f(dataspace_id, hdferr)


        deallocate(u_stacked3,Stacked,Phi,Psi_lambda)

        
        !_______________________________________________________Saving


        print *, ' Storing the data ! '
        print *, '------------------------------------------'
        print *, ' '


        write(filename, '(A)') "Psi"

        dims2(1) = m 
        dims2(2) = m
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        call h5dcreate_f(group_id, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(R), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)

        call FFTW_sub(R,m)


        write(filename, '(A)') "FFT_Psi"
        call h5dcreate_f(file_id, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(abs(R)), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)


        allocate(Freqs(m))
        do i = 1, m
            Freqs(i) = (i-1) / real(m, dp) * fs
        end do

        write(filename, '(A)') "Freq"
        dims1(1) = nt 
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, Freqs, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)
        deallocate(Freqs)

        allocate(Index_m(m))


        call sort_index(W, Index_m, reverse=.true.)


        write(filename, '(A)') "Sigma"
        dims1(1) = m 
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, W, dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)
        deallocate(W)

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

        write(filename, '(A)') "Reconstruction_error"
        dims1(1) = m
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(file_id, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, Re_error, dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)
        deallocate(Re_error)

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

        ! Save FD preprocessing flag
        dims1(1) = 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "fd_preprocessing_flag", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, [fd_preprocessing_flag], dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save dt for time-integration in reconstruction
        dims1(1) = 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "dt", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, [1.0_dp / fs], dims1, hdferr)
        CALL h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save initial position for trapezoidal integration
        if (fd_preprocessing_flag == 1) then
            dims1(1) = n
            call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
            call h5dcreate_f(group_id, "x_initial", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, x_initial, dims1, hdferr)
            CALL h5dclose_f(dataset_id, hdferr)
            call h5sclose_f(dataspace_id, hdferr)
        end if

        deallocate(x_mean, x_slope)
        if (allocated(x_initial)) deallocate(x_initial)
        if (allocated(snap_prev)) deallocate(snap_prev)
        if (allocated(snap_next)) deallocate(snap_next)


        CALL h5fclose_f(file_id, hdferr)
        CALL h5close_f(hdferr)



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

end program pod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! read_snapshot_fd – read one snapshot; optionally convert position → velocity via 2nd-order FD
! fd_flag = 0 : read the position snapshot directly (original behaviour)
! fd_flag = 1 : compute velocity by finite differences of neighbouring position snapshots
! snap       : returned data vector of size n_size (position when fd_flag=0, velocity when fd_flag=1)
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