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
    dt, rho, scale_factor, rho_min_user, rho_max_user                                                   ! Sampling frequency

    integer(HID_T) :: file_id, group_id, group_id2, group_id3, dataset_id, dataspace_id, group_id4, group_id5, &
    group_id6, group_id7, group_id8
    integer ::  hdferr

    integer(HSIZE_T) :: dims2(2), dims1(1)


    integer :: info, lwork, io_status, t, j, i, ii, counter, tt2, &                 ! Dynamic integer variables (loop counters, etc)
    start_c, end_c, start_r, end_r, s_c, e_c, s_c2, e_c2, cc, r_answer, &
    mem_man
    
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
    real(dp), allocatable :: gamval(:), Jval(:), Jval_p(:), rel_error(:)
    real(dp), allocatable :: array_2D(:,:), arr_real(:,:), arr_imag(:,:), array_1D(:), vec_real(:), vec_imag(:)


    real(dp) :: query(1), Ej, F_norm_actual

    real(dp) :: start_time, end_time, start_time2, end_time2, elapsed_time, L_start, L_stop, delta         ! CPU time variables

    character(1) :: jobvl, jobvr

    character(512) :: filename, filename2, filename_results, filename_input        ! Name of the files 
    character(1030) :: data_folder_dir                                             ! Where the velocity fields are stored
    character(512) :: data_folder                                                  ! The folder of the case
    character(512) :: results_folder_dir                                           ! Results folder
    character(512) :: data_folder_dir2, jobz, uplo, data_folder_dir3
    character(256) :: line_buf                                                     ! Buffer for auto-detecting input format


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
    print *, 'Amplitude calculation'
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
    filename_input = trim(results_folder_dir)//"input_"//trim(filename)//".txt"


    open(510, file=TRIM(filename_input), status='OLD', iostat=io_status)
    if (io_status /= 0) then
        print *, 'Error opening file: ', filename_input
        stop
    end if


    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *)
    read(510, *) F_norm_actual



    print *, ' Case :                                       ', filename
    read(2000, *) nt
    print *, ' Number of snapshots :                        ', nt
    print *, '______________________'
    read(2000, *) dt
    print *, ' Sampling Frequency :                         ', 1.0/dt, 'Hz'
    print *, '______________________'
    read(2000, *) time_phase
    write(filename, '(A,I0,A)') 'Stacked_', 1+time_phase, '.csv'
    print *, ' The name of the first snapshot is :          ', trim(filename)
    print *, '______________________'
    read(2000, *) ns
    print *, ' Number of grid points :                      ', ns
    print *, '______________________'
    read(2000, *) nv
    print *, ' 2D or 3D case :                              ', nv,"D"
    print *, '______________________'
    read(2000, *) np
    print *, ' Number of column-wise partition :            ', np
    print *, '______________________'
    read(2000, *) np2
    print *, ' Number of row-wise partition :               ', np2
    print *, '______________________'
    read(2000, *) r_mode
    read(2000, *) r_answer
    if (r_answer==1) then
        print *, ' Truncated mode disabled, User Rank :          ', r_mode
    else
        print *, ' Truncated mode enabled '
    end if
    print *, '______________________'

    read(2000, *) nx

    read(2000, *) ny

    read(2000, *) nz

    read(2000, *) gama_n

    print *, ' Number of Gamma value ', gama_n

    read(2000, *) L_start 

    print *, ' Lower bound ', L_start

    read(2000, *) L_stop 

    print *, ' Upper bound ', L_stop

    read(2000, *) rho

    print *, ' rho ', rho

    ! Auto-detect: does the file have rho_min_user/rho_max_user lines?
    ! Read next line as text. If it contains '.', it's a real (rho_min_user).
    ! Otherwise it's an integer (mem_man).
    rho_min_user = 1.0e-4_dp
    rho_max_user = 1.0e6_dp

    read(2000, '(A)', iostat=io_status) line_buf
    if (io_status == 0) then
        if (index(line_buf, '.') > 0) then
            ! It's a real -> rho_min_user
            read(line_buf, *) rho_min_user
            print *, ' rho_min_user ', rho_min_user
            ! Next line should be rho_max_user
            read(2000, *, iostat=io_status) rho_max_user
            if (io_status /= 0) rho_max_user = 1.0e6_dp
            print *, ' rho_max_user ', rho_max_user
            ! Now read mem_man
            read(2000, *) mem_man
        else
            ! It's an integer -> mem_man directly
            read(line_buf, *) mem_man
        end if
    end if

    print *, ' Memory management (1=on 0=off) ', mem_man

    close(2000)

    print *, ' Is everything fine with the above case (Yes=1 / No=0 )?'
    read(*,*) info  ! User answer

    rank=r_mode

    if (info==1) then

        n=nv*ns
        m=nt

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

        ! Open HDF5 file
        call h5fopen_f(filename_results, H5F_ACC_RDWR_F, file_id, hdferr)

        call h5gopen_f(file_id, "Raw", group_id, hdferr)
    
        call h5gopen_f(file_id, "Standard_Decomposition", group_id5, hdferr)



        allocate(P(rank,rank), G(rank,m-1) ,q(rank), rel_error(gama_n))




        allocate(array_2D(rank,rank),arr_real(rank,rank),arr_imag(rank,rank))

        dims2(1) = rank
        dims2(2) = rank
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
    
        call h5dopen_f(group_id5, "P_Real", dataset_id, hdferr)
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)
    
        arr_real = array_2D(:,:)
    
    
        call h5dopen_f(group_id5, "P_Imag", dataset_id, hdferr)
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)
    
        arr_imag = array_2D(:,:)
    
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)
    
        P = CMPLX(arr_real, arr_imag, kind=dp)
    
        deallocate(array_2D,arr_real,arr_imag)





        allocate(array_2D(rank,m-1))

        dims2(1) = rank
        dims2(2) = m-1
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
    
        call h5dopen_f(group_id5, "G", dataset_id, hdferr)
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)
    
        G = array_2D(:,:)
    
        deallocate(array_2D)





        allocate(array_1D(rank), vec_real(rank), vec_imag(rank))

        dims1(1) = rank
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    
        call h5dopen_f(group_id5, "q_Real", dataset_id, hdferr)
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_1D, dims1, hdferr)
    
        vec_real = array_1D(:)
    
    
        call h5dopen_f(group_id5, "q_Imag", dataset_id, hdferr)
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_1D, dims1, hdferr)
    
        vec_imag = array_1D(:)
    
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)
    
        q = CMPLX(vec_real, vec_imag, kind=dp)
    
        deallocate(array_1D,vec_real,vec_imag)


        allocate(Amp(rank,gama_n), Amp_p(rank,gama_n), Jval(gama_n), Jval_p(gama_n), nzval(gama_n))


        allocate(gamval(gama_n))

        ! ! ! Start and stop values in logarithmic space
        ! L_start = log10(1.0_dp)  ! log10(10^1)
        ! L_stop = log10(10.0_dp)  ! log10(1600)

        ! ! Delta for evenly spaced points in the logarithmic scale
        ! delta = (L_stop - L_start) / (gama_n - 1.0_dp)

        ! do i = 1, gama_n
        !     gamval(i) = 10.0_dp**(L_start + (i - 1.0_dp) * delta)
        ! end do

        delta = (L_stop - L_start) / (gama_n)


        ! Generate the logarithmically spaced vector
        do i = 1, gama_n
            gamval(i) = (L_start + (i - 1) * delta)
        end do

        !gamval = [0.0_dp]

        dims1(1) = gama_n 
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)


        ! Suppress HDF5 error printing for check-if-exists pattern
        call h5eset_auto_f(0, hdferr)

        call h5dopen_f(group_id, "Gamma", dataset_id, hdferr)
        if (hdferr == 0) then
            call h5dclose_f(dataset_id, hdferr)
            call h5ldelete_f(group_id, "Gamma", hdferr)
        end if

        ! Re-enable HDF5 error printing
        call h5eset_auto_f(1, hdferr)

        call h5dcreate_f(group_id, "Gamma", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, gamval, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)


                !__________________________________________________________________________------ Obtaining the matrix of amplitudes
        print *, 'Starting sparse amplitude optimization (', gama_n, ' gamma values, rank=', rank, ')'
        print *, '------------------------------------------'
        call Amplitudes_sp(G, P, q, gamval, rho, Amp, Jval, Amp_p, Jval_p, nzval, &
                          rho_min_user, rho_max_user)


        do i=1, gama_n

            rel_error(i) = (Jval_p(i)/F_norm_actual)*100

        end do

        ! dims1(1) = rank
        ! call h5screate_simple_f(1, dims1, dataspace_id, hdferr)



        ! call h5dopen_f(group_id5, "Rel_error", dataset_id, hdferr)
        ! if (hdferr == 0) then
        !     ! If hdferr is 0, the dataset exists, so close and delete it
        !     call h5dclose_f(dataset_id, hdferr)
        !     call h5ldelete_f(group_id5, "Rel_error", hdferr)   ! Delete the existing dataset
        ! end if

        !         ! Create the dataset within the "/Phi" group with double precision
        ! call h5dcreate_f(group_id5, "Rel_error", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        ! call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, rel_error, dims1, hdferr)
        ! call h5dclose_f(dataset_id, hdferr)



        dims1(1) = gama_n 
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)



        call h5eset_auto_f(0, hdferr)

        call h5dopen_f(group_id5, "nzval", dataset_id, hdferr)
        if (hdferr == 0) then
            call h5dclose_f(dataset_id, hdferr)
            call h5ldelete_f(group_id5, "nzval", hdferr)
        end if

        call h5eset_auto_f(1, hdferr)

                ! Create the dataset within the "/Phi" group with double precision
        call h5dcreate_f(group_id5, "nzval", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, nzval, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        
        
    	print *, ' Non zero values', nzval


        dims2(1) = gama_n
        dims2(2) = rank
        call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
        


        call h5eset_auto_f(0, hdferr)

        call h5dopen_f(group_id5, "B_p_Real", dataset_id, hdferr)
        if (hdferr == 0) then
            call h5dclose_f(dataset_id, hdferr)
            call h5ldelete_f(group_id5, "B_p_Real", hdferr)
        end if

        call h5eset_auto_f(1, hdferr)
        ! Create the dataset within the "/Phi" group with double precision
        call h5dcreate_f(group_id5, "B_p_Real", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(real(Amp_p)), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)



        call h5eset_auto_f(0, hdferr)

        call h5dopen_f(group_id5, "B_p_Imag", dataset_id, hdferr)
        if (hdferr == 0) then
            call h5dclose_f(dataset_id, hdferr)
            call h5ldelete_f(group_id5, "B_p_Imag", hdferr)
        end if

        call h5eset_auto_f(1, hdferr)

        call h5dcreate_f(group_id5, "B_p_Imag", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(aimag(Amp_p)), dims2, hdferr)
        call h5dclose_f(dataset_id, hdferr)


        
        call h5sclose_f(dataspace_id, hdferr)


        print *, ' Saving ! '
        print *, '------------------------------------------'
        print *, ' '



        print *, 'Data saving is complete!'
        print *, '------------------------------------------'


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
