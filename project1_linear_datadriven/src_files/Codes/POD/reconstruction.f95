program reconstruction
    use HDF5
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)     ! Double precision kind.

    real(dp), allocatable :: array_2D(:,:), Phi(:,:), Psi(:,:), velocity(:,:), temp(:,:), Lambda(:)

    real(dp), allocatable :: stacked_X(:), stacked_XPOD(:), Error_r(:,:)
    real(dp), allocatable :: x_mean(:), x_slope(:)

    real(dp) :: summ, L_2normx, elapsed_time, end_time, start_time, t_bar, t_bar_tmp(1)

    integer(HID_T) :: file_id, group_id, group_id2, group_id3, group_id4, group_id5, dataset_id, dataspace_id
    integer ::  hdferr
    integer(HSIZE_T) :: dims2(2), dims1(1)

    integer :: nt, ns, nv, np2, r, r2, n, nc2                     ! Number of snapshots
    integer :: subtract_mean_flag, detrend_flag, fd_preprocessing_flag
    real(dp) :: dt
    real(dp), allocatable :: x_initial(:), vel_temp(:,:)
    real(dp) :: dt_tmp(1)
    integer :: fd_tmp(1)

    integer :: io_status, i, ii, j, k, ij, start_r, end_r, start, stopp

    integer, allocatable :: Rank_vec(:)


    character(512) :: filename, filename2, filename_results                         ! Name of the files 
    character(1030) :: data_folder_dir                                              ! Where the velocity fields are stored
    character(512) :: data_folder                                                   ! The folder of the case
    character(512) :: results_folder_dir                                            ! Results folder
    character(512) :: data_folder_dir2, jobz, uplo, data_folder_dir3


    !******************************************************* Directories
    data_folder_dir2 = '../../database/'                                            ! directory of velocity field
    results_folder_dir = './results/'                                               ! Directory to save the files


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

    filename = trim(data_folder_dir3)//"input_POD.txt"

    open(2000, file=filename, status='OLD')

    read(2000, *)
    read(2000, *)
    read(2000, *)

    read(2000, *) filename

    filename_results = trim(results_folder_dir)//"Results_"//trim(filename)//".h5"


    filename2 = trim(results_folder_dir)//"input_"//trim(filename)//".txt"
    open(2000, file=filename2, status='OLD')

    read(2000, *) nt
    read(2000, *) ns
    read(2000, *) nv
    read(2000, *) 
    read(2000, *) np2
    read(2000, *)  ! nms
    read(2000, *)  ! normal
    read(2000, *)  ! fs
    read(2000, *)  ! nx
    read(2000, *)  ! ny
    read(2000, *)  ! nz
    read(2000, *)  ! scalling_P
    read(2000, *, iostat=io_status) subtract_mean_flag
    if (io_status /= 0) subtract_mean_flag = 0
    read(2000, *, iostat=io_status) detrend_flag
    if (io_status /= 0) detrend_flag = 0
    read(2000, *, iostat=io_status) fd_preprocessing_flag
    if (io_status /= 0) fd_preprocessing_flag = 0
    close(2000)
    print *, ' Mean subtraction (1=on 0=off) ', subtract_mean_flag
    print *, ' Linear detrending (1=on 0=off) ', detrend_flag
    print *, ' FD preprocessing (1=on 0=off)  ', fd_preprocessing_flag

    print *, ' Number of reconstructions?'
    read(*,*) r   ! User answer

    allocate(Rank_vec(r))

    start = 1

    stopp = r

    print *, ' Reconstruction ranks?'
    read(*,*) Rank_vec   ! User answer

    n=nv*ns

    nc2=n/np2

    allocate(Psi(nt, nt), Lambda(nt))

    ! Initialize HDF5
    call h5open_f(hdferr)
    ! Open HDF5 file
    call h5fopen_f(filename_results, H5F_ACC_RDWR_F, file_id, hdferr)
    call h5gopen_f(file_id, "Raw", group_id, hdferr)
    call h5gopen_f(group_id, "Phi_partitioned", group_id2, hdferr)
    call h5gopen_f(group_id, "xPOD_partitioned", group_id3, hdferr)
    call h5gopen_f(group_id, "X_partitioned", group_id5, hdferr)


    allocate(array_2D(nt,nt))
    dims2(1) = nt
    dims2(2) = nt
    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
    call h5dopen_f(group_id, "Psi", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)
    Psi = transpose(array_2D)
    deallocate(array_2D)
    call h5dclose_f(dataset_id, hdferr)  ! Ensure closing the dataset after reading


    dims1 = nt
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id, "Sigma", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, Lambda, dims1, hdferr)
    call h5dclose_f(dataset_id, hdferr)  ! Ensure closing the dataset after reading

    ! Read preprocessing data from HDF5
    allocate(x_mean(n), x_slope(n))
    x_mean = 0.0_dp
    x_slope = 0.0_dp
    t_bar = 0.0_dp

    ! Read x_mean
    dims1(1) = n
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id, "x_mean", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, x_mean, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
    end if
    call h5sclose_f(dataspace_id, hdferr)

    ! Read x_slope
    dims1(1) = n
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id, "x_slope", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, x_slope, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
    end if
    call h5sclose_f(dataspace_id, hdferr)

    ! Read t_bar
    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id, "t_bar", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, t_bar_tmp, dims1, hdferr)
        t_bar = t_bar_tmp(1)
        call h5dclose_f(dataset_id, hdferr)
    end if
    call h5sclose_f(dataspace_id, hdferr)

    ! Read preprocessing flags from HDF5 (overrides text file values)
    ! This ensures reconstruction uses exactly what POD computed
    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id, "subtract_mean_flag", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, fd_tmp, dims1, hdferr)
        subtract_mean_flag = fd_tmp(1)
        call h5dclose_f(dataset_id, hdferr)
    end if
    call h5sclose_f(dataspace_id, hdferr)

    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id, "detrend_flag", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, fd_tmp, dims1, hdferr)
        detrend_flag = fd_tmp(1)
        call h5dclose_f(dataset_id, hdferr)
    end if
    call h5sclose_f(dataspace_id, hdferr)

    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id, "fd_preprocessing_flag", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, fd_tmp, dims1, hdferr)
        fd_preprocessing_flag = fd_tmp(1)
        call h5dclose_f(dataset_id, hdferr)
    end if
    call h5sclose_f(dataspace_id, hdferr)

    dt = 0.0_dp
    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id, "dt", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, dt_tmp, dims1, hdferr)
        dt = dt_tmp(1)
        call h5dclose_f(dataset_id, hdferr)
    end if
    call h5sclose_f(dataspace_id, hdferr)

    allocate(x_initial(n))
    x_initial = 0.0_dp
    if (fd_preprocessing_flag == 1) then
        dims1(1) = n
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dopen_f(group_id, "x_initial", dataset_id, hdferr)
        if (hdferr == 0) then
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, x_initial, dims1, hdferr)
            call h5dclose_f(dataset_id, hdferr)
        end if
        call h5sclose_f(dataspace_id, hdferr)
        print *, 'FD preprocessing: will integrate velocity back to position'
    end if


    dims2(1) = nt
    dims2(2) = nc2
    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

    do ii=1,r

        r2 = Rank_vec(ii)

        write(filename, '(A, I0)') "xPOD_rank",r2
        call h5gcreate_f(group_id3, trim(filename), group_id4, hdferr)

        do i=1,np2

            write(filename, '(A, I0)') "Phi_",i
            allocate(array_2D(nt,nc2))

            call h5dopen_f(group_id2, filename, dataset_id, hdferr)
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)

            allocate(Phi(nc2, nt))

            Phi = transpose(array_2D)
            call h5dclose_f(dataset_id, hdferr)  ! Ensure closing the dataset after reading

            deallocate(array_2D)

            allocate(velocity(nc2,nt), temp(nc2,nt))

            velocity = 0.0_dp
            temp = 0.0_dp

            do ij = nt, nt-(r2-1), -1
                do j=1, nc2
                    do k=1, nt
                        temp(j,k) = sqrt(Lambda((nt-ij)+1)) * Phi(j,ij) * Psi(k,ij)
                    end do
                end do
                velocity = velocity + temp
            end do

            ! Re-add mean and linear trend to reconstructed velocity
            start_r = ((i-1)*nc2) + 1
            end_r   = i*nc2
            if (subtract_mean_flag == 1) then
                do k = 1, nt
                    velocity(:,k) = velocity(:,k) + x_mean(start_r:end_r) &
                        + x_slope(start_r:end_r) * (real(k, dp) - t_bar)
                end do
            end if

            ! Trapezoidal integration: velocity → position when fd_preprocessing was used
            if (fd_preprocessing_flag == 1) then
                allocate(vel_temp(nc2, nt))
                vel_temp = velocity
                velocity(:,1) = x_initial(start_r:end_r)
                do k = 2, nt
                    velocity(:,k) = velocity(:,k-1) &
                        + (vel_temp(:,k-1) + vel_temp(:,k)) / 2.0_dp * dt
                end do
                deallocate(vel_temp)
            end if

            ! Check if the dataset "Velocity_" exists, if it does, delete it
            write(filename, '(A, I0)') "Velocity_",i

            ! Create the dataset with double precision
            call h5dcreate_f(group_id4, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(velocity), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)  ! Close the dataset after writing

            deallocate(velocity,temp,Phi)

            print *, "Partition", i, "Rank", r2, "Reconstrcuted"

        end do

        call h5gclose_f(group_id4, hdferr)


        print *, "Rank", r2

    end do
        
    print *, 'Reconstruction is done'

    deallocate(Psi,Lambda,x_mean,x_slope)
    if (allocated(x_initial)) deallocate(x_initial)

    ! Close HDF5 file
    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

end program reconstruction
