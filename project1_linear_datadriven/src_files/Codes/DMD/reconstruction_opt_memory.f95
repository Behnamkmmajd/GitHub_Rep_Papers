program reconstruction
    use HDF5
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)     ! Double precision kind.

    integer(HID_T) :: file_id, file_id2, group_id, group_id2, group_id3, group_id4, dataset_id, dataspace_id, dataspace_id2, &
    group_id5, group_id6, group_id7, group_id8, group_id9, group_id10
    integer ::  hdferr, D, dif
    integer(HSIZE_T) :: dims2(2), dims1(1), dims2_phi(2)
    integer :: i, j, m, ns, nc2, r1, r3, ii, start, stopp, r4, r5, jj, np2, counter, ij, answerr, answ
    integer :: subtract_mean_flag, detrend_flag, project_eig_flag, fd_preprocessing_flag
    integer, allocatable :: Index_m(:), r2(:), rank(:), flag_arr(:)
    real(dp) :: dt, summ, L_2normx, t_bar_val

    real(dp), allocatable ::  t(:), xdmd(:,:), omega_real(:), omega_imag(:), Phi_real(:,:), &
    Phi_imag(:,:), array_2D(:,:), stacked_X(:), stacked_XDMD(:), Error_r(:,:), &
    vector(:,:), new_col1(:), new_col2(:), new_col3(:), gamma_val(:), x_mean(:), &
    x_slope(:), t_bar_arr(:), x_initial(:), vel_temp(:,:), dt_arr(:)

    complex(dp), allocatable :: omega(:), Phi(:,:), B_p(:,:), C(:,:), B(:,:), PSI(:,:)   ! Eigenvalue of A

    character(512) :: filename, filename2, filename_results, filename_results2      ! Name of the files 
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

    filename = trim(data_folder_dir3)//"input_DMD.txt"

    open(2000, file=filename, status='OLD')

    read(2000, *)
    read(2000, *)
    read(2000, *)

    read(2000, *) filename

    filename_results = trim(results_folder_dir)//"Results_"//trim(filename)//".h5"
    !filename_results2 = trim(results_folder_dir)//"Reconstructed_"//trim(filename)//".h5"



    filename2 = trim(results_folder_dir)//"input_"//trim(filename)//".txt"
    open(2000, file=filename2, status='OLD', action="readwrite")

    read(2000, *) m
    read(2000, *) dt
    read(2000, *) ns
    read(2000, *) D
    read(2000, *) 
    read(2000, *) np2
    read(2000, *) 
    read(2000, *) 
    read(2000, *) 
    read(2000, *) 
    read(2000, *) r1
    read(2000, *) r5

    print *, ' Full reconstruction for all values of Gamma [Y=1 No=0]?'
    read(*,*) answ   ! User answer

    if (answ==1) then

        r3=r5
        ! print *, ' Number of reconstruction ranks?'
        ! read(*,*) r3   ! User answer

        ! print *, ' Delta rank ?'
        ! read(*,*) dif   ! User answer

        dif = 1

        allocate(rank(r3),Index_m(r3))

        ! print *, ' Reconstruction rank vector?'
        ! read(*,*) rank   ! User answer 


        ! ! Generate the logarithmically spaced vector
        do i = 1, r3
            rank(i) = (1 + (i-1)*(dif))
        end do

    else


        print *, ' Number of gamma values, if you enter 3, then 3 indices of gamma should be given based on the recons error ?'
        read(*,*) r3   ! User answer

        allocate(rank(r3),Index_m(r5))

        print *, ' gamma index ?'
        read(*,*) rank   ! User answer

    end if

    ! rank = [5]


    ! Set dimensions
    ns = D*ns

    nc2 = ns/np2

    allocate(xdmd(nc2,m-1),B(r1,r1),gamma_val(r5))


    B = 0.0_dp

    ! Initialize HDF5
    call h5open_f(hdferr)

    ! Open HDF5 file
    call h5fopen_f(filename_results, H5F_ACC_RDWR_F, file_id, hdferr)
    !call h5fcreate_f(filename_results2, H5F_ACC_TRUNC_F, file_id2, hdferr)


    ! Open 'Processed' group
    call h5gopen_f(file_id, "Processed", group_id, hdferr)

    ! Open 'Raw' group
    call h5gopen_f(file_id, "Raw", group_id2, hdferr)
    !call h5gcreate_f(file_id2, "Raw", group_id8, hdferr)

    !call h5gcreate_f(group_id8, "xDMD_E", group_id9, hdferr)
    call h5gopen_f(group_id2, "xDMD_E", group_id10, hdferr)


    call h5gopen_f(file_id, "Standard_Decomposition", group_id5, hdferr)
    call h5gopen_f(group_id5, "A", group_id6, hdferr)
    call h5gopen_f(group_id2, "x", group_id7, hdferr)

    ! Read temporal mean for reconstruction
    allocate(x_mean(ns))
    dims1(1) = ns
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id2, "x_mean", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, x_mean, dims1, hdferr)
    call h5dclose_f(dataset_id, hdferr)
    call h5sclose_f(dataspace_id, hdferr)

    ! Read linear trend slope
    allocate(x_slope(ns))
    dims1(1) = ns
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id2, "x_slope", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, x_slope, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        print *, 'Linear trend slope loaded'
    else
        x_slope = 0.0_dp
        print *, 'No x_slope found, assuming zero (mean-only subtraction)'
    end if
    call h5sclose_f(dataspace_id, hdferr)

    ! Read t_bar
    allocate(t_bar_arr(1))
    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id2, "t_bar", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, t_bar_arr, dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        t_bar_val = t_bar_arr(1)
    else
        t_bar_val = real(m + 1, dp) / 2.0_dp
    end if
    call h5sclose_f(dataspace_id, hdferr)
    deallocate(t_bar_arr)

    ! Read preprocessing flags from HDF5 (with backward-compatible defaults)
    allocate(flag_arr(1))
    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)

    call h5dopen_f(group_id2, "subtract_mean_flag", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, flag_arr, dims1, hdferr)
        subtract_mean_flag = flag_arr(1)
        call h5dclose_f(dataset_id, hdferr)
    else
        subtract_mean_flag = 1   ! default: mean was subtracted
    end if

    call h5dopen_f(group_id2, "detrend_flag", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, flag_arr, dims1, hdferr)
        detrend_flag = flag_arr(1)
        call h5dclose_f(dataset_id, hdferr)
    else
        detrend_flag = 0         ! default: no detrending
    end if

    call h5dopen_f(group_id2, "project_eig_flag", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, flag_arr, dims1, hdferr)
        project_eig_flag = flag_arr(1)
        call h5dclose_f(dataset_id, hdferr)
    else
        project_eig_flag = 0     ! default: no projection
    end if

    call h5sclose_f(dataspace_id, hdferr)
    deallocate(flag_arr)

    ! Read FD preprocessing flag
    allocate(flag_arr(1))
    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id2, "fd_preprocessing_flag", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, flag_arr, dims1, hdferr)
        fd_preprocessing_flag = flag_arr(1)
        call h5dclose_f(dataset_id, hdferr)
    else
        fd_preprocessing_flag = 0  ! default: no FD preprocessing
    end if
    call h5sclose_f(dataspace_id, hdferr)
    deallocate(flag_arr)

    ! Read dt for trapezoidal integration (if FD)
    allocate(dt_arr(1))
    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
    call h5dopen_f(group_id2, "dt", dataset_id, hdferr)
    if (hdferr == 0) then
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, dt_arr, dims1, hdferr)
        dt = dt_arr(1)
        call h5dclose_f(dataset_id, hdferr)
    else
        dt = 0.0_dp
    end if
    call h5sclose_f(dataspace_id, hdferr)
    deallocate(dt_arr)

    ! Read initial position snapshot (if FD preprocessing enabled)
    if (fd_preprocessing_flag == 1) then
        allocate(x_initial(ns))
        dims1(1) = ns
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dopen_f(group_id2, "x_initial", dataset_id, hdferr)
        if (hdferr == 0) then
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, x_initial, dims1, hdferr)
            call h5dclose_f(dataset_id, hdferr)
            print *, 'Initial position snapshot loaded for FD integration'
        else
            print *, 'WARNING: fd_preprocessing_flag=1 but x_initial not found!'
            fd_preprocessing_flag = 0
        end if
        call h5sclose_f(dataspace_id, hdferr)
    end if

    print *, 'Preprocessing flags: subtract_mean=', subtract_mean_flag, &
             ' detrend=', detrend_flag, ' project_eig=', project_eig_flag, &
             ' fd_preprocessing=', fd_preprocessing_flag


    allocate(array_2D(r5,r1),Phi_real(r1,r5),Phi_imag(r1,r5),B_p(r1,r5))

    dims2(1) = r5
    dims2(2) = r1
    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

    call h5dopen_f(group_id5, "B_p_Real", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)

    Phi_real = transpose(array_2D(:,:))


    call h5dopen_f(group_id5, "B_p_Imag", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)

    Phi_imag = transpose(array_2D(:,:))

    call h5dclose_f(dataset_id, hdferr)
    call h5sclose_f(dataspace_id, hdferr)

    B_p = CMPLX(Phi_real, Phi_imag, kind=dp)

    deallocate(array_2D,Phi_imag,Phi_real)


    dims1(1) = r5
    call h5screate_simple_f(1, dims1, dataspace_id2, hdferr)


    write(filename, '(A)') "nzval"

    call h5dopen_f(group_id5, filename, dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, Index_m, dims1, hdferr)
    call h5dclose_f(dataset_id, hdferr)  ! Ensure closing the dataset after reading


    write(filename, '(A)') "Gamma"

    call h5dopen_f(group_id2, filename, dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, gamma_val, dims1, hdferr)
    call h5dclose_f(dataset_id, hdferr)  ! Ensure closing the dataset after reading


    allocate(array_2D(m-1,r1),Phi_real(r1,m-1),Phi_imag(r1,m-1),C(r1,m-1))

    dims2(1) = m-1
    dims2(2) = r1
    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

    call h5dopen_f(group_id5, "C_Real", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)

    Phi_real = transpose(array_2D(:,:))

    call h5dopen_f(group_id5, "C_Imag", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)

    Phi_imag = transpose(array_2D(:,:))

    call h5dclose_f(dataset_id, hdferr)
    call h5sclose_f(dataspace_id, hdferr)

    C = CMPLX(Phi_real, Phi_imag, kind=dp)

    deallocate(array_2D,Phi_imag,Phi_real)


    allocate(array_2D(m-1,nc2),stacked_X(ns))


    dims2(1) = m-1
    dims2(2) = nc2
    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

    ! !_______________________________Calculating_L-2 norm of the X matrix

    ! summ = 0.0_dp

    ! do j=1,m-1

    !     do i=1,np2

    !         write(filename, '(A, I0)') "X_",i

    !         call h5dopen_f(group_id7, filename, dataset_id, hdferr)
    !         call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)
    !         stacked_X(((i-1)*nc2)+1:(i)*nc2) = array_2D(j,:)
    !         call h5dclose_f(dataset_id, hdferr)  ! Ensure closing the dataset after reading

    !     end do

    !     summ = summ + (norm2(stacked_X))**2

    ! end do

    deallocate(array_2D,stacked_X)

    ! L_2normx = sqrt(summ)



    
    allocate(Phi(nc2,r1), Phi_real(nc2,r1), Phi_imag(nc2,r1),array_2D(r1,nc2),PSI(r1,m-1))


    dims2_phi(1) = r1
    dims2_phi(2) = nc2

    dims2(1) = m-1 
    dims2(2) = nc2
    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

    vector = reshape([0.0, 0.0, 0.0], [1, 3])
    new_col1 = [0.0, 0.0, 0.0]
    new_col2 = [0.0, 0.0, 0.0]
    new_col3 = [0.0, 0.0, 0.0]


    do ii=2, r3+1

        answerr = 0

        do ij=1, size(vector,1)

            if (Index_m(ii-1)==vector(ij,2)) then

                answerr = 1

            end if

        end do

        if (answerr==0) then

            r4=rank(ii-1)

            write(filename, '(A, I0)') "Gamma_",Index_m(r4)

            call h5gcreate_f(group_id10, filename, group_id4, hdferr)

            print *, ' Gamma', r4

            do j=1,r1

                B(j,j) = B_p(j,r4)

            end do

            PSI=matmul(B,C)

            do jj=1,np2

                write(filename, '(A, I0)') "A_Real_",jj

                call h5dopen_f(group_id6, filename, dataset_id, hdferr)
                call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2_phi, hdferr)
            
                Phi_real = transpose(array_2D(:,:))

                write(filename, '(A, I0)') "A_Imag_",jj
            
                call h5dopen_f(group_id6, filename, dataset_id, hdferr)
                call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2_phi, hdferr)
            
                Phi_imag = transpose(array_2D(:,:))
            
                call h5dclose_f(dataset_id, hdferr)
            
                Phi = CMPLX(Phi_real, Phi_imag, kind=dp)

                xdmd = matmul(matmul(Phi,B),C)

                ! Add back preprocessing (mean and/or linear trend) if they were subtracted
                if (subtract_mean_flag == 1) then
                    do j = 1, m-1
                        xdmd(:,j) = xdmd(:,j) + x_mean(((jj-1)*nc2)+1:jj*nc2)
                    end do
                end if
                if (detrend_flag == 1) then
                    do j = 1, m-1
                        xdmd(:,j) = xdmd(:,j) + x_slope(((jj-1)*nc2)+1:jj*nc2)*(real(j,dp) - t_bar_val)
                    end do
                end if

                ! FD integration: convert reconstructed velocity back to position
                if (fd_preprocessing_flag == 1) then
                    ! xdmd(:,j) now contains velocity v_j for j=1..m-1
                    ! Apply trapezoidal integration: x_k = x_0 + sum_{j=1}^{k} (v_{j-1}+v_j)/2 * dt
                    ! Allocate temporary to hold velocity before overwriting with position
                    allocate(vel_temp(nc2, m-1))
                    vel_temp = real(xdmd)

                    ! First column: x_1 = x_initial (exact initial condition)
                    xdmd(:,1) = x_initial(((jj-1)*nc2)+1:jj*nc2)

                    ! Trapezoidal integration for columns 2..m-1
                    do j = 2, m-1
                        xdmd(:,j) = xdmd(:,j-1) + (vel_temp(:,j-1) + vel_temp(:,j)) / 2.0_dp * dt
                    end do

                    deallocate(vel_temp)
                end if

                write(filename, '(A, I0)') "xDMD_",jj
                    ! Create the dataset with double precision
                call h5dcreate_f(group_id4, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
                call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(xdmd), dims2, hdferr)
                call h5dclose_f(dataset_id, hdferr)  ! Close the dataset after writing

            end do

            write(filename, '(A)') "PSI_real"

            dims2(1) = m-1 
            dims2(2) = r1
            call h5screate_simple_f(2, dims2, dataspace_id, hdferr)


            ! Create the dataset with double precision
            call h5dcreate_f(group_id4, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(real(PSI)), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)  ! Close the dataset after writing

            write(filename, '(A)') "PSI_imag"

            ! Create the dataset with double precision
            call h5dcreate_f(group_id4, filename, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(aimag(PSI)), dims2, hdferr)
            call h5dclose_f(dataset_id, hdferr)  ! Close the dataset after writing

            dims2(1) = m-1 
            dims2(2) = nc2
            call h5screate_simple_f(2, dims2, dataspace_id, hdferr)


            if (vector(1,1)==0) then

                new_col1    = [rank(ii-1)*1.0_dp]
                new_col2   = [Index_m(ii-1)*1.0_dp]
                new_col3   = [gamma_val(ii-1)*1.0_dp]

                vector = reshape([new_col1, new_col2, new_col3], [size(vector,1), 3])

            else

                new_col1   = [new_col1, rank(ii-1)*1.0_dp]
                new_col2   = [new_col2, Index_m(ii-1)*1.0_dp]
                new_col3   = [new_col3, gamma_val(ii-1)*1.0_dp]

                vector = reshape([new_col1, new_col2, new_col3], [size(vector,1) +1 , 3])

            end if

            ! Append the new row to the matrix

        end if

        B = 0.0_dp


    end do

    deallocate(array_2D,Phi_imag,Phi_real)
    deallocate(x_mean, x_slope)
    if (allocated(x_initial)) deallocate(x_initial)

    ! dims2(1) = size(vector,1)
    ! dims2(2) = size(vector,2)

    ! call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

    ! call h5dcreate_f(group_id5, "Recons_vec_two", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
    ! call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, vector, dims2, hdferr)
    ! call h5dclose_f(dataset_id, hdferr)

    ! Close HDF5 file
    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

    ! write(2000, *) size(vector,1)

    print *, 'Finish'

end program reconstruction
