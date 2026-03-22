program reconstruction
    use HDF5
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)     ! Double precision kind.

    integer(HID_T) :: file_id, file_id2, group_id, group_id2, group_id3, group_id4, dataset_id, dataspace_id, dataspace_id2, &
    group_id5, group_id6, group_id7, group_id8, group_id9
    integer ::  hdferr, D, dif
    integer(HSIZE_T) :: dims2(2), dims1(1), dims2_phi(2)
    integer :: i, j, m, ns, nc2, r1, r3, ii, start, stopp, r4, r5, jj, np2, counter, ij, answerr
    integer, allocatable :: Index_m(:), r2(:), rank(:)
    real(dp) :: dt, summ, L_2normx, F_norm_orgi, F_norm_XDMD

    real(dp), allocatable ::  t(:), xdmd(:,:), omega_real(:), omega_imag(:), Phi_real(:,:), &
    Phi_imag(:,:), array_2D(:,:), stacked_X(:), stacked_XDMD(:), Error_r(:,:), &
    vector(:,:), new_col1(:), new_col2(:), new_col3(:), gamma_val(:), SVD_S(:,:), &
    new_col4(:), G_matrix(:,:)

    complex(dp), allocatable :: omega(:), Phi(:,:), B_p(:,:), C(:,:), B(:,:), Vr(:,:)   ! Eigenvalue of A

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
    filename_results2 = trim(results_folder_dir)//"Reconstructed_"//trim(filename)//".h5"


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

    ! rank = [5]


    ! Set dimensions
    ns = D*ns

    nc2 = ns/np2

    allocate(xdmd(r1,m-1),B(r1,r1),gamma_val(r3))


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

    call h5gopen_f(file_id, "Standard_Decomposition", group_id5, hdferr)
    call h5gopen_f(group_id5, "A", group_id6, hdferr)
    call h5gopen_f(group_id2, "x", group_id7, hdferr)


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


    ! --- Read C (Vandermonde) with dims (m-1, r1) ---
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

    ! --- Read Vr and SVD_S with correct dims (r1, r1) ---
    allocate(array_2D(r1,r1),Phi_real(r1,r1),Phi_imag(r1,r1),Vr(r1,r1),SVD_S(r1,r1))

    dims2(1) = r1
    dims2(2) = r1
    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

    call h5dopen_f(group_id5, "Vr_Real", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)

    Phi_real = array_2D(:,:)

    call h5dopen_f(group_id5, "Vr_Imag", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)

    Phi_imag = array_2D(:,:)

    call h5dclose_f(dataset_id, hdferr)

    Vr = CMPLX(Phi_real, Phi_imag, kind=dp)


    call h5dopen_f(group_id5, "SVD_S", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)

    SVD_S = array_2D(:,:)

    call h5dclose_f(dataset_id, hdferr)
    call h5sclose_f(dataspace_id, hdferr)

    deallocate(array_2D,Phi_imag,Phi_real)


    ! Compute F_norm_orgi from singular values: ||Sigma * VT||_F^2 = ||Sigma||_F^2 = sum(sigma_i^2)
    ! (because V is orthogonal, so ||Sigma * VT||_F = ||Sigma||_F)
    F_norm_orgi = sum(SVD_S**2)

    ! --- Read G matrix (reduced-space data) for true error computation ---
    allocate(G_matrix(r1, m-1), array_2D(r1, m-1))
    dims2(1) = r1
    dims2(2) = m-1
    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)
    call h5dopen_f(group_id5, "G", dataset_id, hdferr)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, array_2D, dims2, hdferr)
    G_matrix = array_2D
    call h5dclose_f(dataset_id, hdferr)
    call h5sclose_f(dataspace_id, hdferr)
    deallocate(array_2D)


    vector = reshape([0.0, 0.0, 0.0, 0.0], [1, 4])
    new_col1 = [0.0, 0.0, 0.0, 0.0]
    new_col2 = [0.0, 0.0, 0.0, 0.0]
    new_col3 = [0.0, 0.0, 0.0, 0.0]
    new_col4 = [0.0, 0.0, 0.0, 0.0]



    do ii=2, r3+1

        answerr = 0

        do ij=1, size(vector,1)

            if (Index_m(ii-1)==vector(ij,2)) then

                answerr = 1

            end if

        end do

        if (answerr==0) then

            r4=rank(ii-1)

            write(filename, '(A, I0)') "Gamma_",Index_m(ii-1)

            !call h5gcreate_f(group_id9, filename, group_id4, hdferr)

            print *, ' Gamma', r4

            do j=1,r1

                B(j,j) = B_p(j,r4)

            end do

            xdmd = real(matmul(matmul(Vr,B),C))

            ! Compute true reconstruction error: ||G - Re(W*diag(b)*C)||_F^2 / ||G||_F^2
            summ = 0.0_dp

            do i=1,r1
                do j=1,m-1
                    summ = summ + (G_matrix(i,j) - xdmd(i,j))**2
                end do
            end do

            F_norm_XDMD = summ

            print *, "Error rank " ,Index_m(ii-1), "_", (F_norm_XDMD/F_norm_orgi)*100

            if (vector(1,1)==0) then

                new_col1    = [rank(ii-1)*1.0_dp]
                new_col2    = [Index_m(ii-1)*1.0_dp]
                new_col3    = [gamma_val(ii-1)*1.0_dp]
                new_col4    = [(F_norm_XDMD/F_norm_orgi)*100]

                vector = reshape([new_col1, new_col2, new_col3, new_col4], [size(vector,1), 4])

            else

                new_col1   = [new_col1, rank(ii-1)*1.0_dp]
                new_col2   = [new_col2, Index_m(ii-1)*1.0_dp]
                new_col3   = [new_col3, gamma_val(ii-1)*1.0_dp]
                new_col4   = [new_col4, (F_norm_XDMD/F_norm_orgi)*100]


                vector = reshape([new_col1, new_col2, new_col3, new_col4], [size(vector,1) +1 , 4])

            end if

        end if

        B = 0.0_dp

    end do



    dims2(1) = size(vector,2)
    dims2(2) = size(vector,1)

    call h5screate_simple_f(2, dims2, dataspace_id, hdferr)

    call h5dopen_f(group_id5, "Recons_vec", dataset_id, hdferr)
    if (hdferr == 0) then
        ! If hdferr is 0, the dataset exists, so close and delete it
        call h5dclose_f(dataset_id, hdferr)
        call h5ldelete_f(group_id5, "Recons_vec", hdferr)   ! Delete the existing dataset
    end if


    call h5dcreate_f(group_id5, "Recons_vec", H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, transpose(vector), dims2, hdferr)
    call h5dclose_f(dataset_id, hdferr)

    ! Close HDF5 file
    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

    write(2000, *) size(vector,1)

    print *, 'Finish'

end program reconstruction
