! This subroutine is part of d77 and writes
! calculated density to a cube file.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE cube_write
  USE global_read_data
  IMPLICIT NONE

CONTAINS

  SUBROUTINE cube_write_header(fid, comment1, comment2)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'cube_write_header'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Input variables
    INTEGER, INTENT(IN)            :: fid
    CHARACTER(LEN=100), INTENT(IN) :: comment1, comment2
    
!   Local variables
    INTEGER :: i_xyz, i_atm
    
    WRITE(fid,'(A)') TRIM(ADJUSTL(comment1))
    WRITE(fid,'(A)') TRIM(ADJUSTL(comment2))
    WRITE(fid,'(I5, 3F12.6)') N_atm,  Xmin,  Ymin,  Zmin
    WRITE(fid,'(I5, 3F12.6)')    Nx,    Dx, 0.0D0, 0.0D0
    WRITE(fid,'(I5, 3F12.6)')    Ny, 0.0D0,    Dy, 0.0D0
    WRITE(fid,'(I5, 3F12.6)')    Nz, 0.0D0, 0.0D0,    Dz
    
    DO i_atm = 1, N_atm
      WRITE(fid, '(I5, 4F12.6)') Atmnum(i_atm), DBLE(Atmnum(i_atm)),&
                                &(Xyznuc(i_xyz, i_atm), i_xyz = 1,3)
    ENDDO ! i_atm = 1, n_atm
    
  END SUBROUTINE cube_write_header

  SUBROUTINE cube_write_density(fid, data)
    
!   Input variables
    INTEGER, INTENT(IN) :: fid
    DOUBLE PRECISION    :: data(:)
    
!   Local variables
    INTEGER :: icol, ncol ! icol = 1, 2, 3, 4, 5, 6; ncol = 6
    
    ncol = UBOUND(data, 1)
    DO icol = 1, ncol
      IF(ABS(data(icol)) < 1.0D-99) data(icol) = 0.0D0
    ENDDO
    
    SELECT CASE(ncol)
      CASE(1)
        WRITE(fid, '(E13.5)') data(1:ncol)
      CASE(2)
        WRITE(fid,'(2E13.5)') data(1:ncol)
      CASE(3)
        WRITE(fid,'(3E13.5)') data(1:ncol)
      CASE(4)
        WRITE(fid,'(4E13.5)') data(1:ncol)
      CASE(5)
        WRITE(fid,'(5E13.5)') data(1:ncol)
      CASE(6)
        WRITE(fid,'(6E13.5)') data(1:ncol)
      CASE DEFAULT
    END SELECT

  END SUBROUTINE cube_write_density

  SUBROUTINE num2char(num09, char09)
    INTEGER, INTENT(IN) :: num09
    CHARACTER(LEN=1), INTENT(OUT) :: char09
 
    char09 = 'N' ! temporary value (not a number)
    SELECT CASE(num09)
      CASE(0); char09 = '0'
      CASE(1); char09 = '1'
      CASE(2); char09 = '2'
      CASE(3); char09 = '3'
      CASE(4); char09 = '4'
      CASE(5); char09 = '5'
      CASE(6); char09 = '6'
      CASE(7); char09 = '7'
      CASE(8); char09 = '8'
      CASE(9); char09 = '9'
      CASE DEFAULT
        WRITE(6,*) num09, ' error: not integer/num2char'
        STOP
    END SELECT
  END SUBROUTINE num2char
  
  SUBROUTINE cube_gen

!   Local variables
    INTEGER                       :: N_atm_1, Nx_1, Ny_1, Nz_1
    DOUBLE PRECISION              :: Xmin_1, Ymin_1, Zmin_1
    DOUBLE PRECISION              :: Dxx_1,  Dxy_1,  Dxz_1, &
                                    &Dyx_1,  Dyy_1,  Dyz_1, &
                                    &Dzx_1,  Dzy_1,  Dzz_1
    INTEGER, ALLOCATABLE          :: Atmnum_1(:)
    DOUBLE PRECISION, ALLOCATABLE :: Xyznuc_1(:,:)
    DOUBLE PRECISION              :: Atmnum_dummy 
    INTEGER :: i_atm, i_xyz
    INTEGER :: fid_1, fid_2, fid_res
    CHARACTER(LEN=80) :: comment1, comment2
    INTEGER           :: ix, iy, ncol, mod_nz_ncol_1, mod_nz_ncol_2, i_line, n_line
    DOUBLE PRECISION  :: data_1(1:6)
    
    ncol = 6
    mod_nz_ncol_1 = 0; mod_nz_ncol_2 = 0

    N_atm_1 = 0; Nx_1 = 0; Ny_1 = 0; Nz_1 = 0
    Xmin_1 = 0.0D0; Ymin_1 = 0.0D0; Zmin_1 = 0.0D0
    Dxx_1  = 0.0D0; Dxy_1  = 0.0D0; Dxz_1  = 0.0D0
    Dyx_1  = 0.0D0; Dyy_1  = 0.0D0; Dyz_1  = 0.0D0
    Dzx_1  = 0.0D0; Dzy_1  = 0.0D0; Dzz_1  = 0.0D0

    comment1 = 'CUBE_FILE_1 = ; CUBE_FILE_2 = '
    comment2 = 'CUBE_FILE_1 - CUBE_FILE_2'

    SELECT CASE(Cube_op)
      CASE('Add')
        comment2 = 'CUBE_FILE_1 + CUBE_FILE_2'
      CASE('Sub')
        comment2 = 'CUBE_FILE_1 - CUBE_FILE_2'
      CASE('Mul')
        comment2 = 'CUBE_FILE_1 * CUBE_FILE_2'
      CASE('Scale')
        comment2 = 'Scaling factor ='
      CASE DEFAULT
    END SELECT    
    
    fid_1   = 10
    fid_2   = 20
    fid_res = 100
    OPEN(fid_1, FILE = 'CUBE_FILE_1')
    SELECT CASE(Cube_op)
      CASE('Add', 'Sub', 'Mul')
        OPEN(fid_2, FILE = 'CUBE_FILE_2')
      CASE DEFAULT
    END SELECT    
    OPEN(fid_res, FILE = 'CUBE_RES')

      READ(fid_1,*) 
      READ(fid_1,*) 
      READ(fid_1,*) N_atm_1, Xmin_1, Ymin_1, Zmin_1
      READ(fid_1,*) Nx_1,    Dxx_1,  Dxy_1,  Dxz_1
      READ(fid_1,*) Ny_1,    Dyx_1,  Dyy_1,  Dyz_1
      READ(fid_1,*) Nz_1,    Dzx_1,  Dzy_1,  Dzz_1

      mod_nz_ncol_1 = MOD(Nz_1, ncol)
      
      ALLOCATE(Atmnum_1(1:N_atm_1)); Atmnum_1 = 0
      ALLOCATE(Xyznuc_1(1:3, 1:N_atm_1)); Xyznuc_1 = 0.0D0
      Atmnum_dummy = 0.0D0
      DO i_atm = 1, N_atm_1
        READ(fid_1, *) Atmnum_1(i_atm), Atmnum_dummy, &
                      &(Xyznuc_1(i_xyz, i_atm), i_xyz = 1,3)
      ENDDO ! i_atm 

      WRITE(fid_res,'(A)') TRIM(ADJUSTL(comment1))
      WRITE(fid_res,'(A)') TRIM(ADJUSTL(comment2))
      WRITE(fid_res,'(I5, 3F12.6)') N_atm_1, Xmin_1, Ymin_1, Zmin_1
      WRITE(fid_res,'(I5, 3F12.6)') Nx_1,    Dxx_1,  Dxy_1,  Dxz_1
      WRITE(fid_res,'(I5, 3F12.6)') Ny_1,    Dyx_1,  Dyy_1,  Dyz_1
      WRITE(fid_res,'(I5, 3F12.6)') Nz_1,    Dzx_1,  Dzy_1,  Dzz_1
      DO i_atm = 1, N_atm_1
        WRITE(fid_res, '(I5, 4F12.6)') Atmnum_1(i_atm), DBLE(Atmnum_1(i_atm)),&
                                      &(Xyznuc_1(i_xyz, i_atm), i_xyz = 1,3)
      ENDDO

      IF(mod_nz_ncol_1 == 0) THEN ! 6 columns
        n_line = nz_1/ncol ! The number of lines for each ix-iy pair
        DO ix = 1, Nx_1
          DO iy = 1, Ny_1
            DO i_line = 1, n_line
              CALL cube_read_density (fid_1,   data_1(1:ncol))
              CALL cube_write_density(fid_res, data_1(1:ncol))
            ENDDO ! i_line = 1, n_line
          ENDDO ! iy = 1, ny
        ENDDO ! ix = 1, nx
      ELSE ! 1-5 columns
        n_line = (Nz_1 - mod_nz_ncol_1)/ncol + 1
        DO ix = 1, Nz_1
          DO iy = 1, Ny_1
            IF(n_line == 1) THEN
              CALL cube_read_density (fid_1,   data_1(1:mod_nz_ncol_1))
              CALL cube_write_density(fid_res, data_1(1:mod_nz_ncol_1))
            ELSE ! n_line > 2
              DO i_line = 1, n_line - 1
                CALL cube_read_density (fid_1,   data_1(1:ncol))
                CALL cube_write_density(fid_res, data_1(1:ncol))
              ENDDO
              CALL cube_read_density (fid_1,   data_1(1:mod_nz_ncol_1))
              CALL cube_write_density(fid_res, data_1(1:mod_nz_ncol_1))
            ENDIF
          ENDDO ! iy = 1, ny
        ENDDO ! ix = 1, nx
      ENDIF

    CLOSE(fid_1)
    SELECT CASE(Cube_op)
      CASE('Add', 'Sub', 'Mul')
        CLOSE(fid_2)
      CASE DEFAULT
    END SELECT
    CLOSE(fid_res)

  END SUBROUTINE cube_gen
  
  SUBROUTINE cube_read_density(fid, data)

!   Input variables
    INTEGER, INTENT(IN)           :: fid
    DOUBLE PRECISION, INTENT(OUT) :: data(:)

!   Local variables
    INTEGER :: ncol 

    ncol = UBOUND(data, 1)
    data = 0.0D0

    SELECT CASE(ncol)
      CASE(1)
        READ(fid,*) data(1:ncol)
      CASE(2)
        READ(fid,*) data(1:ncol)
      CASE(3)
        READ(fid,*) data(1:ncol)
      CASE(4)
        READ(fid,*) data(1:ncol)
      CASE(5)
        READ(fid,*) data(1:ncol)
      CASE(6)
        READ(fid,*) data(1:ncol)
      CASE DEFAULT
    END SELECT

  END SUBROUTINE cube_read_density


END MODULE cube_write
