! This subroutine is part of d77 and writes
! calculated density to a cube file.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE cube
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

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'cube_gen'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables
    INTEGER                       :: N_atm_1, N_atm_2
    INTEGER                       :: N_1(1:3), N_2(1:3)
    DOUBLE PRECISION              :: XYZmin_1(1:3), XYZmin_2(1:3)
    DOUBLE PRECISION              :: D_1(1:3,1:3), D_2(1:3,1:3)
    INTEGER, ALLOCATABLE          :: Atmnum_1(:), Atmnum_2(:)
    DOUBLE PRECISION, ALLOCATABLE :: Xyznuc_1(:,:), Xyznuc_2(:,:)
    DOUBLE PRECISION              :: Atmnum_dummy 
    INTEGER :: i_atm, i_xyz
    INTEGER :: fid_cube_file_name, fid_1, fid_2, fid_res
    CHARACTER(LEN=80) :: cube_file_name_res, cube_file_name_1, cube_file_name_2
    CHARACTER(LEN=80) :: comment1, comment2
    CHARACTER(LEN=13) :: text_scale_cube
    INTEGER           :: ix, iy, ncol, mod_nz_ncol_1, mod_nz_ncol_2, i_line, n_line
    DOUBLE PRECISION  :: data_1(1:6), data_2(1:6), data_res(1:6)
    
    CALL write_messages(2, Text_blank, type_program, name_program)
    
    ncol = 6
    mod_nz_ncol_1 = 0; mod_nz_ncol_2 = 0

    cube_file_name_res = REPEAT(' ', 80)
    cube_file_name_1 = REPEAT(' ', 80); cube_file_name_2  = REPEAT(' ', 80)
    comment1 = REPEAT(' ', 80); comment2  = REPEAT(' ', 80)
    
    N_atm_1 = 0; N_atm_2 = 0

    N_1 = 0; N_2 = 0
    XYZmin_1 = 0.0D0; XYZmin_2 = 0.0D0
    D_1 = 0.0D0; D_2 = 0.0D0
    data_1 = 0.0D0; data_2 = 0.0D0; data_res = 0.0D0

    fid_cube_file_name = 10
    OPEN(fid_cube_file_name, FILE = 'CUBE_FILE_NAME')
      READ(fid_cube_file_name, '(A)') cube_file_name_res
      READ(fid_cube_file_name, '(A)') cube_file_name_1
      SELECT CASE(Cube_op)
        CASE('Add', 'Sub', 'Mul')
          READ(fid_cube_file_name, '(A)') cube_file_name_2
        CASE DEFAULT
      END SELECT
    CLOSE(fid_cube_file_name)

    SELECT CASE(Cube_op)
      CASE('Add')
        comment1 = 'CUBE_1 = '//TRIM(ADJUSTL(cube_file_name_1))//&
                &'; CUBE_2 = '//TRIM(ADJUSTL(cube_file_name_2))
        comment2 = TRIM(ADJUSTL(cube_file_name_res))//' = CUBE_1 + CUBE_2'
      CASE('Sub')
        comment1 = 'CUBE_1 = '//TRIM(ADJUSTL(cube_file_name_1))//&
                &'; CUBE_2 = '//TRIM(ADJUSTL(cube_file_name_2))
        comment2 = TRIM(ADJUSTL(cube_file_name_res))//' = CUBE_1 - CUBE_2'
      CASE('Mul')
        comment1 = 'CUBE_1 = '//TRIM(ADJUSTL(cube_file_name_1))//&
                &'; CUBE_2 = '//TRIM(ADJUSTL(cube_file_name_2))
        comment2 = TRIM(ADJUSTL(cube_file_name_res))//' = CUBE_1 * CUBE_2'
      CASE('Scale')
        WRITE(text_scale_cube, '(E13.5)') Scale_cube
        comment1 = 'CUBE_1 = '//TRIM(ADJUSTL(cube_file_name_1))
        comment2 = TRIM(ADJUSTL(cube_file_name_res))//' = CUBE_1 * '//TRIM(ADJUSTL(text_scale_cube))
      CASE DEFAULT
    END SELECT    
    
    fid_1   = 11
    fid_2   = 12
    fid_res = 13
    OPEN(fid_1, FILE = 'CUBE_1')
    SELECT CASE(Cube_op)
      CASE('Add', 'Sub', 'Mul')
        OPEN(fid_2, FILE = 'CUBE_2')
      CASE DEFAULT
    END SELECT    
    OPEN(fid_res, FILE = 'CUBE_RES')

      WRITE(6,'(X)') 
      WRITE(6,'(X, A)') 'Reading CUBE_1'
      READ(fid_1,*) 
      READ(fid_1,*) 
      READ(fid_1,*) N_atm_1, XYZmin_1(1), XYZmin_1(2), XYZmin_1(3)
      READ(fid_1,*) N_1(1),  D_1(1,1),    D_1(1,2),    D_1(1,3)
      READ(fid_1,*) N_1(2),  D_1(2,1),    D_1(2,2),    D_1(2,3)
      READ(fid_1,*) N_1(3),  D_1(3,1),    D_1(3,2),    D_1(3,3)

      mod_nz_ncol_1 = MOD(N_1(3), ncol)
    
      ALLOCATE(Atmnum_1(1:N_atm_1)); Atmnum_1 = 0
      ALLOCATE(Xyznuc_1(1:3, 1:N_atm_1)); Xyznuc_1 = 0.0D0
      Atmnum_dummy = 0.0D0
      DO i_atm = 1, N_atm_1
        READ(fid_1, *) Atmnum_1(i_atm), Atmnum_dummy, &
                      &(Xyznuc_1(i_xyz, i_atm), i_xyz = 1,3)
      ENDDO ! i_atm 

      SELECT CASE(Cube_op)
        CASE('Add', 'Sub', 'Mul')
          WRITE(6,'(X, A)') 'Reading CUBE_2'
          READ(fid_2,*)
          READ(fid_2,*)
          READ(fid_2,*) N_atm_2, XYZmin_2(1), XYZmin_2(2), XYZmin_2(3)
          READ(fid_2,*) N_2(1),  D_2(1,1),    D_2(1,2),    D_2(1,3)
          READ(fid_2,*) N_2(2),  D_2(2,1),    D_2(2,2),    D_2(2,3)
          READ(fid_2,*) N_2(3),  D_2(3,1),    D_2(3,2),    D_2(3,3)

          IF(N_atm_2 /= N_atm_1) THEN
            WRITE(fid_res,'(A)')     'The number of atoms in CUBE_1 is not equal to that in CUBE_2.'
            WRITE(fid_res,'(X)') 
            WRITE(fid_res,'(A, I0)') 'The number of atoms in CUBE_1: ', N_atm_1
            WRITE(fid_res,'(A, I0)') 'The number of atoms in CUBE_2: ', N_atm_2
            CALL write_messages(-9999, Text_blank, type_program, name_program)
          ELSEIF(N_2(1) /= N_1(1) .OR. N_2(2) /= N_1(2) .OR. N_2(3) /= N_1(3)) THEN
            WRITE(fid_res,'(A)')          'Grid size of CUBE_1 is not equal to that of CUBE_2.'
            WRITE(fid_res,'(X)') 
            WRITE(fid_res,'(A)')          'CUBE_1'
            WRITE(fid_res,'(I5, 3F12.6)') N_atm_1, XYZmin_1(1), XYZmin_1(2), XYZmin_1(3)
            WRITE(fid_res,'(I5, 3F12.6)') N_1(1),  D_1(1,1),    D_1(1,2),    D_1(1,3)
            WRITE(fid_res,'(I5, 3F12.6)') N_1(2),  D_1(2,1),    D_1(2,2),    D_1(2,3)
            WRITE(fid_res,'(I5, 3F12.6)') N_1(3),  D_1(3,1),    D_1(3,2),    D_1(3,3)
            WRITE(fid_res,'(X)') 
            WRITE(fid_res,'(A)')          'CUBE_2'
            WRITE(fid_res,'(I5, 3F12.6)') N_atm_2, XYZmin_2(1), XYZmin_2(2), XYZmin_2(3)
            WRITE(fid_res,'(I5, 3F12.6)') N_2(1),  D_2(1,1),    D_2(1,2),    D_2(1,3)
            WRITE(fid_res,'(I5, 3F12.6)') N_2(2),  D_2(2,1),    D_2(2,2),    D_2(2,3)
            WRITE(fid_res,'(I5, 3F12.6)') N_2(3),  D_2(3,1),    D_2(3,2),    D_2(3,3)
            CALL write_messages(-9999, Text_blank, type_program, name_program)
          ELSE
          ENDIF

          ALLOCATE(Atmnum_2(1:N_atm_2)); Atmnum_2 = 0
          ALLOCATE(Xyznuc_2(1:3, 1:N_atm_2)); Xyznuc_2 = 0.0D0
          Atmnum_dummy = 0.0D0
          DO i_atm = 1, N_atm_2
            READ(fid_2, *) Atmnum_2(i_atm), Atmnum_dummy, &
                          &(Xyznuc_2(i_xyz, i_atm), i_xyz = 1,3)
          ENDDO ! i_atm
        CASE DEFAULT
      END SELECT
      WRITE(6,'(X, A)') 'Done'

      WRITE(fid_res,'(A)') TRIM(ADJUSTL(comment1))
      WRITE(fid_res,'(A)') TRIM(ADJUSTL(comment2))
      WRITE(fid_res,'(I5, 3F12.6)') N_atm_1, XYZmin_1(1), XYZmin_1(2), XYZmin_1(3)
      WRITE(fid_res,'(I5, 3F12.6)') N_1(1),  D_1(1,1),    D_1(1,2),    D_1(1,3)
      WRITE(fid_res,'(I5, 3F12.6)') N_1(2),  D_1(2,1),    D_1(2,2),    D_1(2,3)
      WRITE(fid_res,'(I5, 3F12.6)') N_1(3),  D_1(3,1),    D_1(3,2),    D_1(3,3)
      DO i_atm = 1, N_atm_1
        WRITE(fid_res, '(I5, 4F12.6)') Atmnum_1(i_atm), DBLE(Atmnum_1(i_atm)),&
                                      &(Xyznuc_1(i_xyz, i_atm), i_xyz = 1,3)
      ENDDO

      IF(mod_nz_ncol_1 == 0) THEN ! 6 columns
        n_line = N_1(3)/ncol ! The number of lines for each ix-iy pair
        DO ix = 1, N_1(1)
          DO iy = 1, N_1(2)
            DO i_line = 1, n_line

              CALL cube_read_density (fid_1, data_1(1:ncol))
              SELECT CASE(Cube_op)
                CASE('Add', 'Sub', 'Mul')
                  CALL cube_read_density (fid_2, data_2(1:ncol))
                CASE DEFAULT
              END SELECT
              IF(Cube_op == 'Scale') THEN
                data_res = data_1 * Scale_cube
              ELSEIF(Cube_op == 'Add') THEN
                data_res = data_1 + data_2
              ELSEIF(Cube_op == 'Sub') THEN
                data_res = data_1 - data_2
              ELSEIF(Cube_op == 'Mul') THEN
                data_res = data_1 * data_2
              ELSE
              ENDIF
              CALL cube_write_density(fid_res, data_res(1:ncol))

            ENDDO ! i_line = 1, n_line
          ENDDO ! iy  
        ENDDO ! ix 
      ELSE ! 1-5 columns
        n_line = (N_1(3) - mod_nz_ncol_1)/ncol + 1
        DO ix = 1, N_1(1)
          DO iy = 1, N_1(2)
            IF(n_line == 1) THEN
 
              CALL cube_read_density (fid_1, data_1(1:mod_nz_ncol_1))
              SELECT CASE(Cube_op)
                CASE('Add', 'Sub', 'Mul')
                  CALL cube_read_density (fid_2, data_2(1:mod_nz_ncol_1))
                CASE DEFAULT
              END SELECT
              IF(Cube_op == 'Scale') THEN
                data_res = data_1 * Scale_cube
              ELSEIF(Cube_op == 'Add') THEN
                data_res = data_1 + data_1
              ELSEIF(Cube_op == 'Sub') THEN
                data_res = data_1 - data_2
              ELSEIF(Cube_op == 'Mul') THEN
                data_res = data_1 * data_2
              ELSE
              ENDIF
              CALL cube_write_density(fid_res, data_res(1:mod_nz_ncol_1))

            ELSE ! n_line > 2
              DO i_line = 1, n_line - 1

                CALL cube_read_density (fid_1, data_1(1:ncol))
                SELECT CASE(Cube_op)
                  CASE('Add', 'Sub', 'Mul')
                    CALL cube_read_density (fid_2, data_2(1:ncol))
                  CASE DEFAULT
                END SELECT
                IF(Cube_op == 'Scale') THEN
                  data_res = data_1 * Scale_cube
                ELSEIF(Cube_op == 'Add') THEN
                  data_res = data_1 + data_2
                ELSEIF(Cube_op == 'Sub') THEN
                  data_res = data_1 - data_2
                ELSEIF(Cube_op == 'Mul') THEN
                  data_res = data_1 * data_2
                ELSE
                ENDIF
                CALL cube_write_density(fid_res, data_res(1:ncol))

              ENDDO

              CALL cube_read_density (fid_1,   data_1(1:mod_nz_ncol_1))
              SELECT CASE(Cube_op)
                CASE('Add', 'Sub', 'Mul')
                  CALL cube_read_density (fid_2, data_2(1:mod_nz_ncol_1))
                CASE DEFAULT
              END SELECT
              IF(Cube_op == 'Scale') THEN
                data_res = data_1 * Scale_cube
              ELSEIF(Cube_op == 'Add') THEN
                data_res = data_1 + data_2
              ELSEIF(Cube_op == 'Sub') THEN
                data_res = data_1 - data_2
              ELSEIF(Cube_op == 'Mul') THEN
                data_res = data_1 * data_2
              ELSE
              ENDIF
              CALL cube_write_density(fid_res, data_res(1:mod_nz_ncol_1))

            ENDIF
          ENDDO ! iy 
        ENDDO ! ix 
      ENDIF

    CLOSE(fid_1)
    SELECT CASE(Cube_op)
      CASE('Add', 'Sub', 'Mul')
        CLOSE(fid_2)
      CASE DEFAULT
    END SELECT
    CLOSE(fid_res)

    CALL write_messages(3, Text_blank, type_program, name_program)
  
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


END MODULE cube
