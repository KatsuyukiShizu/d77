! A part of dpp
!
! dpp is a data preprocessing utility for d77.
!
! dpp is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE gen_g16_cicoef_td_unres
  USE global_constants
  USE read_txt
  IMPLICIT NONE

  DOUBLE PRECISION :: threshold_ci_coef 

  INTEGER :: N, K_divided_by_2, n_ex_state
  CHARACTER(LEN=10) mult_ex(1:Max_n_state)

! Dummy variables
  INTEGER :: mo_occ_x_a(1:Max_n_elec_config,1:Max_n_state), &
            &mo_unocc_x_a(1:Max_n_elec_config,1:Max_n_state), &
            &mo_occ_x_b(1:Max_n_elec_config,1:Max_n_state), &
            &mo_unocc_x_b(1:Max_n_elec_config,1:Max_n_state), &
            &mo_occ_y_a(1:Max_n_elec_config,1:Max_n_state), &
            &mo_unocc_y_a(1:Max_n_elec_config,1:Max_n_state), &
            &mo_occ_y_b(1:Max_n_elec_config,1:Max_n_state), &
            &mo_unocc_y_b(1:Max_n_elec_config,1:Max_n_state)
  CHARACTER(LEN=2) :: arrow_x_as_read(1:Max_n_elec_config,1:Max_n_state), &
                      &arrow_y_as_read(1:Max_n_elec_config,1:Max_n_state)
  DOUBLE PRECISION :: x_a(1:Max_n_elec_config,1:Max_n_state), &
                     &x_b(1:Max_n_elec_config,1:Max_n_state), &
                     &y_a(1:Max_n_elec_config,1:Max_n_state), &
                     &y_b(1:Max_n_elec_config,1:Max_n_state), &
                     &yy_a(1:Max_n_elec_config,1:Max_n_state), &
                     &yy_b(1:Max_n_elec_config,1:Max_n_state)
  CHARACTER(LEN=7) :: label_x_a(1:Max_n_elec_config,1:Max_n_state), &
                     &label_x_b(1:Max_n_elec_config,1:Max_n_state), &
                     &label_y_a(1:Max_n_elec_config,1:Max_n_state), &
                     &label_y_b(1:Max_n_elec_config,1:Max_n_state)


  INTEGER, ALLOCATABLE :: so_occ_a(:,:), so_unocc_a(:,:)
  INTEGER, ALLOCATABLE :: so_occ_b(:,:), so_unocc_b(:,:)
  CHARACTER(2), ALLOCATABLE :: arrow(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: x_plus_y_a(:,:), x_minus_y_a(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: x_plus_y_b(:,:), x_minus_y_b(:,:)

  INTEGER :: i_space(1:2)

  INTEGER :: i_elec_config, i_elec_config_x, i_elec_config_y
  INTEGER :: n_elec_config_x_a(1:Max_n_state), &
            &n_elec_config_x_b(1:Max_n_state), &
            &n_elec_config_y_a(1:Max_n_state), &
            &n_elec_config_y_b(1:Max_n_state), &
            &n_elec_config_x_only_a(1:Max_n_state), &
            &n_elec_config_x_only_b(1:Max_n_state), &
            &n_elec_config_y_only_a(1:Max_n_state), &
            &n_elec_config_y_only_b(1:Max_n_state), &
            &n_elec_config_x_and_y_a(1:Max_n_state), &
            &n_elec_config_x_and_y_b(1:Max_n_state)
  INTEGER :: i_ex_state
  INTEGER :: fid
  INTEGER :: counter

  i_space = 0

  N = 0; K_divided_by_2 = 0; n_ex_state = 0
  n_elec_config_x_a = 0; n_elec_config_y_a = 0
  n_elec_config_x_b = 0; n_elec_config_y_b = 0
  n_elec_config_x_only_a = 0; n_elec_config_x_only_b = 0
  n_elec_config_y_only_a = 0; n_elec_config_y_only_b = 0
  n_elec_config_x_and_y_a = 0; n_elec_config_x_and_y_b = 0

  mo_occ_x_a = 0; mo_unocc_x_a = 0
  mo_occ_x_b = 0; mo_unocc_x_b = 0
  mo_occ_y_a = 0; mo_unocc_y_a = 0
  mo_occ_y_b = 0; mo_unocc_y_b = 0
  x_a = 0.0D0; x_b = 0.0D0
  y_a = 0.0D0; y_b = 0.0D0; 
  yy_a = 0.0D0; yy_b = 0.0D0
  arrow_x_as_read = '  '; arrow_y_as_read = '  '
  label_x_a = 'x_only'; label_x_b = 'x_only'
  label_y_a = 'x_only'; label_y_b = 'y_only'

  threshold_ci_coef = threshold_ci_coef_default
  OPEN(FR_threshold_values_ci, FILE = Fname(FR_threshold_values_ci))
    READ(FR_threshold_values_ci, *) threshold_ci_coef
  CLOSE(FR_threshold_values_ci)

  fid = FR_g16_td_log_x
  OPEN(fid, FILE = Fname(fid))
    READ(fid,*) N
    READ(fid,*) K_divided_by_2
    READ(fid,*) n_ex_state
    DO i_ex_state = 1, n_ex_state
      READ(fid,*) mult_ex(i_ex_state)

      READ(fid,*) n_elec_config_x_a(i_ex_state)
      IF(n_elec_config_x_a(i_ex_state) > 0) THEN
        DO i_elec_config = 1, n_elec_config_x_a(i_ex_state)
          READ(fid,*) mo_occ_x_a(i_elec_config, i_ex_state), &
                     &mo_unocc_x_a(i_elec_config, i_ex_state), &
                     &x_a(i_elec_config, i_ex_state)
        ENDDO
      ELSE
      ENDIF

      READ(fid,*) n_elec_config_x_b(i_ex_state)
      IF(n_elec_config_x_b(i_ex_state) > 0) THEN
        DO i_elec_config = 1, n_elec_config_x_b(i_ex_state)
          READ(fid,*) mo_occ_x_b(i_elec_config, i_ex_state), &
                     &mo_unocc_x_b(i_elec_config, i_ex_state), &
                     &x_b(i_elec_config, i_ex_state)
        ENDDO
      ELSE
      ENDIF

    ENDDO
  CLOSE(fid)

  fid = FR_g16_td_log_y
  OPEN(fid, FILE = Fname(fid))
    READ(fid,*) 
    DO i_ex_state = 1, n_ex_state
      READ(fid,*) 
      READ(fid,*) n_elec_config_y_a(i_ex_state)
      IF(n_elec_config_y_a(i_ex_state) > 0) THEN
        DO i_elec_config = 1, n_elec_config_y_a(i_ex_state)
          READ(fid,*) mo_occ_y_a(i_elec_config, i_ex_state),&
                     &mo_unocc_y_a(i_elec_config, i_ex_state),&
                     &y_a(i_elec_config, i_ex_state)
        ENDDO
      ELSE
      ENDIF

      READ(fid,*) n_elec_config_y_b(i_ex_state)
      IF(n_elec_config_y_b(i_ex_state) > 0) THEN
        DO i_elec_config = 1, n_elec_config_y_b(i_ex_state)
          READ(fid,*) mo_occ_y_b(i_elec_config, i_ex_state),&
                     &mo_unocc_y_b(i_elec_config, i_ex_state),&
                     &y_b(i_elec_config, i_ex_state)
        ENDDO
      ELSE
      ENDIF
    ENDDO
  CLOSE(fid)

  DO i_ex_state = 1, n_ex_state
    IF(n_elec_config_x_a(i_ex_state) > 0) THEN
      DO i_elec_config_x = 1, n_elec_config_x_a(i_ex_state)
        IF(n_elec_config_y_a(i_ex_state) > 0) THEN
          DO i_elec_config_y = 1, n_elec_config_y_a(i_ex_state)
            IF(   mo_occ_x_a(i_elec_config_x, i_ex_state) &
              &== mo_occ_y_a(i_elec_config_y, i_ex_state) .AND. &
                 &mo_unocc_x_a(i_elec_config_x, i_ex_state) &
              &== mo_unocc_y_a(i_elec_config_y, i_ex_state)) THEN
    
              label_x_a(i_elec_config_x, i_ex_state) = 'x_and_y'
              label_y_a(i_elec_config_y, i_ex_state) = 'x_and_y'
              yy_a(i_elec_config_x, i_ex_state) &
             &= y_a(i_elec_config_y, i_ex_state)
            ELSE
            ENDIF
          ENDDO
        ELSE  
        ENDIF  
      ENDDO
    ELSE
    ENDIF  

    IF(n_elec_config_x_b(i_ex_state) > 0) THEN
      DO i_elec_config_x = 1, n_elec_config_x_b(i_ex_state)
        IF(n_elec_config_y_b(i_ex_state) > 0) THEN
          DO i_elec_config_y = 1, n_elec_config_y_b(i_ex_state)
            IF(   mo_occ_x_b(i_elec_config_x, i_ex_state) &
              &== mo_occ_y_b(i_elec_config_y, i_ex_state) .AND. &
                 &mo_unocc_x_b(i_elec_config_x, i_ex_state) &
              &== mo_unocc_y_b(i_elec_config_y, i_ex_state)) THEN
          
              label_x_b(i_elec_config_x, i_ex_state) = 'x_and_y'
              label_y_b(i_elec_config_y, i_ex_state) = 'x_and_y'
              yy_b(i_elec_config_x, i_ex_state) &
             &= y_b(i_elec_config_y, i_ex_state)
            ELSE
            ENDIF
          ENDDO
        ELSE
        ENDIF  
      ENDDO
    ELSE
    ENDIF
  ENDDO




  ALLOCATE(so_occ_a(1:2*Max_n_elec_config,1:n_ex_state)); so_occ_a = 0
  ALLOCATE(so_unocc_a(1:2*Max_n_elec_config,1:n_ex_state)); so_unocc_a = 0
  ALLOCATE(so_occ_b(1:2*Max_n_elec_config,1:n_ex_state)); so_occ_b = 0
  ALLOCATE(so_unocc_b(1:2*Max_n_elec_config,1:n_ex_state)); so_unocc_b = 0
  ALLOCATE(arrow(1:2*Max_n_elec_config,1:n_ex_state)); arrow = '  '
  ALLOCATE(x_plus_y_a(1:2*Max_n_elec_config,1:n_ex_state)); x_plus_y_a = 0.0D0
  ALLOCATE(x_minus_y_a(1:2*Max_n_elec_config,1:n_ex_state)); x_minus_y_a= 0.0D0
  ALLOCATE(x_plus_y_b(1:2*Max_n_elec_config,1:n_ex_state)); x_plus_y_b = 0.0D0
  ALLOCATE(x_minus_y_b(1:2*Max_n_elec_config,1:n_ex_state)); x_minus_y_b = 0.0D0

  DO i_ex_state = 1, n_ex_state

!   Alpha electrons
    DO i_elec_config_x = 1, n_elec_config_x_a(i_ex_state)

      so_occ_a(i_elec_config_x, i_ex_state)&
      = 2*mo_occ_x_a(i_elec_config_x, i_ex_state)-1
      so_unocc_a(i_elec_config_x, i_ex_state)&
      = 2*mo_unocc_x_a(i_elec_config_x, i_ex_state)-1
      arrow(i_elec_config_x, i_ex_state) &
     &= arrow_x_as_read(i_elec_config_x, i_ex_state)
      x_plus_y_a(i_elec_config_x, i_ex_state) &
     &= x_a(i_elec_config_x, i_ex_state)
      x_minus_y_a(i_elec_config_x, i_ex_state) &
     &= x_a(i_elec_config_x, i_ex_state)

      SELECT CASE(label_x_a(i_elec_config_x, i_ex_state))
        CASE('x_and_y')
          n_elec_config_x_and_y_a(i_ex_state) &
         &= n_elec_config_x_and_y_a(i_ex_state) + 1 

          x_plus_y_a(i_elec_config_x, i_ex_state) &
         &= x_plus_y_a(i_elec_config_x, i_ex_state) &
         &+ yy_a(i_elec_config_x, i_ex_state)

          x_minus_y_a(i_elec_config_x, i_ex_state) &
         &= x_minus_y_a(i_elec_config_x, i_ex_state) &
         &- yy_a(i_elec_config_x, i_ex_state)

        CASE DEFAULT
      END SELECT  
    ENDDO

!   Beta electrons
    DO i_elec_config_x = 1, n_elec_config_x_b(i_ex_state)

      so_occ_b(i_elec_config_x, i_ex_state)&
      = 2*mo_occ_x_b(i_elec_config_x, i_ex_state)
      so_unocc_b(i_elec_config_x, i_ex_state)&
      = 2*mo_unocc_x_b(i_elec_config_x, i_ex_state)
      arrow(i_elec_config_x, i_ex_state) &
     &= arrow_x_as_read(i_elec_config_x, i_ex_state)
      x_plus_y_b(i_elec_config_x, i_ex_state) &
     &= x_b(i_elec_config_x, i_ex_state)
      x_minus_y_b(i_elec_config_x, i_ex_state) &
     &= x_b(i_elec_config_x, i_ex_state)

      SELECT CASE(label_x_b(i_elec_config_x, i_ex_state))
        CASE('x_and_y')
          n_elec_config_x_and_y_b(i_ex_state) &
         &= n_elec_config_x_and_y_b(i_ex_state) + 1

          x_plus_y_b(i_elec_config_x, i_ex_state) &
         &= x_plus_y_b(i_elec_config_x, i_ex_state) &
         &+ yy_b(i_elec_config_x, i_ex_state)

          x_minus_y_b(i_elec_config_x, i_ex_state) &
         &= x_minus_y_b(i_elec_config_x, i_ex_state) &
         &- yy_b(i_elec_config_x, i_ex_state)

        CASE DEFAULT
      END SELECT
    ENDDO

  ENDDO
  
  n_elec_config_x_only_a = n_elec_config_x_a - n_elec_config_x_and_y_a
  n_elec_config_x_only_b = n_elec_config_x_b - n_elec_config_x_and_y_b

  DO i_ex_state = 1, n_ex_state
!   Alpha electron
    counter = n_elec_config_x_a(i_ex_state) 
    IF(n_elec_config_y_a(i_ex_state) > 0) THEN
      DO i_elec_config_y = 1, n_elec_config_y_a(i_ex_state)
        SELECT CASE(label_y_a(i_elec_config_y, i_ex_state))
          CASE('y_only')
            counter = counter + 1

            so_occ_a(counter, i_ex_state)&
            = 2*mo_occ_y_a(i_elec_config_y, i_ex_state)-1
            so_unocc_a(counter, i_ex_state)&
            = 2*mo_unocc_y_a(i_elec_config_y, i_ex_state)-1
            x_plus_y_a(counter, i_ex_state) &
           &= y_a(i_elec_config_y, i_ex_state)
            x_minus_y_a(counter, i_ex_state) &
           &= - y_a(i_elec_config_y, i_ex_state)
      
        CASE DEFAULT
      END SELECT  
      ENDDO
    n_elec_config_y_only_a(i_ex_state) &
   &= counter - n_elec_config_x_a(i_ex_state)
    ELSE
    ENDIF 

!   Beta electron
    counter = n_elec_config_x_b(i_ex_state)
    IF(n_elec_config_y_b(i_ex_state) > 0) THEN
      DO i_elec_config_y = 1, n_elec_config_y_b(i_ex_state)
        SELECT CASE(label_y_b(i_elec_config_y, i_ex_state))
          CASE('y_only')
            counter = counter + 1

            so_occ_b(counter, i_ex_state)&
            = 2*mo_occ_y_b(i_elec_config_y, i_ex_state)-1
            so_unocc_b(counter, i_ex_state)&
            = 2*mo_unocc_y_b(i_elec_config_y, i_ex_state)-1
            x_plus_y_b(counter, i_ex_state) &
           &= y_b(i_elec_config_y, i_ex_state)
            x_minus_y_b(counter, i_ex_state) &
           &= - y_b(i_elec_config_y, i_ex_state)

        CASE DEFAULT
      END SELECT
      ENDDO
    n_elec_config_y_only_b(i_ex_state) &
   &= counter - n_elec_config_x_b(i_ex_state)
    ELSE
    ENDIF

  ENDDO

  fid = FW_xy_count
  OPEN(fid, FILE = Fname(fid))
    DO i_ex_state = 1, n_ex_state
      WRITE(fid,'(A8, I2)') 'STATE = ', i_ex_state
      WRITE(fid,9001) 'Alpha electron'
      WRITE(fid,'(A9, I2)') 'X only : ', n_elec_config_x_only_a(i_ex_state)
      WRITE(fid,'(A9, I2)') 'Y only : ', n_elec_config_y_only_a(i_ex_state)
      WRITE(fid,'(A9, I2)') 'X and Y: ', n_elec_config_x_and_y_a(i_ex_state)
      WRITE(fid,'(A9, I2)') 'X total: ', n_elec_config_x_a(i_ex_state)
      WRITE(fid,'(A9, I2)') 'y total: ', n_elec_config_y_a(i_ex_state)
      WRITE(fid,9001) 'Beta electron'
      WRITE(fid,'(A9, I2)') 'X only : ', n_elec_config_x_only_b(i_ex_state)
      WRITE(fid,'(A9, I2)') 'Y only : ', n_elec_config_y_only_b(i_ex_state)
      WRITE(fid,'(A9, I2)') 'X and Y: ', n_elec_config_x_and_y_b(i_ex_state)
      WRITE(fid,'(A9, I2)') 'X total: ', n_elec_config_x_b(i_ex_state)
      WRITE(fid,'(A9, I2)') 'y total: ', n_elec_config_y_b(i_ex_state)
      WRITE(fid,9001) 'Beta electron'
      WRITE(fid,'(1X)') 
    ENDDO
  CLOSE(fid)

  fid = FW_xy
  OPEN(fid, FILE = Fname(fid))
    DO i_ex_state = 1, n_ex_state
      WRITE(fid,'(A)') '$EX_STATE'
      WRITE(fid,'(A9, I2)') ' STATE = ', i_ex_state
      WRITE(fid,'(X, A)') mult_ex(i_ex_state)

!     Alpha electron
      IF(n_elec_config_x_a(i_ex_state) > 0) THEN
        DO i_elec_config = 1, n_elec_config_x_a(i_ex_state)
          IF(ABS(x_plus_y_a(i_elec_config, i_ex_state)) &
            &< threshold_ci_coef .AND. &
            &ABS(x_minus_y_a(i_elec_config, i_ex_state)) &
            &< threshold_ci_coef) CYCLE
          WRITE(fid,9011)&
          so_occ_a(i_elec_config, i_ex_state), &
            '->', &
          so_unocc_a(i_elec_config, i_ex_state), &
          x_plus_y_a(i_elec_config, i_ex_state), &
          x_minus_y_a(i_elec_config, i_ex_state)
        ENDDO ! i_elec_config
       
        IF(n_elec_config_y_only_a(i_ex_state) /= 0) THEN
          DO counter = 1, n_elec_config_y_only_a(i_ex_state)
            IF(ABS(x_plus_y_a(counter, i_ex_state)) &
              &< threshold_ci_coef .AND. &
              &ABS(x_minus_y_a(counter, i_ex_state)) &
              &< threshold_ci_coef) CYCLE
            i_elec_config = n_elec_config_x_a(i_ex_state) + counter  
            WRITE(fid,9011)&
            so_occ_a(i_elec_config, i_ex_state), &
            '->', &
            so_unocc_a(i_elec_config, i_ex_state), &
            x_plus_y_a(i_elec_config, i_ex_state), &
            x_minus_y_a(i_elec_config, i_ex_state)
          ENDDO ! i_elec_config
        ELSE
        ENDIF
      ELSE
      ENDIF
!     Beta electron
      IF(n_elec_config_x_b(i_ex_state) > 0) THEN
        DO i_elec_config = 1, n_elec_config_x_b(i_ex_state)
          IF(ABS(x_plus_y_b(i_elec_config, i_ex_state)) &
            &< threshold_ci_coef .AND. &
            &ABS(x_minus_y_b(i_elec_config, i_ex_state)) &
            &< threshold_ci_coef) CYCLE
          WRITE(fid,9011)&
          so_occ_b(i_elec_config, i_ex_state), &
            '->', &
          so_unocc_b(i_elec_config, i_ex_state), &
          x_plus_y_b(i_elec_config, i_ex_state), &
          x_minus_y_b(i_elec_config, i_ex_state)
        ENDDO ! i_elec_config
       
        IF(n_elec_config_y_only_b(i_ex_state) /= 0) THEN
          DO counter = 1, n_elec_config_y_only_b(i_ex_state)
            IF(ABS(x_plus_y_b(counter, i_ex_state)) &
              &< threshold_ci_coef .AND. &
              &ABS(x_minus_y_b(counter, i_ex_state)) &
              &< threshold_ci_coef) CYCLE
            i_elec_config = n_elec_config_x_b(i_ex_state) + counter
            WRITE(fid,9011)&
            so_occ_b(i_elec_config, i_ex_state), &
            '->', &
            so_unocc_b(i_elec_config, i_ex_state), &
            x_plus_y_b(i_elec_config, i_ex_state), &
            x_minus_y_b(i_elec_config, i_ex_state)
          ENDDO ! i_elec_config
        ELSE
        ENDIF
      ELSE
      ENDIF

      WRITE(fid,'(A)') '$END_EX_STATE'
      WRITE(fid,'(A)') ' '
    ENDDO ! i_ex_state
  CLOSE(fid)

  9001 FORMAT(1X, A)
  9011 FORMAT(I5, 2X, A2, I5, 2E20.6)

END SUBROUTINE gen_g16_cicoef_td_unres  
