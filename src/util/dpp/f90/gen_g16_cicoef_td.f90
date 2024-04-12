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

SUBROUTINE gen_g16_cicoef_td
  USE global_constants
  USE read_txt
  IMPLICIT NONE

  DOUBLE PRECISION :: threshold_ci_coef 

  INTEGER :: N, K_divided_by_2, n_ex_state
  CHARACTER(LEN=10) mult_ex(1:Max_n_state)

! Dummy variables
  INTEGER :: mo_occ_x(1:Max_n_elec_config,1:Max_n_state), &
            &mo_unocc_x(1:Max_n_elec_config,1:Max_n_state), &
            &mo_occ_y(1:Max_n_elec_config,1:Max_n_state), &
            &mo_unocc_y(1:Max_n_elec_config,1:Max_n_state)
  CHARACTER(LEN=2) :: arrow_x_as_read(1:Max_n_elec_config,1:Max_n_state), &
                      &arrow_y_as_read(1:Max_n_elec_config,1:Max_n_state)
  DOUBLE PRECISION :: x(1:Max_n_elec_config,1:Max_n_state), &
                     &y(1:Max_n_elec_config,1:Max_n_state), &
                     &yy(1:Max_n_elec_config,1:Max_n_state)
  CHARACTER(LEN=7) :: label_x(1:Max_n_elec_config,1:Max_n_state), &
                     &label_y(1:Max_n_elec_config,1:Max_n_state)


  INTEGER, ALLOCATABLE :: so_occ(:,:), so_unocc(:,:)
  CHARACTER(2), ALLOCATABLE :: arrow(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: x_plus_y(:,:), x_minus_y(:,:)

  INTEGER :: i_space(1:2)

  INTEGER :: i_elec_config, i_elec_config_x, i_elec_config_y
  INTEGER :: n_elec_config_x(1:Max_n_state), &
            &n_elec_config_y(1:Max_n_state), &
            &n_elec_config_x_only(1:Max_n_state), &
            &n_elec_config_y_only(1:Max_n_state), &
            &n_elec_config_x_and_y(1:Max_n_state)
  INTEGER :: i_ex_state
  INTEGER :: fid
  INTEGER :: counter

  i_space = 0

  N = 0; K_divided_by_2 = 0; n_ex_state = 0
  n_elec_config_x = 0; n_elec_config_y = 0
  n_elec_config_x_only = 0; n_elec_config_y_only = 0
  n_elec_config_x_and_y = 0

  mo_occ_x = 0; mo_unocc_x = 0
  mo_occ_y = 0; mo_unocc_y = 0
  x = 0.0D0; y = 0.0D0; yy = 0.0D0
  arrow_x_as_read = '  '; arrow_y_as_read = '  '
  label_x = 'x_only'; label_y = 'y_only'

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
      READ(fid,*) n_elec_config_x(i_ex_state)
      DO i_elec_config = 1, n_elec_config_x(i_ex_state)
        READ(fid,*) mo_occ_x(i_elec_config, i_ex_state),&
                   &mo_unocc_x(i_elec_config, i_ex_state),&
                   &x(i_elec_config, i_ex_state)
      ENDDO
    ENDDO
  CLOSE(fid)

  fid = FR_g16_td_log_y
  OPEN(fid, FILE = Fname(fid))
    READ(fid,*) 
    DO i_ex_state = 1, n_ex_state
      READ(fid,*) 
      READ(fid,*) n_elec_config_y(i_ex_state)
      DO i_elec_config = 1, n_elec_config_y(i_ex_state)
        READ(fid,*) mo_occ_y(i_elec_config, i_ex_state),&
                   &mo_unocc_y(i_elec_config, i_ex_state),&
                   &y(i_elec_config, i_ex_state)
      ENDDO
    ENDDO
  CLOSE(fid)

  DO i_ex_state = 1, n_ex_state
    DO i_elec_config_x = 1, n_elec_config_x(i_ex_state)

      DO i_elec_config_y = 1, n_elec_config_y(i_ex_state)
        IF(   mo_occ_x(i_elec_config_x, i_ex_state) &
          &== mo_occ_y(i_elec_config_y, i_ex_state) .AND. &
             &mo_unocc_x(i_elec_config_x, i_ex_state) &
          &== mo_unocc_y(i_elec_config_y, i_ex_state)) THEN

          label_x(i_elec_config_x, i_ex_state) = 'x_and_y'
          label_y(i_elec_config_y, i_ex_state) = 'x_and_y'
          yy(i_elec_config_x, i_ex_state) &
         &= y(i_elec_config_y, i_ex_state)
        ELSE
        ENDIF
      ENDDO

    ENDDO
  ENDDO


  ALLOCATE(so_occ(1:2*Max_n_elec_config,1:n_ex_state)); so_occ = 0
  ALLOCATE(so_unocc(1:2*Max_n_elec_config,1:n_ex_state)); so_unocc = 0
  ALLOCATE(arrow(1:2*Max_n_elec_config,1:n_ex_state)); arrow = '  '
  ALLOCATE(x_plus_y(1:2*Max_n_elec_config,1:n_ex_state)); x_plus_y = 0.0D0
  ALLOCATE(x_minus_y(1:2*Max_n_elec_config,1:n_ex_state)); x_minus_y= 0.0D0

!  IF(Ms == 0) THEN

    DO i_ex_state = 1, n_ex_state
      DO i_elec_config_x = 1, n_elec_config_x(i_ex_state)
    
!       Alpha electron
        so_occ(2*i_elec_config_x-1, i_ex_state)&
        = 2*mo_occ_x(i_elec_config_x, i_ex_state)-1
        so_unocc(2*i_elec_config_x-1, i_ex_state)&
        = 2*mo_unocc_x(i_elec_config_x, i_ex_state)-1
        arrow(2*i_elec_config_x-1, i_ex_state) &
       &= arrow_x_as_read(i_elec_config_x, i_ex_state)
        x_plus_y(2*i_elec_config_x-1, i_ex_state) &
       &= x(i_elec_config_x, i_ex_state)
        x_minus_y(2*i_elec_config_x-1, i_ex_state) &
       &= x(i_elec_config_x, i_ex_state)
    
!       Beta electron
        so_occ(2*i_elec_config_x, i_ex_state)&
        = 2*mo_occ_x(i_elec_config_x, i_ex_state)
        so_unocc(2*i_elec_config_x, i_ex_state)&
        = 2*mo_unocc_x(i_elec_config_x, i_ex_state)
        arrow(2*i_elec_config_x, i_ex_state)&
        = arrow_x_as_read(i_elec_config_x, i_ex_state)
        x_plus_y(2*i_elec_config_x, i_ex_state) &
       &= x(i_elec_config_x, i_ex_state)
        x_minus_y(2*i_elec_config_x, i_ex_state) &
       &= x(i_elec_config_x, i_ex_state)
    
    
        SELECT CASE(label_x(i_elec_config_x, i_ex_state))
          CASE('x_and_y')
            n_elec_config_x_and_y(i_ex_state) &
           &= n_elec_config_x_and_y(i_ex_state) + 1 
    
!           Alpha electron
            x_plus_y(2*i_elec_config_x-1, i_ex_state) &
           &= x_plus_y(2*i_elec_config_x-1, i_ex_state) &
           &+ yy(i_elec_config_x, i_ex_state)
    
            x_minus_y(2*i_elec_config_x-1, i_ex_state) &
           &= x_minus_y(2*i_elec_config_x-1, i_ex_state) &
           &- yy(i_elec_config_x, i_ex_state)
    
!           Beta electron
            x_plus_y(2*i_elec_config_x, i_ex_state) &
           &= x_plus_y(2*i_elec_config_x, i_ex_state) &
           &+ yy(i_elec_config_x, i_ex_state)
    
            x_minus_y(2*i_elec_config_x, i_ex_state) &
           &= x_minus_y(2*i_elec_config_x, i_ex_state) &
           &- yy(i_elec_config_x, i_ex_state)
          CASE DEFAULT
        END SELECT  

        SELECT CASE(mult_ex(i_ex_state))
          CASE('Triplet')
            x_plus_y(2*i_elec_config_x, i_ex_state) &
           &= - x_plus_y(2*i_elec_config_x, i_ex_state)
            x_minus_y(2*i_elec_config_x, i_ex_state) &
           &= - x_minus_y(2*i_elec_config_x, i_ex_state)
          CASE DEFAULT
        END SELECT
      ENDDO
    ENDDO

  n_elec_config_x_only = n_elec_config_x - n_elec_config_x_and_y

  DO i_ex_state = 1, n_ex_state
    counter = n_elec_config_x(i_ex_state) 
    DO i_elec_config_y = 1, n_elec_config_y(i_ex_state)
      SELECT CASE(label_y(i_elec_config_y, i_ex_state))
        CASE('y_only')
          counter = counter + 1

!         Alpha electron
          so_occ(2*counter-1, i_ex_state)&
          = 2*mo_occ_y(i_elec_config_y, i_ex_state)-1
          so_unocc(2*counter-1, i_ex_state)&
          = 2*mo_unocc_y(i_elec_config_y, i_ex_state)-1
          x_plus_y(2*counter-1, i_ex_state) &
         &= y(i_elec_config_y, i_ex_state)
          x_minus_y(2*counter-1, i_ex_state) &
         &= - y(i_elec_config_y, i_ex_state)
      
!         Beta electron
          so_occ(2*counter, i_ex_state)&
          = 2*mo_occ_y(i_elec_config_y, i_ex_state)
          so_unocc(2*counter, i_ex_state)&
          = 2*mo_unocc_y(i_elec_config_y, i_ex_state)
          x_plus_y(2*counter, i_ex_state) &
         &= y(i_elec_config_y, i_ex_state)
          x_minus_y(2*counter, i_ex_state) &
         &= - y(i_elec_config_y, i_ex_state)
          
      CASE DEFAULT
    END SELECT  
    ENDDO
    n_elec_config_y_only(i_ex_state) &
   &= counter - n_elec_config_x(i_ex_state)
  ENDDO

  fid = FW_xy_count
  OPEN(fid, FILE = Fname(fid))
    DO i_ex_state = 1, n_ex_state
      WRITE(fid,'(A8, I2)') 'STATE = ', i_ex_state
      WRITE(fid,'(A9, I2)') 'X only : ', n_elec_config_x_only(i_ex_state)
      WRITE(fid,'(A9, I2)') 'Y only : ', n_elec_config_y_only(i_ex_state)
      WRITE(fid,'(A9, I2)') 'X and Y: ', n_elec_config_x_and_y(i_ex_state)
      WRITE(fid,'(A9, I2)') 'X total: ', n_elec_config_x(i_ex_state)
      WRITE(fid,'(A9, I2)') 'y total: ', n_elec_config_y(i_ex_state)
      WRITE(fid, '(1X)') 
    ENDDO
  CLOSE(fid)

  fid = FW_xy
  OPEN(fid, FILE = Fname(fid))
    DO i_ex_state = 1, n_ex_state
      WRITE(fid,'(A)') '$EX_STATE'
      WRITE(fid,'(A9, I2)') ' STATE = ', i_ex_state
      WRITE(fid,'(X, A)') mult_ex(i_ex_state)

      DO i_elec_config = 1, 2*n_elec_config_x(i_ex_state)
        IF(ABS(x_plus_y(i_elec_config, i_ex_state)) &
          &< threshold_ci_coef .AND. &
          &ABS(x_minus_y(i_elec_config, i_ex_state)) &
          &< threshold_ci_coef) CYCLE
        WRITE(fid,9011)&
        so_occ(i_elec_config, i_ex_state), &
          '->', &
        !arrow(i_elec_config, i_ex_state), &
        so_unocc(i_elec_config, i_ex_state), &
        x_plus_y(i_elec_config, i_ex_state), &
        x_minus_y(i_elec_config, i_ex_state)
      ENDDO ! i_elec_config

      IF(n_elec_config_y_only(i_ex_state) /= 0) THEN
        DO counter = 1, 2*n_elec_config_y_only(i_ex_state)
          IF(ABS(x_plus_y(counter, i_ex_state)) &
            &< threshold_ci_coef .AND. &
            &ABS(x_minus_y(counter, i_ex_state)) &
            &< threshold_ci_coef) CYCLE
          i_elec_config = 2*n_elec_config_x(i_ex_state) + counter  
          WRITE(fid,9011)&
          so_occ(i_elec_config, i_ex_state), &
          '->', &
          so_unocc(i_elec_config, i_ex_state), &
          x_plus_y(i_elec_config, i_ex_state), &
          x_minus_y(i_elec_config, i_ex_state)
        ENDDO ! i_elec_config
      ELSE
      ENDIF

      WRITE(fid,'(A)') '$END_EX_STATE'
      WRITE(fid,'(A)') ' '
    ENDDO ! i_ex_state
  CLOSE(fid)

  9011 FORMAT(I5, 2X, A2, I5, 2E20.6)

END SUBROUTINE gen_g16_cicoef_td  
