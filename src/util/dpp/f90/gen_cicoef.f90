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

MODULE gen_cicoef
  USE global_constants
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE log2cicoef_ci

    INTEGER N, K_divided_by_2
    CHARACTER(LEN=10) mult_ex(1:Max_n_state)
    DOUBLE PRECISION threshold_ci_coef 
  
    INTEGER so_occ(1:Max_n_elec_config,1:Max_n_state)
    INTEGER so_unocc(1:Max_n_elec_config,1:Max_n_state)
    CHARACTER(2) arrow(1:Max_n_elec_config,1:Max_n_state)
  
    INTEGER so_occ_tdlog(1:Max_n_elec_config,1:Max_n_state)
    INTEGER so_unocc_tdlog(1:Max_n_elec_config,1:Max_n_state)
    CHARACTER(2) arrow_tdlog(1:Max_n_elec_config,1:Max_n_state)
  
    INTEGER i_ex_state, n_ex_state
    INTEGER i_elec_config
    INTEGER n_elec_config(1:Max_n_state)
    INTEGER tuple(1:Max_n_elec_config,1:Max_n_state)
    DOUBLE PRECISION ci_coef(1:Max_n_elec_config,1:Max_n_state),&
               ci_coef_tdlog(1:Max_n_elec_config,1:Max_n_state)
    INTEGER :: Mult
    INTEGER :: fio



    N = 0; K_divided_by_2 = 0
    n_ex_state = 0
    n_elec_config = 0
    tuple = 0

    so_occ = 0; so_unocc = 0; arrow = '  '
    so_occ_tdlog = 0; so_unocc_tdlog = 0; arrow_tdlog = '  '
    ci_coef = 0.0D0; ci_coef_tdlog = 0.0D0
    fio = 0


    OPEN(10, FILE = 'G16_CIS')
      READ(10,*) N
      READ(10,*) K_divided_by_2
      READ(10,*) n_ex_state
      DO i_ex_state = 1, n_ex_state
        READ(10,*) mult_ex(i_ex_state)
        READ(10,*) n_elec_config(i_ex_state)
        DO i_elec_config = 1, n_elec_config(i_ex_state)
          READ(10,*) so_occ_tdlog(i_elec_config, i_ex_state),&
                     arrow_tdlog(i_elec_config, i_ex_state),&
                     so_unocc_tdlog(i_elec_config, i_ex_state),&
                     ci_coef_tdlog(i_elec_config, i_ex_state)
        ENDDO
      ENDDO
    CLOSE(10)
  
    
    Mult = 0
    !Mult = 2
    DO i_ex_state = 1, n_ex_state
      DO i_elec_config = 1, n_elec_config(i_ex_state)
        ! Alpha electron
        tuple(2*i_elec_config-1, i_ex_state) = 1
        so_occ(2*i_elec_config-1, i_ex_state)&
        = 2*so_occ_tdlog(i_elec_config, i_ex_state)-1
        so_unocc(2*i_elec_config-1, i_ex_state)&
        = 2*so_unocc_tdlog(i_elec_config, i_ex_state)-1
        arrow(2*i_elec_config-1, i_ex_state)&
        = arrow_tdlog(i_elec_config, i_ex_state)
        ci_coef(2*i_elec_config-1, i_ex_state)&
        = ci_coef_tdlog(i_elec_config, i_ex_state)
        !IF(Mult == 2) CYCLE
        ! Beta electron
        tuple(2*i_elec_config, i_ex_state) = 1
        so_occ(2*i_elec_config, i_ex_state)&
        = so_occ(2*i_elec_config-1, i_ex_state)+1
        so_unocc(2*i_elec_config, i_ex_state)&
        = so_unocc(2*i_elec_config-1, i_ex_state)+1
        arrow(2*i_elec_config, i_ex_state)&
        = arrow(i_elec_config, i_ex_state)
        ci_coef(2*i_elec_config, i_ex_state)&
        = ci_coef(2*i_elec_config-1, i_ex_state)
      ENDDO
      IF(mult_ex(i_ex_state) == 'Triplet') THEN
        DO i_elec_config = 1, n_elec_config(i_ex_state)
          ci_coef(2*i_elec_config, i_ex_state)&
          = - ci_coef(2*i_elec_config, i_ex_state)
        ENDDO
      ELSE
      ENDIF
    ENDDO

  threshold_ci_coef = Threshold_ci_coef_default
  fio = FR_threshold_values_ci
  OPEN(fio, FILE = Fname(fio))
    READ(fio, *) threshold_ci_coef
  CLOSE(fio)

    IF(Mult == 2) THEN
      OPEN(FW_cicoef, FILE = Fname(FW_cicoef))
        DO i_ex_state = 1, n_ex_state
          WRITE(FW_cicoef,'(A)') '$EX_STATE'
          WRITE(FW_cicoef,'(A9, I2)') ' STATE = ', i_ex_state
          WRITE(FW_cicoef,'(X, A)') mult_ex(i_ex_state)
          DO i_elec_config = 1, 2*n_elec_config(i_ex_state), 2
            IF(ABS(ci_coef(i_elec_config, i_ex_state)) < threshold_ci_coef) CYCLE
            WRITE(FW_cicoef,9011)&
            tuple(i_elec_config, i_ex_state),&
            so_occ(i_elec_config, i_ex_state),&
            arrow(i_elec_config, i_ex_state),&
            so_unocc(i_elec_config, i_ex_state),&
            ci_coef(i_elec_config, i_ex_state)
          ENDDO ! i_elec_config
          WRITE(FW_cicoef,'(A)') '$END_EX_STATE'
          WRITE(FW_cicoef,'(A)') ' '
        ENDDO ! i_ex_state
      CLOSE(FW_cicoef)
    ELSE
      OPEN(FW_cicoef, FILE = Fname(FW_cicoef))
        DO i_ex_state = 1, n_ex_state
          WRITE(FW_cicoef,'(A)') '$EX_STATE'
          WRITE(FW_cicoef,'(A9, I2)') ' STATE = ', i_ex_state
          WRITE(FW_cicoef,'(X, A)') mult_ex(i_ex_state)
          DO i_elec_config = 1, 2*n_elec_config(i_ex_state)
            IF(ABS(ci_coef(i_elec_config, i_ex_state)) < threshold_ci_coef) CYCLE 
            WRITE(FW_cicoef,9011)&
            tuple(i_elec_config, i_ex_state),&
            so_occ(i_elec_config, i_ex_state),&
            arrow(i_elec_config, i_ex_state),&
            so_unocc(i_elec_config, i_ex_state),&
            ci_coef(i_elec_config, i_ex_state)
          ENDDO ! i_elec_config
          WRITE(FW_cicoef,'(A)') '$END_EX_STATE'
          WRITE(FW_cicoef,'(A)') ' '
        ENDDO ! i_ex_state
      CLOSE(FW_cicoef)
    ENDIF
    9011 FORMAT(I2,  I5, 2X, A2,  I5, F12.7)
  
  END SUBROUTINE log2cicoef_ci


  SUBROUTINE out2cicoef_ci

    INTEGER :: N, K_divided_by_2
    DOUBLE PRECISION :: threshold_ci_coef 

    CHARACTER(LEN=10) mult_ex(1:Max_n_state)

    INTEGER so_occ(1:Max_n_elec_config,1:Max_n_state)
    INTEGER so_unocc(1:Max_n_elec_config,1:Max_n_state)
    CHARACTER(2) arrow(1:Max_n_elec_config,1:Max_n_state)

    INTEGER so_occ_tdlog(1:Max_n_elec_config,1:Max_n_state)
    INTEGER so_unocc_tdlog(1:Max_n_elec_config,1:Max_n_state)
    CHARACTER(2) arrow_tdlog(1:Max_n_elec_config,1:Max_n_state)

    INTEGER i_ex_state, n_ex_state
    INTEGER i_elec_config
    INTEGER n_elec_config(1:Max_n_state)
    INTEGER tuple(1:Max_n_elec_config,1:Max_n_state)
    DOUBLE PRECISION ci_coef(1:Max_n_elec_config,1:Max_n_state),&
               ci_coef_tdlog(1:Max_n_elec_config,1:Max_n_state)
    INTEGER :: fio

    N = 0; K_divided_by_2 = 0
    n_ex_state = 0
    n_elec_config = 0
    tuple = 0

    so_occ = 0; so_unocc = 0; arrow = '  '
    so_occ_tdlog = 0; so_unocc_tdlog = 0; arrow_tdlog = '  '
    ci_coef = 0.0D0; ci_coef_tdlog = 0.0D0

  threshold_ci_coef = Threshold_ci_coef_default 
  fio = FR_threshold_values_ci
  OPEN(fio, FILE = Fname(fio))
    READ(fio, *) threshold_ci_coef
  CLOSE(fio)


    OPEN(10, FILE = 'QCHEM_CIS')
      READ(10,*) N
      READ(10,*) K_divided_by_2
      READ(10,*) n_ex_state
      DO i_ex_state = 1, n_ex_state
        READ(10,*) mult_ex(i_ex_state)
        READ(10,*) n_elec_config(i_ex_state)
        DO i_elec_config = 1, n_elec_config(i_ex_state)
          READ(10,*) so_occ_tdlog(i_elec_config, i_ex_state),&
                     arrow_tdlog(i_elec_config, i_ex_state),&
                     so_unocc_tdlog(i_elec_config, i_ex_state),&
                     ci_coef_tdlog(i_elec_config, i_ex_state)
          so_unocc_tdlog(i_elec_config, i_ex_state) &
         &= so_unocc_tdlog(i_elec_config, i_ex_state) + N/2
          ci_coef_tdlog(i_elec_config, i_ex_state) &
          = ci_coef_tdlog(i_elec_config, i_ex_state)/DSQRT(2.0D0) ! QChem specific
        ENDDO
      ENDDO
    CLOSE(10)

    DO i_ex_state = 1, n_ex_state
      DO i_elec_config = 1, n_elec_config(i_ex_state)
        ! Alpha electron
        tuple(2*i_elec_config-1, i_ex_state) = 1
        so_occ(2*i_elec_config-1, i_ex_state)&
        = 2*so_occ_tdlog(i_elec_config, i_ex_state)-1
        so_unocc(2*i_elec_config-1, i_ex_state)&
        = 2*so_unocc_tdlog(i_elec_config, i_ex_state)-1
        arrow(2*i_elec_config-1, i_ex_state)&
        = arrow_tdlog(i_elec_config, i_ex_state)
        ci_coef(2*i_elec_config-1, i_ex_state)&
        = ci_coef_tdlog(i_elec_config, i_ex_state)
        ! Beta electron
        tuple(2*i_elec_config, i_ex_state) = 1
        so_occ(2*i_elec_config, i_ex_state)&
        = so_occ(2*i_elec_config-1, i_ex_state)+1
        so_unocc(2*i_elec_config, i_ex_state)&
        = so_unocc(2*i_elec_config-1, i_ex_state)+1
        arrow(2*i_elec_config, i_ex_state)&
        = arrow(i_elec_config, i_ex_state)
        ci_coef(2*i_elec_config, i_ex_state)&
        = ci_coef(2*i_elec_config-1, i_ex_state)
      ENDDO
      IF(mult_ex(i_ex_state) == 'Triplet') THEN
        DO i_elec_config = 1, n_elec_config(i_ex_state)
          ci_coef(2*i_elec_config, i_ex_state)&
          = - ci_coef(2*i_elec_config, i_ex_state)
        ENDDO
      ELSE
      ENDIF
    ENDDO

    OPEN(FW_cicoef, FILE = Fname(FW_cicoef))
      DO i_ex_state = 1, n_ex_state
        WRITE(FW_cicoef,'(A)') '$EX_STATE'
        WRITE(FW_cicoef,'(A9, I2)') ' STATE = ', i_ex_state
        WRITE(FW_cicoef,'(X, A)') mult_ex(i_ex_state)
        DO i_elec_config = 1, 2*n_elec_config(i_ex_state)
          IF(ABS(ci_coef(i_elec_config, i_ex_state)) < threshold_ci_coef) CYCLE 
          WRITE(FW_cicoef,9011)&
          tuple(i_elec_config, i_ex_state),&
          so_occ(i_elec_config, i_ex_state),&
          arrow(i_elec_config, i_ex_state),&
          so_unocc(i_elec_config, i_ex_state),&
          ci_coef(i_elec_config, i_ex_state)
        ENDDO ! i_elec_config
        WRITE(FW_cicoef,'(A)') '$END_EX_STATE'
        WRITE(FW_cicoef,'(A)') ' '
      ENDDO ! i_ex_state
    CLOSE(FW_cicoef)

    9011 FORMAT(I2,  I5, 2X, A2,  I5, F12.7)

  END SUBROUTINE out2cicoef_ci




  SUBROUTINE log2cicoef_td
  END SUBROUTINE log2cicoef_td


SUBROUTINE get_text_right_symbol&
           (text, i_symbol, text_right)
  CHARACTER(72), INTENT(IN) :: text
  INTEGER, INTENT(IN) :: i_symbol
  CHARACTER(72), INTENT(OUT) :: text_right

  ! Local variables
  INTEGER i

  text_right = REPEAT(' ',72)

  ! text ='   temperature    =    298.0  ...'
  DO i = i_symbol + 1, 72
    text_right(i:i) = text(i:i) ! text_right='... 298.0 ...'
  ENDDO
  text_right = ADJUSTL(text_right) ! text_right='298.0       ...'
  IF(LEN_TRIM(text_right) == 0) THEN
    !CALL write_error(1) ! 'Error termination'
    WRITE(6,'(A)') 'Error termination'
    WRITE(6,'(A)') 'Incorrect input file'
    WRITE(6,'(A)') 'Error detected in ', text
  ENDIF
END SUBROUTINE get_text_right_symbol

END MODULE gen_cicoef
