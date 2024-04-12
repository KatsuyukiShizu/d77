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

SUBROUTINE gen_g16_cicoef_cisd
  USE global_constants
  USE read_txt
  IMPLICIT NONE

  INTEGER :: threshold_tuple 
  DOUBLE PRECISION :: threshold_ci_coef 

  INTEGER :: N, K_divided_by_2, n_state

  CHARACTER(LEN=100) :: text(1:10), text_AB, text_dummy 
  CHARACTER(LEN=4), ALLOCATABLE :: AB_as_read(:), AB(:) 
  INTEGER :: mo_occ(1:2), mo_unocc(1:2)
  INTEGER, ALLOCATABLE :: so_occ(:,:), so_unocc(:,:)

  INTEGER :: i_space(1:2)

  INTEGER :: i_elec_config, n_elec_config_as_read, n_elec_config, counter
  DOUBLE PRECISION :: ci_coef_as_read
  DOUBLE PRECISION, ALLOCATABLE :: ci_coef(:)
  INTEGER :: fio

  i_space = 0

  N = 0; K_divided_by_2 = 0; n_state = 0
  n_elec_config = 0

  mo_occ = 0; mo_unocc = 0


  threshold_tuple = threshold_tuple_default
  threshold_ci_coef = threshold_ci_coef_default
  OPEN(FR_threshold_values_ci, FILE = Fname(FR_threshold_values_ci))
    READ(FR_threshold_values_ci, *) threshold_ci_coef
  CLOSE(FR_threshold_values_ci)

  counter = 0
  fio = FR_g16_cisd_log
  OPEN(fio, FILE = Fname(fio))
    READ(fio,*) N
    READ(fio,*) K_divided_by_2
    READ(fio,*) n_state
    READ(fio,*) n_elec_config_as_read
    ALLOCATE(AB_as_read(1:n_elec_config_as_read)); AB_as_read = '    '
    DO i_elec_config = 1, n_elec_config_as_read       
      READ(fio,'(A)') text(1)
      i_space(1) = INDEX(text(1), '           ')
      CALL get_text_left_symbol&
          &(text(1), i_space(1), text_AB)
      IF(INDEX(text_AB, 'AAAA') /= 0) THEN
        AB_as_read(i_elec_config) = 'AAAA'
        counter = counter + 2
      ELSEIF(INDEX(text_AB, 'ABAB') /= 0) THEN
        AB_as_read(i_elec_config) = 'ABAB'
        counter = counter + 2
      ELSEIF(INDEX(text_AB, 'AA') /= 0) THEN
        AB_as_read(i_elec_config) = 'AA'
        counter = counter + 2
      ELSE
        STOP
      ENDIF
    ENDDO
  CLOSE(fio)

  n_elec_config = counter
  ALLOCATE(AB(1:n_elec_config)); AB = '    '
  ALLOCATE(so_occ(1:2, 1:n_elec_config)); so_occ = 0
  ALLOCATE(so_unocc(1:2, 1:n_elec_config)); so_unocc = 0
  ALLOCATE(ci_coef(1:n_elec_config)); ci_coef = 0.0D0

  counter = 0
  OPEN(fio, FILE = Fname(fio))
    READ(fio,*) 
    READ(fio,*) 
    READ(fio,*) 
    READ(fio,*) 
    DO i_elec_config = 1, n_elec_config_as_read
      IF(AB_as_read(i_elec_config) == 'AAAA') THEN
        READ(fio, *) text_dummy, &
                    &mo_occ(1), mo_occ(2), mo_unocc(1), mo_unocc(2), &
                    &ci_coef_as_read
        ci_coef_as_read = ci_coef_as_read / DSQRT(2.0D0)            

        counter = counter + 1
        AB(counter) = 'AAAA'
        so_occ(1, counter) = 2*mo_occ(1) - 1          
        so_occ(2, counter) = 2*mo_occ(2) - 1          
        so_unocc(1, counter) = 2*mo_unocc(1) - 1          
        so_unocc(2, counter) = 2*mo_unocc(2) - 1        
        ci_coef(counter) = ci_coef_as_read

        counter = counter + 1
        AB(counter) = 'BBBB'
        so_occ(1, counter) = 2*mo_occ(1)          
        so_occ(2, counter) = 2*mo_occ(2)          
        so_unocc(1, counter) = 2*mo_unocc(1)          
        so_unocc(2, counter) = 2*mo_unocc(2)        
        ci_coef(counter) = ci_coef_as_read
      ELSEIF(AB_as_read(i_elec_config) == 'ABAB') THEN
        READ(fio, *) text_dummy, &
                    &mo_occ(1), mo_occ(2), mo_unocc(1), mo_unocc(2), &
                    &ci_coef_as_read
          counter = counter + 1
          AB(counter) = 'ABAB'
          so_occ(1, counter) = 2*mo_occ(1) - 1
          so_occ(2, counter) = 2*mo_occ(2)
          so_unocc(1, counter) = 2*mo_unocc(1) - 1
          so_unocc(2, counter) = 2*mo_unocc(2)
          ci_coef(counter) = ci_coef_as_read
      ELSEIF(AB_as_read(i_elec_config) == 'AA') THEN
        READ(fio, *) text_dummy, &
                    &mo_occ(1), mo_unocc(1), &
                    &ci_coef_as_read
        ci_coef_as_read = ci_coef_as_read / DSQRT(2.0D0)            

        counter = counter + 1
        AB(counter) = 'AA'
        so_occ(1, counter) = 2*mo_occ(1) - 1
        so_unocc(1, counter) = 2*mo_unocc(1) - 1
        ci_coef(counter) = ci_coef_as_read

        counter = counter + 1
        AB(counter) = 'BB'
        so_occ(1, counter) = 2*mo_occ(1)
        so_unocc(1, counter) = 2*mo_unocc(1)
        ci_coef(counter) = ci_coef_as_read
      ELSE
        WRITE(*,*) 'STOP'
        STOP
      ENDIF
    ENDDO
  CLOSE(fio)

  OPEN(Fw_cicoef, FILE = Fname(Fw_cicoef))
    WRITE(Fw_cicoef,'(A)') '$EX_STATE'  
      WRITE(Fw_cicoef,'(A9, I2)') ' STATE = ', 1  
      WRITE(FW_cicoef,'(X, A)') 'Singlet'
    WRITE(Fw_cicoef,9010) 0, 1.0D0
    DO i_elec_config = 1, n_elec_config
      IF(ABS(ci_coef(i_elec_config)) < threshold_ci_coef) CYCLE
      SELECT CASE(AB(i_elec_config))
        CASE('AA', 'BB')
          WRITE(Fw_cicoef,9011) 1,&
          &so_occ(1:1, i_elec_config), '->', so_unocc(1:1, i_elec_config),&
         &ci_coef(i_elec_config)
        CASE('AAAA', 'BBBB', 'ABAB', 'BABA')
          WRITE(Fw_cicoef,9012) 2,&
          &so_occ(1:2, i_elec_config), '->', so_unocc(1:2, i_elec_config),&
          &ci_coef(i_elec_config)
        CASE DEFAULT  
      END SELECT
    ENDDO
    WRITE(Fw_cicoef,'(A)') '$END_EX_STATE'  
  CLOSE(Fw_cicoef)

  9010  FORMAT(I2, E20.6)
  9011  FORMAT(I2,  I5, 2X, A2,  I5, E20.6)
  9012  FORMAT(I2, 2I5, 2X, A2, 2I5, E20.6)

END SUBROUTINE gen_g16_cicoef_cisd  
