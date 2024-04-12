! The main program of dpp
!
! dpp is a data preprocessing utility for d77. 
!
! dpp is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

PROGRAM main
  USE global_constants
  USE fchk2inputs
  USE gen_cicoef
  IMPLICIT NONE

  CHARACTER(LEN=30) :: file_typ, property
  CHARACTER(LEN=10) :: qc_program_elec, qc_program_vib, qc_method
  DOUBLE PRECISION :: threshold_ci_coef 
  INTEGER :: fio

  CALL assign_globals

! Initializing variables

  file_typ        = REPEAT(' ', 20)
  property        = REPEAT(' ', 20)
  qc_program_elec = REPEAT(' ', 10)
  qc_program_vib  = REPEAT(' ', 10)
  qc_method       = REPEAT(' ', 10)
  threshold_ci_coef = Threshold_ci_coef_default 

  fio = Ffile_typ
  OPEN(fio, FILE = Fname(fio))
    READ(fio,*) file_typ
  CLOSE(fio)
  WRITE(*,'(X, A, A)') 'File type = ', file_typ

      fio = FR_control_elec
      OPEN(fio, FILE = Fname(fio))
        READ(fio,*) qc_program_elec
        READ(fio,*) qc_program_vib
        READ(fio,*) qc_method
        READ(fio,*) property
      CLOSE(fio)

! Reading ground-state electronic structure
  SELECT CASE(file_typ)
    CASE('elec_fchk')
      CALL elec_fchk(qc_program_elec)
    CASE('vib_fchk')
      CALL vib_fchk
    CASE('g16_log')
      SELECT CASE(qc_method)
        CASE('cis', 'CIS')
          CALL log2cicoef_ci
        CASE DEFAULT
          WRITE(*,*) 'Incorrect qc_method'
          WRITE(*,*) 'qc_method = ', qc_method
          STOP
      END SELECT ! qc_method
    CASE('g16_cisd_log')
      CALL gen_g16_cicoef_cisd
    CASE('g16_td_log')
      CALL gen_g16_cicoef_td
    CASE('g16_td_log_unres')
      CALL gen_g16_cicoef_td_unres
    CASE DEFAULT
      WRITE(*,*) 'Incorrect file_typ'
      WRITE(*,*) 'file_typ = ', file_typ 
      STOP
  END SELECT ! file_typ 


END PROGRAM main
