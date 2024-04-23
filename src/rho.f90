! This subroutine is part of d77 and computes
! electron density or overlap (transition) density. 

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE rho
  USE ifmod, ONLY: write_messages, &
                  &calc_opr_dmat, &
                  &calc_dmat_cgf, &
                  &cube_rho
  USE global_constants
  USE global_read_input
  USE global_read_data
  IMPLICIT NONE

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'rho'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Local variables
! One-particle reduced density matrix
  INTEGER                       :: left_state, right_state, n_ci_ci
  DOUBLE PRECISION, ALLOCATABLE :: opr_dmat(:,:)
  DOUBLE PRECISION              :: trace_opr_dmat

! Density matrix
  DOUBLE PRECISION, ALLOCATABLE :: dmat_cgf(:,:)

  DOUBLE PRECISION :: int_rho

  CALL write_messages(2, Text_blank, type_program, name_program)

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '--------------'
  WRITE(6,'(1X, A)') 'Density matrix'
  WRITE(6,'(1X, A)') '--------------'

! -------------------------------------------------------
! Calculating one particle reduced density matrix (gamma)
! -------------------------------------------------------

  left_state = 0; right_state = 0
  IF(Bra <= Ket) THEN
    left_state  = Bra
    right_state = Ket
  ELSE
    left_state  = Ket
    right_state = Bra
  ENDIF

! dmat_cgf is calculated from opr_dmat.
  ALLOCATE(dmat_cgf(1:N_cgf, 1:N_cgf)); dmat_cgf = 0.0D0
  ALLOCATE(opr_dmat(1:N_so, 1:N_so));   opr_dmat = 0.0D0
  n_ci_ci = 0
  trace_opr_dmat = 0.0D0
  CALL calc_opr_dmat&
      &(left_state, right_state, n_ci_ci, opr_dmat, trace_opr_dmat)

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'gamma: one-particle reduced density matrix'
  WRITE(6,'(1X, A)') 'between electronic states | i > and | j >'
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Trace of gamma'
  WRITE(6,'(1X, A)') REPEAT('=', 29)
  WRITE(6,'(1X, A3, 1X, A3, 3X, A7, 3X, A)') 'i', 'j', 'n_ci_ci', 'Tr[gamma]'
  WRITE(6,'(1X, A)') REPEAT('-', 29)
  WRITE(6,'(1X, I3, 1X, I3, 3X, I5, 2X, F12.5)') Bra, Ket, n_ci_ci, trace_opr_dmat
  WRITE(6,'(1X, A)') REPEAT('=', 29)

  CALL calc_dmat_cgf(opr_dmat, dmat_cgf)
  int_rho = 0.0D0

  WRITE(6,'(1X)')
  IF(Bra == Ket) THEN
    WRITE(6,'(1X, A)') '----------------'
    WRITE(6,'(1X, A)') 'Electron density'
    WRITE(6,'(1X, A)') '----------------'
  ELSE
    WRITE(6,'(1X, A)') '----------------------------'
    WRITE(6,'(1X, A)') 'Overlap (Transition) density'
    WRITE(6,'(1X, A)') '----------------------------'
  ENDIF
  CALL cube_rho(dmat_cgf, int_rho)
  WRITE(6,'(1X)')
  IF(Bra == Ket) THEN
    WRITE(6,'(1X, A)') '----------------------------------------'
    WRITE(6,'(1X, A, 3F10.5)') 'Integral of electron density: ', int_rho
    WRITE(6,'(1X, A)') '----------------------------------------'
    WRITE(6,'(1X, A, I0, A)') 'This value should be approximately equal to ', Ne, '.'
  ELSE
    WRITE(6,'(1X, A)') '-----------------------------------------------------'
    WRITE(6,'(1X, A, 3F10.5)') 'Integral of overlap (transition) density:  ', int_rho
    WRITE(6,'(1X, A)') '-----------------------------------------------------'
    WRITE(6,'(1X, A)') 'This value should be approximately equal to zero.'
  ENDIF

  DEALLOCATE (Atmnum, Nuccharg, Xyznuc, Elmsym)
  DEALLOCATE (Lmn, Xyzao)
  DEALLOCATE (Coef_so)

  CALL write_messages(3, Text_blank, type_program, name_program)

END SUBROUTINE rho
