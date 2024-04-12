! This subroutine is part of d77 and computes
! density matrices between contracted Gaussian functions.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE calc_dmat_cgf(opr_dmat, dmat_cgf)
  USE ifmod, ONLY: write_messages, calc_opr_dmat
  USE global_constants
  USE global_read_input
  USE global_read_data
  !USE func_int_slater
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'calc_dmat_cgf'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Output variables
  DOUBLE PRECISION, INTENT(IN) :: opr_dmat(:,:)
  DOUBLE PRECISION, INTENT(OUT) :: dmat_cgf(:,:)

! Local variables
  INTEGER i_cgf, j_cgf, i_so, j_so

  dmat_cgf = 0.0D0

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Calculating density matrix in terms of CGFs'
  WRITE(6,'(1X, A)') 'from gamma and CGF coefficients'
! Converting P_so into density matrix in terms of CGFs (P_CGF)
!  !$omp parallel do private(i_so, j_so)
  DO i_cgf = 1, N_cgf
    DO j_cgf = 1, N_cgf
      DO i_so = 1, N_so
        DO j_so = 1, N_so
          dmat_cgf(i_cgf, j_cgf) &
         &= dmat_cgf(i_cgf, j_cgf) &
         &+ opr_dmat(j_so, i_so) &
         &* coef_so(j_cgf, i_so) * coef_so(i_cgf, j_so)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
!  !$omp end parallel do

  WRITE(6,'(1X, A)') 'Done'

END SUBROUTINE calc_dmat_cgf
