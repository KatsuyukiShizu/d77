! This subroutine is part of d77 and prints elements of matrices.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE write_mat(fid, data)
  USE global_constants
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'SUBROUTINE'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'write_mat'

! Input variables
  INTEGER, INTENT(IN) :: fid
  DOUBLE PRECISION, INTENT(IN) :: data(:,:)

! Local variables
  INTEGER :: n_elm, i_elm, j_elm

  n_elm = UBOUND(data, 1)
  OPEN(fid, FILE = Fname(fid))
    DO i_elm = 1, n_elm
      DO j_elm = 1, n_elm
        WRITE(fid,*) data(i_elm, j_elm)
      ENDDO
    ENDDO
  CLOSE(fid)  

END SUBROUTINE write_mat
