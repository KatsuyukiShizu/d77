! This subroutine is part of d77 and computes
! the R_NLM term in Eq (3.13) in the paper
! written by McMurchie and Davidson:
! McMurchie, L. E. & Davidson, E. R.,
! "One- and two-electron integrals over cartesian gaussian functions",
! J. Comput. Phys. 26, 218-231 (1978)
! DOI: https://doi.org/10.1016/0021-9991(78)90092-X

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE calc_rlmn_0&
          &(l, m, n, max_j, &
           &a, b, c, t, &
           &alpha, &
           &rlmn)
  USE ifmod, ONLY: write_messages
  USE func_int_pgf, ONLY: func_fj
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'calc_rlmn_0'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Input variables
  DOUBLE PRECISION, INTENT(IN) :: a, b, c, t
  DOUBLE PRECISION, INTENT(IN) :: alpha
  INTEGER, INTENT(IN)          :: l, m, n, max_j

! Output variables
  DOUBLE PRECISION, INTENT(OUT) :: rlmn

! Local variables
  DOUBLE PRECISION, DIMENSION(0:max_j)                :: fj
  DOUBLE PRECISION, DIMENSION(0:l, 0:m, 0:n, 0:max_j) :: rlmnj
  INTEGER :: j, i_l, i_m, i_n
  

  fj = 0.0D0
  rlmnj = 0.0D0; rlmn = 0.0D0

  fj(0:max_j) = func_fj(max_j, t)
  DO j = 0, max_j
    rlmnj(0, 0, 0, j) = ((-2.0D0*alpha)**DBLE(j))*fj(j) 
  ENDDO

! Calculating rlmn
  IF(max_j == 0) THEN
  ELSE
    IF(n >= 1) THEN
      DO i_n = 0, n - 1
        DO j = 0, max_j - (i_n + 1)
          IF(i_n == 0) THEN
            rlmnj(0, 0, i_n+1, j) = rlmnj(0, 0, i_n, j+1) * c
          ELSEIF(i_n >= 1) THEN
            rlmnj(0, 0, i_n+1, j) = rlmnj(0, 0, i_n,   j+1) * c &
                                 &+ rlmnj(0, 0, i_n-1, j+1) * DBLE(i_n)
          ELSE
          ENDIF
        ENDDO
      ENDDO
    ELSE
    ENDIF
  
    IF(m >= 1) THEN
      DO i_m = 0, m-1
        DO j = 0, max_j - n - (i_m + 1)
          IF(i_m == 0) THEN
            rlmnj(0, i_m+1, n, j) = rlmnj(0, i_m,   n, j+1) * b
          ELSEIF(i_m >= 1) THEN
            rlmnj(0, i_m+1, n, j) = rlmnj(0, i_m,   n, j+1) * b &
                                 &+ rlmnj(0, i_m-1, n, j+1) * DBLE(i_m)
          ELSE
          ENDIF
        ENDDO
      ENDDO
    ELSE
    ENDIF
  
    IF(l >= 1) THEN
      DO i_l = 0, l-1
        DO j = 0, max_j - n - m - (i_l+1)
          IF(i_l == 0) THEN
            rlmnj(i_l+1, m, n, j) = rlmnj(i_l,   m, n, j+1) * a
          ELSEIF(i_l >= 1) THEN
            rlmnj(i_l+1, m, n, j) = rlmnj(i_l,   m, n, j+1) * a &
                                 &+ rlmnj(i_l-1, m, n, j+1) * DBLE(i_l)
          ELSE
          ENDIF
        ENDDO
      ENDDO
    ELSE
    ENDIF
  ENDIF  
  
  rlmn = rlmnj(l, m, n, 0)*((-1.0D0)**DBLE(max_j))

END SUBROUTINE calc_rlmn_0
