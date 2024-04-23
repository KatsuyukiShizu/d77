! This subroutine is part of d77 and computes 
! the d*e*f term in Eq (2.27) in the paper 
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

SUBROUTINE calc_def&
          &(xyz_a, xyz_b, xyz_p, &
           &lmn_a, lmn_b, cntexp_p, &
           &d, e, f)
  IMPLICIT NONE

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'calc_def'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Input variables
  DOUBLE PRECISION, INTENT(IN) :: xyz_a(1:3), xyz_b(1:3), xyz_p(1:3)
  INTEGER, INTENT(IN)          :: lmn_a(1:3), lmn_b(1:3)
  DOUBLE PRECISION, INTENT(IN) :: cntexp_p

! Output variables
  DOUBLE PRECISION, INTENT(OUT) :: d(0:lmn_a(1)+lmn_b(1), 0:lmn_a(1), 0:lmn_b(1)) 
  DOUBLE PRECISION, INTENT(OUT) :: e(0:lmn_a(2)+lmn_b(2), 0:lmn_a(2), 0:lmn_b(2)) 
  DOUBLE PRECISION, INTENT(OUT) :: f(0:lmn_a(3)+lmn_b(3), 0:lmn_a(3), 0:lmn_b(3)) 

! Local variables
  DOUBLE PRECISION :: const_p
  INTEGER          :: l, m, n, i_a, i_b
  DOUBLE PRECISION :: d_temp(-1:lmn_a(1)+lmn_b(1)+2, 0:lmn_a(1), 0:lmn_b(1)), & 
                     &e_temp(-1:lmn_a(2)+lmn_b(2)+2, 0:lmn_a(2), 0:lmn_b(2)), & 
                     &f_temp(-1:lmn_a(3)+lmn_b(3)+2, 0:lmn_a(3), 0:lmn_b(3)) 
  
  const_p = 0.5D0/cntexp_p

  d = 0.0D0; e = 0.0D0; f = 0.0D0
  d_temp = 0.0D0; e_temp = 0.0D0; f_temp = 0.0D0

! Calculating d(L; lmn_a(1), lmn_b(1)), e(M; lmn_a(2), lmn_b(2)), or f(N; lmn_a(3), lmn_b(3))
  d(0, 0, 0) = 1.0D0; e(0, 0, 0) = 1.0D0; f(0, 0, 0) = 1.0D0
  d_temp(0, 0, 0) = 1.0D0; e_temp(0, 0, 0) = 1.0D0; f_temp(0, 0, 0) = 1.0D0

! i_b = 0
  IF(lmn_a(1) >= 1) THEN
    DO i_a = 0, lmn_a(1)-1
      DO l = 0, i_a+1
        d_temp(l, i_a+1, 0) = d_temp(l+1, i_a, 0) * DBLE(l+1) &
                           &+ d_temp(l,   i_a, 0) * (xyz_p(1) - xyz_a(1)) &
                           &+ d_temp(l-1, i_a, 0) * const_p
      ENDDO ! l
    ENDDO ! i_a
  ELSE ! lmn_a(1) = 0
  ENDIF

  IF(lmn_b(1) >= 1) THEN
    DO i_b = 0, lmn_b(1)-1
      DO i_a = 0, lmn_a(1)
        DO l = 0, i_a+i_b+1
          d_temp(l, i_a, i_b+1) = d_temp(l+1, i_a, i_b) * DBLE(l+1) &
                               &+ d_temp(l,   i_a, i_b) * (xyz_p(1) - xyz_b(1)) &
                               &+ d_temp(l-1, i_a, i_b) * const_p
        ENDDO ! l
      ENDDO ! i_a
    ENDDO ! i_b
  ELSE ! lmn_b(1) = 0
  ENDIF

! Calculating e
  ! i_b = 0
  IF(lmn_a(2) >= 1) THEN
    DO i_a = 0, lmn_a(2)-1
      DO m = 0, i_a+1
        e_temp(m, i_a+1, 0) = e_temp(m+1, i_a, 0) * DBLE(m+1) &
                           &+ e_temp(m,   i_a, 0) * (xyz_p(2) - xyz_a(2)) &
                           &+ e_temp(m-1, i_a, 0) * const_p
      ENDDO ! m
    ENDDO ! i_a
  ELSE ! lmn_a(2) = 0
  ENDIF

  IF(lmn_b(2) >= 1) THEN
    DO i_b = 0, lmn_b(2)-1
      DO i_a = 0, lmn_a(2)
        DO m = 0, i_a+i_b+1
          e_temp(m, i_a, i_b+1) = e_temp(m+1, i_a, i_b) * DBLE(m+1) &
                               &+ e_temp(m,   i_a, i_b) * (xyz_p(2) - xyz_b(2)) &
                               &+ e_temp(m-1, i_a, i_b) * const_p
        ENDDO ! m
      ENDDO ! i_a
    ENDDO ! i_b
  ELSE ! lmn_b(2) = 0
  ENDIF


  ! i_b = 0
  IF(lmn_a(3) >= 1) THEN
    DO i_a = 0, lmn_a(3)-1
      DO n = 0, i_a+1
        f_temp(n, i_a+1, 0) = f_temp(n+1, i_a, 0) * DBLE(n+1) &
                           &+ f_temp(n,   i_a, 0) * (xyz_p(3) - xyz_a(3)) &
                           &+ f_temp(n-1, i_a, 0) * const_p
      ENDDO ! n
    ENDDO ! i_a
  ELSE ! lmn_a(3) = 0
  ENDIF

  IF(lmn_b(3) >= 1) THEN
    DO i_b = 0, lmn_b(3)-1
      DO i_a = 0, lmn_a(3)
        DO n = 0, i_a+i_b+1
          f_temp(n, i_a, i_b+1) = f_temp(n+1, i_a, i_b) * DBLE(n+1) &
                               &+ f_temp(n,   i_a, i_b) * (xyz_p(3) - xyz_b(3)) &
                               &+ f_temp(n-1, i_a, i_b) * const_p
        ENDDO ! n
      ENDDO ! i_a
    ENDDO ! i_b
  ELSE ! lmn_b(3) = 0
  ENDIF

  d(0:lmn_a(1)+lmn_b(1), 0:lmn_a(1), 0:lmn_b(1)) = d_temp(0:lmn_a(1)+lmn_b(1), 0:lmn_a(1), 0:lmn_b(1))
  e(0:lmn_a(2)+lmn_b(2), 0:lmn_a(2), 0:lmn_b(2)) = e_temp(0:lmn_a(2)+lmn_b(2), 0:lmn_a(2), 0:lmn_b(2)) 
  f(0:lmn_a(3)+lmn_b(3), 0:lmn_a(3), 0:lmn_b(3)) = f_temp(0:lmn_a(3)+lmn_b(3), 0:lmn_a(3), 0:lmn_b(3))

END SUBROUTINE calc_def
