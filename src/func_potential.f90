! This subroutine is part of d77 and computes potential.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE func_potential
  USE global_read_data
  IMPLICIT NONE
  CONTAINS

  FUNCTION func_vne(xyz) RESULT(f)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'func_vne'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'FUNCTION'
    
!   Input variable
    DOUBLE PRECISION, INTENT(IN) :: xyz(1:3)
    
!   Output variable
    DOUBLE PRECISION :: f

!   Local variables
    INTEGER i_atm
    DOUBLE PRECISION :: dxyz(1:3)
    DOUBLE PRECISION :: norm_dxyz
    
    dxyz = 0.0D0; norm_dxyz = 0.0D0
    f = 0.0D0

    DO i_atm = 1, N_atm
      dxyz(1:3) = xyz(1:3) - Xyznuc(1:3,i_atm)
      norm_dxyz = DSQRT(DOT_PRODUCT(dxyz, dxyz))
      f = f - nuccharg(i_atm) / norm_dxyz
    ENDDO
  END FUNCTION func_vne

  FUNCTION func_dvne(xyz, num_mode) RESULT(f)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'func_dvne'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'FUNCTION'
    
!   Input variable
    DOUBLE PRECISION :: xyz(:)
    INTEGER          :: num_mode
    
!   Output variable
    DOUBLE PRECISION f(1:num_mode)
    
!   Local variables
    INTEGER i_atm, i_mode
    DOUBLE PRECISION, DIMENSION(1:3) :: dxyz
    DOUBLE PRECISION :: norm_dxyz
    
    dxyz = 0.0D0; norm_dxyz = 0.0D0
    f = 0.0D0

    DO i_mode = 1, num_mode
      DO i_atm = 1, N_atm
        dxyz(:) = xyz(:) - Xyznuc(:,i_atm)
        norm_dxyz = DSQRT(DOT_PRODUCT(dxyz, dxyz))
        f(i_mode)&
        = f(i_mode)&
        - Nuccharg(i_atm) / Sqrt_atmwt(i_atm)&
        * DOT_PRODUCT(Vibmode_calc(:,i_atm, i_mode), dxyz)&
        / (norm_dxyz**3)
      ENDDO
    ENDDO
  END FUNCTION func_dvne

END MODULE func_potential
