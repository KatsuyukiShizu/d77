! This subroutine is part of d77 and writes
! calculated density to a cube file.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE cube_write
  USE global_read_data
  IMPLICIT NONE

CONTAINS
  SUBROUTINE cube_write_header(fid, comment1, comment2)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'cube_write_header'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Input variables
    INTEGER, INTENT(IN)            :: fid
    CHARACTER(LEN=100), INTENT(IN) :: comment1, comment2
    
!   Local variables
    INTEGER i_xyz, i_atm
    
    WRITE(fid,'(A)') TRIM(ADJUSTL(comment1))
    WRITE(fid,'(A)') TRIM(ADJUSTL(comment2))
    WRITE(fid,'(I5, 3F12.6)') N_atm,  Xmin,  Ymin,  Zmin
    WRITE(fid,'(I5, 3F12.6)')    Nx,    Dx, 0.0D0, 0.0D0
    WRITE(fid,'(I5, 3F12.6)')    Ny, 0.0D0,    Dy, 0.0D0
    WRITE(fid,'(I5, 3F12.6)')    Nz, 0.0D0, 0.0D0,    Dz
    
    DO i_atm = 1, N_atm
      WRITE(fid, '(I5, 4F12.6)') Atmnum(i_atm), DBLE(Atmnum(i_atm)),&
                                &(Xyznuc(i_xyz, i_atm), i_xyz = 1,3)
    ENDDO ! i_atm = 1, n_atm
    
  END SUBROUTINE cube_write_header

  SUBROUTINE cube_write_density(fid, data)
    
!   Input variables
    INTEGER, INTENT(IN) :: fid
    DOUBLE PRECISION    :: data(:)
    
!   Local variables
    INTEGER :: icol, ncol ! icol = 1, 2, 3, 4, 5, 6; ncol = 6
    
    ncol = UBOUND(data, 1)
    DO icol = 1, ncol
      IF(ABS(data(icol)) < 1.0D-99) data(icol) = 0.0D0
    ENDDO
    
    SELECT CASE(ncol)
      CASE(1)
        WRITE(fid, '(E13.5)') data(1:ncol)
      CASE(2)
        WRITE(fid,'(2E13.5)') data(1:ncol)
      CASE(3)
        WRITE(fid,'(3E13.5)') data(1:ncol)
      CASE(4)
        WRITE(fid,'(4E13.5)') data(1:ncol)
      CASE(5)
        WRITE(fid,'(5E13.5)') data(1:ncol)
      CASE(6)
        WRITE(fid,'(6E13.5)') data(1:ncol)
      CASE DEFAULT
    END SELECT

  END SUBROUTINE cube_write_density

  SUBROUTINE num2char(num09, char09)
    INTEGER, INTENT(IN) :: num09
    CHARACTER(LEN=1), INTENT(OUT) :: char09
 
    char09 = 'N' ! temporary value (not a number)
    SELECT CASE(num09)
      CASE(0); char09 = '0'
      CASE(1); char09 = '1'
      CASE(2); char09 = '2'
      CASE(3); char09 = '3'
      CASE(4); char09 = '4'
      CASE(5); char09 = '5'
      CASE(6); char09 = '6'
      CASE(7); char09 = '7'
      CASE(8); char09 = '8'
      CASE(9); char09 = '9'
      CASE DEFAULT
        WRITE(6,*) num09, ' error: not integer/num2char'
        STOP
    END SELECT
  END SUBROUTINE num2char

END MODULE cube_write
