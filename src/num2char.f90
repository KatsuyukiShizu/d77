! This subroutine is part of d77 and converts integers to characters.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

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

