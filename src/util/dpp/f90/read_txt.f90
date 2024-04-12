! This module is part of dpp and reads texts in an input file.
!
! dpp is a data preprocessing utility for d77.
!
! dpp is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE read_txt
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE get_text_right_symbol&
             (text, i_symbol, text_right)
    CHARACTER(100), INTENT(IN) :: text
    INTEGER, INTENT(IN) :: i_symbol
    CHARACTER(100), INTENT(OUT) :: text_right

    ! Local variables
    INTEGER i

    text_right = REPEAT(' ',100)

    DO i = i_symbol + 1, 100
      text_right(i:i) = text(i:i) 
    ENDDO
    text_right = ADJUSTL(text_right) 
    IF(LEN_TRIM(text_right) == 0) THEN
      WRITE(6,'(A)') 'Error termination'
      WRITE(6,'(A)') 'Incorrect input file'
      WRITE(6,'(A)') 'Error detected in ', text
    ENDIF
  END SUBROUTINE get_text_right_symbol

  SUBROUTINE get_text_left_symbol&
             (text, i_symbol, text_left)
!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'get_text_left_symbol'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

    CHARACTER(LEN=100), INTENT(IN) :: text
    INTEGER, INTENT(IN) :: i_symbol
    CHARACTER(LEN=100), INTENT(OUT) :: text_left

!   Local variables
    INTEGER i

!   Initializing variables
    text_left = REPEAT(' ',100)

    DO i = 1, i_symbol - 1
      text_left(i:i) = text(i:i) 
    ENDDO
    text_left = ADJUSTL(text_left) 
  END SUBROUTINE get_text_left_symbol


END MODULE read_txt
