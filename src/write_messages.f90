! This subroutine is part of d77 and prints messages.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE write_messages(message_number, text, type_program, name_program)
  USE global_constants
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Input variables
  INTEGER, INTENT(IN) :: message_number
  CHARACTER(LEN=100), INTENT(IN) :: text, type_program, name_program

! Local variables
  CHARACTER(LEN=8)  :: date
  CHARACTER(LEN=10) :: time
  CHARACTER(LEN=5)  :: zone
  INTEGER :: date_time(1:8)

! Writing start time
  IF(message_number == 0) THEN
    WRITE(6,'(1X, A)') '****************************************************'
    WRITE(6,'(1X, A)') '**                                                **'
    WRITE(6,'(1X, A)') '**             d77: version '//TRIM(Program_ver)//'              **'
    WRITE(6,'(1X, A)') '**                                                **'
    WRITE(6,'(1X, A)') '**  Visual understanding of molecular properties  **'
    WRITE(6,'(1X, A)') '**                                                **'
    WRITE(6,'(1X, A)') '**        Coded by SHIZU Katsuyuki (Japan)        **'
    WRITE(6,'(1X, A)') '**                                                **'
    WRITE(6,'(1X, A)') '****************************************************'
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'd77 is free software and can be redistributed and/or modified'
    WRITE(6,'(1X, A)') 'under the terms of the GNU General Public License v3.0'
    WRITE(6,'(1X, A)') 'as published by the Free Software Foundation.'
    WRITE(6,'(1X, A)') 'https://www.gnu.org/licenses/gpl-3.0.html'
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp'
    WRITE(6,'(1X)')
    WRITE(6,'(1X)')
    WRITE(6,'(1X)')
  ELSEIF(message_number == 1 ) THEN
    CALL DATE_AND_TIME(date, time, zone, date_time)
    WRITE(6,'(1X, A, A, I4, A, I0, A, I0, A, I0, A, I0, A, I0, A, A, A)') ADJUSTL(TRIM(text)),&
   &' on ', date_time(1), '/', date_time(2), '/', date_time(3), &
   &' at ', date_time(5), ':', date_time(6), ':', date_time(7), &
   &' (time zone = UTC', zone, ')'
  ELSEIF(message_number == 2) THEN
    WRITE(6,'(1X)')
    WRITE(6,'(1X)')
    WRITE(6,'(1X, 2A, 1X, A)') '--- Entering ', TRIM(type_program), TRIM(name_program)
    !CALL FLUSH(6)
  ELSEIF(message_number == 3) THEN
    WRITE(6,'(1X)')
    WRITE(6,'(1X, 2A, 1X, A)') '--- Leaving ',  TRIM(type_program), TRIM(name_program)
    !CALL FLUSH(6)
  ELSEIF(message_number == 12) THEN
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'gamma_ij: one-particle reduced density matrix between'
    WRITE(6,'(1X, A)') 'electronic states |i> and |j> in terms of spin orbitals'
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Trace of gamma_ij'

  ELSEIF(message_number == 9999) THEN
    WRITE(6,'(1X)')
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '****************************'
    WRITE(6,'(1X, A)') '**                        **'
    WRITE(6,'(1X, A)') '**   Normal termination   **'
    WRITE(6,'(1X, A)') '**                        **'
    WRITE(6,'(1X, A)') '****************************'
    WRITE(6,'(1X)')
    WRITE(6,'(1X)')
    CALL CPU_TIME(Time_end)
    WRITE(6,'(1X, A, F8.2, A)') 'Total CPU time = ', Time_end-Time_start, ' seconds'
    CALL DATE_AND_TIME(date, time, zone, date_time)
    WRITE(6,'(1X, A, A, I4, A, I0, A, I0, A, I0, A, I0, A, I0, A, A, A)') &
   &ADJUSTL(TRIM(text)), &
   &' on ', date_time(1), '/', date_time(2), '/', date_time(3), &
   &' at ', date_time(5), ':', date_time(6), ':', date_time(7), &
   &' (time zone = UTC', zone, ')'
    WRITE(6,'(1X)')

! ------------------------------
! Messages for error termination
! ------------------------------

  ELSEIF(message_number == -9999) THEN 
    IF(text /= Text_blank) THEN
      WRITE(6,'(1X, A)') TRIM(ADJUSTL(text))
    ELSE
    ENDIF
    WRITE(6,'(1X)')
    WRITE(6,'(1X, 2A, 1X, A, A)') 'd77 was terminated in ', TRIM(type_program), TRIM(name_program)
    WRITE(6,'(1X)')
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '***************************'
    WRITE(6,'(1X, A)') '**                       **'
    WRITE(6,'(1X, A)') '**   Error termination   **'
    WRITE(6,'(1X, A)') '**                       **'
    WRITE(6,'(1X, A)') '***************************'
    WRITE(6,'(1X)')

    STOP
  ELSE
  ENDIF

END SUBROUTINE write_messages

