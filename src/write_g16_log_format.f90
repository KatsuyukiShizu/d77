! This subroutine is part of d77 and prints data in Gaussian 16 formats.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE write_g16_log_format(fid, data)
  USE global_constants
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Input variables
  INTEGER, INTENT(IN) :: fid
  DOUBLE PRECISION, INTENT(IN) :: data(:,:)

! Local variables
  INTEGER :: n_elm
  INTEGER :: n_col, n_block, n_col_final
  INTEGER :: i_elm, i_block, counter

  n_block = 0; n_col_final = 0
  counter = 0
  n_col = 5

  n_elm = UBOUND(data, 1)
  n_col_final = MOD(n_elm, n_col)
  OPEN(fid, FILE = Fname(fid))

    IF(n_col_final == 0) THEN
      n_block = n_elm/n_col
      DO i_block = 1, n_block
        counter = (i_block-1) * n_col
        WRITE(fid,9005) counter+1, counter+2, counter+3, counter+4, counter+5
        WRITE(fid,9011) counter+1, data(counter+1, counter+1)
        WRITE(fid,9012) counter+2, data(counter+2, counter+1:counter+2)
        WRITE(fid,9013) counter+3, data(counter+3, counter+1:counter+3)
        WRITE(fid,9014) counter+4, data(counter+4, counter+1:counter+4)
        WRITE(fid,9015) counter+5, data(counter+5, counter+1:counter+5)
        DO i_elm = counter+6, n_elm
          WRITE(fid,9015) i_elm, data(i_elm, counter+1:counter+5)
        ENDDO
      ENDDO
    ELSE
      IF(n_elm < n_col) THEN
        SELECT CASE(n_col_final)
          CASE(1)
            WRITE(fid,9005) 1
            WRITE(fid,9011) 1, data(1, 1)
          CASE(2)
            WRITE(fid,9005) 1, 2
            WRITE(fid,9011) 1, data(1, 1)
            WRITE(fid,9012) 2, data(2, 1:2)
          CASE(3)
            WRITE(fid,9005) 1, 2, 3
            WRITE(fid,9011) 1, data(1, 1)
            WRITE(fid,9012) 2, data(2, 1:2)
            WRITE(fid,9013) 3, data(3, 1:3)
          CASE(4)
            WRITE(fid,9005) 1, 2, 3, 4
            WRITE(fid,9011) 1, data(1, 1)
            WRITE(fid,9012) 2, data(2, 1:2)
            WRITE(fid,9013) 3, data(3, 1:3)
            WRITE(fid,9014) 4, data(4, 1:4)
          CASE DEFAULT
            STOP
        END SELECT
      ELSE
        n_block = (n_elm - MOD(n_elm, n_col))/n_col + 1
        DO i_block = 1, n_block - 1
          counter = (i_block-1) * n_col
          WRITE(fid,9005) counter+1, counter+2, counter+3, counter+4, counter+5
          WRITE(fid,9011) counter+1, data(counter+1, counter+1)
          WRITE(fid,9012) counter+2, data(counter+2, counter+1:counter+2)
          WRITE(fid,9013) counter+3, data(counter+3, counter+1:counter+3)
          WRITE(fid,9014) counter+4, data(counter+4, counter+1:counter+4)
          WRITE(fid,9015) counter+5, data(counter+5, counter+1:counter+5)
          DO i_elm = counter+6, n_elm
            WRITE(fid,9015) i_elm, data(i_elm, counter+1:counter+5)
          ENDDO
        ENDDO
  
        counter = (n_block-1) * n_col
        SELECT CASE(n_col_final)
          CASE(1)
            WRITE(fid,9005) counter+1
            WRITE(fid,9011) counter+1, data(counter+1, counter+1)
          CASE(2)
            WRITE(fid,9005) counter+1, counter+2
            WRITE(fid,9011) counter+1, data(counter+1, counter+1)
            WRITE(fid,9012) counter+2, data(counter+2, counter+1:counter+2)
          CASE(3)
            WRITE(fid,9005) counter+1, counter+2, counter+3
            WRITE(fid,9011) counter+1, data(counter+1, counter+1)
            WRITE(fid,9012) counter+2, data(counter+2, counter+1:counter+2)
            WRITE(fid,9013) counter+3, data(counter+3, counter+1:counter+3)
          CASE(4)
            WRITE(fid,9005) counter+1, counter+2, counter+3, counter+4
            WRITE(fid,9011) counter+1, data(counter+1, counter+1)
            WRITE(fid,9012) counter+2, data(counter+2, counter+1:counter+2)
            WRITE(fid,9013) counter+3, data(counter+3, counter+1:counter+3)
            WRITE(fid,9014) counter+4, data(counter+4, counter+1:counter+4)
          CASE DEFAULT
            STOP
        END SELECT
      ENDIF
    ENDIF

  CLOSE(fid)

  9005 FORMAT(3X, 5I14)
  9011 FORMAT(I7,  E14.6)
  9012 FORMAT(I7, 2E14.6)
  9013 FORMAT(I7, 3E14.6)
  9014 FORMAT(I7, 4E14.6)
  9015 FORMAT(I7, 5E14.6)
END SUBROUTINE write_g16_log_format
