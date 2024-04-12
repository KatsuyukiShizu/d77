! This subroutine is part of d77 and computes
! density matrices between contracted Gaussian functions
! required for calculating spin-orbit couplings.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU Lesser General Public License v3.0
! as published by the Free Software Foundation.
! http://www.gnu.org/licenses
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE calc_dmat_soc_cgf(opr_dmat, dmat_soc_cgf)
  USE ifmod, ONLY: write_messages
  USE global_constants
  USE global_read_input
  USE global_read_data
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'calc_dmat_soc_cgf'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Output variables
  DOUBLE PRECISION, INTENT(IN) :: opr_dmat(:,:)
  DOUBLE PRECISION, INTENT(OUT) :: dmat_soc_cgf(:,:,:)

! Local variables
  INTEGER :: mu_cgf, nu_cgf, i_odd_so, j_odd_so, i_even_so, j_even_so
  INTEGER :: fid

  IF(Debug == 'Yes') &
 &CALL write_messages(2, Text_blank, type_program, name_program)

  dmat_soc_cgf = 0.0D0

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Calculating density matrix in terms of CGFs'
  IF(Ms == 0) THEN
    DO nu_cgf = 1, N_cgf
      DO mu_cgf = 1, N_cgf
    
!       Calculating z component
        DO i_odd_so = 1, N_so, 2
          DO j_odd_so = 1, N_so, 2
            dmat_soc_cgf(nu_cgf, mu_cgf, 3) &
           &= dmat_soc_cgf(nu_cgf, mu_cgf, 3) &
           &+ opr_dmat(j_odd_so, i_odd_so) &
           &* coef_so(mu_cgf, i_odd_so) * coef_so(nu_cgf, j_odd_so)
          ENDDO
        ENDDO
!  
        DO i_even_so = 2, N_so, 2
          DO j_even_so = 2, N_so, 2
            dmat_soc_cgf(nu_cgf, mu_cgf, 3) &
           &= dmat_soc_cgf(nu_cgf, mu_cgf, 3) &
           &- opr_dmat(j_even_so, i_even_so) &
           &* coef_so(mu_cgf, i_even_so) * coef_so(nu_cgf, j_even_so)
          ENDDO
        ENDDO
    
      ENDDO
    ENDDO

    IF(Save_dmat_soc_cgf == 'Yes') THEN
      fid = fdmat_soc_cgf_z
      OPEN(fid, FILE = Fname(fid))
        DO mu_cgf = 1, N_cgf
          DO nu_cgf = mu_cgf, N_cgf
            WRITE(fid,*) mu_cgf, nu_cgf, dmat_soc_cgf(mu_cgf, nu_cgf, 3)
          ENDDO
        ENDDO
      CLOSE(fid)  
    ELSE
    ENDIF

  ELSEIF(Ms == 1 .OR. Ms == -1) THEN
    DO nu_cgf = 1, N_cgf
      DO mu_cgf = 1, N_cgf

!       Calculating x component
        DO i_odd_so = 1, N_so, 2
          DO j_even_so = 2, N_so, 2
            dmat_soc_cgf(nu_cgf, mu_cgf, 1) &
           &= dmat_soc_cgf(nu_cgf, mu_cgf, 1) &
           &+ opr_dmat(j_even_so, i_odd_so) &
           &* coef_so(mu_cgf, i_odd_so) * coef_so(nu_cgf, j_even_so)
          ENDDO
        ENDDO
!     
        DO j_odd_so = 1, N_so, 2
          DO i_even_so = 2, N_so, 2
            dmat_soc_cgf(nu_cgf, mu_cgf, 1) &
           &= dmat_soc_cgf(nu_cgf, mu_cgf, 1) &
           &+ opr_dmat(j_odd_so, i_even_so) &
           &* coef_so(mu_cgf, i_even_so) * coef_so(nu_cgf, j_odd_so)
          ENDDO
        ENDDO

!       Calculating y component
        DO i_odd_so = 1, N_so, 2
          DO j_even_so = 2, N_so, 2
            dmat_soc_cgf(nu_cgf, mu_cgf, 2) &
           &= dmat_soc_cgf(nu_cgf, mu_cgf, 2) &
           &+ opr_dmat(j_even_so, i_odd_so) &
           &* coef_so(mu_cgf, i_odd_so) * coef_so(nu_cgf, j_even_so)
          ENDDO
        ENDDO
!  
        DO j_odd_so = 1, N_so, 2
          DO i_even_so = 2, N_so, 2
            dmat_soc_cgf(nu_cgf, mu_cgf, 2) &
           &= dmat_soc_cgf(nu_cgf, mu_cgf, 2) &
           &- opr_dmat(j_odd_so, i_even_so) &
           &* coef_so(mu_cgf, i_even_so) * coef_so(nu_cgf, j_odd_so)
          ENDDO
        ENDDO

      ENDDO
    ENDDO

    IF(Save_dmat_soc_cgf == 'Yes') THEN
      fid = fdmat_soc_cgf_x
      OPEN(fid, FILE = Fname(fid))
        DO mu_cgf = 1, N_cgf
          DO nu_cgf = mu_cgf, N_cgf
            WRITE(fid,*) mu_cgf, nu_cgf, dmat_soc_cgf(mu_cgf, nu_cgf, 1)
          ENDDO
        ENDDO
      CLOSE(fid)
  
      fid = fdmat_soc_cgf_y
      OPEN(fid, FILE = Fname(fid))
        DO mu_cgf = 1, N_cgf
          DO nu_cgf = mu_cgf, N_cgf
            WRITE(fid,*) mu_cgf, nu_cgf, dmat_soc_cgf(mu_cgf, nu_cgf, 2)
          ENDDO
        ENDDO
      CLOSE(fid)
    ELSE
    ENDIF

  ELSE  
  ENDIF

  WRITE(6,'(1X, A)') 'Done'

  IF(Debug == 'Yes') &
 &CALL write_messages(3, Text_blank, type_program, name_program)
END SUBROUTINE calc_dmat_soc_cgf
