! This subroutine is part of d77 and computes
! spin-orbit couplings between singlet and triplet states.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE soc_int_pgf
  USE ifmod, ONLY: write_messages, &
                  &calc_opr_dmat, &
                  &calc_dmat_soc_cgf, &
                  &write_mat
  USE global_constants
  USE global_read_input
  USE global_read_data
  USE func_int_pgf, ONLY: func_int_pgf_r_pgf, func_int_pgf_soc_pgf
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'soc_int_pgf'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! One-particle reduced density matrix
  INTEGER :: left_state, right_state, n_ci_ci
  DOUBLE PRECISION, ALLOCATABLE :: opr_dmat(:,:)
  DOUBLE PRECISION :: trace_opr_dmat
  INTEGER :: i_so, j_so

  INTEGER :: i_cgf, j_cgf, mu_cgf, nu_cgf 
  !INTEGER :: i_cgf, j_cgf, mu_cgf, nu_cgf, i_pgf, j_pgf

! Spin-orbit coupling
  DOUBLE PRECISION, ALLOCATABLE :: dmat_soc_cgf(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: soc_cgf(:,:,:)
  !DOUBLE PRECISION :: int_soc_pgf_atm(1:3)
  DOUBLE PRECISION :: int_soc_pgf(1:3)
  DOUBLE PRECISION :: soc_atm(1:3, 1:N_atm), soc(1:3)

! Device number
  INTEGER :: fid


  CALL write_messages(2, Text_blank, type_program, name_program)

  left_state = 0; right_state = 0; n_ci_ci = 0
  trace_opr_dmat = 0.0D0

! Spin-orbit coupling
  soc_atm = 0.0D0; soc = 0.0D0
  int_soc_pgf = 0.0D0

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '--------------------------'
  WRITE(6,'(1X, A)') 'Calculating density matrix'
  WRITE(6,'(1X, A)') '--------------------------'

! -------------------------------------------------------
! Calculating one particle reduced density matrix (gamma)
! -------------------------------------------------------

  IF(Bra <= Ket) THEN
    left_state  = Bra
    right_state = Ket
  ELSE
    left_state  = Ket
    right_state = Bra
  ENDIF

  ALLOCATE(dmat_soc_cgf(1:N_cgf, 1:N_cgf, 1:3)); dmat_soc_cgf = 0.0D0
  IF(Gen_dmat_soc_cgf == 'Read') THEN
    IF(Ms == 0) THEN
      fid = Fdmat_soc_cgf_z
      WRITE(6,'(1X)')
      WRITE(6,'(1X, A)') 'Reading dmat_soc_cgf from ', Fname(fid)
      OPEN(fid, FILE = Fname(fid))
        DO i_cgf = 1, N_cgf
          DO j_cgf = 1, N_cgf
            READ(fid,*) dmat_soc_cgf(i_cgf, j_cgf, 3)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)
      WRITE(6,'(1X, A)') 'Done'
      WRITE(6,'(1X)')
    ELSEIF(Ms == 1 .OR. Ms == -1) THEN
      fid = Fdmat_soc_cgf_x
      WRITE(6,'(1X, A)') 'Reading dmat_soc_cgf from ', Fname(fid)
      OPEN(fid, FILE=Fname(fid))
        DO i_cgf = 1, N_cgf
          DO j_cgf = 1, N_cgf
            READ(fid,*) dmat_soc_cgf(i_cgf, j_cgf, 1)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)
      fid = Fdmat_soc_cgf_y
      WRITE(6,'(1X, A)') 'Reading dmat_soc_cgf from ', Fname(fid)
      OPEN(fid, FILE=Fname(fid))
        DO i_cgf = 1, N_cgf
          DO j_cgf = 1, N_cgf
            READ(fid,*) dmat_soc_cgf(i_cgf, j_cgf, 2)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)
      WRITE(6,'(1X, A)') 'Done'
      WRITE(6,'(1X)')
    ELSE
    ENDIF
  ELSEIF(Gen_dmat_soc_cgf == 'Calc') THEN
    ALLOCATE(opr_dmat(1:N_so, 1:N_so)); opr_dmat = 0.0D0
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'gamma: one-particle reduced density matrix between'
    WRITE(6,'(1X, A)') 'electronic states |i> and |j> in terms of spin orbitals'
    WRITE(6,'(1X)')
    IF(Gen_opr_dmat == 'Read') THEN
      fid = Fopr_dmat
      WRITE(6,'(1X, A)') 'Reading opr_dmat from ', Fname(fid)
      OPEN(fid, FILE = Fname(fid))
        DO i_so = 1, N_so
          DO j_so = 1, N_so
            READ(fid,*) opr_dmat(i_so, j_so)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)
      WRITE(6,'(1X, A)') 'Done'
      WRITE(6,'(1X)')
    ELSEIF(Gen_opr_dmat == 'Calc') THEN
      CALL calc_opr_dmat&
     &(left_state, right_state, n_ci_ci, opr_dmat, trace_opr_dmat)
      IF(Save_opr_dmat == 'Yes') THEN
        fid = Fopr_dmat
        WRITE(6,'(1X, A)') 'Saving opr_dmat to ', Fname(fid)
        CALL write_mat(fid, opr_dmat)
      ELSE
      ENDIF
      WRITE(6,'(1X, A)') 'Done'
      WRITE(6,'(1X)')
      WRITE(6,'(1X, A)') 'Trace of gamma'
      WRITE(6,'(1X, A)') REPEAT('=', 29)
      WRITE(6,9006) 'i', 'j', 'n_ci_ci', 'Tr[gamma]'
      WRITE(6,'(1X, A)') REPEAT('-', 29)
      WRITE(6,9007) Bra, Ket, n_ci_ci, trace_opr_dmat
      WRITE(6,'(1X, A)') REPEAT('=', 29)
      WRITE(6,'(1X)')
      9006 FORMAT(1X, A3, 1X, A3, 3X, A7, 3X, A)
      9007 FORMAT(1X, I3, 1X, I3, 3X, I5, 2X, F12.5)
    ENDIF

    CALL calc_dmat_soc_cgf(opr_dmat, dmat_soc_cgf) ! For Ms = 0, +1, and -1
    DEALLOCATE(opr_dmat)
  ENDIF

  ALLOCATE(soc_cgf(1:N_cgf, 1:N_cgf, 1:3)); soc_cgf = 0.0D0
  WRITE(6,'(1X)')
  WRITE(6,'(1X)')
  IF(Gen_soc_cgf == 'Read') THEN
!   soc_cgf can be calculated with Gaussian 16
    IF(Ms == 0) THEN
      fid = Fsoc_cgf_z
      WRITE(6,'(1X, A)') 'Reading z components of soc_cgf from ', Fname(fid)
      OPEN(fid, FILE=Fname(fid))
        DO j_cgf = 1, n_cgf
          DO i_cgf = j_cgf, n_cgf
            READ(fid,*) soc_cgf(i_cgf, j_cgf, 3)
            soc_cgf(i_cgf, j_cgf, 3) = - soc_cgf(i_cgf, j_cgf, 3)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)

      DO j_cgf = 1, n_cgf-1
        DO i_cgf = j_cgf+1, n_cgf
          soc_cgf(j_cgf, i_cgf, 3) &
          &= - soc_cgf(i_cgf, j_cgf, 3)
        ENDDO
      ENDDO
  
      DO mu_cgf = 1, N_cgf
        DO nu_cgf = 1, N_cgf
          soc(3) = soc(3) &
                &+ dmat_soc_cgf(nu_cgf, mu_cgf, 3) &
                &* soc_cgf(mu_cgf, nu_cgf, 3)
        ENDDO
      ENDDO

    ELSEIF(Ms == 1 .OR. Ms == -1) THEN
      fid = Fsoc_cgf_x
      WRITE(6,'(1X, A)') 'Reading x components of soc_cgf from ', Fname(fid)
      OPEN(fid, FILE=Fname(fid))
        DO j_cgf = 1, n_cgf
          DO i_cgf = j_cgf, n_cgf
            READ(fid,*) soc_cgf(i_cgf, j_cgf, 1)
            soc_cgf(i_cgf, j_cgf, 1) = - soc_cgf(i_cgf, j_cgf, 1)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)

      fid = Fsoc_cgf_y
      WRITE(6,'(1X, A)') 'Reading y components of soc_cgf from ', Fname(fid)
      OPEN(fid, FILE=Fname(fid))
        DO j_cgf = 1, n_cgf
          DO i_cgf = j_cgf, n_cgf
            READ(fid,*) soc_cgf(i_cgf, j_cgf, 2)
            soc_cgf(i_cgf, j_cgf, 2) = - soc_cgf(i_cgf, j_cgf, 2)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)

    ELSE
    ENDIF

    DO j_cgf = 1, n_cgf-1
      DO i_cgf = j_cgf+1, n_cgf
        soc_cgf(j_cgf, i_cgf, 1:2) &
        &= - soc_cgf(i_cgf, j_cgf, 1:2)
      ENDDO
    ENDDO

    DO mu_cgf = 1, N_cgf
      DO nu_cgf = 1, N_cgf
        soc(1:2) = soc(1:2) &
              &+ dmat_soc_cgf(nu_cgf, mu_cgf, 1:2) &
              &* soc_cgf(mu_cgf, nu_cgf, 1:2)
      ENDDO
    ENDDO

  ELSE
  ENDIF


  WRITE(6,'(1X)')
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '----------------------------------'
  WRITE(6,'(1X, A)') 'Spin-orbit coupling constant (SOC)'
  WRITE(6,'(1X, A)') '----------------------------------'
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Alpha: fine-structure constant'
  WRITE(6,'(1X, A, F13.6)') 'Alpha   = ', Fsc
  WRITE(6,'(1X, A, F12.8)') '1/Alpha = ', Inv_fsc
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A, I0, 3A, I0, A)') '< ', Bra, ' | ', TRIM(ADJUSTL(property)), ' | ', ket, ' > (a.u.)'
  WRITE(6,'(1X, A)') '========================================'
  WRITE(6,'(1X, A)') '      x             y             z      '
  WRITE(6,'(1X, A)') '----------------------------------------'
  WRITE(6,'(E13.5, 1X, E13.5, 1X, E13.5)') soc(1:3)
  WRITE(6,'(1X, A)') '========================================'
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A, I0, 3A, I0, A)') '< ', Bra, ' | ', TRIM(ADJUSTL(property)), ' | ', ket, ' > * Alpha**2/2 (a.u.)'
  WRITE(6,'(1X, A)') '========================================'
  WRITE(6,'(1X, A)') '      x             y             z      '
  WRITE(6,'(1X, A)') '----------------------------------------'
  WRITE(6,'(E13.5, 1X, E13.5, 1X, E13.5)') soc(1:3)*(fsc**2.0D0)*0.25D0
  WRITE(6,'(1X, A)') '========================================'

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A, I0, 3A, I0, A)') '< ', Bra, ' | ', TRIM(ADJUSTL(property)), ' | ', ket, ' > * Alpha**2/2 (cm-1)'
  WRITE(6,'(1X, A)') '========================================'
  WRITE(6,'(1X, A)') '      x             y             z      '
  WRITE(6,'(1X, A)') '----------------------------------------'
  WRITE(6,'(E13.5, 1X, E13.5, 1X, E13.5)') soc(1:3)*(fsc**2.0D0)*0.25D0*Au_to_wavenumber_soc
  WRITE(6,'(1X, A)') '========================================'


  fid = Fsoc_au
  OPEN(fid, FILE = Fname(fid))
    WRITE(fid,*) soc(1:3)
  CLOSE(fid)



  CALL write_messages(3, Text_blank, type_program, name_program)
END SUBROUTINE soc_int_pgf

