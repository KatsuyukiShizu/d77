! This subroutine is part of d77 and computes
! transition/permanent electric dipole moments.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

!   CASSCF, RAS-2SF
!     state = -1 : HF determinant
!     state =  0 : ground state      <- CICOEF
!     state =  1 : 1st excited state <- CICOEF
!     state =  n : nth excited state <- CICOEF
!   CIS, TD
!     state =  0 : ground state (HF determinant)
!     state =  1 : 1st excited state <- CICOEF
!     state =  n : nth excited state <- CICOEF

SUBROUTINE dipole_int_pgf
  USE ifmod, ONLY: write_messages, calc_opr_dmat, calc_dmat_cgf
  USE global_constants
  USE global_read_input
  USE global_read_data
  USE func_int_pgf, ONLY: func_int_pgf_r_pgf
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'dipole_int_pgf'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! ---------------
! Local variables
! ---------------

  INTEGER :: left_state, right_state, n_ci_ci
  DOUBLE PRECISION, ALLOCATABLE :: opr_dmat(:,:)
  DOUBLE PRECISION :: trace_opr_dmat
  DOUBLE PRECISION, ALLOCATABLE :: dmat_cgf(:,:)

  INTEGER i_atm, i_cgf, j_cgf, i_pgf, j_pgf, i_so, j_so

  DOUBLE PRECISION :: xyz_c(1:3)
  DOUBLE PRECISION, ALLOCATABLE :: dipole_elec_cgf(:,:,:)
  DOUBLE PRECISION :: dipole_elec(1:3), dipole_nuc(1:3), dipole(1:3), dipole_norm

  INTEGER :: fid
  INTEGER :: n_repeat

  CALL write_messages(2, Text_blank, type_program, name_program)

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '--------------------------'
  WRITE(6,'(1X, A)') 'Calculating density matrix'
  WRITE(6,'(1X, A)') '--------------------------'

! -------------------------------------------------------
! Calculating one particle reduced density matrix (gamma)
! -------------------------------------------------------

  left_state = 0; right_state = 0
  IF(Bra <= Ket) THEN
    left_state  = Bra
    right_state = Ket
  ELSE
    left_state  = Ket
    right_state = Bra
  ENDIF

  ALLOCATE(dmat_cgf(1:N_cgf, 1:N_cgf)); dmat_cgf = 0.0D0
  IF(Gen_dmat_cgf == 'Read') THEN
    fid = Fdmat_cgf
    WRITE(6,'(1X)')
    WRITE(6,'(1X, 2A)') 'Reading dmat_cgf from ', Fname(fid)
    OPEN(fid, FILE = Fname(fid))
      DO i_cgf = 1, N_cgf
        DO j_cgf = 1, N_cgf
          READ(fid,*) dmat_cgf(i_cgf, j_cgf)
        ENDDO ! j_cgf
      ENDDO ! i_cgf
    CLOSE(fid)
    WRITE(6,'(1X, A)') 'Done'
    WRITE(6,'(1X)')
  ELSEIF(Gen_dmat_cgf == 'Calc') THEN
    ALLOCATE(opr_dmat(1:N_so, 1:N_so)); opr_dmat = 0.0D0
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'gamma: one-particle reduced density matrix between'
    WRITE(6,'(1X, A)') 'electronic states |i> and |j> in terms of spin orbitals'
    IF(Gen_opr_dmat == 'Read') THEN
      fid = Fopr_dmat
      WRITE(6,'(1X, 2A)') 'Reading opr_dmat from ', Fname(fid)
      OPEN(fid, FILE = Fname(fid))
        DO i_so = 1, N_so
          DO j_so = 1, N_so
            READ(fid,*) opr_dmat(i_so, j_so)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)
      WRITE(6,'(1X, A)') 'Done'
    ELSEIF(Gen_opr_dmat == 'Calc') THEN
      n_ci_ci = 0
      trace_opr_dmat = 0.0D0
      CALL calc_opr_dmat&
          &(left_state, right_state, n_ci_ci, opr_dmat, trace_opr_dmat)
      IF(Save_opr_dmat == 'Yes') THEN
        fid = Fopr_dmat
        WRITE(6,'(1X, 2A)') 'Saving opr_dmat to ', Fname(fid)
        CALL write_mat(fid, opr_dmat)
      ELSE
      ENDIF
      WRITE(6,'(1X, A)') 'Done'
      WRITE(6,'(1X)')
      WRITE(6,'(1X, A)') 'Trace of gamma'
      WRITE(6,'(1X, A)') REPEAT('=', 29)
      WRITE(6,'(1X, A3, 1X, A3, 3X, A7, 3X, A)') 'i', 'j', 'n_ci_ci', 'Tr[gamma]'
      WRITE(6,'(1X, A)') REPEAT('-', 29)
      WRITE(6,'(1X, I3, 1X, I3, 3X, I5, 2X, F12.5)') Bra, Ket, n_ci_ci, trace_opr_dmat
      WRITE(6,'(1X, A)') REPEAT('=', 29)
    ENDIF

    CALL calc_dmat_cgf(opr_dmat, dmat_cgf)
    DEALLOCATE(opr_dmat)
    IF(Save_dmat_cgf == 'Yes') THEN
      fid = Fdmat_cgf
      WRITE(6,'(1X, 2A)') 'Saving dmat_cgf to ', Fname(fid)
      OPEN(fid, FILE=Fname(fid))
        DO i_cgf = 1, N_cgf
          DO j_cgf = 1, N_cgf
            WRITE(fid,*) dmat_cgf(i_cgf, j_cgf)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)
      WRITE(6,'(1X, A)') 'Done'
    ELSEIF(Save_dmat_cgf == 'No') THEN
    ELSE
    ENDIF
  ENDIF


  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') '-------------'
  WRITE(6,'(1X, A)') 'Dipole moment'
  WRITE(6,'(1X, A)') '-------------'
  
! -----------------------------------------
! Calculating dipole_cgf from int_pgf_r_pgf
! -----------------------------------------

  xyz_c = 0.0D0
  dipole_elec = 0.0D0; dipole_nuc = 0.0D0; dipole = 0.0D0
  dipole_norm = 0.0D0
  
  ALLOCATE(dipole_elec_cgf(1:N_cgf, 1:N_cgf, 1:3)); dipole_elec_cgf = 0.0D0
  DO i_cgf = 1, N_cgf
    DO j_cgf = 1, N_cgf
      DO i_pgf = 1, N_pgf(i_cgf)
        DO j_pgf = 1, N_pgf(j_cgf)
          dipole_elec_cgf(i_cgf, j_cgf, :) &
         &= dipole_elec_cgf(i_cgf, j_cgf, :) &
         &- cntcoef(i_pgf, i_cgf) * cntcoef(j_pgf, j_cgf) &
         &* cntnormfac(i_pgf, i_cgf) * cntnormfac(j_pgf, j_cgf) &
         &* func_int_pgf_r_pgf &
           &(xyzao(:, i_cgf), xyzao(:, j_cgf), &
           &cntexp(i_pgf, i_cgf), cntexp(j_pgf, j_cgf), &
           &Lmn(:, i_cgf), Lmn(:, j_cgf), &
           &xyz_c)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO i_cgf = 1, N_cgf
    DO j_cgf = 1, N_cgf
      dipole_elec(1:3) = dipole_elec(1:3)&
                      &+ dmat_cgf(i_cgf, j_cgf) * dipole_elec_cgf(i_cgf, j_cgf, 1:3)
    ENDDO
  ENDDO
  DEALLOCATE (dipole_elec_cgf)

  IF(Bra == Ket) THEN
    DO i_atm = 1, N_atm
      dipole_nuc(1:3) = dipole_nuc(1:3) &
                     &+ Nuccharg(i_atm) * Xyznuc(1:3, i_atm)
    ENDDO
  ELSE
  ENDIF
  
  dipole(1:3) = dipole_elec(1:3) + dipole_nuc(1:3)
  dipole_norm = DSQRT(DOT_PRODUCT(dipole, dipole))

  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') 'Dipole moments in atomic units'
  n_repeat = 52
  WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
  WRITE(6,'(1X, A3, 1X, A3, 7X, A)') 'i', 'j', 'x          y          z        Total'
  WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
  WRITE(6,'(1X, I3, 1X, I3, 4F11.5)') bra, ket, dipole(1:3), dipole_norm
  WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') 'Dipole moments in debye'
  WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
  WRITE(6,'(1X, A3, 1X, A3, 7X, A)') 'i', 'j', 'x          y          z        Total'
  WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
  WRITE(6,'(1X, I3, 1X, I3, 4F11.5)') bra, ket, dipole(1:3)*Debye, dipole_norm*Debye
  WRITE(6,'(1X, A)') REPEAT('=', n_repeat)

  fid = Fdipole_au
  OPEN(fid, FILE = Fname(fid))
    WRITE(fid,*) dipole(1:3)
  CLOSE(fid)

  CALL write_messages(3, Text_blank, type_program, name_program)
END SUBROUTINE dipole_int_pgf



