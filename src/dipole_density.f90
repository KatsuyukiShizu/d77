! This subroutine is part of d77 and computes
! transition dipole density.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE dipole_density
  USE ifmod, ONLY: write_messages, &
                  &calc_opr_dmat, &
                  &calc_dmat_cgf, &
                  &cube_dipole_density, &
                  &write_mat
  USE global_constants
  USE global_read_input
  USE global_read_data
  IMPLICIT NONE

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'dipole_density'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! One-particle reduced density matrix
  INTEGER                       :: left_state, right_state, n_ci_ci
  DOUBLE PRECISION, ALLOCATABLE :: opr_dmat(:,:)
  DOUBLE PRECISION              :: trace_opr_dmat

! Density matrix
  DOUBLE PRECISION, ALLOCATABLE :: dmat_cgf(:,:)
  DOUBLE PRECISION              :: int_rho

! Dipole moment 
  DOUBLE PRECISION              :: dipole(1:4), dipole_norm

! DO Loop variables
  INTEGER :: i_cgf, j_cgf, i_so, j_so

  INTEGER :: fid, n_repeat

  left_state = 0; right_state = 0; n_ci_ci = 0
  trace_opr_dmat = 0.0D0

  int_rho = 0.0D0
  dipole(1:3) = 0.0D0; dipole_norm = 0.0D0

  CALL write_messages(2, Text_blank, type_program, name_program)
    
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '--------------'
  WRITE(6,'(1X, A)') 'Density matrix'
  WRITE(6,'(1X, A)') '--------------'

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

  ALLOCATE(dmat_cgf(1:N_cgf, 1:N_cgf)); dmat_cgf = 0.0D0
  IF(Gen_dmat_cgf == 'Read') THEN
    fid = Fdmat_cgf
    WRITE(6,'(1X)')
    WRITE(6,'(1X, 2A)') 'Reading dmat_cgf from INPUTS/', Fname(fid)
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
      WRITE(6,'(1X, 2A)') 'Reading opr_dmat from INPUTS/', Fname(fid)
      OPEN(fid, FILE = Fname(fid))
        DO i_so = 1, N_so
          DO j_so = 1, N_so
            READ(fid,*) opr_dmat(i_so, j_so)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      CLOSE(fid)
      WRITE(6,'(1X, A)') 'Done'
    ELSEIF(Gen_opr_dmat == 'Calc') THEN
      CALL calc_opr_dmat&
          &(left_state, right_state, n_ci_ci, opr_dmat, trace_opr_dmat)
      IF(Save_opr_dmat == 'Yes') THEN
        fid = Fopr_dmat
        WRITE(6,'(1X, 2A)') 'Saving opr_dmat to ', Fname(fid)
        CALL write_mat(fid, opr_dmat)
      ELSE
      ENDIF
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
  WRITE(6,'(1X, A)') '---------------------'
  WRITE(6,'(1X, A)') 'Dipole moment density'
  WRITE(6,'(1X, A)') '---------------------'
  CALL cube_dipole_density(dmat_cgf, int_rho, dipole)
  dipole_norm = DSQRT(DOT_PRODUCT(dipole, dipole))

  WRITE(6,'(1X)')
  IF(Bra == Ket) THEN
    WRITE(6,'(1X, A)') '----------------------------------------'
    WRITE(6,'(1X, A, 3F10.5)') 'Integral of electron density: ', int_rho
    WRITE(6,'(1X, A)') '----------------------------------------'
    WRITE(6,'(1X, A, I0, A)') 'This value should be approximately equal to ', Ne, '.'
  ELSE
    WRITE(6,'(1X, A)') '-----------------------------------------------------'
    WRITE(6,'(1X, A, 3F10.5)') 'Integral of overlap (transition) density:  ', int_rho
    WRITE(6,'(1X, A)') '-----------------------------------------------------'
    WRITE(6,'(1X, A)') 'This value should be approximately equal to zero.'
  ENDIF

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '-------------'
  WRITE(6,'(1X, A)') 'Dipole moment'
  WRITE(6,'(1X, A)') '-------------'
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Calculating x, y, and z components of dipole moment'
  WRITE(6,'(1X, A)') 'by integrating dipole moment densities.'
  n_repeat = 52
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Dipole moment in atomic units'
  WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
  WRITE(6,'(1X, A3, 1X, A3, 6X, A)') 'i', 'j', 'mu_x       mu_y       mu_z       |mu|'
  WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
  WRITE(6,'(1X, I3, 1X, I3, 4F11.5)') bra, ket, dipole(1:3), dipole_norm
  WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Dipole moment in debye'
  WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
  WRITE(6,'(1X, A3, 1X, A3, 6X, A)') 'i', 'j', 'mu_x       mu_y       mu_z       |mu|'
  WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
  WRITE(6,'(1X, I3, 1X, I3, 4F11.5)') bra, ket, dipole(1:3)*Debye, dipole_norm*Debye
  WRITE(6,'(1X, A)') REPEAT('=', n_repeat)


  DEALLOCATE (Atmnum, Nuccharg, Xyznuc, Elmsym)
  DEALLOCATE (Lmn, Xyzao)
  DEALLOCATE (Coef_so)

  CALL write_messages(3, Text_blank, type_program, name_program)

END SUBROUTINE dipole_density 
