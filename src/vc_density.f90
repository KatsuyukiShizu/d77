! This subroutine is part of d77 and computes
! vibronic coupling density.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE vc_density
  USE ifmod, ONLY: write_messages, &
                  &calc_opr_dmat, &
                  &calc_dmat_cgf, &
                  &write_mat, &
                  &cube_vc_density
  USE global_constants
  USE global_read_input
  USE global_read_data
  IMPLICIT NONE

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'vc_density'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! One-particle reduced density matrix
  INTEGER                       :: left_state, right_state, n_ci_ci
  DOUBLE PRECISION, ALLOCATABLE :: opr_dmat(:,:)
  DOUBLE PRECISION              :: trace_opr_dmat

! Density matrix
  DOUBLE PRECISION, ALLOCATABLE :: dmat_cgf(:,:)

! Vibronic coupling density
  DOUBLE PRECISION :: int_rho
  DOUBLE PRECISION :: vcc_en(1:N_mode_calc), &
                     &vcc_nn(1:N_mode_calc), &
                     &vcc(1:N_mode_calc)

! VCCnn: nuclear repulsion potential
  DOUBLE PRECISION :: dxyz(1:3), norm_dxyz

! DO Loop variables
  INTEGER :: i_atm, j_atm, i_xyz, i_mode_calc, i_cgf, j_cgf, i_so, j_so

  INTEGER :: fid

  left_state = 0; right_state = 0; n_ci_ci = 0
  trace_opr_dmat = 0.0D0

  int_rho = 0.0D0
  vcc_en = 0.0D0; vcc_nn = 0.0D0; vcc = 0.0D0

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
    WRITE(6,'(1X, 2A)') 'Reading dmat_cgf ', Fname(fid)
    OPEN(fid, FILE = Fname(fid))
      DO i_cgf = 1, N_cgf
        DO j_cgf = 1, N_cgf
          READ(fid,*) dmat_cgf(i_cgf, j_cgf)
        ENDDO ! j_cgf
      ENDDO ! i_cgf
    CLOSE(fid)
    WRITE(6,'(1X, A)') 'Done'
    WRITE(6,'(1X)')
  ELSEIF(Gen_dmat_cgf == 'Calc') THEN ! Default
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
  WRITE(6,'(1X, A)') '-------------------------'
  WRITE(6,'(1X, A)') 'Vibronic coupling density'
  WRITE(6,'(1X, A)') '-------------------------'
  CALL cube_vc_density(dmat_cgf, int_rho, vcc_en)

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

  IF(Bra == Ket) THEN
    DO i_atm = 1, N_atm
      DO j_atm = 1, N_atm
        dxyz(:) = Xyznuc(:, j_atm) - Xyznuc(:,i_atm)
        norm_dxyz = DSQRT(DOT_PRODUCT(dxyz, dxyz))
        IF(j_atm /= i_atm) THEN
          DO i_xyz = 1, 3
            vcc_nn(:)&
            = vcc_nn(:)&
            + Nuccharg(i_atm) * Nuccharg(j_atm)&
            * (Xyznuc(i_xyz, j_atm) - Xyznuc(i_xyz, i_atm))&
            * Vibmode_calc(i_xyz, i_atm, :)&
            / Sqrt_atmwt(i_atm)&
            / (norm_dxyz**3)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  ELSE
  ENDIF

  vcc = vcc_en + vcc_nn
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '--------------------------------'
  WRITE(6,'(1X, A)') 'Vibronic coupling constant (VCC)'
  WRITE(6,'(1X, A)') '--------------------------------'
  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') 'Calculating VCCs by integrating vibronic coupling densities.'
  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Freq: frequency in cm-1'
  WRITE(6,'(1X, A)') 'VCC_en: electronic part of VCC in a.u.'
  WRITE(6,'(1X, A)') 'VCC_nn: nucelar part of VCC in a.u.'
  WRITE(6,'(1X, A)') 'VCC = VCC_en + VCC_nn (a.u.)'
  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') '==============================================='
  WRITE(6,'(A5, A10, 3A11)') 'Mode', 'Freq', 'VCC_en', 'VCC_nn', 'VCC'
  WRITE(6,'(1X, A)') '-----------------------------------------------'
  DO i_mode_calc = 1, n_mode_calc
    WRITE(6,'(I5, F10.2, 3F11.6)') mode_calc(i_mode_calc),&
                                   freq_calc(i_mode_calc),&
                                   vcc_en(i_mode_calc),&
                                   vcc_nn(i_mode_calc),&
                                   vcc(i_mode_calc)
  ENDDO ! i_mode_calc = 1, n_mode_calc
  WRITE(6,'(1X, A)') '==============================================='

  DEALLOCATE (Atmnum, Nuccharg, Xyznuc, Elmsym)
  DEALLOCATE (Lmn, Xyzao)
  DEALLOCATE (Coef_so)
  DEALLOCATE (Freq, Vibmode, Atmwt, Sqrt_atmwt)
  DEALLOCATE (Freq_calc, Vibmode_calc)

  CALL write_messages(3, Text_blank, type_program, name_program)
END SUBROUTINE vc_density 
