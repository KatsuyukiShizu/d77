! This subroutine is part of d77 and computes
! diagonal/off-diagonal vibronic coupling constants.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE vc_int_pgf
  USE ifmod, ONLY: write_messages, &
                  &calc_opr_dmat, &
                  &calc_dmat_cgf, &
                  &write_mat
  USE global_constants
  USE global_read_input
  USE global_read_data
  USE func_int_pgf, ONLY: func_int_pgf_r_pgf, func_int_pgf_elfld_pgf
  IMPLICIT NONE

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'vc_int_pgf'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Local variables
! One-particle reduced density matrix
  INTEGER :: left_state, right_state, n_ci_ci
  DOUBLE PRECISION, ALLOCATABLE :: opr_dmat(:,:)
  DOUBLE PRECISION :: trace_opr_dmat

! Density matrix
  DOUBLE PRECISION, ALLOCATABLE :: dmat_cgf(:,:)

! Dipole moment
  DOUBLE PRECISION :: xyz_c(1:3)
  DOUBLE PRECISION, ALLOCATABLE :: dipole_elec_cgf(:,:,:)
  DOUBLE PRECISION :: dipole_elec(1:3), dipole_nuc(1:3), dipole(1:3), dipole_norm

! VCC
  DOUBLE PRECISION :: vcc_en(1:N_mode_calc), &
                     &vcc_nn(1:N_mode_calc), &
                     &vcc(1:N_mode_calc), &
                     &vcc_en_atm(1:N_atm, 1:N_mode_calc), &
                     &vcc_nn_atm(1:N_atm, 1:N_mode_calc), &
                     &vcc_atm(1:N_atm, 1:N_mode_calc)

! DO loop variables
  INTEGER :: i_mode, i_atm, j_atm, i_xyz, i_cgf, j_cgf, i_pgf, j_pgf, i_so, j_so

! Electric field integrals
  DOUBLE PRECISION :: int_elfld_gf(1:3)
  DOUBLE PRECISION, ALLOCATABLE :: elfld_cgf_atm(:,:,:,:)

! VCCnn: nuclear repulsion potential
  DOUBLE PRECISION :: dxyz(1:3), norm_dxyz

  INTEGER :: fid
  INTEGER :: n_repeat

  CALL write_messages(2, Text_blank, type_program, name_program)


  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '--------------'
  WRITE(6,'(1X, A)') 'Density matrix'
  WRITE(6,'(1X, A)') '--------------'

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
!   dmat_cgf is calculated from opr_dmat.
    ALLOCATE(opr_dmat(1:N_so, 1:N_so)); opr_dmat = 0.0D0
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'gamma: one-particle reduced density matrix'
    WRITE(6,'(1X, A)') 'between electronic states | i > and | j >'
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
      WRITE(6,'(1X)')
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
                      &+ dmat_cgf(i_cgf, j_cgf)&
                      &* dipole_elec_cgf(i_cgf, j_cgf, 1:3)
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
  n_repeat = 52
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


  fid = Fdipole_au
  OPEN(fid, FILE = Fname(fid))
    WRITE(fid,*) dipole(1:3)
  CLOSE(fid)

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') '-----------------'
  WRITE(6,'(1X, A)') 'Vibronic coupling'
  WRITE(6,'(1X, A)') '-----------------'
  WRITE(6,'(1X)')

  ALLOCATE(elfld_cgf_atm(1:N_cgf, 1:N_cgf, 1:3, 1:N_atm)); elfld_cgf_atm = 0.0D0
  IF(Gen_elfld_cgf_atm == 'Read') THEN
    fid = Felfld_cgf_atm
    WRITE(6,'(1X, 2A)') 'Reading elfld_cgf_atm from ', Fname(fid)
    OPEN(fid, FILE = Fname(fid))
      DO i_atm = 1, n_atm
        DO i_xyz = 1, 3 

          DO i_cgf = 1, N_cgf
            DO j_cgf = i_cgf, n_cgf
              READ(fid,*) elfld_cgf_atm(i_cgf, j_cgf, i_xyz, i_atm)
            ENDDO ! j_cgf
          ENDDO ! i_cgf

        ENDDO ! i_xyz
      ENDDO ! i_atm
    CLOSE(fid)
    DO i_cgf = 1, n_cgf-1
      DO j_cgf = i_cgf+1, n_cgf
        elfld_cgf_atm(j_cgf, i_cgf, :, :) &
     &= elfld_cgf_atm(i_cgf, j_cgf, :, :)
      ENDDO
    ENDDO
    WRITE(6,'(1X, A)') 'Done'

  ELSEIF(Gen_elfld_cgf_atm == 'Calc') THEN

    WRITE(6,'(1X, A)') 'Calculating electric field integrals between CGFs'
    int_elfld_gf(1:3) = 0.0D0
    DO i_atm = 1, n_atm
      DO i_cgf = 1, n_cgf
        DO j_cgf = i_cgf, n_cgf
          DO i_pgf = 1, n_pgf(i_cgf)
            DO j_pgf = 1, n_pgf(j_cgf)
              int_elfld_gf(1:3) &
             &= func_int_pgf_elfld_pgf &
               &(Xyzao(:, i_cgf), Xyzao(:, j_cgf), &
                &Cntexp(i_pgf, i_cgf), Cntexp(j_pgf, j_cgf), &
                &Lmn(:, i_cgf), Lmn(:, j_cgf), &
                &Xyznuc(:, i_atm))
              elfld_cgf_atm(i_cgf, j_cgf, 1:3, i_atm)&
              = elfld_cgf_atm(i_cgf, j_cgf, 1:3, i_atm)&
              + Cntcoef(i_pgf, i_cgf) * Cntcoef(j_pgf, j_cgf)&
              * Cntnormfac(i_pgf, i_cgf) * Cntnormfac(j_pgf, j_cgf)&
              * int_elfld_gf(1:3)
            ENDDO ! j_pgf
          ENDDO ! i_pgf
        ENDDO ! j_cgf
      ENDDO ! i_cgf
    ENDDO ! i_atm

    DO i_cgf = 1, n_cgf-1
      DO j_cgf = i_cgf+1, n_cgf
        elfld_cgf_atm(j_cgf, i_cgf, :, :) &
        &= elfld_cgf_atm(i_cgf, j_cgf, :, :)
      ENDDO
    ENDDO

  ELSE
  ENDIF ! Gen_elfld_cgf_atm == 

  vcc_en = 0.0D0; vcc_nn = 0.0D0; vcc = 0.0D0
  vcc_en_atm = 0.0D0; vcc_nn_atm = 0.0D0; vcc_atm = 0.0D0
  DO i_atm = 1, n_atm

    DO i_mode = 1, n_mode_calc
      DO i_xyz = 1, 3
        DO i_cgf = 1, N_cgf
          DO j_cgf = 1, N_cgf
            vcc_en_atm(i_atm, i_mode)&
           &= vcc_en_atm(i_atm, i_mode)&
           &- dmat_cgf(i_cgf, j_cgf)&
           &* nuccharg(i_atm)&
           &* elfld_cgf_atm(i_cgf, j_cgf, i_xyz, i_atm)&
           &* vibmode_calc(i_xyz, i_atm, i_mode)&
           &/ sqrt_atmwt(i_atm)
          ENDDO
        ENDDO
      ENDDO
    ENDDO ! i_mode = 1, n_mode_calc

    vcc_en(:) = vcc_en(:) + vcc_en_atm(i_atm, :)
  ENDDO
  DEALLOCATE(dmat_cgf)

  IF(Bra == Ket) THEN
    DO i_atm = 1, N_atm
      DO j_atm = 1, N_atm
        dxyz(:) = Xyznuc(:, j_atm) - Xyznuc(:,i_atm)
        norm_dxyz = DSQRT(DOT_PRODUCT(dxyz, dxyz))
        IF(j_atm /= i_atm) THEN
          !vcc_nn_atm = 0.0D0
          DO i_xyz = 1, 3
            vcc_nn_atm(i_atm, :)&
            = vcc_nn_atm(i_atm, :)&
            + Nuccharg(i_atm) * Nuccharg(j_atm)&
            * (Xyznuc(i_xyz, j_atm) - Xyznuc(i_xyz, i_atm))&
            * Vibmode_calc(i_xyz, i_atm, :)&
            / Sqrt_atmwt(i_atm)&
            / (norm_dxyz**3)
          ENDDO
        ENDIF
      ENDDO
      vcc_nn(:) = vcc_nn(:) + vcc_nn_atm(i_atm, :)
    ENDDO
  ELSE
  ENDIF


  vcc_atm = vcc_en_atm + vcc_nn_atm
  vcc = vcc_en + vcc_nn
  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') '--------------------------------'
  WRITE(6,'(1X, A)') 'Vibronic coupling constant (VCC)'
  WRITE(6,'(1X, A)') '--------------------------------'
  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') 'Freq   : frequency (cm-1)'
  WRITE(6,'(1X, A)') 'VCC_en : electronic part of VCC (a.u.)'
  WRITE(6,'(1X, A)') 'VCC_nn : nucelar part of VCC (a.u.)'
  WRITE(6,'(1X, A)') 'VCC = VCC_en + VCC_nn'
  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') '===================================================='
  WRITE(6,'(A5, A8, A12, A13, A10)') 'Mode', 'Freq', 'VCC_en', 'VCC_nn', 'VCC'
  WRITE(6,'(1X, A)') '----------------------------------------------------'
  DO i_mode = 1, n_mode_calc
    WRITE(6,'(I4, F10.2, 3F13.8)') mode_calc(i_mode),&
                                   freq_calc(i_mode),&
                                   vcc_en(i_mode),&
                                   vcc_nn(i_mode),&
                                   vcc(i_mode)
  ENDDO ! i_mode = 1, n_mode_calc
  WRITE(6,'(1X, A)') '===================================================='

  fid = Fvcc
  OPEN(fid, FILE = Fname(fid))
    DO i_mode = 1, n_mode_calc
    WRITE(fid,'(I5, 4E14.5)') mode_calc(i_mode),&
                              freq_calc(i_mode),&
                              vcc_en(i_mode),&
                              vcc_nn(i_mode),&
                              vcc(i_mode)
    ENDDO ! i_mode = 1, n_mode_calc
  CLOSE(fid)

! -------------------
! Writing atomic VCCs
! -------------------

  IF(Save_avcc == 'Yes') THEN
    fid = Favcc
    OPEN(fid, FILE = Fname(fid))
      DO i_mode = 1, n_mode_calc
        WRITE(fid,'(1X, A, I5)') 'Mode', i_mode
        WRITE(fid,'(1X, A, F8.2, A)') 'Wavenumber ', Freq(i_mode), ' cm-1'
        WRITE(fid,'(1X, A)') '==============================================='
        WRITE(fid,'(1X, A5, A10, A14, A14)') 'Atom', 'VCC_en', 'VCC_nn', 'VCC'
        WRITE(fid,'(1X, A)') '-----------------------------------------------'
        DO i_atm = 1, n_atm
          WRITE(fid,'(1X, I5, 3E14.5)') i_atm, &
                                       &vcc_en_atm(i_atm, i_mode), &
                                       &vcc_nn_atm(i_atm, i_mode), &
                                       &vcc_atm(i_atm, i_mode)
        ENDDO ! i_atm = 1, n_atm
        WRITE(fid,'(1X, A)') '-----------------------------------------------'
        WRITE(fid,'(1X, A5, 3E14.5)') 'Total', &
                                     &vcc_en(i_mode),&
                                     &vcc_nn(i_mode),&
                                     &vcc(i_mode)
        WRITE(fid,'(1X, A)') '==============================================='
        WRITE(6,'(1X)')
      ENDDO ! i_mode = 1, n_mode_calc
    CLOSE(fid)
  ELSE
  ENDIF


  DEALLOCATE (atmnum)
  DEALLOCATE (nuccharg)
  DEALLOCATE (xyznuc)
  DEALLOCATE (elmsym)
  DEALLOCATE (ao2atm)
  DEALLOCATE (n_pgf)
  DEALLOCATE (orbtyp)
  DEALLOCATE (lmn)
  DEALLOCATE (xyzao)
  DEALLOCATE (cntexp)
  DEALLOCATE (cntcoef)
  DEALLOCATE (cntnormfac)
  DEALLOCATE (coef_so)
  DEALLOCATE (ex_energy)
  DEALLOCATE (freq)
  DEALLOCATE (vibmode)
  DEALLOCATE (atmwt, sqrt_atmwt)
  DEALLOCATE (freq_calc)
  DEALLOCATE (vibmode_calc)

 CALL write_messages(3, Text_blank, type_program, name_program)
END SUBROUTINE vc_int_pgf

