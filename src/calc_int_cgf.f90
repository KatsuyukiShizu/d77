! This subroutine is part of d77 and computes
! one-electron integrals between contracted Gaussian functions.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE calc_int_cgf
  USE ifmod, ONLY: write_messages, write_g16_log_format
  USE global_constants
  USE global_read_input
  USE global_read_data
  USE func_int_pgf
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'calc_int_cgf'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'


! ---------------
! Local variables
! ---------------

! DO loop variables
  INTEGER :: i_atm, i_cgf, j_cgf, i_pgf, j_pgf, i_xyz

! Kinetic energy integral < pgf |  | pgf >
  DOUBLE PRECISION :: int_ke_gf(1:3) ! (1,2,3) = (x,y,z)
  DOUBLE PRECISION, ALLOCATABLE :: ke_cgf(:,:), ke_cgf_xyz(:,:,:)

! Potential energy integral: < pgf | 1/r | pgf >
  DOUBLE PRECISION :: int_pe_gf ! < pgf | 1/r | pgf >
  DOUBLE PRECISION, ALLOCATABLE :: pe_cgf_atm(:,:,:), pe_cgf(:,:)

! Electric field integral: < pgf | x/r^3 | pgf >, < pgf | y/r^3 | pgf >, < pgf | z/r^3 | pgf >
  DOUBLE PRECISION :: int_elfld_gf(1:3)
  DOUBLE PRECISION, ALLOCATABLE :: elfld_cgf_atm(:,:,:,:)

! Spin-orbit coupling
  DOUBLE PRECISION :: int_soc_pgf(1:3)
  DOUBLE PRECISION :: int_soc_pgf_atm(1:3)
  DOUBLE PRECISION, ALLOCATABLE :: soc_cgf_atm(:,:,:,:), soc_cgf(:,:,:)

  INTEGER :: n_repeat

! Device number
  INTEGER :: fid

  CALL write_messages(2, Text_blank, type_program, name_program)

  n_repeat = 0


  IF(Property == 'Ke') THEN
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '--------------'
    WRITE(6,'(1X, A)') 'Kinetic energy'
    WRITE(6,'(1X, A)') '--------------'
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Calculating kinetic energies between CGFs and'
    ALLOCATE(ke_cgf(1:N_cgf, 1:N_cgf)); ke_cgf = 0.0D0
    ALLOCATE(ke_cgf_xyz(1:N_cgf, 1:N_cgf, 1:3)); ke_cgf_xyz = 0.0D0
    int_ke_gf = 0.0D0
    DO i_cgf = 1, N_cgf
      DO j_cgf = i_cgf, N_cgf
        DO i_pgf = 1, N_pgf(i_cgf)
          DO j_pgf = 1, N_pgf(j_cgf)
            int_ke_gf &
           &= func_int_pgf_ke_pgf& ! < pgf | ke | pgf >
             &(Xyzao(:, i_cgf), Xyzao(:, j_cgf), &
              &Cntexp(i_pgf, i_cgf), Cntexp(j_pgf, j_cgf), &
              &Lmn(:, i_cgf), Lmn(:, j_cgf))
            ke_cgf_xyz(i_cgf, j_cgf, 1:3)&
           &= ke_cgf_xyz(i_cgf, j_cgf, 1:3)&
           &+ Cntcoef(i_pgf, i_cgf) * Cntcoef(j_pgf, j_cgf)&
           &* Cntnormfac(i_pgf, i_cgf) * Cntnormfac(j_pgf, j_cgf)&
           &* int_ke_gf
          ENDDO ! j_pgf
        ENDDO ! i_pgf
      ENDDO ! j_cgf
    ENDDO ! i_cgf
  
    DO i_cgf = 1, N_cgf-1
      DO j_cgf = i_cgf+1, N_cgf
        ke_cgf_xyz(j_cgf, i_cgf, 1:3) = ke_cgf_xyz(i_cgf, j_cgf, 1:3)
      ENDDO ! j_cgf
    ENDDO ! i_cgf
    ke_cgf_xyz = 0.5D0 * ke_cgf_xyz
  
    DO i_cgf = 1, N_cgf
      DO j_cgf = 1, N_cgf
        ke_cgf(i_cgf, j_cgf) = SUM(ke_cgf_xyz(i_cgf, j_cgf, 1:3))
      ENDDO ! j_cgf
    ENDDO ! i_cgf
  
    IF(Save_ke_cgf_g16_format == 'Yes') THEN
      fid = Fke_cgf_g16_format
      WRITE(6,'(1X, A)') 'saving them in the IOP(3/33=1) format of Gaussian 16'
      CALL write_g16_log_format(fid, ke_cgf)
      WRITE(6,'(1X, A)') 'Done'
    ELSE
    ENDIF  
  ELSE
  ENDIF

  IF(Property == 'Pe') THEN
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '----------------'
    WRITE(6,'(1X, A)') 'Potential energy'
    WRITE(6,'(1X, A)') '----------------'
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Calculating potential energies between CGFs and'
    ALLOCATE(pe_cgf_atm(1:N_cgf, 1:N_cgf, 1:N_atm)); pe_cgf_atm = 0.0D0
    ALLOCATE(pe_cgf(1:N_cgf, 1:N_cgf)); pe_cgf = 0.0D0
    int_pe_gf = 0.0D0
    DO i_atm = 1, N_atm
      DO i_cgf = 1, N_cgf
        DO j_cgf = i_cgf, N_cgf
          DO i_pgf = 1, N_pgf(i_cgf)
            DO j_pgf = 1, N_pgf(j_cgf)
              int_pe_gf &
             &= func_int_pgf_inv_r_pgf &
               &(Xyzao(:, i_cgf), Xyzao(:, j_cgf), &
                &Cntexp(i_pgf, i_cgf), Cntexp(j_pgf, j_cgf), &
                &Lmn(:, i_cgf), Lmn(:, j_cgf), &
                &Xyznuc(:, i_atm))
              pe_cgf_atm(i_cgf, j_cgf, i_atm)&
             &= pe_cgf_atm(i_cgf, j_cgf, i_atm)&
             &- Cntcoef(i_pgf, i_cgf) * Cntcoef(j_pgf, j_cgf)&
             &* Cntnormfac(i_pgf, i_cgf) * Cntnormfac(j_pgf, j_cgf)&
             &* int_pe_gf
            ENDDO ! j_pgf
          ENDDO ! i_pgf
          pe_cgf_atm(i_cgf, j_cgf, i_atm)&
          &= pe_cgf_atm(i_cgf, j_cgf, i_atm)&
          &* Nuccharg(i_atm)
        ENDDO ! j_cgf
      ENDDO ! i_cgf
    ENDDO ! i_atm
  
    DO i_cgf = 1, N_cgf-1
      DO j_cgf = i_cgf+1, N_cgf
        pe_cgf_atm(j_cgf, i_cgf, 1:N_atm) = pe_cgf_atm(i_cgf, j_cgf, 1:N_atm)
      ENDDO ! j_cgf
    ENDDO ! i_cgf
  
    DO i_cgf = 1, N_cgf
      DO j_cgf = 1, N_cgf
        pe_cgf(i_cgf, j_cgf) = SUM(pe_cgf_atm(i_cgf, j_cgf, :))
      ENDDO
    ENDDO
  
    IF(Save_pe_cgf_g16_format == 'Yes') THEN
      fid = Fpe_cgf_g16_format
      WRITE(6,'(1X, A)') 'saving them in the IOP(3/33=1) format of Gaussian 16'
      CALL write_g16_log_format(Fid, -pe_cgf)
      WRITE(6,'(1X, A)') 'Done'
    ELSE
    ENDIF  

    fid = Fpe_cgf_atm ! Data file for d77
    OPEN(fid, FILE=Fname(fid))
      DO i_atm = 1, N_atm
        DO i_cgf = 1, N_cgf
          DO j_cgf = i_cgf, N_cgf
            WRITE(fid, *) pe_cgf_atm(i_cgf, j_cgf, i_atm)
          ENDDO ! j_cgf
        ENDDO ! i_cgf
      ENDDO ! i_atm
    CLOSE(fid)

  ELSE
  ENDIF

  IF(Property == 'Vc') THEN
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '-----------------------'
    WRITE(6,'(1X, A)') 'Electric field integral'
    WRITE(6,'(1X, A)') '-----------------------'
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Calculating electric field integrals between CGFs and'
    int_elfld_gf(1:3) = 0.0D0
    ALLOCATE(elfld_cgf_atm(1:N_cgf, 1:N_cgf, 1:3, 1:N_atm)); elfld_cgf_atm = 0.0D0
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

    fid = Felfld_cgf_atm
    WRITE(6,'(1X, 2A)') 'saving them to a file ', TRIM(ADJUSTL(Fname(fid)))
    OPEN(fid, FILE=Fname(fid))
      DO i_atm = 1, n_atm
        DO i_xyz = 1, 3
          DO i_cgf = 1, n_cgf
            DO j_cgf = i_cgf, n_cgf
              WRITE(fid, *) elfld_cgf_atm(i_cgf, j_cgf, i_xyz, i_atm)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    CLOSE(fid)
    WRITE(6,'(1X, A)') 'Done'
  ELSE
  ENDIF

  IF(Property == 'Soc') THEN
    WRITE(6,'(1X)')
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '-------------------'
    WRITE(6,'(1X, A)') 'Spin-orbit coupling'
    WRITE(6,'(1X, A)') '-------------------'
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Calculating x, y, and z components of'
    WRITE(6,'(1X, A)') 'spin-orbit couplings between CGFs.'
    WRITE(6,'(1X)') 
    ALLOCATE(soc_cgf_atm(1:N_cgf, 1:N_cgf, 1:N_atm, 1:3)); soc_cgf_atm = 0.0D0
    ALLOCATE(soc_cgf(1:N_cgf, 1:N_cgf, 1:3)); soc_cgf = 0.0D0
    
!   ------------
!   z components
!   ------------
    
    DO i_atm = 1, N_atm
    
      DO i_cgf = 1, N_cgf
        DO j_cgf = i_cgf, N_cgf
    
          DO i_pgf = 1, N_pgf(i_cgf)
            DO j_pgf = 1, N_pgf(j_cgf)
    
              int_soc_pgf(1:3) = 0.0D0
    
              int_soc_pgf_atm(1:3) &
             &= func_int_pgf_soc_pgf &
               &(Xyzao(1:3, i_cgf), Xyzao(1:3, j_cgf), &
                &Cntexp(i_pgf, i_cgf), Cntexp(j_pgf, j_cgf), &
                &Lmn(1:3, i_cgf), Lmn(1:3, j_cgf), &
                &Xyznuc(1:3, i_atm))

              soc_cgf_atm(i_cgf, j_cgf, i_atm, 1:3) &
             &= soc_cgf_atm(i_cgf, j_cgf, i_atm, 1:3) &
             &- Cntcoef(i_pgf, i_cgf) * Cntcoef(j_pgf, j_cgf) &
             &* Cntnormfac(i_pgf, i_cgf) * Cntnormfac(j_pgf, j_cgf) &
             &* int_soc_pgf_atm(1:3)
    
            ENDDO ! j_pgf
          ENDDO ! i_pgf
    
          soc_cgf_atm(i_cgf, j_cgf, i_atm, 1:3) &
         &= soc_cgf_atm(i_cgf, j_cgf, i_atm, 1:3) &
         &* Nuccharg_soc(i_atm)
        ENDDO ! j_cgf
      ENDDO ! i_cgf
    
    ENDDO ! i_atm
    
    
    DO i_cgf = 1, N_cgf-1
      DO j_cgf = i_cgf+1, N_cgf
        soc_cgf_atm(j_cgf, i_cgf, 1:N_atm, 1:3) &
        &= - soc_cgf_atm(i_cgf, j_cgf, 1:N_atm, 1:3)
      ENDDO
    ENDDO
    
    DO i_cgf = 1, N_cgf
      DO j_cgf = 1, N_cgf
        DO i_atm = 1, N_atm
          soc_cgf(i_cgf, j_cgf, 1:3) &
          &= soc_cgf(i_cgf, j_cgf, 1:3) &
          &+ soc_cgf_atm(i_cgf, j_cgf, i_atm, 1:3)
        ENDDO
      ENDDO
    ENDDO

!   Saving soc_cgf in the d77 format 
    fid = Fsoc_cgf_x
    WRITE(6,'(1X, 2A)') 'Saving the x components to a file ', TRIM(ADJUSTL(Fname(fid)))
    OPEN(fid, FILE=Fname(fid))
      DO i_cgf = 1, n_cgf
        DO j_cgf = i_cgf, n_cgf
          WRITE(fid, *) soc_cgf(i_cgf, j_cgf, 1)
        ENDDO
      ENDDO
    CLOSE(fid)
    WRITE(6,'(1X, A)') 'Done'

    fid = Fsoc_cgf_y
    WRITE(6,'(1X, 2A)') 'Saving the y components to a file ', TRIM(ADJUSTL(Fname(fid)))
    OPEN(fid, FILE=Fname(fid))
      DO i_cgf = 1, n_cgf
        DO j_cgf = i_cgf, n_cgf
          WRITE(fid, *) soc_cgf(i_cgf, j_cgf, 2)
        ENDDO
      ENDDO
    CLOSE(fid)
    WRITE(6,'(1X, A)') 'Done'

    fid = Fsoc_cgf_z
    WRITE(6,'(1X, 2A)') 'Saving the z components to a file ', TRIM(ADJUSTL(Fname(fid)))
    OPEN(fid, FILE=Fname(fid))
      DO i_cgf = 1, n_cgf
        DO j_cgf = i_cgf, n_cgf
          WRITE(fid, *) soc_cgf(i_cgf, j_cgf, 3)
        ENDDO
      ENDDO
    CLOSE(fid)
    WRITE(6,'(1X, A)') 'Done'    
    WRITE(6,'(1X)')     
    
!   Saving soc_cgf in the IOP(3/33=1) format of Gaussian 16
    fid = Fsoc_cgf_x_g16_format
    WRITE(6,'(1X, A)') 'saving them in the IOP(3/33=1) format of Gaussian 16.'
    WRITE(6,'(1X, 2A)') 'Saving the x components to a file ', &
                 &TRIM(ADJUSTL(Fname(fid)))
    CALL write_g16_log_format(fid, -soc_cgf(:,:,1))
    !CALL write_mat_atm(fid, soc_cgf_atm(:,:,:,1)) ! 1 means the x components
    fid = Fsoc_cgf_y_g16_format
    WRITE(6,'(1X, 2A)') 'Saving the y components to a file ', &
                 &TRIM(ADJUSTL(Fname(fid)))
    CALL write_g16_log_format(fid, -soc_cgf(:,:,2))
    !CALL write_mat_atm(fid, soc_cgf_atm(:,:,:,2)) ! 2 means the y components
    fid = Fsoc_cgf_z_g16_format
    WRITE(6,'(1X, 2A)') 'Saving the z components to a file ', &
                 &TRIM(ADJUSTL(Fname(fid)))
    CALL write_g16_log_format(fid, -soc_cgf(:,:,3))
    !CALL write_mat_atm(fid, soc_cgf_atm(:,:,:,3)) ! 3 means the z components
    WRITE(6,'(1X, A)') 'Done'
  ELSE
  ENDIF

 CALL write_messages(3, Text_blank, type_program, name_program)
END SUBROUTINE calc_int_cgf
