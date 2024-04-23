! This subroutine is part of d77 and computes
! vibronic coupling density.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE cube_vc_density(dmat_cgf, int_rho, int_vcd)
  USE ifmod, ONLY: write_messages
  USE global_constants
  USE global_read_input
  USE global_read_data
  USE func_rho, ONLY: func_phi_phi_cgf, func_rho_0
  USE func_potential, ONLY: func_dvne
  USE cube_write
  IMPLICIT NONE

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'cube_vc_density'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Input variable
  DOUBLE PRECISION, INTENT(IN) :: dmat_cgf(:,:)

! Output variables
  DOUBLE PRECISION, INTENT(OUT) :: int_rho
  DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT) :: int_vcd(:)
  
! Local variables
  INTEGER :: fid_rho
  CHARACTER(LEN=15) :: fname_rho
  INTEGER, DIMENSION(1:N_mode_calc) :: fid_dvne, fid_vcd
  CHARACTER(LEN=15), DIMENSION(1:N_mode_calc) :: fname_dvne, fname_vcd
  CHARACTER(LEN=100) :: comment1_rho, comment2_rho
  CHARACTER(LEN=2) :: text_bra, text_ket
  CHARACTER(LEN=100), DIMENSION(1:N_mode_calc) :: &
 &comment1_dvne, comment2_dvne
  CHARACTER(LEN=100), DIMENSION(1:N_mode_calc) :: &
 &comment1_vcd, comment2_vcd
  INTEGER :: ix, iy, iz, icol, ncol, mod_nz_ncol, i_line, n_line
  DOUBLE PRECISION, DIMENSION(1:3) :: xyz, dxyz, dxyz_i, dxyz_j
  DOUBLE PRECISION :: dxyz2, dxyz_i2
  DOUBLE PRECISION, DIMENSION(1:6) :: rho

  INTEGER :: i_mode_calc
  INTEGER :: n, n1, n2, n3, n4
  CHARACTER(LEN=1) :: n1_char, n2_char, n3_char, n4_char
  DOUBLE PRECISION, DIMENSION(1:N_mode_calc, 1:6) :: dvne, vcd
  DOUBLE PRECISION, DIMENSION(1:N_cgf, 1:N_cgf) :: phi_phi_cgf 

  IF(Debug == 'Yes') &
 &CALL write_messages(2, Text_blank, type_program, name_program)

! Local variables
  xyz = 0.0D0; dxyz = 0.0D0; dxyz2 = 0.0D0
  dxyz_i = 0.0D0; dxyz_j = 0.0D0; dxyz_i2 = 0.0D0; 
  phi_phi_cgf = 0.0D0
  comment1_rho = REPEAT(' ', 100); comment2_rho = REPEAT(' ', 100)
  text_bra = REPEAT(' ', 2); text_ket = REPEAT(' ', 2)

  int_vcd = 0.0D0

  CALL cube_grid
  CALL gen_lattice_point

  ncol = 6 ! The number of column of cube files
  mod_nz_ncol = MOD(nz, ncol)
  fid_rho = 1000
  fname_rho = 'RHO.cube'
  IF(Bra == Ket) THEN
    comment1_rho = 'Electron density'
  ELSE
    comment1_rho = 'Overlap density'
  ENDIF
  WRITE(text_bra, '(I2)') Bra
  WRITE(text_ket, '(I2)') Ket
  comment2_rho = '< '//TRIM(ADJUSTL(text_bra))//&
                &' | density operator | '//TRIM(ADJUSTL(text_ket))//' >'


  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') 'Generating and opening cube files'
  OPEN(fid_rho, FILE = fname_rho)
  CALL cube_write_header(fid_rho, comment1_rho, comment2_rho)

  DO i_mode_calc = 1, n_mode_calc
    fid_dvne(i_mode_calc) = i_mode_calc + 1000
    fid_vcd(i_mode_calc)  = i_mode_calc + 2000
    n  = mode_calc(i_mode_calc)
    n4 = (n - MOD(n,1000))/1000
    n  = n - n4*1000
    n3 = (n - MOD(n,100))/100
    n  = n - n3*100
    n2 = (n - MOD(n,10))/10
    n1 = n - n2*10
    CALL num2char(n4, n4_char)
    CALL num2char(n3, n3_char)
    CALL num2char(n2, n2_char)
    CALL num2char(n1, n1_char)
    fname_dvne(i_mode_calc) = 'DVNE'//n4_char//n3_char//&
                             n2_char//n1_char//'.cube'
    fname_vcd(i_mode_calc) = 'VCD'//n4_char//n3_char//&
                            n2_char//n1_char//'.cube'
    comment1_dvne(i_mode_calc) =&
    'Derivative of nuclear-electronic potential'
    comment2_dvne(i_mode_calc) = 'Mode = '//n4_char//n3_char//&
                              n2_char//n1_char
    comment1_vcd(i_mode_calc) =&
    'Vibronic coupling density'
    comment2_vcd(i_mode_calc) = 'Mode = '//n4_char//n3_char//&
                             n2_char//n1_char
    OPEN(fid_dvne(i_mode_calc), FILE = fname_dvne(i_mode_calc))
    CALL cube_write_header&
         (fid_dvne(i_mode_calc), &
          comment1_dvne(i_mode_calc), comment2_dvne(i_mode_calc))

    OPEN(fid_vcd(i_mode_calc),  FILE = fname_vcd(i_mode_calc))
    CALL cube_write_header&
         (fid_vcd(i_mode_calc), &
          comment1_vcd(i_mode_calc), comment2_vcd(i_mode_calc))

  ENDDO ! i_mode_calc
  WRITE(6,'(1X, A)') 'Done'

  WRITE(6,'(1X)') 
  IF(bra == ket) THEN
    WRITE(6,'(1X, A)') 'Calculating and writing electron density,'
  ELSE
    WRITE(6,'(1X, A)') 'Calculating and writing overlap density,'
  ENDIF
  WRITE(6,'(1X, A)') 'derivative of nuclear-electronic potential,'
  WRITE(6,'(1X, A)') 'and vibronic coupling density'
  
  int_rho = 0.0D0; int_vcd = 0.0D0
  IF(mod_nz_ncol == 0) THEN ! 6 columns
    n_line = nz/ncol ! The number of lines for each ix-iy pair
    DO ix = 1, nx
      xyz(1) = xyz_lattice_point(1, ix, 1, 1)  ! x coordinate
      DO iy = 1, ny
        xyz(2) = xyz_lattice_point(2, ix, iy, 1)  ! y coordinate
        DO i_line = 1, n_line
          rho = 0.0D0
          vcd = 0.0D0
          DO icol = 1, ncol
            iz = (i_line - 1)*ncol + icol
            xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
            phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
            rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
            dvne(:,icol) = func_dvne(xyz, N_mode_calc)
            vcd(:, icol) = rho(icol) * dvne(:, icol)

            int_rho = int_rho + rho(icol)
            int_vcd(:) = int_vcd(:) + vcd(:, icol)
          ENDDO ! icol = 1, ncol

          CALL cube_write_density(fid_rho, rho)
          DO i_mode_calc = 1, n_mode_calc 
            CALL cube_write_density(fid_dvne(i_mode_calc),&
                 dvne(i_mode_calc, 1:6))
            CALL cube_write_density(fid_vcd(i_mode_calc),&
                 vcd(i_mode_calc, 1:6))
          ENDDO
        ENDDO ! i_line = 1, n_line
      ENDDO ! iy = 1, ny
    ENDDO ! ix = 1, nx
  ELSE ! 1-5 columns
    n_line = (nz - mod_nz_ncol)/ncol + 1
    DO ix = 1, nx
      xyz(1) = xyz_lattice_point(1, ix, 1, 1)  ! x coordinate
      DO iy = 1, ny
        xyz(2) = xyz_lattice_point(2, ix, iy, 1)  ! y coordinate
        rho = 0.0D0
        vcd = 0.0D0
        IF(n_line == 1) THEN
          DO icol = 1, mod_nz_ncol
            iz = icol
            xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
            phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
            rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
 
            dvne(:,icol) = func_dvne(xyz, N_mode_calc)
            vcd(:, icol) = rho(icol) * dvne(:, icol)

            int_rho = int_rho + rho(icol)
            int_vcd(:) = int_vcd(:) + vcd(:, icol)
          ENDDO ! icol = 1, ncol

          CALL cube_write_density(fid_rho, rho(1:mod_nz_ncol))
          DO i_mode_calc = 1, n_mode_calc 
            CALL cube_write_density(fid_dvne(i_mode_calc),&
                 dvne(i_mode_calc, 1:mod_nz_ncol))
            CALL cube_write_density(fid_vcd(i_mode_calc),&
                 vcd(i_mode_calc, 1:mod_nz_ncol))
          ENDDO
        ELSE ! n_line > 2
          DO i_line = 1, n_line - 1
            rho = 0.0D0
            vcd = 0.0D0
            DO icol = 1, ncol
              iz = (i_line - 1)*ncol + icol
              xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
              phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
              rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)

              dvne(:,icol) = func_dvne(xyz, N_mode_calc)
              vcd(:, icol) = rho(icol) * dvne(:, icol)

              int_rho = int_rho + rho(icol)
              int_vcd(:) = int_vcd(:) + vcd(:, icol)
            ENDDO ! icol = 1, ncol
            CALL cube_write_density(fid_rho, rho)
            DO i_mode_calc = 1, n_mode_calc
              CALL cube_write_density(fid_dvne(i_mode_calc),&
                   dvne(i_mode_calc, 1:6))
              CALL cube_write_density(fid_vcd(i_mode_calc),&
                   vcd(i_mode_calc, 1:6))
            ENDDO

          ENDDO ! i_line = 1, n_line - 1
          ! The last line
          rho = 0.0D0
          vcd = 0.0D0
          DO icol = 1, mod_nz_ncol
            iz = (n_line - 1)*ncol + icol
            xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
            phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
            rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)

            dvne(:,icol) = func_dvne(xyz, N_mode_calc)
            vcd(:, icol) = rho(icol) * dvne(:, icol)

            int_rho = int_rho + rho(icol)
            int_vcd(:) = int_vcd(:) + vcd(:, icol)
          ENDDO ! icol = 1, ncol

          CALL cube_write_density(fid_rho, rho(1:mod_nz_ncol))
          DO i_mode_calc = 1, n_mode_calc
            CALL cube_write_density(fid_dvne(i_mode_calc),&
                 dvne(i_mode_calc, 1:mod_nz_ncol))
            CALL cube_write_density(fid_vcd(i_mode_calc),&
                 vcd(i_mode_calc, 1:mod_nz_ncol))
          ENDDO
        ENDIF
      ENDDO ! iy = 1, ny
    ENDDO ! ix = 1, nx
  ENDIF

  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') 'Closing cube files'
  DO i_mode_calc = 1, n_mode_calc
    CLOSE(fid_dvne(i_mode_calc))
    CLOSE(fid_vcd(i_mode_calc))
  ENDDO ! i_mode_calc = 1, n_mode_calc
  CLOSE(fid_rho)
  WRITE(6,'(1X, A)') 'Done'


  int_rho = int_rho*Dtau
  int_vcd = int_vcd*Dtau

  IF(Debug == 'Yes') &
 &CALL write_messages(3, Text_blank, type_program, name_program)
END SUBROUTINE cube_vc_density 
