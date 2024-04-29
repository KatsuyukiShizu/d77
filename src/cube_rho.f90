! This subroutine is part of d77 and writes 
! electron density or overlap (transition) density
! to a cube file, RHO.cube.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE cube_rho(dmat_cgf, int_rho)
  USE ifmod, ONLY: write_messages
  USE global_constants
  USE global_read_input
  USE global_read_data
  USE func_rho, ONLY: func_phi_phi_cgf, func_rho_0
  USE cube
  IMPLICIT NONE

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'cube_rho'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Input variable
  DOUBLE PRECISION, INTENT(IN)  :: dmat_cgf(:,:) 

! Output variable
  DOUBLE PRECISION, INTENT(OUT) :: int_rho

! Local variable
  INTEGER            :: fid_density
  CHARACTER(15)      :: fname_density
  CHARACTER(LEN=100) :: comment1, comment2
  CHARACTER(LEN=2)   :: text_bra, text_ket
  INTEGER            :: ix, iy, iz, icol, ncol, mod_nz_ncol, i_line, n_line
  DOUBLE PRECISION   :: xyz(1:3)
  DOUBLE PRECISION   :: rho(1:6)

  DOUBLE PRECISION :: phi_phi_cgf(1:N_cgf, 1:N_cgf)

  IF(Debug == 'Yes') &
 &CALL write_messages(2, Text_blank, type_program, name_program)

  xyz = 0.0D0
  phi_phi_cgf = 0.0D0
  comment1 = REPEAT(' ', 100); comment2 = REPEAT(' ', 100)
  text_bra = REPEAT(' ', 2); text_ket = REPEAT(' ', 2)

  CALL cube_grid ! global_read_data.f90
  ncol = 6 ! The number of column of cube files
  mod_nz_ncol = MOD(nz, ncol)

  CALL gen_lattice_point ! global_read_data.f90

  fid_density = 1000
  fname_density = 'RHO.cube'
  WRITE(text_bra, '(I2)') Bra
  WRITE(text_ket, '(I2)') Ket
  IF(Bra == Ket) THEN
    comment1 = 'Electron density of | '//TRIM(ADJUSTL(text_ket))//' >'
  ELSE
    comment1 = 'Overlap (transition) density between | '//TRIM(ADJUSTL(text_bra))//&
               ' > and | '//TRIM(ADJUSTL(text_ket))//' >'
  ENDIF
  comment2 = 'rho = < '//TRIM(ADJUSTL(text_bra))//&
            &' | density operator | '//TRIM(ADJUSTL(text_ket))//' >'

! Opening cube files          
  WRITE(6,'(1X)') 
  WRITE(6,'(1X, A)') 'Generating and opening cube files'
  OPEN(fid_density, FILE = fname_density)
  CALL cube_write_header(fid_density, comment1, comment2)
  WRITE(6,'(1X, A)') 'Done'

  WRITE(6,'(1X)')
  IF(bra == ket) THEN
    WRITE(6,'(1X, A)') 'Calculating and writing electron density'
  ELSE
    WRITE(6,'(1X, A)') 'Calculating and writing overlap density'
  ENDIF

  int_rho = 0.0D0
  IF(mod_nz_ncol == 0) THEN ! 6 columns
    n_line = nz/ncol ! The number of lines for each ix-iy pair
    DO ix = 1, nx
      xyz(1) = xyz_lattice_point(1, ix, 1, 1)  ! x coordinate
      DO iy = 1, ny
        xyz(2) = xyz_lattice_point(2, ix, iy, 1)  ! y coordinate
        DO i_line = 1, n_line
          rho = 0.0D0
          DO icol = 1, ncol
            iz = (i_line - 1)*ncol + icol
            xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
            phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
            rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
            int_rho = int_rho + rho(icol)
          ENDDO ! icol = 1, ncol
          CALL cube_write_density(fid_density, rho(1:ncol))
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
        IF(n_line == 1) THEN
          DO icol = 1, mod_nz_ncol
            iz = icol
            xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
            phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
            rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
            int_rho = int_rho + rho(icol)
          ENDDO ! icol = 1, ncol
          CALL cube_write_density(fid_density, rho(1:mod_nz_ncol))
        ELSE ! n_line > 2
          DO i_line = 1, n_line - 1
            rho = 0.0D0
            DO icol = 1, ncol
              iz = (i_line - 1)*ncol + icol
              xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
              phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
              rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
              int_rho = int_rho + rho(icol)
            ENDDO ! icol = 1, ncol
            CALL cube_write_density(fid_density, rho(1:ncol))
          ENDDO ! i_line = 1, n_line - 1
          ! The last line
          rho = 0.0D0
          DO icol = 1, mod_nz_ncol
            iz = (n_line - 1)*ncol + icol
            xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
            phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
            rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
            int_rho = int_rho + rho(icol)
          ENDDO ! icol = 1, ncol
          CALL cube_write_density(fid_density, rho(1:mod_nz_ncol))
        ENDIF
      ENDDO ! iy = 1, ny
    ENDDO ! ix = 1, nx
  ENDIF

  int_rho = int_rho*Dtau

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Closing cube files'
  CLOSE(fid_density)
  WRITE(6,'(1X, A)') 'Done'

  IF(Debug == 'Yes') &
 &CALL write_messages(3, Text_blank, type_program, name_program)

END SUBROUTINE cube_rho

