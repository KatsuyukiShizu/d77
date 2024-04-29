! This subroutine is part of d77 and computes
! transition dipole density.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE cube_dipole_density(projection, dmat_cgf, int_rho, int_dd, int_dd_p)
  USE ifmod, ONLY: write_messages
  USE global_constants
  USE global_read_input
  USE global_read_data
  USE func_rho, ONLY: func_phi_phi_cgf, func_rho_0
  USE cube
  IMPLICIT NONE

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'cube_dipole_density'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Input variable
  CHARACTER(LEN=4), INTENT(IN) :: projection
  DOUBLE PRECISION, INTENT(IN) :: dmat_cgf(:,:)

! Input and output variables
  DOUBLE PRECISION, INTENT(INOUT) :: int_rho
  DOUBLE PRECISION, INTENT(INOUT) :: int_dd(1:3), int_dd_p
  
! Local variables
  INTEGER :: i_xyz

  INTEGER :: fid_rho
  CHARACTER(LEN=15) :: fname_rho
  INTEGER           :: fid_er(1:3)
  CHARACTER(LEN=15) :: fname_er(1:3)
  INTEGER           :: fid_dd(1:4)
  CHARACTER(LEN=15) :: fname_dd(1:4)

  CHARACTER(LEN=100) :: comment1_rho, comment2_rho
  CHARACTER(LEN=2)   :: text_bra, text_ket
  CHARACTER(LEN=100) :: comment1_er(1:3), comment2_er(1:3)
  CHARACTER(LEN=100) :: comment1_dd(1:4), comment2_dd(1:4)

  INTEGER          :: ix, iy, iz, icol, ncol, mod_nz_ncol, i_line, n_line
  DOUBLE PRECISION :: xyz(1:3), dxyz(1:3), dxyz_i(1:3), dxyz_j(1:3)
  DOUBLE PRECISION :: dxyz2, dxyz_i2
  DOUBLE PRECISION :: rho(1:6)

  DOUBLE PRECISION :: phi_phi_cgf(1:N_cgf, 1:N_cgf)
  DOUBLE PRECISION :: minus_er(1:3, 1:6)
  DOUBLE PRECISION :: dd(1:3, 1:6), dd_p(1:6)
  DOUBLE PRECISION :: dipole_unit_vector(1:3)

  IF(Debug == 'Yes') &
 &CALL write_messages(2, Text_blank, type_program, name_program)

! Local variables
  xyz = 0.0D0; dxyz = 0.0D0; dxyz2 = 0.0D0
  dxyz_i = 0.0D0; dxyz_j = 0.0D0; dxyz_i2 = 0.0D0; 
  phi_phi_cgf = 0.0D0
  comment1_rho = REPEAT(' ', 100); comment2_rho = REPEAT(' ', 100)
  text_bra = REPEAT(' ', 2); text_ket = REPEAT(' ', 2)

  ncol = 6 ! The number of column of cube files

  IF(projection == 'No') THEN
    CALL cube_grid
    CALL gen_lattice_point
  ELSE
  ENDIF  
  mod_nz_ncol = MOD(nz, ncol)

  IF(projection == 'No') THEN
    
    int_dd = 0.0D0
  
    fid_er(1) = 1001; fid_er(2) = 1002; fid_er(3) = 1003 
    fid_dd(1) = 2001; fid_dd(2) = 2002; fid_dd(3) = 2003
    fname_er(1) = 'ER_X.cube'
    fname_er(2) = 'ER_Y.cube'
    fname_er(3) = 'ER_Z.cube'
    IF(Bra == Ket) THEN
      fname_dd(1) = 'PDMD_X.cube'
      fname_dd(2) = 'PDMD_Y.cube'
      fname_dd(3) = 'PDMD_Z.cube'
      comment1_dd(1:3) = 'Permanent dipole moment density in a.u.'
      comment2_dd(1) = 'x component'
      comment2_dd(2) = 'y component'
      comment2_dd(3) = 'z component'
    ELSE
      fname_dd(1) = 'TDMD_X.cube'
      fname_dd(2) = 'TDMD_Y.cube'
      fname_dd(3) = 'TDMD_Z.cube'
      comment1_dd(1:3) = 'Transition dipole moment density in a.u.'
      comment2_dd(1) = 'x component'
      comment2_dd(2) = 'y component'
      comment2_dd(3) = 'z component'
    ENDIF
  
    comment1_er(1) = 'Spatial distribution of - e * x in a.u.'
    comment1_er(2) = 'Spatial distribution of - e * y in a.u.'
    comment1_er(3) = 'Spatial distribution of - e * z in a.u.'
    comment2_er(1) = '- e * x is independent of electronic states'
    comment2_er(2) = '- e * y is independent of electronic states'
    comment2_er(3) = '- e * z is independent of electronic states'

    DO i_xyz = 1, 3
      OPEN(fid_er(i_xyz), FILE = fname_er(i_xyz))
      CALL cube_write_header&
           (fid_er(i_xyz), &
            comment1_er(i_xyz), comment2_er(i_xyz))
      OPEN(fid_dd(i_xyz), FILE = fname_dd(i_xyz))
      CALL cube_write_header&
           (fid_dd(i_xyz), &
            comment1_dd(i_xyz), comment2_dd(i_xyz))
    ENDDO ! i_xyz
   
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Generating and opening cube files.'
    WRITE(6,'(1X)')
    IF(bra == ket) THEN
      WRITE(6,'(1X, A)') 'Calculating electron density and permanent dipole moment densities.'
    ELSE
      WRITE(6,'(1X, A)') 'Calculating overlap density and transition dipole moment densities.'
    ENDIF
    WRITE(6,'(1X)')
  
    int_rho = 0.0D0; int_dd = 0.0D0
    IF(Save_rho_cube == 'Yes') THEN
      WRITE(6,'(1X, A)') 'Saving them to cube files.'
      fid_rho = 1000
      fname_rho = 'RHO.cube'
      WRITE(text_bra, '(I2)') Bra
      WRITE(text_ket, '(I2)') Ket
      IF(Bra == Ket) THEN
        comment1_rho = 'Electron density of | '//TRIM(ADJUSTL(text_ket))//' >'
      ELSE
        comment1_rho = 'Overlap (transition) density between | '//TRIM(ADJUSTL(text_bra))//&
                      &' > and | '//TRIM(ADJUSTL(text_ket))//' >'
      ENDIF
      comment2_rho = 'rho = < '//TRIM(ADJUSTL(text_bra))//&
                    &' | density operator | '//TRIM(ADJUSTL(text_ket))//' >'
  
      OPEN(fid_rho, FILE = fname_rho)
      CALL cube_write_header(fid_rho, comment1_rho, comment2_rho)

      IF(mod_nz_ncol == 0) THEN ! 6 columns
        n_line = nz/ncol ! The number of lines for each ix-iy pair
        DO ix = 1, nx
          xyz(1) = xyz_lattice_point(1, ix, 1, 1)  ! x coordinate
          DO iy = 1, ny
            xyz(2) = xyz_lattice_point(2, ix, iy, 1)  ! y coordinate
            DO i_line = 1, n_line
              rho = 0.0D0
              minus_er = 0.0D0
              dd = 0.0D0
              DO icol = 1, ncol
                iz = (i_line - 1)*ncol + icol
                xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
                minus_er(1:3,icol) = - xyz(1:3)
                dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)
    
                int_rho = int_rho + rho(icol)
                int_dd(1:3) = int_dd(1:3) + dd(1:3, icol)
              ENDDO ! icol = 1, ncol
    
              CALL cube_write_density(fid_rho, rho)
              DO i_xyz = 1, 3 
                CALL cube_write_density(fid_er(i_xyz),&
                     minus_er(i_xyz, 1:6))
                CALL cube_write_density(fid_dd(i_xyz),&
                     dd(i_xyz, 1:6))
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
            minus_er = 0.0D0
            dd = 0.0D0
            IF(n_line == 1) THEN
              DO icol = 1, mod_nz_ncol
                iz = icol
                xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
     
                minus_er(1:3,icol) = - xyz(1:3)
                dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)
    
                int_rho = int_rho + rho(icol)
                int_dd(1:3) = int_dd(1:3) + dd(1:3, icol)
              ENDDO ! icol = 1, ncol
    
              CALL cube_write_density(fid_rho, rho(1:mod_nz_ncol))
              DO i_xyz = 1, 3 
                CALL cube_write_density(fid_er(i_xyz),&
                     minus_er(i_xyz, 1:mod_nz_ncol))
                CALL cube_write_density(fid_dd(i_xyz),&
                     dd(i_xyz, 1:mod_nz_ncol))
              ENDDO
            ELSE ! n_line > 2
              DO i_line = 1, n_line - 1
                rho = 0.0D0
                minus_er = 0.0D0
                dd = 0.0D0
                DO icol = 1, ncol
                  iz = (i_line - 1)*ncol + icol
                  xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                  phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                  rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
    
                  minus_er(1:3,icol) = - xyz(1:3)
                  dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)
             
                  int_rho = int_rho + rho(icol)
                  int_dd(1:3) = int_dd(1:3) + dd(1:3, icol)
                ENDDO ! icol = 1, ncol
    
                CALL cube_write_density(fid_rho, rho)
                DO i_xyz = 1, 3 
                  CALL cube_write_density(fid_er(i_xyz),&
                       minus_er(i_xyz, 1:6))
                  CALL cube_write_density(fid_dd(i_xyz),&
                       dd(i_xyz, 1:6))
                ENDDO
    
              ENDDO ! i_line = 1, n_line - 1
              ! The last line
              rho = 0.0D0
              minus_er = 0.0D0
              dd = 0.0D0
              DO icol = 1, mod_nz_ncol
                iz = (n_line - 1)*ncol + icol
                xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
    
                minus_er(1:3,icol) = - xyz(1:3)
                dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)
             
                int_rho = int_rho + rho(icol)
                int_dd(1:3) = int_dd(1:3) + dd(1:3, icol)
              ENDDO ! icol = 1, ncol
    
              CALL cube_write_density(fid_rho, rho(1:mod_nz_ncol))
              DO i_xyz = 1, 3 
                CALL cube_write_density(fid_er(i_xyz),&
                     minus_er(i_xyz, 1:mod_nz_ncol))
                CALL cube_write_density(fid_dd(i_xyz),&
                     dd(i_xyz, 1:mod_nz_ncol))
              ENDDO
            ENDIF
          ENDDO ! iy = 1, ny
        ENDDO ! ix = 1, nx
      ENDIF
      CLOSE(fid_rho)
    ELSE ! Save_rho_cube /= Yes
      IF(bra == ket) THEN
        WRITE(6,'(1X, A)') 'Saving permanent dipole moment densities to cube files.'
      ELSE
        WRITE(6,'(1X, A)') 'Saving transition dipole moment densities to cube files.'
      ENDIF
      IF(mod_nz_ncol == 0) THEN ! 6 columns
        n_line = nz/ncol ! The number of lines for each ix-iy pair
        DO ix = 1, nx
          xyz(1) = xyz_lattice_point(1, ix, 1, 1)  ! x coordinate
          DO iy = 1, ny
            xyz(2) = xyz_lattice_point(2, ix, iy, 1)  ! y coordinate
            DO i_line = 1, n_line
              rho = 0.0D0
              minus_er = 0.0D0
              dd = 0.0D0
              DO icol = 1, ncol
                iz = (i_line - 1)*ncol + icol
                xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
                minus_er(1:3,icol) = - xyz(1:3)
                dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)

                int_rho = int_rho + rho(icol)
                int_dd(1:3) = int_dd(1:3) + dd(1:3, icol)
              ENDDO ! icol = 1, ncol

              DO i_xyz = 1, 3
                CALL cube_write_density(fid_er(i_xyz),&
                     minus_er(i_xyz, 1:6))
                CALL cube_write_density(fid_dd(i_xyz),&
                     dd(i_xyz, 1:6))
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
            minus_er = 0.0D0
            dd = 0.0D0
            IF(n_line == 1) THEN
              DO icol = 1, mod_nz_ncol
                iz = icol
                xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)

                minus_er(1:3,icol) = - xyz(1:3)
                dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)

                int_rho = int_rho + rho(icol)
                int_dd(1:3) = int_dd(1:3) + dd(1:3, icol)
              ENDDO ! icol = 1, ncol

              DO i_xyz = 1, 3
                CALL cube_write_density(fid_er(i_xyz),&
                     minus_er(i_xyz, 1:mod_nz_ncol))
                CALL cube_write_density(fid_dd(i_xyz),&
                     dd(i_xyz, 1:mod_nz_ncol))
              ENDDO
            ELSE ! n_line > 2
              DO i_line = 1, n_line - 1
                rho = 0.0D0
                minus_er = 0.0D0
                dd = 0.0D0
                DO icol = 1, ncol
                  iz = (i_line - 1)*ncol + icol
                  xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                  phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                  rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)

                  minus_er(1:3,icol) = - xyz(1:3)
                  dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)

                  int_rho = int_rho + rho(icol)
                  int_dd(1:3) = int_dd(1:3) + dd(1:3, icol)
                ENDDO ! icol = 1, ncol

                DO i_xyz = 1, 3
                  CALL cube_write_density(fid_er(i_xyz),&
                       minus_er(i_xyz, 1:6))
                  CALL cube_write_density(fid_dd(i_xyz),&
                       dd(i_xyz, 1:6))
                ENDDO

              ENDDO ! i_line = 1, n_line - 1
              ! The last line
              rho = 0.0D0
              minus_er = 0.0D0
              dd = 0.0D0
              DO icol = 1, mod_nz_ncol
                iz = (n_line - 1)*ncol + icol
                xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)

                minus_er(1:3,icol) = - xyz(1:3)
                dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)

                int_rho = int_rho + rho(icol)
                int_dd(1:3) = int_dd(1:3) + dd(1:3, icol)
              ENDDO ! icol = 1, ncol

              DO i_xyz = 1, 3
                CALL cube_write_density(fid_er(i_xyz),&
                     minus_er(i_xyz, 1:mod_nz_ncol))
                CALL cube_write_density(fid_dd(i_xyz),&
                     dd(i_xyz, 1:mod_nz_ncol))
              ENDDO
            ENDIF
          ENDDO ! iy = 1, ny
        ENDDO ! ix = 1, nx
      ENDIF ! mod_nz_ncol = 0 or /= 0
    ENDIF ! Save_rho_cube = Yes or /= Yes

    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Closing cube files.'
    DO i_xyz = 1, 3
      CLOSE(fid_er(i_xyz))
      CLOSE(fid_dd(i_xyz))
    ENDDO ! i_xyz = 1, 3
  
    int_rho = int_rho*Dtau
    int_dd = int_dd*Dtau
  
  ELSEIF(projection == 'Yes') THEN
    
    dipole_unit_vector = 0.0D0
    dipole_unit_vector = int_dd/DSQRT(DOT_PRODUCT(int_dd, int_dd))

    fid_dd(4) = 2004
    fname_dd(4) = 'TDMD_MU.cube' ! Projected dipole moment
    comment1_dd(4) = 'Transition dipole moment density in a.u.'
    comment2_dd(4) = 'Projected to the direction of dipole moment vector'
    OPEN(fid_dd(4), FILE = fname_dd(4))
    CALL cube_write_header&
         (fid_dd(4), &
          comment1_dd(4), comment2_dd(4))

    int_dd_p = 0.0D0
    IF(mod_nz_ncol == 0) THEN ! 6 columns
      n_line = nz/ncol ! The number of lines for each ix-iy pair
      DO ix = 1, nx
        xyz(1) = xyz_lattice_point(1, ix, 1, 1)  ! x coordinate
        DO iy = 1, ny
          xyz(2) = xyz_lattice_point(2, ix, iy, 1)  ! y coordinate
          DO i_line = 1, n_line
            rho = 0.0D0
            minus_er = 0.0D0
            dd = 0.0D0
            dd_p = 0.0D0
            DO icol = 1, ncol
              iz = (i_line - 1)*ncol + icol
              xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
              phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
              rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)
              minus_er(1:3,icol) = - xyz(1:3)
              dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)
              dd_p(icol) =  DOT_PRODUCT(dd(1:3, icol), dipole_unit_vector)
              int_dd_p = int_dd_p + dd_p(icol)
            ENDDO ! icol = 1, ncol
            CALL cube_write_density(fid_dd(4), dd_p(1:6))
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
          minus_er = 0.0D0
          dd = 0.0D0
          dd_p = 0.0D0
          IF(n_line == 1) THEN
            DO icol = 1, mod_nz_ncol
              iz = icol
              xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
              phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
              rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)

              minus_er(1:3,icol) = - xyz(1:3)
              dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)
              dd_p(icol) =  DOT_PRODUCT(dd(1:3, icol), dipole_unit_vector)
              int_dd_p = int_dd_p + dd_p(icol)
            ENDDO ! icol = 1, ncol
            CALL cube_write_density(fid_dd(4), dd_p(1:mod_nz_ncol))
          ELSE ! n_line > 2
            DO i_line = 1, n_line - 1
              rho = 0.0D0
              minus_er = 0.0D0
              dd = 0.0D0
              dd_p = 0.0D0
              DO icol = 1, ncol
                iz = (i_line - 1)*ncol + icol
                xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
                phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
                rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)

                minus_er(1:3,icol) = - xyz(1:3)
                dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)
                dd_p(icol) =  DOT_PRODUCT(dd(1:3, icol), dipole_unit_vector)
                int_dd_p = int_dd_p + dd_p(icol)
              ENDDO ! icol = 1, ncol
              CALL cube_write_density(fid_dd(4), dd_p(1:6))
            ENDDO ! i_line = 1, n_line - 1
            ! The last line
            rho = 0.0D0
            minus_er = 0.0D0
            dd = 0.0D0
            dd_p = 0.0D0
            DO icol = 1, mod_nz_ncol
              iz = (n_line - 1)*ncol + icol
              xyz(3) = xyz_lattice_point(3, ix, iy, iz)  ! z coordinate
              phi_phi_cgf = func_phi_phi_cgf(xyz, N_cgf)
              rho(icol) = func_rho_0(dmat_cgf, phi_phi_cgf)

              minus_er(1:3,icol) = - xyz(1:3)
              dd(1:3, icol) = rho(icol) * minus_er(1:3, icol)
              dd_p(icol) =  DOT_PRODUCT(dd(1:3, icol), dipole_unit_vector)
              int_dd_p = int_dd_p + dd_p(icol)
            ENDDO ! icol = 1, ncol
            CALL cube_write_density(fid_dd(4), dd_p(1:mod_nz_ncol))
          ENDIF
        ENDDO ! iy = 1, ny
      ENDDO ! ix = 1, nx
    ENDIF
    
    CLOSE(fid_dd(4))
    int_dd_p = int_dd_p*Dtau
  
  ELSE
    CALL write_messages(-9999, Text_blank, type_program, name_program)
  ENDIF

  IF(Debug == 'Yes') &
 &CALL write_messages(3, Text_blank, type_program, name_program)

END SUBROUTINE cube_dipole_density 
