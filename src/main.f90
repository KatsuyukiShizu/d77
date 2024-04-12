! The main program of d77
!
! d77 provides visual understanding of molecular properties.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

! Global constants are defined in MODULE 'global_constants.f90'.
PROGRAM main 
  USE ifmod
  USE global_constants
  USE global_read_input
  USE global_read_data
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'main'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'PROGRAM'

! CPU time
  CHARACTER(LEN=100) :: text_date_time


! ---------------
! Local variables
! ---------------

! DO loop variables  
  INTEGER :: i_atm

! Assigning global variables
  CALL assign_globals ! SUBROUTINE 'assign_globals' in MODULE 'global_constants.f90'.

  CALL write_messages(0, text_date_time, type_program, name_program)

! Writing start time
  Time_start = 0.0D0; Time_end = 0.0D0
  CALL CPU_TIME(Time_start)
  text_date_time = '--- Entering d77.exe'
  CALL write_messages(1, text_date_time, type_program, name_program)
  
  CALL read_input_control ! SUBROUTINE 'read_input_control' in MODULE 'global_read_input.f90'

! ---------------------------------------
! Reading general information of a system
! ---------------------------------------

! N_atm: Total number of atoms 
! Charg: Total charge
! Mult:  Total spin multiplicity 
! Ne:    Toal number of electrons
! N_cgf: Total number of contracted Gaussian functions (CGFs) 

  CALL read_data_elec_0 ! MODULE global_read_data.f90

  SELECT CASE(Runtyp)
    CASE('Int_pgf', 'Density')
      SELECT CASE(Property)
        CASE('Vc')
          CALL read_data_vib_0
          CALL read_data_vib_1
          CALL read_mode_calc
        CASE DEFAULT
      END SELECT
    CASE DEFAULT
  END SELECT

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)')     'Information of the system'
  WRITE(6,'(1X, A)')     '============================='
  WRITE(6,'(1X, A, I0)') 'Number of atoms:         ', N_atm
  WRITE(6,'(1X, A, I0)') 'Charge:                  ', Charg
  WRITE(6,'(1X, A, I0)') 'Spin multiplicity:       ', Mult
  WRITE(6,'(1X, A, I0)') 'Number of electrons:     ', Ne
  WRITE(6,'(1X, A, I0)') 'Number of CGFs:          ', N_cgf
  WRITE(6,'(1X, A, I0)') 'Number of spin orbitals: ', N_so
  IF(N_mode_calc == 0) THEN
  ELSE
    WRITE(6,'(1X, A, I0)') 'Number of vib. modes:    ', N_mode_calc
  ENDIF
  WRITE(6,'(1X, A)')     '============================='


! -----------------------------------
! Read information of nuclei and PGFs
! -----------------------------------

! --- Information of nuclei ---
! Atmnum:   atomic numbers
! Nuccharg: nuclear charges
! Xyznuc:   xyz coordinates of nuclei
! Elmsym:   element symbols
! --- Information of PGFs ---
! Xyzao:      xyz cooedinates of PGF centers
! Orbtyp:     orbital types
! Lmn:        (l, m, n) sets
! Cntexp:     contraction exponents
! Cntcoef:    expantion coefficients for CGFs
! Cntnormfac: normalization factors

  CALL read_data_elec_1 ! MODULE global_read_data.f90
  CALL so_num2so_name   ! MODULE global_read_data.f90

  WRITE(6,'(1X)')
  WRITE(6,'(1X, A)') 'Nuclear coordinates (bohr)'
  WRITE(6,'(1X, A)') '=============================================='
  WRITE(6,'(1X, A)') 'Atom        x             y             z'
  WRITE(6,'(1X, A)') '----------------------------------------------'
  DO i_atm = 1, N_atm
    WRITE(6,'(1X, I3, F15.8, 2F14.8)') i_atm, Xyznuc(1:3, i_atm)
  ENDDO
  WRITE(6,'(1X, A)') '=============================================='
  WRITE(6,'(1X)')
  IF(N_mode_calc /= 0) THEN
    WRITE(6,'(1X, A)') 'Nuclear charges (elementary charge) and'
    WRITE(6,'(1X, A)') 'atomic weights (atomic mass unit)'
    WRITE(6,'(1X, A)') '============================================='
    WRITE(6,'(1X, A)') 'Atom   Element   Atomic   Nuclear     Atomic '
    WRITE(6,'(1X, A)') '       symbol    number   charge      weight '
    WRITE(6,'(1X, A)') '---------------------------------------------'
    DO i_atm = 1, N_atm
      WRITE(6,'(1X, I3, 7X, A2, 6X, I3, 5X, F5.1, 3X, F11.5)') &
     &i_atm, Elmsym(i_atm), Atmnum(i_atm),&
      Nuccharg(i_atm), Atmwt(i_atm)*Me/Amc
    ENDDO
    WRITE(6,'(1X, A)') '============================================='
  ELSE
    WRITE(6,'(1X, A)') 'Nuclear charges (elementary charge)'
    WRITE(6,'(1X, A)') '================================='
    WRITE(6,'(1X, A)') 'Atom   Element   Atomic   Nuclear'
    WRITE(6,'(1X, A)') '       symbol    number   charge '
    WRITE(6,'(1X, A)') '---------------------------------'
    DO i_atm = 1, N_atm
      WRITE(6,'(1X, I3, 7X, A2, 6X, I3, 5X, F5.1)') &
     &i_atm, Elmsym(i_atm), Atmnum(i_atm), Nuccharg(i_atm)
    ENDDO
    WRITE(6,'(1X, A)') '================================='
  ENDIF
  SELECT CASE(Property)
    CASE('Soc')
      WRITE(6,'(1X)')
      WRITE(6,'(1X, A)') 'Zeff for SOC calculations (elementary charge)'
      WRITE(6,'(1X, A)') '================================='
      WRITE(6,'(1X, A)') 'Atom   Element   Atomic   Zeff   '
      WRITE(6,'(1X, A)') '       symbol    number          '
      WRITE(6,'(1X, A)') '---------------------------------'
      DO i_atm = 1, N_atm
        WRITE(6,'(1X, I3, 7X, A2, 6X, I3, 5X, F5.1)') &
       &i_atm, Elmsym(i_atm), Atmnum(i_atm), Nuccharg_soc(i_atm)
      ENDDO
      WRITE(6,'(1X, A)') '================================='
    CASE DEFAULT
  END SELECT


! -----------------------------------------------
! Reading CI coefficients and excitation energies
! -----------------------------------------------

! Cntcoef:   CI coefficients of slater determinants
! Ex_energy: Excitation energies of electronic states

  IF(Runtyp == 'Calc_int_cgf' .OR. Property == 'Dvne_only') THEN
    Read_ci_coef = 'No'
  ELSE
  ENDIF

  IF(Read_ci_coef == 'Yes') THEN ! Default
    CALL read_data_ex_energy ! global_read_data.f90
    SELECT CASE(Method)
      CASE('Td')
        CALL read_data_ci_coef_td ! global_read_data.f90
      CASE DEFAULT
    END SELECT
  ELSE
  ENDIF

  SELECT CASE(Runtyp)

    CASE('Calc_int_cgf') ! Runtyp = Calc_int_cgf
      CALL calc_int_cgf 

    CASE('Int_pgf') ! Runtyp = Int_pgf
      SELECT CASE(Property)
        CASE('Dipole')
          CALL dipole_int_pgf
        CASE('Vc')
          CALL vc_int_pgf
        CASE('Soc')
          CALL soc_int_pgf
        CASE DEFAULT ! Invalid Property
      END SELECT  

    !CASE('Density') ! Runtyp = Density
    !  SELECT CASE(Property)
    !    CASE('Vc')
    !      CALL vc_density
    !    CASE('Soc')
    !      CALL soc_density
    !    CASE DEFAULT ! Invalid Property
    !  END SELECT  

    CASE DEFAULT ! Invalid Runtyp
  END SELECT  





  !SELECT CASE(Property)
    !CASE('Calc_int_cgf')
    !  CALL calc_int_cgf 
    !CASE('Dvne_only')
    !  CALL cube_mode_only 
    !CASE('Opr_dmat')
    !  CALL opr_dmat_only
    !CASE('Electron_density', 'Rho')
    !  CALL electron_density
    !CASE('Pe')
    !  SELECT CASE(Runtyp)
    !    CASE('Density') ! Runtyp = Density
    !      CALL pe_density
    !    CASE('Int_pgf') ! Runtyp = Int_pgf
    !      CALL pe_int_pgf
    !    CASE DEFAULT
    !  END SELECT
    !CASE('Kinetic_energy', 'Ke')
    !  SELECT CASE(Runtyp)
    !    CASE('Density') ! Runtyp = Density
    !      CALL ke_density
    !    CASE('Int_pgf') ! Runtyp = Int_pgf
    !      CALL ke_int_pgf
    !    CASE DEFAULT
    !  END SELECT
    !CASE('Hcore', 'Kepe')
    !  SELECT CASE(Runtyp)
    !    CASE('Density') ! Runtyp = Density
    !      CALL hcore_density
    !    CASE('Int_pgf') ! Runtyp = Int_pgf
    !      CALL hcore_int_pgf
    !    CASE DEFAULT
    !  END SELECT
    !CASE('Dipole')
    !  SELECT CASE(Runtyp)
    !    CASE('Density') ! Runtyp = Density
    !      CALL dipole_density
    !    CASE('Int_pgf') ! Runtyp = Int_pgf
    !      CALL dipole_int_pgf
    !    CASE DEFAULT
    !  END SELECT
    !CASE('Vc')
    !  SELECT CASE(Runtyp)
    !    !CASE('Density') ! Runtyp = Density
    !    !  CALL vc_density
    !    CASE('Int_pgf') ! Runtyp = Int_pgf
    !      CALL vc_int_pgf
    !    CASE DEFAULT
    !  END SELECT
    !CASE('Soc')
    !  SELECT CASE(runtyp)
    !    !CASE('density') ! runtyp = density
    !    !  CALL soc_density
    !    CASE('Int_pgf') ! runtyp = int_pgf
    !      IF(Asoc == 'Yes') THEN
    !        !CALL asoc_int_pgf
    !      ELSE
    !        CALL soc_int_pgf
    !      ENDIF
    !    CASE DEFAULT
    !  END SELECT
    !CASE DEFAULT ! Property = ?
  !END SELECT  

! Writing end time
  text_date_time = '--- Leaving d77.exe'
  CALL write_messages(9999, text_date_time, type_program, name_program)

END PROGRAM main
