! This module is part of d77 and reads global variables.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE global_read_data

  USE ifmod, ONLY: write_messages
  USE global_constants
  USE global_read_input
  USE func_int_pgf, ONLY: func_int_pgf_pgf
  IMPLICIT NONE

! SUBROUTINE read_data_elec_0

! N_atm  : Total number of atoms
! Charg  : Total charge
! Mult   : Spin multiplicity
! Ne     : Toal number of electrons
! N_cgf  : Total number of independent contracted Gaussian functions (CGFs)
! N_so_a : Total number of alpha spin orbitals (N_cgf <= N_so_a: typically N_cgf = N_so_a)
! N_so   : Total number of spin orbitals (N_so = 2 * N_so_a)
  INTEGER, SAVE :: N_atm, Charg, Mult, Ne, N_cgf, N_so_a, N_so

! SUBROUTINE read_data_elec_1

! --- Information of nuclei ---
! Atmnum   : Set of atomic numbers
! Nuccharg : Set of nuclear charges (eg, 1.0 for H, 6.0 for C)
! Xyznuc   : Set of xyz coordinates of nuclei
! Elmsym   : Set of element symbols (H, He, Li ...)
  INTEGER, ALLOCATABLE, SAVE :: Atmnum(:) ! Atmnum(1:N_atm)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Nuccharg(:) ! Nuccharg(1:N_atm): 
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Xyznuc(:,:) ! Xyznuc(1:3, 1:N_atm): (1,2,3) = (x,y,z)
  CHARACTER(LEN=2), ALLOCATABLE, SAVE :: Elmsym(:) ! Elmsym(1:N_atm)

! --- Information of atomic orbitals (AOs) and primitive Gaussian functions (PGFs) ---
! Xyzao      : xyz coordinate of each AO center
! Orbtyp     : Orbital type (s, px, py, pz, dxx ...) of each AO
! Ao2atm     : Atom label to which each AO belongs
! N_pgf      : The number of PGFs of each AO
! Lmn        : (l, m, n) set of each AO
! Cntexp     : Contraction exponents for each AO
! Cntcoef    : Expantion coefficients for each AO
! Cntnormfac : The normalization factor for each AO
  INTEGER, ALLOCATABLE, SAVE :: Ao2atm(:), N_pgf(:) ! Ao2atm(1:N_cgf), N_pgf(1:N_cgf)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Xyzao(:,:) ! Xyzao(1:3, 1:N_cgf): (1,2,3) = (x,y,z)
  CHARACTER(LEN=5), ALLOCATABLE, SAVE :: Orbtyp(:) ! Orbtyp(1:N_cgf)
  INTEGER, ALLOCATABLE, SAVE :: Lmn(:,:) ! Lmn(1:3, 1:N_cgf): (1,2,3) = (x,y,z)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Cntexp(:,:) ! Cntexp(1:Max_n_pgf, 1:N_cgf)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Cntcoef(:,:) ! Cntcoef(1:Max_n_pgf, 1:N_cgf)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Cntnormfac(:,:) ! Cntnormfac(1:Max_n_pgf, 1:N_cgf)

! --- Information of molecular orbitals (MOs) ---
! Coef_so   : Expansion coefficients (MO coefficients) of spin orbitals
! Indentity : Identity matrix (size = N_so X N_so)  
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Coef_so(:,:) ! Coef_so(1:N_cgf, 1:N_so)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Identity(:,:) ! Identity(1:N_so, 1:N_so)

! so_num2so_name
  CHARACTER(LEN=21), ALLOCATABLE, SAVE :: So_name(:)  

! read_data_ci_coef
  INTEGER, SAVE :: N_state
  INTEGER, ALLOCATABLE, SAVE :: N_elec_config(:)
  INTEGER, ALLOCATABLE, SAVE :: Tuple(:,:)
  INTEGER, ALLOCATABLE, SAVE :: A_ras(:,:,:), R_ras(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Ci_coef(:,:)

! read_data_ci_coef_td
  !INTEGER, ALLOCATABLE, SAVE :: A(:,:), R(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: X_plus_y(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: X_minus_y(:,:)

! $GRID keywords
  DOUBLE PRECISION, SAVE :: Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
  DOUBLE PRECISION, SAVE :: Dx, Dy, Dz, Dtau
  DOUBLE PRECISION, SAVE :: Delta
  INTEGER, SAVE :: Nx, Ny, Nz

! gen_lattice_point
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Xyz_lattice_point(:,:,:,:)

! read_data_vib_0
  INTEGER, SAVE :: N_mode

! read_data_vib_1
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Freq(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vibmode(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Atmwt(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Sqrt_atmwt(:)

! read_mode_calc
  INTEGER, SAVE :: N_mode_calc
  INTEGER, ALLOCATABLE, SAVE :: Mode_calc(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Freq_calc(:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vibmode_calc(:,:,:)

! read_data_ex_energy
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Ex_energy(:)

! read_data_effcharg_soc
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: Nuccharg_soc(:)

  CONTAINS

!***************************************************************************************************
  SUBROUTINE read_data_elec_0

! read_data_elec_0 
! (1) reads the number of atoms, electrons, and CGFs. 
! (2) reads charge and spin multiplicity of the system. 

!***************************************************************************************************

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_data_elec_0'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Reading general information of the system from '//TRIM(Fname(Fcontrol_elec))
    OPEN(Fcontrol_elec, FILE = Fname(Fcontrol_elec))
      READ(Fcontrol_elec,*) N_atm  ! Total number of atoms                                             
      READ(Fcontrol_elec,*) Charg  ! Total charge
      READ(Fcontrol_elec,*) Mult   ! Spin multiplicity
      READ(Fcontrol_elec,*) Ne     ! Toal number of electrons
      READ(Fcontrol_elec,*) N_cgf  ! Total number of independent contracted Gaussian functions (CGFs)
      READ(Fcontrol_elec,*) N_so_a ! Total number of alpha spin orbitals (N_cgf <= N_so_a: typically N_cgf = N_so_a)
    CLOSE(Fcontrol_elec)
    N_so = 2*N_so_a ! Total number of spin orbitals 

!   Writing error messages
    IF(N_atm > Max_n_atm) THEN ! N_atm exceeds the upper limit.
      WRITE(6,'(1X, A)') 'Too many atoms'
      WRITE(6,'(1X, A, I0)') 'The number of atoms &
     &must be less than or equal to ', Max_n_atm
      CALL write_messages(-9999, Text_blank, type_program, name_program)
    ELSEIF(Ne > Max_ne) THEN ! Ne exceeds the upper limit.
      WRITE(6,'(1X, A)') 'Too many electrons'
      WRITE(6,'(1X, A, I0)') 'The number of electrons &
     &must be less than or equal to ', Max_ne
      CALL write_messages(-9999, Text_blank, type_program, name_program)
    ELSEIF(N_cgf > Max_n_cgf) THEN ! N_cgf exceeds the upper limit.
      WRITE(6,'(1X, A)') 'Too many CGFs'
      WRITE(6,'(1X, A, I0)') 'The number of contracted Gaussian functions (CGFs) &
     &must be less than or equal to ', Max_n_cgf
      CALL write_messages(-9999, Text_blank, type_program, name_program)
    ELSE
    ENDIF
    WRITE(6,'(1X, A)') 'Done'

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE read_data_elec_0


!***************************************************************************************************
  SUBROUTINE read_data_elec_1

! read_data_elec_1  
! (1) reads Cartesian nuclear coordinates, atomic numbers, nuclear chages of atoms.
! (2) determines element symbols of the atoms from their atomic numbers.
! (3) reads information of PGFs: 
!     orbital types, contraction exponents, and contraction coefficients.
! (4) determines (l, m, n) and normalization factor of each PGF
!     from its orbital type and contraction exponents. 
! (5) assigns the atom to which each PGF blongs and 
!     determines the coordinate of the center of the PGF.

!***************************************************************************************************

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_data_elec_1'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables
    INTEGER :: i_atm, i_xyz, i_cgf, j_cgf, i_pgf, j_pgf, i_so, j_so

!   normfac : Normalization factor of CGF
!   s_cgf   : Ovrlap integral between CGFs
    DOUBLE PRECISION, ALLOCATABLE :: normfac(:,:) ! normfac(1:Max_n_pgf, 1:N_cgf)
    DOUBLE PRECISION, ALLOCATABLE :: s_cgf(:,:)   ! s_cgf(1:N_cgf, 1:N_cgf)
    DOUBLE PRECISION :: orthonorm

    CHARACTER(LEN=100) :: text_error
    CHARACTER(LEN=3)   :: text_i_atm, text_atmnum
    CHARACTER(LEN=4)   :: text_i_cgf

    CHARACTER(LEN=3) :: large_overlap

!   Effective nuclear charges for SOC calculation 
    DOUBLE PRECISION :: fm

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   Initializing global variables
    ALLOCATE(Xyznuc(1:3, 1:N_atm)); Xyznuc  = 0.0D0
    ALLOCATE(Atmnum(1:N_atm)); Atmnum  = 0
    ALLOCATE(Nuccharg(1:N_atm)); Nuccharg = 0.0D0
    ALLOCATE(Elmsym(1:N_atm)); Elmsym = REPEAT(' ', 2)
    ALLOCATE(Ao2atm(1:N_cgf)); Ao2atm = 0
    ALLOCATE(N_pgf(1:N_cgf)); N_pgf = 0
    ALLOCATE(Orbtyp(1:N_cgf)); Orbtyp  = '    '
    ALLOCATE(Lmn(1:3, 1:N_cgf)); Lmn = 0
    ALLOCATE(Xyzao(1:3, 1:N_cgf)); Xyzao   = 0.0D0
    ALLOCATE(Cntexp(1:Max_n_pgf, 1:N_cgf)); Cntexp  = 0.0D0
    ALLOCATE(Cntcoef(1:Max_n_pgf, 1:N_cgf)); Cntcoef = 0.0D0
    ALLOCATE(Cntnormfac(1:Max_n_pgf, 1:N_cgf)); Cntnormfac = 0.0D0
    ALLOCATE(Coef_so(1:N_cgf, 1:N_so)); Coef_so = 0.0D0
    ALLOCATE(Identity(1:N_so, 1:N_so)); Identity = 0.0D0

!   Initializing local variables
    text_error = REPEAT(' ', 100)
    text_i_atm = REPEAT(' ', 3); text_atmnum = REPEAT(' ', 3)
    text_i_cgf = REPEAT(' ', 4)
    ALLOCATE(normfac(1:Max_n_pgf, 1:N_cgf)); normfac = 0.0D0
    ALLOCATE(s_cgf(1:N_cgf, 1:N_cgf)); s_cgf = 0.0D0

    DO i_so = 1, N_so
      Identity(i_so, i_so) = 1.0D0
    ENDDO

!   Reading the xyz coordinates (in bohr) of the atoms
    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Reading nuclear coordinates from '//TRIM(Fname(Fcoord))
    OPEN(Fcoord, FILE = Fname(Fcoord))
      DO i_atm = 1, N_atm
        DO i_xyz = 1, 3
          READ(Fcoord,*) Xyznuc(i_xyz, i_atm)
        ENDDO
      ENDDO
    CLOSE(Fcoord)

!   Reading the atomic numbers and nuclear charges of the atoms
    WRITE(6,'(1X, A)') 'Reading atomic numbers from '//TRIM(Fname(Fatmnum))
    WRITE(6,'(1X, A)') 'Reading nuclear charges from '//TRIM(Fname(Fnuccharg))
    OPEN(Fatmnum, FILE = Fname(Fatmnum))
    OPEN(Fnuccharg, FILE = Fname(Fnuccharg))
      DO i_atm = 1, N_atm
        READ(Fatmnum,*) Atmnum(i_atm)
        READ(Fnuccharg,*) Nuccharg(i_atm)
      ENDDO
    CLOSE(Fatmnum)
    CLOSE(Fnuccharg)
    WRITE(6,'(1X, A)') 'Done'

!   Assigning effective nuclear charges for SOC calculation
    IF(Property == 'Soc') THEN
      ALLOCATE(Nuccharg_soc(1:N_atm)); Nuccharg_soc(1:N_atm) = 0.0D0
      SELECT CASE(Effcharg_soc)
        CASE('Read') ! Default
          OPEN(Fnuccharg_soc, FILE = Fname(Fnuccharg_soc))
            DO i_atm = 1, N_atm
              READ(Fnuccharg_soc,*) Nuccharg_soc(i_atm)
            ENDDO
          CLOSE(Fnuccharg_soc)
          Nuccharg_soc = ABS(Nuccharg_soc)
        CASE('Atmnum')
          Nuccharg_soc(1:N_atm) = DBLE(Atmnum(1:N_atm))
        CASE('Nuccharg')
          Nuccharg_soc(1:N_atm) = Nuccharg(1:N_atm)
        CASE('Koseki')
          DO i_atm = 1, N_atm
            IF(2 <= Atmnum(i_atm) .AND. Atmnum(i_atm) <= 10) THEN
              fm = 0.40D0 + 0.05D0 * DBLE((Atmnum(i_atm) - 2))
              Nuccharg_soc(i_atm) = fm * Nuccharg(i_atm)
            ELSEIF(11 <= Atmnum(i_atm) .AND. Atmnum(i_atm) <= 18) THEN
              fm = 0.925D0 - 0.0125D0 * DBLE((Atmnum(i_atm) - 10))
              Nuccharg_soc(i_atm) = fm * Nuccharg(i_atm)
            ELSEIF(19 <= Atmnum(i_atm) .AND. Atmnum(i_atm) <= 20) THEN
              fm = 0.925D0 - 0.0125D0 * DBLE((Atmnum(i_atm) - 18))
              Nuccharg_soc(i_atm) = fm * Nuccharg(i_atm)
            ELSEIF(Atmnum(i_atm) == 53) THEN
              Nuccharg_soc(i_atm) = 65.72
            ELSE
            ENDIF
          ENDDO
        CASE DEFAULT
      END SELECT
    ELSE
    ENDIF

!   Determining element symbols from atomic numbers
    DO i_atm = 1, N_atm
      SELECT CASE(Atmnum(i_atm))
!       First row elements      
        CASE(1);  Elmsym(i_atm) = 'H'
        CASE(2);  Elmsym(i_atm) = 'He'

!       Second row elements       
        CASE(3);  Elmsym(i_atm) = 'Li'
        CASE(4);  Elmsym(i_atm) = 'Be'

        CASE(5);  Elmsym(i_atm) = 'B'
        CASE(6);  Elmsym(i_atm) = 'C'
        CASE(7);  Elmsym(i_atm) = 'N'
        CASE(8);  Elmsym(i_atm) = 'O'
        CASE(9);  Elmsym(i_atm) = 'F'
        CASE(10); Elmsym(i_atm) = 'Ne'

!       Third row elements               
        CASE(11); Elmsym(i_atm) = 'Na'
        CASE(12); Elmsym(i_atm) = 'Mg'

        CASE(13); Elmsym(i_atm) = 'Al'
        CASE(14); Elmsym(i_atm) = 'Si'
        CASE(15); Elmsym(i_atm) = 'P'
        CASE(16); Elmsym(i_atm) = 'S'
        CASE(17); Elmsym(i_atm) = 'Cl'
        CASE(18); Elmsym(i_atm) = 'Ar'

!       Fourth row elements               
        CASE(19); Elmsym(i_atm) = 'K'
        CASE(20); Elmsym(i_atm) = 'Ca'

        CASE(21); Elmsym(i_atm) = 'Sc'
        CASE(22); Elmsym(i_atm) = 'Ti'
        CASE(23); Elmsym(i_atm) = 'V'
        CASE(24); Elmsym(i_atm) = 'Cr'
        CASE(25); Elmsym(i_atm) = 'Mn'
        CASE(26); Elmsym(i_atm) = 'Fe'
        CASE(27); Elmsym(i_atm) = 'Co'
        CASE(28); Elmsym(i_atm) = 'Ni'
        CASE(29); Elmsym(i_atm) = 'Cu'
        CASE(30); Elmsym(i_atm) = 'Zn'

        CASE(31); Elmsym(i_atm) = 'Ga'
        CASE(32); Elmsym(i_atm) = 'Ge'
        CASE(33); Elmsym(i_atm) = 'As'
        CASE(34); Elmsym(i_atm) = 'Se'
        CASE(35); Elmsym(i_atm) = 'Br'
        CASE(36); Elmsym(i_atm) = 'Kr'

!       Fifth row elements               
        CASE(37); Elmsym(i_atm) = 'Rb'
        CASE(38); Elmsym(i_atm) = 'Sr'

        CASE(39); Elmsym(i_atm) = 'Y'
        CASE(40); Elmsym(i_atm) = 'Zr'
        CASE(41); Elmsym(i_atm) = 'Nb'
        CASE(42); Elmsym(i_atm) = 'Mo'
        CASE(43); Elmsym(i_atm) = 'Tc'
        CASE(44); Elmsym(i_atm) = 'Ru'
        CASE(45); Elmsym(i_atm) = 'Rh'
        CASE(46); Elmsym(i_atm) = 'Pd'
        CASE(47); Elmsym(i_atm) = 'Ag'
        CASE(48); Elmsym(i_atm) = 'Cd'

        CASE(49); Elmsym(i_atm) = 'In'
        CASE(50); Elmsym(i_atm) = 'Sn'
        CASE(51); Elmsym(i_atm) = 'Sb'
        CASE(52); Elmsym(i_atm) = 'Te'
        CASE(53); Elmsym(i_atm) = 'I'
        CASE(54); Elmsym(i_atm) = 'Xe'

!       Sixth row elements              
        CASE(55); Elmsym(i_atm) = 'Cs' 
        CASE(56); Elmsym(i_atm) = 'Ba'

        CASE(57); Elmsym(i_atm) = 'La' ! ^
        CASE(58); Elmsym(i_atm) = 'Ce' ! |
        CASE(59); Elmsym(i_atm) = 'Pr' ! |
        CASE(60); Elmsym(i_atm) = 'Nd' ! |
        CASE(61); Elmsym(i_atm) = 'Pm' ! |
        CASE(62); Elmsym(i_atm) = 'Sm' ! |
        CASE(63); Elmsym(i_atm) = 'Eu' ! |  Lanthanoid
        CASE(64); Elmsym(i_atm) = 'Gd' ! |
        CASE(65); Elmsym(i_atm) = 'Tb' ! | 
        CASE(66); Elmsym(i_atm) = 'Dy' ! |
        CASE(67); Elmsym(i_atm) = 'Ho' ! |
        CASE(68); Elmsym(i_atm) = 'Er' ! |
        CASE(69); Elmsym(i_atm) = 'Tm' ! |
        CASE(70); Elmsym(i_atm) = 'Yb' ! |
        CASE(71); Elmsym(i_atm) = 'Lu' ! V 

        CASE(72); Elmsym(i_atm) = 'Hf'
        CASE(73); Elmsym(i_atm) = 'Ta'
        CASE(74); Elmsym(i_atm) = 'W'
        CASE(75); Elmsym(i_atm) = 'Re'
        CASE(76); Elmsym(i_atm) = 'Os'
        CASE(77); Elmsym(i_atm) = 'Ir'
        CASE(78); Elmsym(i_atm) = 'Pt'
        CASE(79); Elmsym(i_atm) = 'Au'
        CASE(80); Elmsym(i_atm) = 'Hg'

        CASE(81); Elmsym(i_atm) = 'Tl'
        CASE(82); Elmsym(i_atm) = 'Pb'
        CASE(83); Elmsym(i_atm) = 'Bi'
        CASE(84); Elmsym(i_atm) = 'Po'
        CASE(85); Elmsym(i_atm) = 'At'
        CASE(86); Elmsym(i_atm) = 'Rn'

!       Seventh row elements              
        CASE(87); Elmsym(i_atm) = 'Fr' 
        CASE(88); Elmsym(i_atm) = 'Ra'

        CASE(89);  Elmsym(i_atm) = 'Ac' ! ^
        CASE(90);  Elmsym(i_atm) = 'Th' ! |
        CASE(91);  Elmsym(i_atm) = 'Pa' ! |
        CASE(92);  Elmsym(i_atm) = 'U'  ! |
        CASE(93);  Elmsym(i_atm) = 'Np' ! |
        CASE(94);  Elmsym(i_atm) = 'Pu' ! |
        CASE(95);  Elmsym(i_atm) = 'Am' ! |  Actinoid
        CASE(96);  Elmsym(i_atm) = 'Cm' ! |
        CASE(97);  Elmsym(i_atm) = 'Bk' ! |
        CASE(98);  Elmsym(i_atm) = 'Cf' ! |
        CASE(99);  Elmsym(i_atm) = 'Es' ! |
        CASE(100); Elmsym(i_atm) = 'Fm' ! |
        CASE(101); Elmsym(i_atm) = 'Md' ! |
        CASE(102); Elmsym(i_atm) = 'No' ! |
        CASE(103); Elmsym(i_atm) = 'Lr' ! V

        CASE(104); Elmsym(i_atm) = 'Rf'
        CASE(105); Elmsym(i_atm) = 'Db'
        CASE(106); Elmsym(i_atm) = 'Sg'
        CASE(107); Elmsym(i_atm) = 'Bh'
        CASE(108); Elmsym(i_atm) = 'Hs'
        CASE(109); Elmsym(i_atm) = 'Mt'
        CASE(110); Elmsym(i_atm) = 'Ds'
        CASE(111); Elmsym(i_atm) = 'Rg'
        CASE(112); Elmsym(i_atm) = 'Cn'

        CASE(113); Elmsym(i_atm) = 'Nh'
        CASE(114); Elmsym(i_atm) = 'Fl'
        CASE(115); Elmsym(i_atm) = 'Mc'
        CASE(116); Elmsym(i_atm) = 'Lv'
        CASE(117); Elmsym(i_atm) = 'Ts'
        CASE(118); Elmsym(i_atm) = 'Og'

!       Writing error messages
        CASE DEFAULT ! Atmnum(i_atm) > 119 or Atmnum(i_atm) < 0
          WRITE(text_i_atm, '(I3)') i_atm 
          WRITE(text_atmnum, '(I3)') Atmnum(i_atm) 
          SELECT CASE(i_atm)
            CASE(1)
              text_error = 'Invalid atomic number '//TRIM(ADJUSTL(text_atmnum))//' for 1st atom.' 
            CASE(2)
              text_error = 'Invalid atomic number '//TRIM(ADJUSTL(text_atmnum))//' for 2nd atom.' 
            CASE(3)
              text_error = 'Invalid atomic number '//TRIM(ADJUSTL(text_atmnum))//' for 3rd atom.' 
            CASE DEFAULT
              text_error = 'Invalid atomic number '//TRIM(ADJUSTL(text_atmnum))//& 
                          &' for '//TRIM(ADJUSTL(text_i_atm))//'th atom.' 
          END SELECT
          WRITE(6,'(1X, A)') 'Invalid description in your ATMNUM file'
          CALL write_messages(-9999, Text_blank, type_program, name_program)
      END SELECT
    ENDDO

    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Reading primitive Gaussian functions (PGFs) from '//TRIM(Fname(Fao))
    OPEN(Fao, FILE = Fname(Fao)) ! Fname(Fao) = 'ATMORB'
      DO i_cgf = 1, N_cgf
!       Reading the atom to which the 'i_cgf'th CGF belongs  
        READ(Fao,*)
        READ(Fao,*) Ao2atm(i_cgf)
!       Determining the xyz coordinate of the 'i_cgf'th CGF center
        xyzao(1:3, i_cgf) = Xyznuc(1:3, Ao2atm(i_cgf)) ! The CGF center is equal to the position of the nucleus.
   
!       Reading the orbital type and the number of PGFs of the 'i_cgf'th CGF
        READ(Fao,*) Orbtyp(i_cgf), N_pgf(i_cgf)
   
        IF(N_pgf(i_cgf) > Max_n_pgf) THEN ! ! N_pgf exceeds the upper limit.
          WRITE(6,'(1X, A)') 'Too many PGFs'
          WRITE(6,'(1X, A, I0)') 'The number of PGFs &
          &must be less than or equal to ', Max_n_pgf
          CALL write_messages(-9999, Text_blank, type_program, name_program)
        ELSE
        ENDIF
   
!       Determining the (l, m, n) set of the 'i_cgf'th CGF
        SELECT CASE(Orbtyp(i_cgf))
!         s orbitals
          CASE('s')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 0
!         p orbitals
          CASE('px')
            Lmn(1, i_cgf) = 1; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 0
          CASE('py')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 1; Lmn(3, i_cgf) = 0
          CASE('pz')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 1
!         d orbitals
          CASE('dxx')
            Lmn(1, i_cgf) = 2; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 0
          CASE('dyy')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 2; Lmn(3, i_cgf) = 0
          CASE('dzz')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 2
          CASE('dxy')
            Lmn(1, i_cgf) = 1; Lmn(2, i_cgf) = 1; Lmn(3, i_cgf) = 0
          CASE('dxz')
            Lmn(1, i_cgf) = 1; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 1
          CASE('dyz')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 1; Lmn(3, i_cgf) = 1
!         f orbitals
          CASE('fxxx')
            Lmn(1, i_cgf) = 3; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 0
          CASE('fyyy')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 3; Lmn(3, i_cgf) = 0
          CASE('fzzz')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 3
          CASE('fxyy')
            Lmn(1, i_cgf) = 1; Lmn(2, i_cgf) = 2; Lmn(3, i_cgf) = 0
          CASE('fxxy')
            Lmn(1, i_cgf) = 2; Lmn(2, i_cgf) = 1; Lmn(3, i_cgf) = 0
          CASE('fxxz')
            Lmn(1, i_cgf) = 2; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 1
          CASE('fxzz')
            Lmn(1, i_cgf) = 1; Lmn(2, i_cgf) = 0; Lmn(3, i_cgf) = 2
          CASE('fyzz')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 1; Lmn(3, i_cgf) = 2
          CASE('fyyz')
            Lmn(1, i_cgf) = 0; Lmn(2, i_cgf) = 2; Lmn(3, i_cgf) = 1
          CASE('fxyz')
            Lmn(1, i_cgf) = 1; Lmn(2, i_cgf) = 1; Lmn(3, i_cgf) = 1

!         Writing error messages
          CASE DEFAULT ! Orbtyp is neither 's', 'p', 'd', nor 'f'.
            WRITE(text_i_cgf, '(I4)') i_cgf 
            SELECT CASE(i_cgf)
              CASE(1)
                text_error = 'Invalid orbital type '//TRIM(ADJUSTL(Orbtyp(i_cgf)))//' for 1st CGF.'
              CASE(2)
                text_error = 'Invalid orbital type '//TRIM(ADJUSTL(Orbtyp(i_cgf)))//' for 2nd CGF.'
              CASE(3)
                text_error = 'Invalid orbital type '//TRIM(ADJUSTL(Orbtyp(i_cgf)))//' for 3rd CGF.'
              CASE DEFAULT
                text_error = 'Invalid orbital type '//TRIM(ADJUSTL(Orbtyp(i_cgf)))//&
                            &' for '//TRIM(ADJUSTL(text_i_cgf))//'th CGF.'
                CALL write_messages(-9999, Text_blank, type_program, name_program)
          END SELECT
          WRITE(6,'(1X, A)') 'Invalid description in your ATMORB file' 
        END SELECT
!       Reading contraction exponents and contraction coefficients of PGFs
        DO i_pgf = 1, N_pgf(i_cgf)
          READ(Fao,*) cntexp(i_pgf, i_cgf), cntcoef(i_pgf, i_cgf)
        ENDDO ! i_pgf = 1, N_pgf(i_cgf)
      ENDDO ! i_cgf = 1, N_cgf
    CLOSE(Fao)


!   Calculating the normalization factors of the PGFs
    DO i_cgf = 1, N_cgf
      SELECT CASE(SUM(Lmn(:, i_cgf)))
        CASE(0) ! s orbital
          DO i_pgf = 1, N_pgf(i_cgf)
            cntnormfac(i_pgf, i_cgf)&
            = (8.0D0*(cntexp(i_pgf, i_cgf)**3)/(Pi**3))**0.25D0
          ENDDO
        CASE(1) ! p orbital
          DO i_pgf = 1, N_pgf(i_cgf)
            cntnormfac(i_pgf, i_cgf)&
            = (128.0D0*(cntexp(i_pgf, i_cgf)**5)/(Pi**3))**0.25D0
          ENDDO
        CASE(2) ! d orbital (6D)
          SELECT CASE(MAXVAL(Lmn(:, i_cgf)))
            CASE(2) ! dxx, dyy, dzz
              DO i_pgf = 1, N_pgf(i_cgf)
                cntnormfac(i_pgf, i_cgf)&
                = (2048.0D0*(cntexp(i_pgf, i_cgf)**7)&
                /(Pi**3))**0.25D0&
                /DSQRT(3.0D0)
              ENDDO
            CASE(1) ! dxy, dxz, dyz
              DO i_pgf = 1, N_pgf(i_cgf)
                cntnormfac(i_pgf, i_cgf)&
                = (2048.0D0*(cntexp(i_pgf, i_cgf)**7)&
                /(Pi**3))**0.25D0
              ENDDO
            CASE DEFAULT
              STOP
          END SELECT
        CASE(3) ! f orbitals (10F) 
          SELECT CASE(MAXVAL(Lmn(:, i_cgf)))
            CASE(3) ! fxxx, fyyy, fzzz
              DO i_pgf = 1, N_pgf(i_cgf)
                cntnormfac(i_pgf, i_cgf)&
                = (32768.0D0*(cntexp(i_pgf, i_cgf)**9)&
                /(Pi**3))**0.25D0&
                /DSQRT(15.0D0)
              ENDDO
            CASE(2) ! fxyy, fxxy, fxxz, fxzz, fyzz, fyyz
              DO i_pgf = 1, N_pgf(i_cgf)
                cntnormfac(i_pgf, i_cgf)&
                = (32768.0D0*(cntexp(i_pgf, i_cgf)**9)&
                /(Pi**3))**0.25D0&
                /DSQRT(3.0D0)
              ENDDO
            CASE(1) ! fxyz
              DO i_pgf = 1, N_pgf(i_cgf)
                cntnormfac(i_pgf, i_cgf)&
                = (32768.0D0*(cntexp(i_pgf, i_cgf)**9)&
                /(Pi**3))**0.25D0
              ENDDO
            CASE DEFAULT
              STOP
          END SELECT
        CASE DEFAULT
          STOP
      END SELECT
    ENDDO  

!   Calculating the overlap integrals between the CGFs
    DO i_cgf = 1, N_cgf
      DO j_cgf = 1, N_cgf
 
        DO i_pgf = 1, N_pgf(i_cgf)
          DO j_pgf = 1, N_pgf(j_cgf)
            s_cgf(i_cgf, j_cgf)&
            = s_cgf(i_cgf, j_cgf)&
            + cntcoef(i_pgf, i_cgf) * cntcoef(j_pgf, j_cgf)&
            * cntnormfac(i_pgf, i_cgf) * cntnormfac(j_pgf, j_cgf)&
            * func_int_pgf_pgf&
             &(xyzao(:, i_cgf), xyzao(:, j_cgf),&
              &cntexp(i_pgf, i_cgf), cntexp(j_pgf, j_cgf),&
              &Lmn(:, i_cgf), Lmn(:, j_cgf))
          ENDDO
        ENDDO

      ENDDO
    ENDDO
  
!   Reading expantion (MO) coefficients for CGFs
    WRITE(6,'(1X, A)') 'Reading expantion coefficients for CGFs (MO coefficients) from'
    WRITE(6,'(1X, A)') TRIM(Fname(Fcoef_soa_ref))//' and '//TRIM(Fname(Fcoef_sob_ref))
  
!   Reading alpha MO coefficients
    OPEN(Fcoef_soa_ref, FILE = Fname(Fcoef_soa_ref))
      DO i_so = 1, n_so, 2
        DO i_cgf = 1, N_cgf
          READ(Fcoef_soa_ref,*) coef_so(i_cgf, i_so)
        ENDDO
      ENDDO
    CLOSE(Fcoef_soa_ref)
!   Reading beta MO coefficients
    OPEN(Fcoef_sob_ref, FILE = Fname(Fcoef_sob_ref))
      DO i_so = 2, n_so, 2
        DO i_cgf = 1, N_cgf
          READ(Fcoef_sob_ref,*) coef_so(i_cgf, i_so)
        ENDDO
      ENDDO
    CLOSE(Fcoef_sob_ref)
    WRITE(6,'(1X, A)') 'Constructing contracted Gaussian functions (CGFs) from the PGFs'

!   Checking orthonormalization condition  
    IF(Check_orthonorm_s_cgf == 'Yes') THEN
      WRITE(6,'(1X, A)') 'Checking the overlap integrals between the CGFs'
      large_overlap = 'No'

!     Calculating overlap integrals between alpha spin orbitals
      DO i_so = 1, n_so, 2
        DO j_so = 1, n_so, 2

          orthonorm = 0.0D0
          DO i_cgf = 1, N_cgf
            DO j_cgf = 1, N_cgf
              orthonorm&
              = orthonorm&
              + coef_so(i_cgf, i_so)&
              * coef_so(j_cgf, j_so)&
              * s_cgf(i_cgf, j_cgf)
            ENDDO
          ENDDO
          IF(ABS(orthonorm - Identity(i_so, j_so)) > Threshold_s_cgf) THEN
            large_overlap = 'Yes' 
          ELSE
          ENDIF

        ENDDO
      ENDDO

!     Calculating overlap integrals between beta spin orbitals
      DO i_so = 2, n_so, 2
        DO j_so = 2, n_so, 2
          
          orthonorm = 0.0D0
          DO i_cgf = 1, N_cgf
            DO j_cgf = 1, N_cgf
              orthonorm&
              = orthonorm&
              + coef_so(i_cgf, i_so)&
              * coef_so(j_cgf, j_so)&
              * s_cgf(i_cgf, j_cgf)
            ENDDO
          ENDDO
          
          IF(ABS(orthonorm - Identity(i_so, j_so)) > Threshold_s_cgf) THEN
            large_overlap = 'Yes' 
          ELSE
          ENDIF
        ENDDO
      ENDDO

!     Writing error messages
      IF(large_overlap == 'Yes') THEN
        WRITE(6,'(1X)') 
        WRITE(6,'(1X, A)') '********************************'
        WRITE(6,'(1X, A)') '*           Warning            *'
        WRITE(6,'(1X, A)') '* Inaccurate Overlap integrals *'
        WRITE(6,'(1X, A)') '********************************'
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ELSE
        WRITE(6,'(1X, A)') 'Done'
      ENDIF
    ELSE
    ENDIF

    DEALLOCATE(normfac, s_cgf)

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE read_data_elec_1

!***************************************************************************************************
  SUBROUTINE read_data_ci_coef

! read_data_ci_coef reads CI coefficeints and electronic excitations from CICOEF.

!***************************************************************************************************

!   ------------------------
!   Declaration of variables
!   ------------------------
   
!   ------------
!   Program name
!   ------------

    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_data_ci_coef'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   ---------------
!   Local variables
!   ---------------

    INTEGER :: i_line, n_line
    CHARACTER(LEN=100), ALLOCATABLE :: text(:)
    INTEGER :: i, i_state, i_elec_config
    INTEGER :: n_state_as_read
    INTEGER, ALLOCATABLE :: n_line_ex_state(:)
    INTEGER, ALLOCATABLE :: n_line_end_ex_state(:)
    INTEGER, ALLOCATABLE :: n_elec_config_as_read(:)
    INTEGER :: fid

!   Spin multuplicity of electronic states
    CHARACTER(LEN=10), ALLOCATABLE :: Mult_ex(:)
    CHARACTER(LEN=10), ALLOCATABLE :: mult_ex_as_read(:)
    CHARACTER(LEN=100) :: mult_text

!   CI coefficeints
    INTEGER, ALLOCATABLE :: tuple_as_read(:,:)
    INTEGER, ALLOCATABLE :: a_ras_as_read(:,:,:), r_ras_as_read(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: ci_coef_as_read(:,:)
    CHARACTER(LEN=100) :: tuple_text
    CHARACTER(LEN=100) :: value_text 
    CHARACTER(LEN=100) :: a_r_ci_coef_text 
    CHARACTER(LEN=100) :: a_text(1:Max_n_tuple)
    CHARACTER(LEN=100) :: r_ci_coef_text(1:Max_n_tuple)
    CHARACTER(LEN=100) :: ci_coef_text 
    INTEGER :: i_symbol 
    INTEGER :: value_int 
    DOUBLE PRECISION :: value_real 
    CHARACTER(LEN=100) arrow_left_text 
    CHARACTER(LEN=100) arrow_right_text 
    INTEGER tuple_temp, tuple_left, tuple_right
    INTEGER a_ras_temp(1:Max_n_tuple), r_ras_temp(1:Max_n_tuple)
    CHARACTER(LEN=21) :: so_occ_name(1:Max_n_tuple), so_unocc_name(1:Max_n_tuple)
    DOUBLE PRECISION :: ci_coef_temp 
    DOUBLE PRECISION :: contribution_temp 
    CHARACTER(2) :: arrow = '->'
    DOUBLE PRECISION sum_ci2_temp 
    DOUBLE PRECISION, ALLOCATABLE :: sum_ci2_as_read(:)
    DOUBLE PRECISION, ALLOCATABLE :: sum_ci2(:)
    DOUBLE PRECISION, ALLOCATABLE :: contribution(:,:), contribution_sort(:,:)

    DOUBLE PRECISION :: contribution_maxval
    INTEGER :: contribution_maxloc
    INTEGER, ALLOCATABLE :: n_contribution(:)
    INTEGER, ALLOCATABLE :: tuple_sort(:,:)
    INTEGER, ALLOCATABLE :: a_ras_sort(:,:,:), r_ras_sort(:,:,:)

    INTEGER :: i_mocube, j_mocube, n_mocube, i_tuple
    INTEGER, ALLOCATABLE :: num_mocube(:), num_mocube_written(:)
    CHARACTER(LEN=21), ALLOCATABLE, SAVE :: so_name_mocube(:)
    CHARACTER(LEN=3) :: write_num_mocube

    INTEGER n_repeat
    CHARACTER(LEN=100) :: text_error 
    CHARACTER(LEN=2) :: text_n_state 

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   Initializing local variables
    mult_text = REPEAT(' ', 10)
    tuple_text = REPEAT(' ', 100)
    value_text = REPEAT(' ', 100)
    a_r_ci_coef_text = REPEAT(' ', 100)
    ci_coef_text =  REPEAT(' ', 100)
    a_text =  REPEAT(' ', 100)
    r_ci_coef_text =  REPEAT(' ', 100)
    i_symbol = 0
    value_int = 0
    value_real = 0.0D0
    arrow_left_text  = REPEAT(' ', 100)
    arrow_right_text = REPEAT(' ', 100)
    tuple_temp = 0; tuple_left = 0; tuple_right = 0
    a_ras_temp = 0; r_ras_temp = 0
    ci_coef_temp = 0.0D0
    n_repeat = 0
   

!   Reading the number of lines in CICOEF 
    n_line = 0
    OPEN(Fcicoef, FILE = Fname(Fcicoef))
      DO
        READ(Fcicoef, '()', END=100)
        n_line = n_line + 1
      ENDDO
100 CLOSE(Fcicoef)


!   Reading the number of electronic states
    ALLOCATE(text(1:n_line)); text = REPEAT(' ', 100)
    n_state_as_read = 0
    OPEN(Fcicoef, FILE = Fname(Fcicoef))
      DO i_line = 1, n_line
        READ(Fcicoef, '(A)') text(i_line)
        IF(INDEX(text(i_line),'$EX_STATE')/=0) THEN
          n_state_as_read = n_state_as_read + 1
        ELSE
        ENDIF
      ENDDO
    CLOSE(Fcicoef)
    IF(n_state_as_read > Max_n_state) n_state_as_read = Max_n_state

    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Reading CI coefficients'
    ALLOCATE(n_line_ex_state(1:n_state_as_read)); n_line_ex_state = 0
    ALLOCATE(n_line_end_ex_state(1:n_state_as_read)); n_line_end_ex_state = 0
    i_state = 0
    DO i_line = 1, n_line
      IF(INDEX(text(i_line),'$EX_STATE')/=0) THEN
        i_state = i_state + 1
        n_line_ex_state(i_state) = i_line
      ELSEIF(INDEX(text(i_line),'$END_EX_STATE')/=0) THEN
        n_line_end_ex_state(i_state) = i_line
      ELSE
      ENDIF
    ENDDO
    WRITE(6,'(1X, A)') 'Done'


!   *******************************************
!   *                                         *
!   *  N_state: the number of excited states  *
!   *                                         *
!   *******************************************
!
!   RAS-2SF
!     N_state = (the number of states in CICOEF) - 1
!     i_state = -1 : HF determinant
!     i_state =  0 : ground state      <- state 1 in CICOEF
!     i_state =  1 : 1st excited state <- state 2 in CICIEF
!     i_state =  n : nth excited state <- state n in CICOEF
!   CIS
!     N_state = the number of states in CICOEF
!     i_state =  0 : ground state (= HF determinant)
!     i_state =  1 : 1st excited state <- state 1 in CICOEF
!     i_state =  n : nth excited state <- state n in CICIEF
    !WRITE(*,*) n_state_as_read
    !STOP                                      
    IF(n_state_as_read > Max_n_state) THEN
      WRITE(6,'(1X, A)') 'Too many excited states in your CICOEF file'
      WRITE(6,'(1X, A, I0)') 'The number of excited states &
     &must be less than or equal to ', Max_n_state
    ELSE
    ENDIF

!   Spin multiplicity of electronic states
    ALLOCATE(mult_ex_as_read(1:n_state_as_read)) 
    mult_ex_as_read = REPEAT(' ', 10)

!   Reading CI coefficients
!   CI coefficients
    ALLOCATE(n_elec_config_as_read(1:n_state_as_read))
    n_elec_config_as_read = 0
    ALLOCATE(tuple_as_read(1:Max_n_elec_config, 1:n_state_as_read))
    tuple_as_read = 0
    ALLOCATE(a_ras_as_read(1:4,1:Max_n_elec_config,1:n_state_as_read))
    a_ras_as_read = 0
    ALLOCATE(r_ras_as_read(1:4,1:Max_n_elec_config,1:n_state_as_read))
    r_ras_as_read = 0
    ALLOCATE(ci_coef_as_read(1:Max_n_elec_config, 1:n_state_as_read))
    ci_coef_as_read = 0.0D0
    ALLOCATE(sum_ci2_as_read(1:n_state_as_read))
    sum_ci2_as_read = 0.0D0
  
    DO i_state = 1, n_state_as_read
      n_elec_config_as_read(i_state)&
      = n_line_end_ex_state(i_state)&
      - n_line_ex_state(i_state) - 3

      IF(n_elec_config_as_read(i_state) > Max_n_elec_config) THEN
        WRITE(6,'(1X, A)') 'Too many electronic configurations'
        WRITE(6,'(1X, A, I0)') 'The number of electronic configurations'
        WRITE(6,'(1X, A, I0)') 'must be less than or equal to ', Max_n_elec_config
        WRITE(6,'(1X)')
        WRITE(6,'(1X, A)') 'Reduce electronic configurations'
        WRITE(6,'(1X, A)') 'in your CICOEF file'
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ELSE
      ENDIF

!     Assigning spin multiplicity of electronic states
      i_line = n_line_ex_state(i_state) + 2
      mult_text = text(i_line)
      Mult_ex_as_read(i_state) = ADJUSTL(TRIM(mult_text))

      sum_ci2_temp = 0.0D0
      DO i_elec_config = 1, n_elec_config_as_read(i_state)
        i_line = n_line_ex_state(i_state) + 2 + i_elec_config

        tuple_text = ADJUSTL(text(i_line))
        i_symbol = INDEX(tuple_text, ' ')
        CALL get_text_left_symbol(tuple_text, i_symbol, value_text)
        CALL get_int(value_text, value_int)
        tuple_as_read(i_elec_config, i_state) = value_int

        CALL get_text_right_symbol&
             (tuple_text, i_symbol, a_r_ci_coef_text)
        SELECT CASE(tuple_as_read(i_elec_config, i_state))
          CASE(0)
            value_text = TRIM(a_r_ci_coef_text)
            CALL get_real(value_text, value_real)
            ci_coef_as_read(i_elec_config, i_state) = value_real
          CASE(1)
            i_symbol = INDEX(a_r_ci_coef_text, '->')
            CALL get_text_left_symbol&
                 (a_r_ci_coef_text, i_symbol, a_text(1))
            CALL get_text_right_symbol&
                 (a_r_ci_coef_text, i_symbol+2, r_ci_coef_text(1))

            ! Read 'a'
            i_symbol = INDEX(a_text(1), ' ')
            CALL get_text_left_symbol&
                 (a_text(1), i_symbol, value_text)
            CALL get_int(value_text, value_int)
            a_ras_as_read(1, i_elec_config, i_state) = value_int

            ! Read 'r'
            i_symbol = INDEX(r_ci_coef_text(1), ' ')
            CALL get_text_left_symbol&
                 (r_ci_coef_text(1), i_symbol, value_text)
            CALL get_int(value_text, value_int)
            r_ras_as_read(1, i_elec_config, i_state) = value_int

            CALL get_text_right_symbol(r_ci_coef_text(1),&
                 i_symbol, ci_coef_text)
            CALL get_real(ci_coef_text, value_real)
            ci_coef_as_read(i_elec_config, i_state) = value_real
          CASE DEFAULT
            i_symbol = INDEX(a_r_ci_coef_text, arrow)
            CALL get_text_left_symbol&
                 (a_r_ci_coef_text, i_symbol, a_text(1))
            CALL get_text_right_symbol&
                 (a_r_ci_coef_text, i_symbol+2, r_ci_coef_text(1))

            DO i = 1, tuple_as_read(i_elec_config, i_state)
              i_symbol = INDEX(a_text(i), ' ')
              ! Read 'a'
              CALL get_text_left_symbol&
                   (a_text(i), i_symbol, value_text)
              CALL get_int(value_text, value_int)
              a_ras_as_read(i, i_elec_config, i_state) = value_int

              ! Read 'r'
              i_symbol = INDEX(r_ci_coef_text(i), ' ')
              CALL get_text_left_symbol&
                   (r_ci_coef_text(i), i_symbol, value_text)
              CALL get_int(value_text, value_int)
              r_ras_as_read(i, i_elec_config, i_state) = value_int

              IF(i == tuple_as_read(i_elec_config, i_state)) EXIT
              CALL get_text_right_symbol&
                   (a_text(i), i_symbol, a_text(i+1))
              CALL get_text_right_symbol&
                   (r_ci_coef_text(i), i_symbol,&
                    r_ci_coef_text(i+1))
            ENDDO

            ! Read CI coefficients
            CALL get_text_right_symbol(r_ci_coef_text&
                 (tuple_as_read(i_elec_config, i_state)),&
                 i_symbol, ci_coef_text)
            CALL get_real(ci_coef_text, value_real)
            ci_coef_as_read(i_elec_config, i_state) = value_real
        END SELECT
        sum_ci2_temp&!(i_state)&
        = sum_ci2_temp&!(i_state)&
        + ci_coef_as_read(i_elec_config, i_state)**2  
      ENDDO ! i_elec_config = 1, n_elec_config(i_state)
      sum_ci2_as_read(i_state) = sum_ci2_temp
    ENDDO
    DEALLOCATE(text, n_line_ex_state, n_line_end_ex_state)

    SELECT CASE(method)
      CASE('cisd')
        WRITE(*,*) 'Norm(A)', DSQRT(sum_ci2_as_read(1))
        ci_coef_as_read(:, 1) &
       &= ci_coef_as_read(:, 1)/DSQRT(sum_ci2_as_read(1))
        sum_ci2_as_read(1) = 0.0D0
        sum_ci2_as_read(1) = SUM(ci_coef_as_read(:, 1)**2.0D0)
      CASE DEFAULT
    END SELECT

!   *******************************************
!   *                                         *
!   *  N_state: the number of excited states  *
!   *                                         *
!   *******************************************
!
!   RAS-2SF
!     N_state = (the number of states in CICOEF) - 1
!     i_state = -1 : HF determinant
!     i_state =  0 : ground state      <- state 1 in CICOEF
!     i_state =  1 : 1st excited state <- state 2 in CICIEF
!     i_state =  n : nth excited state <- state n in CICOEF
!   CIS
!     N_state = the number of states in CICOEF
!     i_state =  0 : ground state (= HF determinant)
!     i_state =  1 : 1st excited state <- state 1 in CICOEF
!     i_state =  n : nth excited state <- state n in CICIEF
                                          

    SELECT CASE(method)
      CASE('Cis')
        N_state = n_state_as_read ! N_state = (the number of states in CICOEF) 
        WRITE(text_n_state, '(I2)') N_state 
        text_error = 'the number of excited states, '//text_n_state
        IF(Bra > N_state) THEN
          WRITE(6,'(1X, A, A)') 'Too large Bra value: ', Bra
          WRITE(6,'(1X, A, I0)') 'Bra value is less than or euqal to', N_state
          CALL write_messages(-9999, Text_blank, type_program, name_program)
        ELSEIF(Ket > N_state) THEN
          WRITE(6,'(1X, A, A)') 'Too large Ket value: ', Ket
          WRITE(6,'(1X, A, I0)') 'Ket value is less than or euqal to ', N_state
          CALL write_messages(-9999, Text_blank, type_program, name_program)
        ELSE
        ENDIF

        ALLOCATE(Mult_ex(0:N_state)); Mult_ex = REPEAT(' ', 10)
        ALLOCATE(N_elec_config(0:N_state)); N_elec_config = 0
        ALLOCATE(Tuple(1:Max_n_elec_config, 0:N_state)); Tuple = 0
        ALLOCATE(A_ras(1:4,1:Max_n_elec_config,0:N_state)); A_ras = 0
        ALLOCATE(R_ras(1:4,1:Max_n_elec_config,0:N_state)); R_ras = 0
        ALLOCATE(Ci_coef(1:Max_n_elec_config, 0:N_state)); Ci_coef = 0.0D0
        ALLOCATE(sum_ci2(0:n_state)); sum_ci2 = 0.0D0
        Mult_ex(1:N_state) = mult_ex_as_read(1:n_state_as_read)
        N_elec_config(1:N_state) = n_elec_config_as_read(1:n_state_as_read)
        Tuple(:, 1:N_state) = tuple_as_read(:, 1:n_state_as_read)
        A_ras(:, :, 1:N_state) = a_ras_as_read(:, :, 1:n_state_as_read)
        R_ras(:, :, 1:N_state) = r_ras_as_read(:, :, 1:n_state_as_read)
        Ci_coef(:, 1:N_state) = ci_coef_as_read(:, 1:n_state_as_read)
        sum_ci2(1:N_state) = sum_ci2_as_read(1:n_state_as_read)
      CASE DEFAULT
    END SELECT  
    DEALLOCATE(mult_ex_as_read)
    DEALLOCATE(n_elec_config_as_read)
    DEALLOCATE(a_ras_as_read)
    DEALLOCATE(r_ras_as_read)
    DEALLOCATE(ci_coef_as_read)
    DEALLOCATE(sum_ci2_as_read)
!   End of reading CI coefficients

!   Calculating contributions from individual electronic configurations
!   to electronic states
    SELECT CASE(method)
      CASE('Cis')
        ALLOCATE(contribution(1:Max_n_elec_config, 0:N_state)); contribution = 0.0D0
        ALLOCATE(contribution_sort(1:Max_n_elec_config, 0:N_state)); contribution_sort = 0.0D0
        ALLOCATE(n_contribution(0:N_state)); n_contribution = 0
        ALLOCATE(tuple_sort(1:Max_n_elec_config, 0:N_state)); tuple_sort = 0 
        ALLOCATE(a_ras_sort(1:1,1:Max_n_elec_config,0:N_state)); a_ras_sort = 0
        ALLOCATE(r_ras_sort(1:1,1:Max_n_elec_config,0:N_state)); r_ras_sort = 0
        DO i_state = 1, N_state
          SELECT CASE(Mult_ex(i_state))
            CASE('Singlet', 'Unknown')
              DO i = 1, n_elec_config(i_state)/2
                contribution(i, i_state) = 200.0D0*Ci_coef(2*i-1, i_state)**2/sum_ci2(i_state)
              ENDDO
            CASE('Triplet')
              IF(Ms == 0) THEN
                DO i = 1, n_elec_config(i_state)/2
                  contribution(i, i_state) = 200.0D0*Ci_coef(2*i-1, i_state)**2/sum_ci2(i_state)
                ENDDO
              ELSEIF(Ms ==  1 .OR. Ms == -1) THEN
                DO i = 1, n_elec_config(i_state)
                  contribution(i, i_state) = 100.0D0*Ci_coef(i, i_state)**2/sum_ci2(i_state)
                ENDDO
              ELSE
              ENDIF
            CASE DEFAULT
          END SELECT
        ENDDO
    
        DO i_state = 1, N_state
          SELECT CASE(Mult_ex(i_state))
            CASE('Singlet', 'Unknown')
              n_contribution(i_state) = n_elec_config(i_state)/2
              DO i = 1, n_elec_config(i_state)/2
                contribution_maxval = MAXVAL(contribution(:, i_state))
                IF(ABS(contribution_maxval) < Threshold_contribution) THEN
                  n_contribution(i_state) = i - 1
                  EXIT
                ELSE
                ENDIF    
                contribution_maxloc = MAXLOC(contribution(:, i_state), 1)
                contribution_sort(i, i_state) = contribution_maxval
                a_ras_sort(1, i, i_state) = A_ras(1, 2*contribution_maxloc-1, i_state)
                r_ras_sort(1, i, i_state) = R_ras(1, 2*contribution_maxloc-1, i_state)
                contribution(contribution_maxloc, i_state) = 0.0D0
              ENDDO
            CASE('Triplet')
              IF(Ms == 0) THEN
                n_contribution(i_state) = n_elec_config(i_state)/2
                DO i = 1, n_elec_config(i_state)/2
                  contribution_maxval = MAXVAL(contribution(:, i_state))
                  IF(ABS(contribution_maxval) < Threshold_contribution) THEN
                    n_contribution(i_state) = i - 1
                    EXIT
                  ELSE
                  ENDIF
                  contribution_maxloc = MAXLOC(contribution(:, i_state), 1)
                  contribution_sort(i, i_state) = contribution_maxval
                  a_ras_sort(1, i, i_state) = A_ras(1, 2*contribution_maxloc-1, i_state)
                  r_ras_sort(1, i, i_state) = R_ras(1, 2*contribution_maxloc-1, i_state)
                  contribution(contribution_maxloc, i_state) = 0.0D0
                ENDDO
              ELSEIF(Ms ==  1 .OR. Ms == -1) THEN
                n_contribution(i_state) = n_elec_config(i_state)
                DO i = 1, n_elec_config(i_state)
                  contribution_maxval = MAXVAL(contribution(:, i_state))
                  IF(ABS(contribution_maxval) < Threshold_contribution) THEN
                    n_contribution(i_state) = i - 1
                    EXIT
                  ELSE
                  ENDIF
                  contribution_maxloc = MAXLOC(contribution(:, i_state), 1)
                  contribution_sort(i, i_state) = contribution_maxval
                  a_ras_sort(1, i, i_state) = A_ras(1, contribution_maxloc, i_state)
                  r_ras_sort(1, i, i_state) = R_ras(1, contribution_maxloc, i_state)
                  contribution(contribution_maxloc, i_state) = 0.0D0
                ENDDO
              ELSE
              ENDIF
            CASE DEFAULT
          END SELECT
        ENDDO
      CASE DEFAULT
    END SELECT  


!   Writing CI coefficients in terms of so numbers    
    WRITE(6,'(1X)')
    n_repeat = 45
    WRITE(6,'(1X, A)') 'CI coefficients and electronic configurations'
    WRITE(6,'(1X, A)') 'in terms of spin-orbital numbers'
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    SELECT CASE(method)
      CASE('Cis')
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
        DO i_state = 1, n_state
          WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
          WRITE(6,'(1X, A, I0, A)') '| ', i_state, ' >'
          WRITE(6,'(1X, A, I3)') 'CIS excited state ', i_state
          WRITE(6,'(1X, A)') mult_ex(i_state)
          WRITE(6,'(1X)')
          DO i_elec_config = 1, n_elec_config(i_state)
            tuple_temp = Tuple(i_elec_config, i_state)
            a_ras_temp = A_ras(:, i_elec_config, i_state)
            r_ras_temp = R_ras(:, i_elec_config, i_state)
            ci_coef_temp = Ci_coef(i_elec_config, i_state)
            SELECT CASE(tuple_temp)
              CASE(0)
                WRITE(6,'(1X, F11.7, 1X, A)') &
               &ci_coef_temp, 'HF determinant'
              CASE(1)
                WRITE(6,'(1X, F11.7, F13.7, I4, A3, I4)') &
               &ci_coef_temp, &
               &(a_ras_temp(i), i = 1, tuple_temp), arrow, &
               &(r_ras_temp(i), i = 1, tuple_temp)
              CASE DEFAULT
            END SELECT
          ENDDO
        ENDDO
        WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
        WRITE(6,'(1X)')
      CASE DEFAULT
    END SELECT ! method
!   End of writing CI coefficients in terms of so numbers    

!   Writing CI coefficients in terms of so names
    WRITE(6,'(1X)')
    n_repeat = 45
    WRITE(6,'(1X, A)') 'CI coefficients and electronic configurations'
    WRITE(6,'(1X, A)') 'in terms of spin-orbital names'
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    SELECT CASE(method)
      CASE('Cis')
        IF(Property == 'Check' .OR. Property == 'Cicoef_only') THEN
          WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
          WRITE(6,'(1X, A)') 'HF ground state'
          WRITE(6,'(1X, A)') 'Singlet'
          WRITE(6,'(1X)')
          WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
          DO i_state = 1, n_state
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', i_state, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', i_state
            WRITE(6,'(1X, A)') mult_ex(i_state)
            WRITE(6,'(1X)')
            DO i_elec_config = 1, n_elec_config(i_state)
              so_occ_name(1) = So_name(A_ras(1, i_elec_config, i_state))
              so_unocc_name(1) = So_name(R_ras(1, i_elec_config, i_state))
              WRITE(6,'(1X, F11.7, 1X, A9, 1X, A9)') ci_coef(i_elec_config, i_state), &
                           &so_unocc_name(1)
              WRITE(6,'(13X,           A9, 1X, A9)') so_occ_name(1) 
            ENDDO
          ENDDO
        ELSE ! Property /= 'Check'
      
          IF(Bra == 0) THEN
            WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
            WRITE(6,'(1X, A)') 'HF ground state'
            WRITE(6,'(1X, A)') 'Singlet'
            WRITE(6,'(1X)')
            WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
          ELSE ! Bra /= 0
            WRITE(6,'(1X, A, I0, A)') '| ', Bra, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', Bra
            WRITE(6,'(1X, A)') mult_ex(Bra)
            WRITE(6,'(1X)')
            DO i_elec_config = 1, n_elec_config(Bra)
              so_occ_name(1) = So_name(A_ras(1, i_elec_config, Bra))
              so_unocc_name(1) = So_name(R_ras(1, i_elec_config, Bra))
              WRITE(6,'(1X, F11.7, 1X, A9, 1X, A9)') ci_coef(i_elec_config, Bra), &
                           &so_unocc_name(1)
            ENDDO               
            WRITE(6,'(13X,           A9, 1X, A9)') so_occ_name(1)
          ENDIF
      
          IF(Ket == Bra) THEN
          ELSEIF(Ket == 0) THEN
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
            WRITE(6,'(1X, A)') 'HF ground state'
            WRITE(6,'(1X, A)') 'Singlet'
            WRITE(6,'(1X)')
            WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
          ELSE ! Ket /= 0
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', Ket, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', Ket
            WRITE(6,'(1X, A)') mult_ex(Ket)
            WRITE(6,'(1X)')
            DO i_elec_config = 1, n_elec_config(Ket)
              so_occ_name(1) = So_name(A_ras(1, i_elec_config, Ket))
              so_unocc_name(1) = So_name(R_ras(1, i_elec_config, Ket))
              WRITE(6,'(1X, F11.7, 1X, A9, 1X, A9)') ci_coef(i_elec_config, Ket), &
                           &so_unocc_name(1)
              WRITE(6,'(13X,           A9, 1X, A9)') so_occ_name(1)
            ENDDO  
          ENDIF
      
        ENDIF
        WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
        WRITE(6,'(1X)')

      CASE DEFAULT
    END SELECT ! method
!   End of writing CI coefficients in terms of so numbers


    WRITE(6,'(1X, A)') 'Sum of squared CI coefficients of'
    WRITE(6,'(1X, A)') 'electronic state |i>'
    n_repeat = 16
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    WRITE(6,'(1X, A)') '  i   Sum[ci**2]'
    WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
    SELECT CASE(method)
      CASE('Cis')
        IF(Property == 'Check' .OR. Property == 'Cicoef_only') THEN
          WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
          WRITE(6,'(1X, A)') 'HF ground state'
          WRITE(6,'(1X, A)') 'Singlet'
          WRITE(6,'(1X)')
          WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
    
          DO i_state = 1, N_state
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', i_state, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', i_state
            WRITE(6,'(1X, A)') mult_ex(i_state)
            WRITE(6,'(1X)')
            DO i_elec_config = 1, n_contribution(i_state)
              tuple_temp = tuple_sort(i_elec_config, i_state)
              a_ras_temp(1) = a_ras_sort(1, i_elec_config, i_state)
              r_ras_temp(1) = r_ras_sort(1, i_elec_config, i_state)
              contribution_temp = contribution_sort(i_elec_config, i_state)
              tuple_temp = 1
              WRITE(6,'(1X, F11.7, F13.7, I4, A3, I4)')&
              contribution_temp,&
              (a_ras_temp(i), i = 1, tuple_temp), arrow,&
              (r_ras_temp(i), i = 1, tuple_temp)
            ENDDO
          ENDDO
        ELSE ! Property /= 'Check'
          IF(Bra == 0) THEN
            WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
            WRITE(6,'(1X, A)') 'HF ground state'
            WRITE(6,'(1X, A)') 'Singlet'
            WRITE(6,'(1X)')
            WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
          ELSE ! Bra /= 0
            WRITE(6,'(1X, A, I0, A)') '| ', Bra, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', Bra
            WRITE(6,'(1X, A)') mult_ex(Bra)
            WRITE(6,'(1X)')
            DO i_elec_config = 1, n_contribution(Bra)
              tuple_temp = tuple_sort(i_elec_config, Bra)
              a_ras_temp = a_ras_sort(:, i_elec_config, Bra)
              r_ras_temp = r_ras_sort(:, i_elec_config, Bra)
              contribution_temp = contribution_sort(i_elec_config, Bra)
              tuple_temp = 1
              WRITE(6,'(1X, F11.7, F13.7, I4, A3, I4)')&
              contribution_temp,&
              (a_ras_temp(i), i = 1, tuple_temp), arrow,&
              (r_ras_temp(i), i = 1, tuple_temp)
            ENDDO
          ENDIF
    
          IF(Ket == Bra) THEN
          ELSEIF(Ket == 0) THEN
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
            WRITE(6,'(1X, A)') 'HF ground state'
            WRITE(6,'(1X, A)') 'Singlet'
            WRITE(6,'(1X)')
            WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
          ELSE ! Ket /= 0
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', Ket, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', Ket
            WRITE(6,'(1X, A)') mult_ex(Ket)
            WRITE(6,'(1X)')
            DO i_elec_config = 1, n_contribution(Ket)
              tuple_temp = tuple_sort(i_elec_config, Ket)
              a_ras_temp = a_ras_sort(:, i_elec_config, Ket)
              r_ras_temp = r_ras_sort(:, i_elec_config, Ket)
              contribution_temp = contribution_sort(i_elec_config, Ket)
              tuple_temp = 1
              WRITE(6,'(1X, F11.7, F13.7, I4, A3, I4)')&
              contribution_temp,&
              (a_ras_temp(i), i = 1, tuple_temp), arrow,&
              (r_ras_temp(i), i = 1, tuple_temp)
            ENDDO
          ENDIF
        ENDIF
        WRITE(6,'(1X, A)') REPEAT('=', n_repeat)

      CASE DEFAULT
    END SELECT

    WRITE(6,'(1X)')
    n_repeat = 46
    WRITE(6,'(1X, A)') 'Contributions from individual electronic configurations'
    WRITE(6,'(1X, A)') 'to electronic states in terms of spin-orbital names'
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    SELECT CASE(method)
      CASE('Cis')

        IF(Property == 'Check' .OR. Property == 'Cicoef_only') THEN
          WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
          WRITE(6,'(1X, A)') 'HF ground state'
          WRITE(6,'(1X, A)') 'Singlet'
          WRITE(6,'(1X)')
          WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'

          DO i_state = 1, N_state
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', i_state, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', i_state
            WRITE(6,'(1X, A)') mult_ex(i_state)
            WRITE(6,'(1X)')
            DO i = 1, n_contribution(i_state)
              WRITE(6, '(1X, F6.2, 2X, A9)') contribution_sort(i, i_state), &
                           &So_name(r_ras_sort(1, i, i_state))
              WRITE(6, '(9X, A9)') So_name(a_ras_sort(1, i, i_state))
            ENDDO
          ENDDO
        ELSE
          IF(Bra == 0) THEN
            WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
            WRITE(6,'(1X, A)') 'HF ground state'
            WRITE(6,'(1X, A)') 'Singlet'
            WRITE(6,'(1X)')
            WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
          ELSE ! Bra /= 0
            WRITE(6,'(1X, A, I0, A)') '| ', Bra, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', Bra
            WRITE(6,'(1X, A)') mult_ex(Bra)
            WRITE(6,'(1X)')
            DO i = 1, n_contribution(Bra)
              WRITE(6, '(1X, F6.2, 2X, A9)') contribution_sort(i, Bra), &
                           &So_name(r_ras_sort(1, i, Bra))
              WRITE(6, '(9X, A9)') So_name(a_ras_sort(1, i, Bra))
            ENDDO
          ENDIF

          IF(Ket == Bra) THEN
          ELSEIF(Ket == 0) THEN
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
            WRITE(6,'(1X, A)') 'HF ground state'
            WRITE(6,'(1X, A)') 'Singlet'
            WRITE(6,'(1X)')
            WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
          ELSE ! Ket /= 0
            WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
            WRITE(6,'(1X, A, I0, A)') '| ', Ket, ' >'
            WRITE(6,'(1X, A, I3)') 'CIS excited state ', Ket
            WRITE(6,'(1X, A)') mult_ex(Ket)
            WRITE(6,'(1X)')
            DO i = 1, n_contribution(Ket)
              WRITE(6, '(1X, F6.2, 2X, A9)') contribution_sort(i, Ket), &
                           &So_name(r_ras_sort(1, i, Ket))
              WRITE(6, '(9X, A9)') So_name(a_ras_sort(1, i, Ket))
            ENDDO
          ENDIF
        ENDIF
        WRITE(6,'(1X, A)') REPEAT('=', n_repeat)

      CASE DEFAULT
    END SELECT

!   2021/01/28
    IF(Property == 'Cicoef_only') THEN
      n_mocube = 0
      DO i_state = 1, N_state
        DO i = 1, n_contribution(i_state)
          IF(contribution_sort(i, i_state) >= 1.0D0) THEN
            n_mocube = n_mocube + 1
          ELSE
          ENDIF
        ENDDO
      ENDDO
      n_mocube = n_mocube * 2 * Max_n_tuple

      i_mocube = 0
      ALLOCATE(num_mocube(1:n_mocube)); num_mocube = 0
      ALLOCATE(so_name_mocube(1:n_mocube)); so_name_mocube = REPEAT(' ', 21)

      DO i_state = 1, N_state
        DO i_elec_config = 1, n_contribution(i_state)
          IF(contribution_sort(i_elec_config, i_state) >= 1.0D0) THEN
            DO i_tuple = 1, tuple_sort(i_elec_config, i_state)
              i_mocube = i_mocube + 1
              IF(MOD(a_ras_sort(i_tuple, i_elec_config, i_state), 2) == 1) THEN
                num_mocube(i_mocube) = (a_ras_sort(i_tuple, i_elec_config, i_state)+1)/2
              ELSEIF(MOD(a_ras_sort(i_tuple, i_elec_config, i_state), 2) == 0) THEN
                num_mocube(i_mocube) = a_ras_sort(i_tuple, i_elec_config, i_state)/2
              ELSE
              ENDIF
              so_name_mocube(i_mocube) = So_name(a_ras_sort(i_tuple, i_elec_config, i_state))
              
              i_mocube = i_mocube + 1
              IF(MOD(r_ras_sort(i_tuple, i_elec_config, i_state), 2) == 1) THEN
                num_mocube(i_mocube) = (r_ras_sort(i_tuple, i_elec_config, i_state)+1)/2
              ELSEIF(MOD(r_ras_sort(i_tuple, i_elec_config, i_state), 2) == 0) THEN
                num_mocube(i_mocube) = r_ras_sort(i_tuple, i_elec_config, i_state)/2
              ELSE
              ENDIF
              so_name_mocube(i_mocube) = So_name(r_ras_sort(i_tuple, i_elec_config, i_state))
            ENDDO 
          ELSE
          ENDIF
        ENDDO
      ENDDO


!     Generate sh file ! 2021/01/28
      ALLOCATE(num_mocube_written(1:n_mocube)); num_mocube_written = 0
      fid = Fmocube
      OPEN(fid, FILE = Fname(fid))

      WRITE(fid, '(A)') 'TD_FCHK=../TD_PBE0.fchk' 
      WRITE(fid, '(A)') 'if [ -e $TD_FCHK ] ; then' 

      i_mocube = 1
      WRITE(fid, '(A,I0,A,I0,A)') &
     &'  cubegen 0 MO=', num_mocube(1), &
     &' '//TRIM(Fchk_mocube)//' ', num_mocube(1), '_'//TRIM(So_name_mocube(1))//'.cube 0 h'

      num_mocube_written(1) = num_mocube(1)
      DO i_mocube = 2, n_mocube
        write_num_mocube = 'Yes'
        DO j_mocube = 1, i_mocube - 1
          IF(num_mocube(i_mocube) == num_mocube_written(j_mocube)) THEN
            write_num_mocube = 'No'
          ELSE
          ENDIF
        ENDDO
        IF(write_num_mocube == 'Yes') THEN
          WRITE(fid, '(A,I0,A,I0,A)') &
         &'  cubegen 0 MO=', num_mocube(i_mocube), &
         &' '//TRIM(Fchk_mocube)//' ', num_mocube(i_mocube), '_'//TRIM(So_name_mocube(i_mocube))//'.cube 0 h'
          num_mocube_written(i_mocube) = num_mocube(i_mocube)
        ELSE
        ENDIF
      ENDDO
      DEALLOCATE(num_mocube, num_mocube_written)
      DEALLOCATE(So_name_mocube)
    
      WRITE(fid, '(A)') 'else' 
      WRITE(fid, '(A)') "  echo 'Missing '$TD_FCHK" 
      WRITE(fid, '(A)') '  exit' 
      WRITE(fid, '(A)') 'fi' 
    ELSE
    ENDIF ! IF(Property == 'Cicoef_only') THEN
      
    CLOSE(fid)

    DEALLOCATE(contribution, contribution_sort, n_contribution)
    DEALLOCATE(tuple_sort)
    DEALLOCATE(a_ras_sort, r_ras_sort)
    DEALLOCATE(sum_ci2)

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE read_data_ci_coef


!***************************************************************************************************
  SUBROUTINE read_data_ci_coef_td

! read_data_ci_coef_td reads CI coefficeints and electronic excitations from XPLUS_Y and X_MINUS_Y.

!***************************************************************************************************

!   ------------------------
!   Declaration of variables
!   ------------------------

!   ------------
!   Program name
!   ------------

    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_data_ci_coef_td'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   ---------------
!   Local variables
!   ---------------

    INTEGER :: i_line, n_line
    CHARACTER(LEN=100), ALLOCATABLE :: text(:)
    INTEGER :: i_state, i_elec_config, i
    INTEGER, ALLOCATABLE :: n_line_ex_state(:)
    INTEGER, ALLOCATABLE :: n_line_end_ex_state(:)

!   Spin multuplicity of electronic states
    CHARACTER(LEN=10), ALLOCATABLE :: Mult_ex(:)
    CHARACTER(LEN=100) :: mult_text

!   CI coefficeints
    INTEGER :: i_symbol
    CHARACTER(2) :: arrow 
    DOUBLE PRECISION, ALLOCATABLE :: sum_ci2(:)
    INTEGER n_repeat
    INTEGER :: fid

    INTEGER, ALLOCATABLE :: a_ras_as_read(:,:), r_ras_as_read(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: X_plus_y_as_read(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: X_minus_y_as_read(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: contribution(:,:), contribution_sort(:,:)
    DOUBLE PRECISION :: contribution_maxval
    INTEGER :: contribution_maxloc
    INTEGER, ALLOCATABLE :: n_contribution(:)
    INTEGER, ALLOCATABLE :: n_elec_config_as_read(:)
    CHARACTER(LEN=21) :: so_occ_name(1:1), so_unocc_name(1:1)
    INTEGER, ALLOCATABLE :: a_ras_sort(:,:,:), r_ras_sort(:,:,:)

    INTEGER :: i_mocube, j_mocube, n_mocube
    INTEGER, ALLOCATABLE :: num_mocube(:), num_mocube_written(:)
    CHARACTER(LEN=21), ALLOCATABLE, SAVE :: so_name_mocube(:)
    CHARACTER(LEN=3) :: write_num_mocube

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   Initializing local variables
    mult_text = REPEAT(' ', 10)
    i_symbol = 0
    n_repeat = 0
    contribution_maxval = 0.0D0; contribution_maxloc = 0


!   Reading the number of lines in CICOEF
    n_line = 0
    fid = Fxy
    OPEN(fid, FILE = Fname(fid))
      DO
        READ(fid, '()', END=100)
        n_line = n_line + 1
      ENDDO
100 CLOSE(fid)


!   Reading the number of electronic states
    ALLOCATE(text(1:n_line)); text = REPEAT(' ', 100)
    OPEN(fid, FILE = Fname(fid))
      DO i_line = 1, n_line
        READ(fid, '(A)') text(i_line)
        IF(INDEX(text(i_line),'$EX_STATE')/=0) THEN
          N_state = N_state + 1
        ELSE
        ENDIF
      ENDDO
    CLOSE(fid)
    IF(N_state > Max_n_state) N_state = Max_n_state

!   *******************************************
!   *                                         *
!   *  N_state: the number of excited states  *
!   *                                         *
!   *******************************************
!
!   TD-DFT
!     N_state = the number of states in CICOEF
!     i_state =  0 : ground state (= HF determinant)
!     i_state =  1 : 1st excited state <- state 1 in CICOEF
!     i_state =  n : nth excited state <- state n in CICOEF


    ALLOCATE(n_line_ex_state(1:N_state)); n_line_ex_state = 0
    ALLOCATE(n_line_end_ex_state(1:N_state)); n_line_end_ex_state = 0
    i_state = 0
    DO i_line = 1, n_line
      IF(INDEX(text(i_line),'$EX_STATE')/=0) THEN
        i_state = i_state + 1
        n_line_ex_state(i_state) = i_line
      ELSEIF(INDEX(text(i_line),'$END_EX_STATE')/=0) THEN
        n_line_end_ex_state(i_state) = i_line
      ELSE
      ENDIF
      IF(i_state > N_state) EXIT
    ENDDO

    ALLOCATE(Mult_ex(0:N_state)); Mult_ex = REPEAT(' ', 10)
    ALLOCATE(n_elec_config_as_read(0:N_state)); n_elec_config_as_read = 0
    ALLOCATE(a_ras_as_read(1:Max_n_elec_config,0:N_state)); a_ras_as_read = 0
    ALLOCATE(r_ras_as_read(1:Max_n_elec_config,0:N_state)); r_ras_as_read = 0
    ALLOCATE(x_plus_y_as_read(1:Max_n_elec_config, 0:N_state)); X_plus_y_as_read = 0.0D0
    ALLOCATE(x_minus_y_as_read(1:Max_n_elec_config, 0:N_state)); X_minus_y_as_read = 0.0D0

    ALLOCATE(N_elec_config(0:N_state)); N_elec_config = 0
    ALLOCATE(A_ras(1:1,1:Max_n_elec_config,0:N_state)); A_ras = 0
    ALLOCATE(R_ras(1:1,1:Max_n_elec_config,0:N_state)); R_ras = 0
    ALLOCATE(X_plus_y(1:Max_n_elec_config, 0:N_state)); X_plus_y = 0.0D0
    ALLOCATE(X_minus_y(1:Max_n_elec_config, 0:N_state)); X_minus_y = 0.0D0
    ALLOCATE(sum_ci2(0:N_state)); sum_ci2 = 0.0D0

    DO i_state = 1, N_state
      n_elec_config_as_read(i_state)&
      = n_line_end_ex_state(i_state)&
      - n_line_ex_state(i_state) - 3

      IF(n_elec_config_as_read(i_state) > Max_n_elec_config) THEN
          WRITE(6,'(1X, A)') 'Too many electronic configurations'
          WRITE(6,'(1X, A, I0)') 'The number of electronic configurations'
          WRITE(6,'(1X, A, I0)') 'must be less than or equal to ', Max_n_elec_config
          WRITE(6,'(1X)')
          WRITE(6,'(1X, A)') 'Reduce electronic configurations'
          WRITE(6,'(1X, A)') 'in your CICOEF file'
      ELSE
      ENDIF
    ENDDO

    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Reading X+Y and X-Y from '//TRIM(Fname(fid))

!   fid = Fxy
    OPEN(fid, FILE = Fname(fid))
      DO i_state = 1, N_state
        READ(fid, *) ! Skip reading '$EX_STATE' 
        READ(fid, *) ! Skip reading 'STATE = i_state' 
        READ(fid, *) Mult_ex(i_state) 
        DO i_elec_config = 1, n_elec_config_as_read(i_state)
          READ(fid,*) a_ras_as_read(i_elec_config, i_state), arrow, &
                     &r_ras_as_read(i_elec_config, i_state), &
                     &x_plus_y_as_read(i_elec_config, i_state), &
                     &x_minus_y_as_read(i_elec_config, i_state)
        ENDDO
        READ(fid, *) ! Skip reading '$END_EX_STATE' 
        READ(fid, *) ! Skip reading blank line 
      ENDDO
    CLOSE(fid)
    DEALLOCATE(text)

    DO i_state = 1, N_state

      SELECT CASE(Mult_ex(i_state))
        CASE('Singlet', 'Unknown')
          n_elec_config(i_state) = n_elec_config_as_read(i_state)
          A_ras(1,1:n_elec_config(i_state), i_state) &
         &= a_ras_as_read(1:n_elec_config(i_state), i_state)
          R_ras(1,1:n_elec_config(i_state), i_state) &
         &= r_ras_as_read(1:n_elec_config(i_state), i_state)
          X_plus_y(1:n_elec_config(i_state), i_state) &
         &= x_plus_y_as_read(1:n_elec_config(i_state), i_state)
          X_minus_y(1:n_elec_config(i_state), i_state) &
         &= x_minus_y_as_read(1:n_elec_config(i_state), i_state)

        CASE('Triplet')
          IF(Ms == 0) THEN
            n_elec_config(i_state) = n_elec_config_as_read(i_state)
            A_ras(1,1:n_elec_config(i_state), i_state) &
           &= a_ras_as_read(1:n_elec_config(i_state), i_state)
            R_ras(1,1:n_elec_config(i_state), i_state) &
           &= r_ras_as_read(1:n_elec_config(i_state), i_state)
            X_plus_y(1:n_elec_config(i_state), i_state) &
           &= x_plus_y_as_read(1:n_elec_config(i_state), i_state)
            X_minus_y(1:n_elec_config(i_state), i_state) &
           &= x_minus_y_as_read(1:n_elec_config(i_state), i_state)
          ELSEIF(Ms ==  1) THEN
            n_elec_config(i_state) = n_elec_config_as_read(i_state)/2
            DO i_elec_config = 1, n_elec_config(i_state)
              A_ras(1,i_elec_config, i_state) &
             &= a_ras_as_read(2*i_elec_config-1, i_state) + 1
              R_ras(1,i_elec_config, i_state) &
             &= r_ras_as_read(2*i_elec_config-1, i_state) 
              X_plus_y(i_elec_config, i_state) &
             &= Sqrt2*x_plus_y_as_read(2*i_elec_config-1, i_state)
              X_minus_y(i_elec_config, i_state) &
             &= Sqrt2*x_minus_y_as_read(2*i_elec_config-1, i_state)
            ENDDO  

          ELSEIF(Ms == -1) THEN
            n_elec_config(i_state) = n_elec_config_as_read(i_state)/2
            DO i_elec_config = 1, n_elec_config(i_state)
              A_ras(1,i_elec_config, i_state) &
             &= a_ras_as_read(2*i_elec_config-1, i_state) 
              R_ras(1,i_elec_config, i_state) &
             &= r_ras_as_read(2*i_elec_config-1, i_state) + 1
              X_plus_y(i_elec_config, i_state) &
             &= Sqrt2*x_plus_y_as_read(2*i_elec_config-1, i_state)
              X_minus_y(i_elec_config, i_state) &
             &= Sqrt2*x_minus_y_as_read(2*i_elec_config-1, i_state)
            ENDDO

          ELSE
            STOP
          ENDIF
        CASE DEFAULT
      END SELECT

      DO i_elec_config = 1, n_elec_config(i_state)
        sum_ci2(i_state)&
        = sum_ci2(i_state)&
       &+ x_plus_y(i_elec_config, i_state) &
       &* x_minus_y(i_elec_config, i_state)
      ENDDO

    ENDDO
    DEALLOCATE(n_elec_config_as_read)
    DEALLOCATE(a_ras_as_read, r_ras_as_read)
    DEALLOCATE(x_plus_y_as_read, x_minus_y_as_read)
    WRITE(6,'(1X, A)') 'Done'


!   Calculating contributions from individual electronic configurations
!   to electronic states
    ALLOCATE(contribution(1:Max_n_elec_config, 0:N_state)); contribution = 0.0D0
    ALLOCATE(contribution_sort(1:Max_n_elec_config/2, 0:N_state)); contribution_sort = 0.0D0
    ALLOCATE(n_contribution(0:N_state)); n_contribution = 0
    ALLOCATE(a_ras_sort(1:1,1:Max_n_elec_config/2,0:N_state)); a_ras_sort = 0
    ALLOCATE(r_ras_sort(1:1,1:Max_n_elec_config/2,0:N_state)); r_ras_sort = 0
    DO i_state = 1, N_state
      SELECT CASE(Mult_ex(i_state))
        CASE('Singlet', 'Unknown')
          DO i = 1, n_elec_config(i_state)/2
            contribution(i, i_state) &
           &= 200.0D0 &
           &* X_plus_y(2*i-1, i_state) &
           &* X_minus_y(2*i-1, i_state) &
           &/ sum_ci2(i_state)
          ENDDO
        CASE('Triplet')
          IF(Ms == 0) THEN
            DO i = 1, n_elec_config(i_state)/2
              contribution(i, i_state) &
             &= 200.0D0 &
             &* X_plus_y(2*i-1, i_state) &
             &* X_minus_y(2*i-1, i_state) &
             &/ sum_ci2(i_state)
            ENDDO
          ELSEIF(Ms ==  1 .OR. Ms == -1) THEN
            DO i = 1, n_elec_config(i_state)
              contribution(i, i_state) &
             &= 100.0D0 &
             &* X_plus_y(i, i_state) &
             &* X_minus_y(i, i_state) &
             &/ sum_ci2(i_state)
            ENDDO
          ELSE
          ENDIF
        CASE DEFAULT
      END SELECT
    ENDDO

    DO i_state = 1, N_state
      SELECT CASE(Mult_ex(i_state))
        CASE('Singlet', 'Unknown')
          n_contribution(i_state) = n_elec_config(i_state)/2
          DO i = 1, n_elec_config(i_state)/2
            contribution_maxval = MAXVAL(contribution(:, i_state))
            IF(ABS(contribution_maxval) < Threshold_contribution) THEN
              n_contribution(i_state) = i - 1
              EXIT
            ELSE
            ENDIF    
            contribution_maxloc = MAXLOC(contribution(:, i_state), 1)
            contribution_sort(i, i_state) = contribution_maxval
            a_ras_sort(1, i, i_state) = A_ras(1, 2*contribution_maxloc-1, i_state)
            r_ras_sort(1, i, i_state) = R_ras(1, 2*contribution_maxloc-1, i_state)
            contribution(contribution_maxloc, i_state) = 0.0D0
          ENDDO
        CASE('Triplet')
          IF(Ms == 0) THEN
            n_contribution(i_state) = n_elec_config(i_state)/2
            DO i = 1, n_elec_config(i_state)/2
              contribution_maxval = MAXVAL(contribution(:, i_state))
              IF(ABS(contribution_maxval) < Threshold_contribution) THEN
                n_contribution(i_state) = i - 1
                EXIT
              ELSE
              ENDIF
              contribution_maxloc = MAXLOC(contribution(:, i_state), 1)
              contribution_sort(i, i_state) = contribution_maxval
              a_ras_sort(1, i, i_state) = A_ras(1, 2*contribution_maxloc-1, i_state)
              r_ras_sort(1, i, i_state) = R_ras(1, 2*contribution_maxloc-1, i_state)
              contribution(contribution_maxloc, i_state) = 0.0D0
            ENDDO
          ELSEIF(Ms ==  1 .OR. Ms == -1) THEN
            n_contribution(i_state) = n_elec_config(i_state)
            DO i = 1, n_elec_config(i_state)
              contribution_maxval = MAXVAL(contribution(:, i_state))
              IF(ABS(contribution_maxval) < Threshold_contribution) THEN
                n_contribution(i_state) = i - 1
                EXIT
              ELSE
              ENDIF
              contribution_maxloc = MAXLOC(contribution(:, i_state), 1)
              contribution_sort(i, i_state) = contribution_maxval
              a_ras_sort(1, i, i_state) = A_ras(1, contribution_maxloc, i_state)
              r_ras_sort(1, i, i_state) = R_ras(1, contribution_maxloc, i_state)
              contribution(contribution_maxloc, i_state) = 0.0D0
            ENDDO
          ELSE
          ENDIF
        CASE DEFAULT
      END SELECT
    ENDDO


!   Writing CI coefficients in terms of spin-orbital numbers
    WRITE(6,'(1X)')
    n_repeat = 45
    WRITE(6,'(1X, A)') 'CI coefficients and electronic configurations'
    WRITE(6,'(1X, A)') 'in terms of spin-orbital numbers'
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    IF(Property == 'Check' .OR. Property == 'Cicoef_only') THEN
      WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
      WRITE(6,'(1X, A)') 'HF ground state'
      WRITE(6,'(1X, A)') 'Singlet'
      WRITE(6,'(1X)')
      WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
      DO i_state = 1, n_state
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', i_state, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', i_state
        WRITE(6,'(1X, A)') mult_ex(i_state)
        WRITE(6,'(1X)')
        DO i_elec_config = 1, n_elec_config(i_state)
          WRITE(6,'(1X, F11.7, F13.7, I4, A3, I4)') &
         &X_plus_y(i_elec_config, i_state), X_minus_y(i_elec_config, i_state), &
         &A_ras(1, i_elec_config, i_state), arrow, &
         &R_ras(1, i_elec_config, i_state)
        ENDDO
      ENDDO
    ELSE
      IF(Bra == 0) THEN
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
      ELSE ! Bra /= 0  
        WRITE(6,'(1X, A, I0, A)') '| ', Bra, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', Bra
        WRITE(6,'(1X, A)') mult_ex(Bra)
        WRITE(6,'(1X)')
        DO i_elec_config = 1, n_elec_config(Bra)
          WRITE(6,'(1X, F11.7, F13.7, I4, A3, I4)') &
         &X_plus_y(i_elec_config, Bra), X_minus_y(i_elec_config, Bra), &
         &A_ras(1, i_elec_config, Bra), arrow, &
         &R_ras(1, i_elec_config, Bra)
        ENDDO
      ENDIF

      IF(Ket == Bra) THEN
      ELSEIF(Ket == 0) THEN
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
      ELSE ! Ket /= 0
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', Ket, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', Ket
        WRITE(6,'(1X, A)') mult_ex(Ket)
        WRITE(6,'(1X)')
        DO i_elec_config = 1, n_elec_config(Ket)
          WRITE(6,'(1X, F11.7, F13.7, I4, A3, I4)') &
         &X_plus_y(i_elec_config, Ket), X_minus_y(i_elec_config, Ket), &
         &A_ras(1, i_elec_config, Ket), arrow, &
         &R_ras(1, i_elec_config, Ket)
        ENDDO
      ENDIF

    ENDIF 
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    WRITE(6,'(1X)')

!   Writing CI coefficients in terms of spin-orbital names
    WRITE(6,'(1X)')
    n_repeat = 45
    WRITE(6,'(1X, A)') 'CI coefficients and electronic configurations'
    WRITE(6,'(1X, A)') 'in terms of spin-orbital names'
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    IF(Property == 'Check' .OR. Property == 'Cicoef_only') THEN
      WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
      WRITE(6,'(1X, A)') 'HF ground state'
      WRITE(6,'(1X, A)') 'Singlet'
      WRITE(6,'(1X)')
      WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
      DO i_state = 1, n_state
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', i_state, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', i_state
        WRITE(6,'(1X, A)') mult_ex(i_state)
        WRITE(6,'(1X)')
        DO i_elec_config = 1, n_elec_config(i_state)
          so_occ_name(1) = So_name(A_ras(1, i_elec_config, i_state))
          so_unocc_name(1) = So_name(R_ras(1, i_elec_config, i_state))
          WRITE(6,'(1X, F11.7, F13.7, 1X, A9)') &
         &X_plus_y(i_elec_config, i_state), &
         &X_minus_y(i_elec_config, i_state), &
         &so_unocc_name(1)
          WRITE(6,'(27X, A9)') so_occ_name(1) 
        ENDDO
      ENDDO
    
    ELSE ! Property /= 'Check'

      IF(Bra == 0) THEN
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
      ELSE ! Bra /= 0
        WRITE(6,'(1X, A, I0, A)') '| ', Bra, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', Bra
        WRITE(6,'(1X, A)') mult_ex(Bra)
        WRITE(6,'(1X)')
        DO i_elec_config = 1, n_elec_config(Bra)
          so_occ_name(1) = So_name(A_ras(1, i_elec_config, Bra))
          so_unocc_name(1) = So_name(R_ras(1, i_elec_config, Bra))
          WRITE(6,'(1X, F11.7, F13.7, 2X, A9)') X_plus_y(i_elec_config, Bra), &
                       &X_minus_y(i_elec_config, Bra), &
                       &so_unocc_name(1)
        ENDDO               
        WRITE(6,'(27X, A9)') so_occ_name(1)
      ENDIF

      IF(Ket == Bra) THEN
      ELSEIF(Ket == 0) THEN
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6,'(1X, F11.7, 1X, A)') 1.0D0, 'HF determinant'
      ELSE ! Ket /= 0
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', Ket, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', Ket
        WRITE(6,'(1X, A)') mult_ex(Ket)
        WRITE(6,'(1X)')
        DO i_elec_config = 1, n_elec_config(Ket)
          so_occ_name(1) = So_name(A_ras(1, i_elec_config, Ket))
          so_unocc_name(1) = So_name(R_ras(1, i_elec_config, Ket))
          WRITE(6,'(1X, F11.7, F13.7, 2X, A9)') X_plus_y(i_elec_config, Ket), &
                       &X_minus_y(i_elec_config, Ket), &
                       &so_unocc_name(1)
          WRITE(6,'(27X, A9)') so_occ_name(1)
        ENDDO  
      ENDIF

    ENDIF
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    WRITE(6,'(1X)')
!   End of writing CI coefficients in terms of so numbers


    WRITE(6,'(1X, A)') 'Sum of squared CI coefficients of'
    WRITE(6,'(1X, A)') 'electronic state |i>'
    n_repeat = 16
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    WRITE(6,'(1X, A)') '  i   Sum[ci**2]'
    WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
    IF(Property == 'Check' .OR. Property == 'Cicoef_only') THEN
      WRITE(6,'(1X, I3, F13.8)') 0, 1.0D0
      DO i_state = 1, N_state
        WRITE(6,'(1X, I3, F13.8)') i_state, sum_ci2(i_state)
      ENDDO
    ELSE ! Property /= 'Check'
      IF(Bra == 0) THEN
        WRITE(6,'(1X, I3, F13.8)') 0, 1.0D0
      ELSE ! Bra /= 0
        WRITE(6,'(1X, I3, F13.8)') Bra, sum_ci2(Bra)
      ENDIF

      IF(Ket == Bra) THEN
      ELSEIF(Ket == 0) THEN
        WRITE(6,'(1X, I3, F13.8)') 0, 1.0D0
      ELSE ! Ket /= 0
        WRITE(6,'(1X, I3, F13.8)') Ket, sum_ci2(Ket)
      ENDIF
    ENDIF
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)

!   Writing contributions from individual electronic configurations
!   to electronic states
    WRITE(6,'(1X)')
    n_repeat = 45
    WRITE(6,'(1X, A)') 'Contributions from individual electronic configurations'
    WRITE(6,'(1X, A)') 'to electronic states in terms of spin-orbital numbers'
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    IF(Property == 'Check' .OR. Property == 'Cicoef_only') THEN
      WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
      WRITE(6,'(1X, A)') 'HF ground state'
      WRITE(6,'(1X, A)') 'Singlet'
      WRITE(6,'(1X)')
      WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'

      DO i_state = 1, N_state
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', i_state, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', i_state
        WRITE(6,'(1X, A)') mult_ex(i_state)
        WRITE(6,'(1X)')
        DO i = 1, n_contribution(i_state)
          WRITE(6,'(1X, F6.2, I4, A3,  I4)') contribution_sort(i, i_state), &
                       &(a_ras_sort(1, i, i_state)+1)/2, arrow, &
                       &(r_ras_sort(1, i, i_state)+1)/2
        ENDDO
      ENDDO
    ELSE ! Property /= 'Check'
      IF(Bra == 0) THEN
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
      ELSE ! Bra /= 0
        WRITE(6,'(1X, A, I0, A)') '| ', Bra, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', Bra
        WRITE(6,'(1X, A)') mult_ex(Bra)
        WRITE(6,'(1X)')
        DO i = 1, n_contribution(Bra)
          WRITE(6,'(1X, F6.2, I4, A3,  I4)') contribution_sort(i, Bra), &
                       &(a_ras_sort(1, i, Bra)+1)/2, arrow, &
                       &(r_ras_sort(1, i, Bra)+1)/2
        ENDDO
      ENDIF

      IF(Ket == Bra) THEN
      ELSEIF(Ket == 0) THEN
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
      ELSE ! Ket /= 0
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', Ket, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', Ket
        WRITE(6,'(1X, A)') mult_ex(Ket)
        WRITE(6,'(1X)')
        DO i = 1, n_contribution(Ket)
          WRITE(6,'(1X, F6.2, I4, A3,  I4)') contribution_sort(i, Ket), &
                       &(a_ras_sort(1, i, Ket)+1)/2, arrow, &
                       &(r_ras_sort(1, i, Ket)+1)/2
        ENDDO
      ENDIF
    ENDIF
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    WRITE(6,'(1X)')
    n_repeat = 45
    WRITE(6,'(1X, A)') 'Contributions from individual electronic configurations'
    WRITE(6,'(1X, A)') 'to electronic states in terms of spin-orbital names'
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)
    IF(Property == 'Check' .OR. Property == 'Cicoef_only') THEN
      WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
      WRITE(6,'(1X, A)') 'HF ground state'
      WRITE(6,'(1X, A)') 'Singlet'
      WRITE(6,'(1X)')
      WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'

      DO i_state = 1, N_state
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', i_state, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', i_state
        WRITE(6,'(1X, A)') mult_ex(i_state)
        WRITE(6,'(1X)')
        DO i = 1, n_contribution(i_state)
          WRITE(6, '(1X, F6.2, 2X, A9)') contribution_sort(i, i_state), &
                       &So_name(r_ras_sort(1, i, i_state))
          WRITE(6, '(9X, A9)') So_name(a_ras_sort(1, i, i_state)) 
        ENDDO
      ENDDO
    ELSE
      IF(Bra == 0) THEN
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
      ELSE ! Bra /= 0
        WRITE(6,'(1X, A, I0, A)') '| ', Bra, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', Bra
        WRITE(6,'(1X, A)') mult_ex(Bra)
        WRITE(6,'(1X)')
        DO i = 1, n_contribution(Bra)
          WRITE(6, '(1X, F6.2, 2X, A9)') contribution_sort(i, Bra), &
                       &So_name(r_ras_sort(1, i, Bra))
          WRITE(6, '(9X, A9)') So_name(a_ras_sort(1, i, Bra)) 
        ENDDO
      ENDIF

      IF(Ket == Bra) THEN
      ELSEIF(Ket == 0) THEN
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', 0, ' >'
        WRITE(6,'(1X, A)') 'HF ground state'
        WRITE(6,'(1X, A)') 'Singlet'
        WRITE(6,'(1X)')
        WRITE(6, '(1X, F6.2, 2X, A)') 100.0D0,  'HF determinant'
      ELSE ! Ket /= 0
        WRITE(6,'(1X, A)') REPEAT('-', n_repeat)
        WRITE(6,'(1X, A, I0, A)') '| ', Ket, ' >'
        WRITE(6,'(1X, A, I3)') 'TD-HF/DFT excited state ', Ket
        WRITE(6,'(1X, A)') mult_ex(Ket)
        WRITE(6,'(1X)')
        DO i = 1, n_contribution(Ket)
          WRITE(6, '(1X, F6.2, 2X, A9)') contribution_sort(i, Ket), &
                       &So_name(r_ras_sort(1, i, Ket))
          WRITE(6, '(9X, A9)') So_name(a_ras_sort(1, i, Ket)) 
        ENDDO
      ENDIF
    ENDIF
    WRITE(6,'(1X, A)') REPEAT('=', n_repeat)

!   2021/01/06 
    IF(Property == 'Cicoef_only') THEN
      n_mocube = 0
      DO i_state = 1, N_state
        DO i = 1, n_contribution(i_state)
          IF(contribution_sort(i, i_state) >= 1.0D0) THEN
            n_mocube = n_mocube + 1
          ELSE
          ENDIF
        ENDDO
      ENDDO
      n_mocube = n_mocube*2
   
      i_mocube = 0
      ALLOCATE(num_mocube(1:n_mocube)); num_mocube = 0
      ALLOCATE(so_name_mocube(1:n_mocube)); so_name_mocube = REPEAT(' ', 21)
   
      DO i_state = 1, N_state
        DO i = 1, n_contribution(i_state)
          IF(contribution_sort(i, i_state) >= 1.0D0) THEN
            i_mocube = i_mocube + 1
            num_mocube(i_mocube) = (a_ras_sort(1, i, i_state)+1)/2
            so_name_mocube(i_mocube) = So_name(a_ras_sort(1, i, i_state))
            i_mocube = i_mocube + 1
            num_mocube(i_mocube) = (r_ras_sort(1, i, i_state)+1)/2
            so_name_mocube(i_mocube) = So_name(r_ras_sort(1, i, i_state))
          ELSE
          ENDIF
        ENDDO
      ENDDO
      !DO i_mocube = 1, n_mocube
      !  WRITE(*,*) i_mocube, num_mocube(i_mocube)
      !ENDDO
   
      ALLOCATE(num_mocube_written(1:n_mocube)); num_mocube_written = 0
      fid = Fmocube
      OPEN(fid, FILE = Fname(fid)) 

      WRITE(fid, '(A)') 'TD_FCHK=../TD_PBE0.fchk' 
      WRITE(fid, '(A)') 'if [ -e $TD_FCHK ] ; then' 
      
      i_mocube = 1
      WRITE(fid, '(A,I0,A,I0,A)') &
     &'  cubegen 0 MO=', num_mocube(1), &
     &' '//TRIM(Fchk_mocube)//' ', num_mocube(1), '_'//TRIM(So_name_mocube(1))//'.cube 0 h'
   
      num_mocube_written(1) = num_mocube(1)
      DO i_mocube = 2, n_mocube
        write_num_mocube = 'Yes'
        DO j_mocube = 1, i_mocube - 1
          IF(num_mocube(i_mocube) == num_mocube_written(j_mocube)) THEN
            write_num_mocube = 'No'
          ELSE
          ENDIF
        ENDDO 
        IF(write_num_mocube == 'Yes') THEN
          WRITE(fid, '(A,I0,A,I0,A)') &
         &'  cubegen 0 MO=', num_mocube(i_mocube), &
         &' '//TRIM(Fchk_mocube)//' ', num_mocube(i_mocube), &
         &'_'//TRIM(So_name_mocube(i_mocube))//'.cube 0 h'
          num_mocube_written(i_mocube) = num_mocube(i_mocube)
        ELSE
        ENDIF
      ENDDO
      DEALLOCATE(num_mocube, num_mocube_written)
      DEALLOCATE(So_name_mocube)

      WRITE(fid, '(A)') 'else'
      WRITE(fid, '(A)') "  echo 'Missing '$TD_FCHK"
      WRITE(fid, '(A)') '  exit'
      WRITE(fid, '(A)') 'fi'
    ELSE
    ENDIF

    CLOSE(fid) 

    DEALLOCATE(contribution, contribution_sort, n_contribution)
    DEALLOCATE(a_ras_sort, r_ras_sort)
    DEALLOCATE(sum_ci2)

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE read_data_ci_coef_td


  SUBROUTINE so_num2so_name

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Constants
    CHARACTER(LEN=100), PARAMETER :: name_program = 'so_num2so_name'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables

    INTEGER :: i_cgf
    CHARACTER(LEN=1), ALLOCATABLE :: spin_name(:)
    INTEGER :: int_mo_name
    CHARACTER(LEN=1) :: char_int_mo_name_1
    CHARACTER(LEN=2) :: char_int_mo_name_2
    CHARACTER(LEN=3) :: char_int_mo_name_3
    CHARACTER(LEN=4) :: char_int_mo_name_4
    CHARACTER(LEN=5) :: char_int_mo_name_5
  

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   Converting so number to so name

    ALLOCATE(So_name(1:N_so)); So_name = REPEAT(' ', 21)
    ALLOCATE(spin_name(1:N_so)); spin_name = ' '
    DO i_cgf = 1, N_so_a
      spin_name(2*i_cgf-1) = 'a'
      spin_name(2*i_cgf)   = 'b'
      IF(1 <= i_cgf .AND. i_cgf <= Ne/2) THEN ! i_cgf = 1, 2, ..., HOMO
        int_mo_name = i_cgf - Ne/2
        IF(-9999 <= int_mo_name .AND. int_mo_name <= -1000) THEN
          WRITE(char_int_mo_name_5, '(I5)') int_mo_name
          So_name(2*i_cgf-1) = 'HOMO'//char_int_mo_name_5//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'HOMO'//char_int_mo_name_5//spin_name(2*i_cgf)
        ELSEIF(-999 <= int_mo_name .AND. int_mo_name <= -100) THEN
          WRITE(char_int_mo_name_4, '(I4)') int_mo_name
          So_name(2*i_cgf-1) = 'HOMO'//char_int_mo_name_4//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'HOMO'//char_int_mo_name_4//spin_name(2*i_cgf)
        ELSEIF(-99 <= int_mo_name .AND. int_mo_name <= -10) THEN
          WRITE(char_int_mo_name_3, '(I3)') int_mo_name
          So_name(2*i_cgf-1) = 'HOMO'//char_int_mo_name_3//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'HOMO'//char_int_mo_name_3//spin_name(2*i_cgf)
        ELSEIF(-9 <= int_mo_name .AND. int_mo_name <= -1) THEN
          WRITE(char_int_mo_name_2, '(I2)') int_mo_name
        !WRITE(*,*) char_int_mo_name_2
          So_name(2*i_cgf-1) = 'HOMO'//char_int_mo_name_2//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'HOMO'//char_int_mo_name_2//spin_name(2*i_cgf)
        ELSEIF(int_mo_name == 0) THEN
          So_name(2*i_cgf-1) = 'HOMO'//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'HOMO'//spin_name(2*i_cgf)
        ELSE
          STOP
        ENDIF

      ELSEIF(Ne/2 + 1 <= i_cgf) THEN ! i_cgf = LUMO, LUMO + 1, ...
        int_mo_name = i_cgf - Ne/2 - 1
        IF(int_mo_name == 0) THEN
          So_name(2*i_cgf-1) = 'LUMO'//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'LUMO'//spin_name(2*i_cgf)
        ELSEIF(1 <= int_mo_name .AND. int_mo_name <= 9) THEN
          WRITE(char_int_mo_name_1, '(I1)') int_mo_name
          So_name(2*i_cgf-1) = 'LUMO+'//char_int_mo_name_1//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'LUMO+'//char_int_mo_name_1//spin_name(2*i_cgf)
        ELSEIF(10 <= int_mo_name .AND. int_mo_name <= 99) THEN
          WRITE(char_int_mo_name_2, '(I2)') int_mo_name
          So_name(2*i_cgf-1) = 'LUMO+'//char_int_mo_name_2//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'LUMO+'//char_int_mo_name_2//spin_name(2*i_cgf)
        ELSEIF(100 <= int_mo_name .AND. int_mo_name <= 999) THEN
          WRITE(char_int_mo_name_3, '(I3)') int_mo_name
          So_name(2*i_cgf-1) = 'LUMO+'//char_int_mo_name_3//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'LUMO+'//char_int_mo_name_3//spin_name(2*i_cgf)
        ELSEIF(1000 <= int_mo_name .AND. int_mo_name <= 9999) THEN
          WRITE(char_int_mo_name_4, '(I4)') int_mo_name
          So_name(2*i_cgf-1) = 'LUMO+'//char_int_mo_name_4//spin_name(2*i_cgf-1)
          So_name(2*i_cgf) = 'LUMO+'//char_int_mo_name_4//spin_name(2*i_cgf)
        ELSE
          STOP
        ENDIF
      ELSE
        STOP
      ENDIF
    ENDDO

    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Names of spin orbitals'
    WRITE(6,'(1X, A)') REPEAT('=', 33)
      DO i_cgf = N_so/2, 1, -1
        WRITE(*, '(1X, 2I5, 3X, A9, 2X, A9)') &
       &2*i_cgf-1, 2*i_cgf, So_name(2*i_cgf-1), So_name(2*i_cgf)
        IF(i_cgf == Ne/2+1) WRITE(6,'(1X, A)') REPEAT('-', 33)
      ENDDO  
    WRITE(6,'(1X, A)') REPEAT('=', 33)

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE so_num2so_name


  SUBROUTINE gen_lattice_point

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Constants
    CHARACTER(LEN=100), PARAMETER :: name_program = 'gen_lattice_point'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables
    INTEGER mx, my, mz
    INTEGER ix, iy, iz
    DOUBLE PRECISION xyz_center(1:3)

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

    mx = 0; my = 0; mz = 0
    xyz_center = 0.0D0

    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Calculating xyz coordinates of lattice points'

    ALLOCATE(xyz_lattice_point(1:3, 1:Nx, 1:Ny, 1:Nz)); xyz_lattice_point = 0.0D0

    IF(MOD(Nx,2)==1 .OR. MOD(Ny,2)==1 .OR. MOD(Nz,2)==1) THEN
      WRITE(6,'(1X, A)') 'Only single lattice point'
      WRITE(6,'(1X, A)') 'Unable to generate cube files'
      CALL write_messages(-9999, Text_blank, type_program, name_program)
    ELSE
      mx = nx/2
      my = ny/2
      mz = nz/2
    ENDIF

    xyz_center(1) = 0.5D0*(xmax + xmin)
    xyz_center(2) = 0.5D0*(ymax + ymin)
    xyz_center(3) = 0.5D0*(zmax + zmin)

    ! x-coordibate
    DO ix = 1, mx
      xyz_lattice_point(1, mx+ix, :, :)&
      = xyz_center(1) + (-0.5D0 + DBLE(ix)) * dx
      xyz_lattice_point(1, mx+1-ix, :, :)&
      = xyz_center(1) + (0.5D0 - DBLE(ix)) * dx
    ENDDO
    ! y-coordibate
    DO iy = 1, my
      xyz_lattice_point(2, :, my+iy, :)&
      = xyz_center(2) + (-0.5D0 + DBLE(iy)) * dy
      xyz_lattice_point(2, :, my+1-iy, :)&
      = xyz_center(2) + (0.5D0 - DBLE(iy)) * dy
    ENDDO
    ! z-coordibate
    DO iz = 1, mz
      xyz_lattice_point(3, :, :, mz+iz)&
      = xyz_center(3) + (-0.5D0 + DBLE(iz)) * dz
      xyz_lattice_point(3, :, :, mz+1-iz)&
      = xyz_center(3) + (0.5D0 - DBLE(iz)) * dz
    ENDDO

    xmin = xyz_lattice_point(1,  1,  1,  1)
    xmax = xyz_lattice_point(1, nx,  1,  1)
    ymin = xyz_lattice_point(2,  1,  1,  1)
    ymax = xyz_lattice_point(2,  1, ny,  1)
    zmin = xyz_lattice_point(3,  1,  1,  1)
    zmax = xyz_lattice_point(3,  1,  1, nz)

    WRITE(6,'(1X, A)') 'Done'

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE gen_lattice_point


  SUBROUTINE cube_grid

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Constants
    CHARACTER(LEN=100), PARAMETER :: name_program = 'cube_grid'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables 
    INTEGER :: i_line
    CHARACTER(LEN=100) text_temp
    INTEGER :: i_symbol
    CHARACTER(LEN=100) :: text_left
    CHARACTER(LEN=100) :: text_right
    CHARACTER(LEN=4) :: real_int_char
    DOUBLE PRECISION :: value_real
    INTEGER :: value_int

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   ----------------------
!   Initializing variables
!   ----------------------

!   Global variables
    Dtau = 0.0D0
    Dx = 0.25D0; Dy = 0.25D0; Dz = 0.25D0
    Delta = 2.0D0
    Nx=0; Ny=0; Nz=0

    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Calculating grid size for cube files'

!   Default grids
    Xmin = MINVAL(Xyznuc(1,:)) - Delta; Xmax = MAXVAL(Xyznuc(1,:)) + Delta
    Ymin = MINVAL(Xyznuc(2,:)) - Delta; Ymax = MAXVAL(Xyznuc(2,:)) + Delta
    Zmin = MINVAL(Xyznuc(3,:)) - Delta; Zmax = MAXVAL(Xyznuc(3,:)) + Delta

!   Local variables
    text_temp = REPEAT(' ', 100)
    text_left = REPEAT(' ', 100)
    text_right = REPEAT(' ', 100)
    real_int_char = REPEAT(' ', 4)
    value_real = 0.0D0
    value_int = 0

    IF(n_line_grid /= 0 .AND.&
      n_line_grid + 1 < n_line_end_grid) THEN
      DO i_line = n_line_grid + 1, n_line_end_grid - 1
        text_temp = Text_input(i_line)
        i_symbol = INDEX(text_temp, '=')
        CALL get_text_left_symbol(text_temp, i_symbol, text_left)
        CALL get_text_right_symbol(text_temp, i_symbol, text_right)
        CALL get_typ(text_left, real_int_char)
        SELECT CASE(real_int_char)
          CASE('REAL'); CALL get_real(text_right, value_real)
          CASE('INT');  CALL get_int(text_right, value_int)
          CASE DEFAULT
        END SELECT
        SELECT CASE(text_left)
          CASE('Xmin'); Xmin = value_real
          CASE('Xmax'); Xmax = value_real
          CASE('Ymin'); Ymin = value_real
          CASE('Ymax'); Ymax = value_real
          CASE('Zmin'); Zmin = value_real
          CASE('Zmax'); Zmax = value_real
          CASE('Dx'); Dx = value_real
          CASE('Dy'); Dy = value_real
          CASE('Dz'); Dz = value_real
          CASE DEFAULT
            WRITE(6,'(1X, A)') 'Incorrect GRID keyward'
            CALL write_messages(-9999, Text_blank, type_program, name_program)
        END SELECT
      ENDDO
    ELSE
    ENDIF
    Dtau = Dx*Dy*Dz

    DO
      Nx = INT((Xmax - Xmin)/Dx) + 1
      IF(Nx < 0) THEN
        WRITE(6,'(1X, A)') 'The number of x points (Nx) is negative'
        WRITE(6,'(1X, A)') 'Nx must be positive'
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ELSEIF(Nx > Max_Nx) THEN
        WRITE(6,'(1X, A)')     'The number of x points (Nx) is too large'
        WRITE(6,'(1X, A, I0)') 'Nx must be less than or equal to ', Max_nx
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ELSE
      ENDIF
      SELECT CASE(MOD(Nx,2))
        CASE(0)
          EXIT
        CASE(1)
          Xmin = Xmin - 0.05D0*Dx
          Xmax = Xmax + 0.05D0*Dx
          CYCLE
      END SELECT
    ENDDO

    DO
      Ny = INT((Ymax - Ymin)/Dy) + 1
      IF(Ny < 0) THEN
        WRITE(6,'(1X, A)') 'The number of y points (Ny) is negative'
        WRITE(6,'(1X, A)') 'Ny must be positive'
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ELSEIF(Ny > Max_ny) THEN
        WRITE(6,'(1X, A)')     'The number of y points (Ny) is too large'
        WRITE(6,'(1X, A, I0)') 'Nz must be less than or equal to ', Max_ny
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ELSE
      ENDIF
      SELECT CASE(MOD(Ny,2))
        CASE(0)
          EXIT
        CASE(1)
          Ymin = Ymin - 0.05D0*Dy
          Ymax = Ymax + 0.05D0*Dy
          CYCLE
      END SELECT
    ENDDO

    DO
      Nz = INT((Zmax - Zmin)/Dz) + 1
      IF(Nz < 0) THEN
        WRITE(6,'(1X, A)') 'The number of z points (Nz) is negative'
        WRITE(6,'(1X, A)') 'Nz must be positive'
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ELSEIF(Nz > 500) THEN
        WRITE(6,'(1X, A)')     'The number of z points (Nz) is too large'
        WRITE(6,'(1X, A, I0)') 'Ny must be less than or equal to ', Max_nz
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ELSE
      ENDIF
      SELECT CASE(MOD(nz,2))
        CASE(0)
          EXIT
        CASE(1)
          Zmin = Zmin - 0.05D0*Dz
          Zmax = Zmax + 0.05D0*Dz
          CYCLE
      END SELECT
    ENDDO
    WRITE(6,'(1X, A)') 'Done'

    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '--------------------------'
    WRITE(6,'(1X, A)') 'Information for cube grids'
    WRITE(6,'(1X, A)') '--------------------------'
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Minimum and maximum values of'
    WRITE(6,'(1X, A)') 'x, y, and z coordinates (a.u.)'
    WRITE(6,'(1X, A)') '========================='
    WRITE(6,'(1X, A)') '        Min         Max  '
    WRITE(6,'(1X, A)') '-------------------------'
    WRITE(6,'(1X, A, 2F12.5)') 'x', Xmin, Xmax
    WRITE(6,'(1X, A, 2F12.5)') 'y', Ymin, Ymax
    WRITE(6,'(1X, A, 2F12.5)') 'z', Zmin, Zmax
    WRITE(6,'(1X, A)') '========================='

    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') 'Grid size (a.u.) and the number of'
    WRITE(6,'(1X, A)') 'lattice ponts of x, y, and z coordinates'
    WRITE(6,'(1X, A)') '===================='
    WRITE(6,'(1X, A, F8.4, 4X, A, I4)') 'Dx', Dx, 'Nx', Nx
    WRITE(6,'(1X, A, F8.4, 4X, A, I4)') 'Dy', Dy, 'Ny', Ny
    WRITE(6,'(1X, A, F8.4, 4X, A, I4)') 'Dz', Dz, 'Nz', Nz
    WRITE(6,'(1X, A)') '===================='

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE cube_grid


! -----------------
! Vibrational modes
! -----------------

  SUBROUTINE read_data_vib_0

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_data_vib_0'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   Initializing global variable
    N_mode = 0
    OPEN(Fcontrol_vib, FILE = Fname(Fcontrol_vib))
      READ(Fcontrol_vib,*) N_mode ! Number of vibrational modes
    CLOSE(Fcontrol_vib)

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE read_data_vib_0

  SUBROUTINE read_data_vib_1

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_data_vib_1'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables
    DOUBLE PRECISION, ALLOCATABLE :: vibmode_cart(:,:,:)
    DOUBLE PRECISION :: normfac_vibmode
    INTEGER :: i_atm, i_xyz, i_mode, j_mode

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   ----------------------
!   Initializing variables
!   ----------------------

!   Global variables
    ALLOCATE(Atmwt(1:N_atm)); Atmwt = 0.0D0
    ALLOCATE(Sqrt_atmwt(1:N_atm)); Sqrt_atmwt = 0.0D0
    ALLOCATE(Vibmode(1:3, 1:N_atm, 1:N_mode)); Vibmode = 0.0D0
    ALLOCATE(Freq(1:N_mode)); Freq = 0.0D0

!   Local variables
    ALLOCATE(Vibmode_cart(1:3, 1:N_atm, 1:N_mode)); Vibmode_cart = 0.0D0
    normfac_vibmode = 0.0D0

    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Reading frequencies from '//TRIM(Fname(Ffreq))
    WRITE(6,'(1X, A)') 'Reading vibrational modes from '//TRIM(Fname(Fvibmode))
    OPEN(Ffreq, FILE = Fname(Ffreq))
    OPEN(Fvibmode, FILE = Fname(Fvibmode))
       DO i_mode = 1, n_mode
         READ(Ffreq,*) freq(i_mode)
         DO i_atm = 1, n_atm
           DO i_xyz = 1, 3
             READ(Fvibmode,*) vibmode_cart(i_xyz, i_atm, i_mode)
             ! vibmode_cart [u -1/2]
           ENDDO
         ENDDO
       ENDDO
    CLOSE(Ffreq)
    CLOSE(Fvibmode)
    WRITE(6,'(1X, A)') 'Done'

    OPEN(Fatmwt, FILE = Fname(Fatmwt))
      DO i_atm = 1, n_atm
        READ(Fatmwt,*) atmwt(i_atm) ! Atom weight in [u]
      ENDDO
    CLOSE(Fatmwt)
    ! atmwt [u]
    ! atmwt*Amc [kg]
    ! atmwt*Amc/Me [a.u.]
    atmwt = atmwt*Amc/Me ! [a.u.]
    sqrt_atmwt = DSQRT(atmwt) ! [a.u.]

    ! Normalization factor for vibrational modes
    DO i_atm = 1, n_atm
      vibmode(:,i_atm,:) = sqrt_atmwt(i_atm)&
                           * vibmode_cart(:,i_atm,:)
    ENDDO
    DO i_mode = 1, n_mode
      normfac_vibmode = 0.0D0
      DO i_atm = 1, n_atm
        DO i_xyz = 1, 3
          normfac_vibmode&
          = normfac_vibmode&
          + vibmode(i_xyz, i_atm, i_mode)**2.0D0
        ENDDO
      ENDDO
      DO i_atm = 1, n_atm
        vibmode(:, i_atm, i_mode)&
        = vibmode(:, i_atm, i_mode)&
        / DSQRT(normfac_vibmode)
      ENDDO
    ENDDO ! i_mode = 1, n_mode

    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Orthonormalization of vibrational vectors'
    DO i_mode = 1, n_mode
      DO j_mode = 1, n_mode
        normfac_vibmode = 0.0D0
        DO i_atm = 1, n_atm
          DO i_xyz = 1, 3
            normfac_vibmode&
            = normfac_vibmode&
            + vibmode(i_xyz, i_atm, i_mode)&
            * vibmode(i_xyz, i_atm, j_mode)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    WRITE(6,'(1X, A)') 'Done'

  END SUBROUTINE read_data_vib_1


  SUBROUTINE read_mode_calc

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Constants
    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_inp_mode_calc'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables
    INTEGER :: i_line, i_mode_calc
    INTEGER :: value_int
    CHARACTER(LEN=100) :: text_temp
    INTEGER :: counter

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

    value_int = 0
    text_temp = REPEAT(' ', 100)

    N_mode_calc = 0
    IF(n_line_mode_calc /= 0 .AND.&
      n_line_mode_calc + 1 < n_line_end_mode_calc) THEN
      N_mode_calc = n_line_end_mode_calc - n_line_mode_calc - 1
    ELSE
      N_mode_calc = N_mode
    ENDIF

!   ----------------------
!   Initializing variables
!   ----------------------

!   Global variables
    ALLOCATE(Mode_calc(1:N_mode_calc)); Mode_calc = 0
    ALLOCATE(Freq_calc(1:N_mode_calc)); Freq_calc = 0.0D0
    ALLOCATE(Vibmode_calc(1:3, 1:N_atm, 1:N_mode_calc)); Vibmode_calc = 0.0D0

    IF(n_line_mode_calc /= 0 .AND.&
      n_line_mode_calc + 1 < n_line_end_mode_calc) THEN
      counter = 0
      DO i_line = n_line_mode_calc + 1, n_line_end_mode_calc - 1
        text_temp = Text_input(i_line)
        counter = counter + 1
        CALL get_int(text_temp, value_int)
        Mode_calc(counter) = value_int
      ENDDO
    ELSE
      DO i_mode_calc = 1, N_mode_calc
        Mode_calc(i_mode_calc) = i_mode_calc
      ENDDO
    ENDIF
    DO i_mode_calc = 1, N_mode_calc
      Freq_calc(i_mode_calc) = Freq(Mode_calc(i_mode_calc))
      Vibmode_calc(1:3, 1:N_atm, i_mode_calc)&
      = Vibmode(1:3, 1:N_atm, Mode_calc(i_mode_calc))
    ENDDO
    
    Freq_calc = Freq_calc * Scale_fac_freq

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE read_mode_calc


! -----------------
! Excitation energy
! -----------------

  SUBROUTINE read_data_ex_energy

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_data_ex_energy'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables
    INTEGER :: i_state

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   ----------------------
!   Initializing variables
!   ----------------------

!   Global variables
    ALLOCATE(Ex_energy(0:N_state)); Ex_energy = 0.0D0

    OPEN(Fex_energy, FILE = Fname(Fex_energy))
      SELECT CASE(Method)
        CASE('Cis', 'Td')
          Ex_energy(0) = 0.0D0
          DO i_state = 1, N_state
            READ(Fex_energy, *) Ex_energy(i_state)
          ENDDO
        CASE DEFAULT
          WRITE(6,'(1X, A, A)') 'Invalid method: ', Method
          CALL write_messages(-9999, Text_blank, type_program, name_program)
      END SELECT      
    CLOSE(Fex_energy)

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  END SUBROUTINE read_data_ex_energy

END MODULE global_read_data
