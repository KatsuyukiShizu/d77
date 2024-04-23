! This module is part of d77 and defines global constants.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE global_constants
  IMPLICIT NONE

  CHARACTER(LEN=8),   SAVE :: Program_ver 
  CHARACTER(LEN=100), SAVE :: Text_blank 

! --------
! CPU time
! --------

! Time_start : Start time
! Time_end   : End time
  DOUBLE PRECISION, SAVE :: Time_start, Time_end

! -------------------
! General information
! -------------------
!
! Max_ne            : Maximum number of electrons
! Max_n_atm         : Maximum number of atoms
! Max_n_mode        : Maximum number of vibrational modes
! Max_n_cgf         : Maximum number of contracted Gaussian functions (CGFs)
! Max_n_pgf         : Maximum number of primitive Gaussian functions (PGFs) for each CGF
! Max_n_state       : Maximum number of electronic states 
! Max_n_elec_config : Maximum number of electronic configurations
! Max_n_ci_ci       : Maximum number of CI pairs
! Max_n_tuple       : Maximum number of electronic excitations
  INTEGER, SAVE :: Max_ne, Max_n_atm, Max_n_mode
  INTEGER, SAVE :: Max_n_cgf, Max_n_pgf
  INTEGER, SAVE :: Max_n_state 
  INTEGER, SAVE :: Max_n_elec_config, Max_n_ci_ci, Max_n_tuple

! ----------------
! Threshold values
! ----------------

! Threshold values used in global_read_input.f90
! Threshold_s_cgf_default: Default value of the allowable error (Threshold_s_cgf)
!                          for the overlap integrals between CGFs (s_cgf).
!                          Matrix elements of s_cgf should be larger than Threshold_s_cgf.
  DOUBLE PRECISION, SAVE :: Threshold_s_cgf_default  
  DOUBLE PRECISION, SAVE :: Threshold_e_ab_default 
  DOUBLE PRECISION, SAVE :: Threshold_ci_ci_default  
  DOUBLE PRECISION, SAVE :: Threshold_contribution_default  
  DOUBLE PRECISION, SAVE :: Threshold_distance_default  

! -----
! Grids 
! -----

  DOUBLE PRECISION, SAVE :: Max_nx, Max_ny, Max_nz

! ----------------------
! Mathematical constants
! ----------------------

  DOUBLE PRECISION, SAVE :: Pi, Sqrt2

! ------------------
! Physical constants
! ------------------

! Physical constants obtained from
! The National Institute of Standards and Technology (NIST)
! https://physics.nist.gov/cuu/index.html

! Values of Fundamental Physical Constants
! https://physics.nist.gov/cuu/Constants/index.html

! -------------------
! Universal constants
! -------------------

! Planck : Planck constant (h) in [J s]
! Dirac  : h/2Pi in [J s]
! C_m    : Speed of light in vacuum in [m s-1]
! C_cm   : Speed of light in vacuum in [cm s-1]
! Eps0   : Electric constant (permittivity of vacuum) in [F m-1]
  DOUBLE PRECISION, SAVE :: Planck, Dirac, C_m, C_cm, Eps0

! ----------------------------
! Atomic and nuclear constants
! ----------------------------

! Fsc     : Fine-structure constant
! Inv_fsc : Inverse fine-structure constant
  DOUBLE PRECISION, SAVE :: Fsc, Inv_fsc

! -------------------------
! Electromagnetic constants
! -------------------------

! Echarg : Atomic unit of charge (elementary charge)
! MuB    : Bohr magneton
  DOUBLE PRECISION, SAVE :: Echarg, Joule_to_ev, MuB

! -----------------------------------
! Atomic units and conversion factors
! -----------------------------------

! Me        : Mass of electron (atomic unit of mass)
! Bohr      : Bohr radius (atomic unit of length)
! Hartree   : Atomic unit of energy
! Time_au   : Atomic unit of time (Dirac/Hartree)
! StatC     : statcoulomb (electrostatic system of units)
! Debye     : 1.0D-18 [statcoulomb cm]
! Debye_si  : debye in international system of units (SI units)
! Dipole_au : Unit of dipole moment in atomic units
  DOUBLE PRECISION, SAVE :: Me, Bohr, Hartree, Time_au
  DOUBLE PRECISION, SAVE :: Meter_to_bohr, Ang_to_bohr
  DOUBLE PRECISION, SAVE :: Hartree_to_ev, Ev_to_hartree
  DOUBLE PRECISION, SAVE :: StatC, Debye, Debye_si, Dipole_au
  DOUBLE PRECISION, SAVE :: Joule_to_wavenumber, Hartree_to_wavenumber

! --------------------------
! Physico-chemical constants
! --------------------------

! Amc      : Atomic mass constant
! Avogadro : Avogadro constant
! Kb       : Boltzmann constant
  DOUBLE PRECISION, SAVE :: Amc, Avogadro, Kb

! ---------
! I/O files 
! ---------
!
! Max_n_file: Maximum number of I/O files 
! Fname: File names
  INTEGER, SAVE :: Max_n_file 
  CHARACTER(LEN=30), ALLOCATABLE, SAVE :: Fname(:)
!
! Unit numbers 
!
  INTEGER, SAVE :: Finp
  INTEGER, SAVE :: Fcontrol_elec, Fcontrol_vib

  INTEGER, SAVE :: Fatmnum, Fnuccharg, Fatmwt, Fcoord, Fnuccharg_soc

  INTEGER, SAVE :: Fshtyp, Fnpgf_sh, Fsh2atm

  INTEGER, SAVE :: Fao
  INTEGER, SAVE :: Fcoef_soa_ref, Fcoef_sob_ref
  INTEGER, SAVE :: Fcoef_soa_sys, Fcoef_sob_sys

  INTEGER, SAVE :: Ffreq, Fvibmode

  INTEGER, SAVE :: Fexstate
  INTEGER, SAVE :: Fxy
  INTEGER, SAVE :: Fcicoef
  INTEGER, SAVE :: Fex_energy

  INTEGER, SAVE :: Fmocube

  INTEGER, SAVE :: Fopr_dmat
  INTEGER, SAVE :: Fdmat_cgf
  INTEGER, SAVE :: Fdmat_soc_cgf_x, Fdmat_soc_cgf_y, Fdmat_soc_cgf_z
  INTEGER, SAVE :: Fhcore
  INTEGER, SAVE :: Fke_cgf, Fpe_cgf_atm, Felfld_cgf_atm
  INTEGER, SAVE :: Fke_cgf_g16_format, Fpe_cgf_g16_format, Felfld_cgf_g16_format
  INTEGER, SAVE :: Fsoc_cgf_x, Fsoc_cgf_y, Fsoc_cgf_z
  INTEGER, SAVE :: Fsoc_cgf_atm_x, Fsoc_cgf_atm_y, Fsoc_cgf_atm_z
  INTEGER, SAVE :: Fsoc_cgf_x_g16_format, Fsoc_cgf_y_g16_format, Fsoc_cgf_z_g16_format
  INTEGER, SAVE :: Fvcc, Favcc, Fnac, Fkic, Fpm
  INTEGER, SAVE :: Fdipole_au, Fsoc_au

CONTAINS

  SUBROUTINE assign_globals
    Program_ver = '20240423'
    Text_blank = REPEAT(' ',100)
    Max_n_atm = 30; Max_n_mode = 100
    Max_ne = 100
    Max_n_cgf = 1000
    Max_n_pgf = 30
    Max_n_state = 20 
    Max_n_elec_config = 1000
    Max_n_ci_ci = Max_n_elec_config**2
    Max_n_tuple = 4 
    Threshold_s_cgf_default = 1.0D-5
    Threshold_e_ab_default = 1.0D-4 
    Threshold_ci_ci_default = 1.0D-8 
    Threshold_contribution_default = 1.0D-1 
    Threshold_distance_default = 4.0D0 ! [A]   
    Max_nx = 500; Max_ny = 500; Max_nz = 500

!   -----------------------------------
!   Mathematical and Physical constants
!   -----------------------------------

!   Physical constants obtained from
!   The National Institute of Standards and Technology (NIST)
!   https://physics.nist.gov/cuu/index.html
!   https://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=universal_in!

!   ----------------------
!   Mathematical constants
!   ----------------------

    Pi = 4.0D0*DATAN(1.0D0)
    Sqrt2 = DSQRT(2.0D0)

!   -------------------
!   Universal constants
!   -------------------

    Planck = 6.62607015D-34  ! [J s]
    Dirac  = 1.054571817D-34 ! [J s]
    C_m    = 2.99792458D+8   ! [m s-1]
    C_cm   = 2.99792458D+10  ! [cm s-1]

!   Electric constant (epsilon0)
    Eps0 = 8.8541878128D-12 ! [F m-1]
  
!   -------------------------
!   Electromagnetic constants
!   -------------------------

    Echarg = 1.6021766208D-19 ! [m]: 1 [e] = Echarg [C]
    Joule_to_ev = 1.0/Echarg ! [eV]: 1 [J] = Joule_to_ev [eV]

    MuB = 9.2740100783D-24 ! [J T-1]

!   ------------
!   Atomic units
!   ------------

    Me = 9.1093837015D-31 ! 1 [me] = Me [kg]
    Bohr = 4.0* Pi * Eps0 * Dirac**2.0D0 / (Me * Echarg**2.0D0) ! 1 [bohr] = Bohr [m]
!   Bohr = 5.2917721897717314E-011 (0.529177210903D-10, NIST)
    Meter_to_bohr = 1/Bohr ! 1 [m] = Meter_to_bohr [bohr]
    Ang_to_bohr = 1.0D-10*Meter_to_bohr ! 1 [A] = Ang_to_bohr [bohr] 
    Hartree = Echarg**2.0D0 / (4.0D0*Pi*Eps0*Bohr) ! 1 [hartree] = Hartree [J]
!   Hartree = 4.3597445838389197E-018 (4.3597447222071D-18, NIST)     
    Hartree_to_ev = Hartree * Joule_to_ev ! 1 [hartree] = Hartree_to_ev [eV]
!   Hartree_to_ev = 27.211386668405989     
    Ev_to_hartree = 1/Hartree_to_ev ! 1 [eV] = Ev_to_hartre [hartree]
!   Ev_to_hartree = 3.6749321605174148E-002

!   Atomic unit of time 
    Time_au = Dirac/Hartree ! 1 [atomic unit of time] = Time_au [s]
!   Time_au = 2.4188843074825784E-017

    StatC = 1.0D-1/C_m ! 1 [StatC] = StatC [C]
!   StatC = 3.3356409519815207E-010 
    Debye_si = 1.0D-18 * StatC * 1.0D-2 ! 1 [debye in SI unit] = Debye_si [C m]
!   Debye_si = 3.3356409519815207E-030
    Dipole_au = 1.0D0 * Echarg * Bohr ! atomic unit for dipole moment
!   Dipole_au = 8.4783535556893754E-030
    Debye = Dipole_au / Debye_si
!   Debye = 2.5417464522531574

!   ----------------------------
!   Atomic and nuclear constants
!   ----------------------------

!   Fine-structure constant
    Fsc = (Echarg**2.0D0)/(4.0D0*Pi*Eps0*Dirac*C_m)
!   Fsc = 7.2973524535065231E-003 (7.2973525693D-3, NIST)
    Inv_fsc = 1.0D0/Fsc
!   Inv_fsc = 137.03600125816592 (137.035999084D+0, NIST)

    Joule_to_wavenumber = 1.0D0/(Planck*C_cm) ! 1 [J] = Joule_to_wavenumber [cm-1]
!   Joule_to_wavenumber = 5.0341165675427087E+022
    Hartree_to_wavenumber = Hartree * Joule_to_wavenumber  ! 1 [J] = Joule_to_wavenumber [cm-1]
!   Hartree_to_wavenumber = 219474.62439758098

!   --------------------------
!   Physico-chemical constants
!   --------------------------

    Amc = 1.66053906660D-27 ! [kg]
    Avogadro = 6.02214076D+23 ! [mol-1]
    Kb = 1.380649D-23 ! [J K-1]

!   Unit numbers for input files
    Finp          =  7
    Fcontrol_elec =  8
    Fcontrol_vib  =  9
 
    Fatmnum       = 11
    Fnuccharg     = 12
    Fatmwt        = 13
    Fcoord        = 14
    Fnuccharg_soc = 15
 
    Fshtyp        = 21
    Fnpgf_sh      = 22
    Fsh2atm       = 23
 
    Fao           = 31
    Fcoef_soa_ref = 32
    Fcoef_sob_ref = 33
    Fcoef_soa_sys = 34
    Fcoef_sob_sys = 35
 
    Ffreq         = 51
    Fvibmode      = 52
 
    Fexstate      = 61
    Fxy           = 62
    Fcicoef       = 64
    Fex_energy    = 65

!   Unit numbers for output files
    Fke_cgf               = 101
    Fke_cgf_g16_format    = 102

    Fpe_cgf_atm           = 103
    Fpe_cgf_g16_format    = 104
    
    Felfld_cgf_atm        = 105
    Felfld_cgf_g16_format = 106
    
    Fhcore                = 107
    
    Fdipole_au            = 108
    Fsoc_au               = 109

    Favcc                 = 112
    Fvcc                  = 113
    Fnac                  = 116
    Fpm                   = 114
    Fkic                  = 115

    Fopr_dmat             = 201
    Fdmat_cgf             = 202
    Fdmat_soc_cgf_x       = 211
    Fdmat_soc_cgf_y       = 212
    Fdmat_soc_cgf_z       = 213
    Fke_cgf_g16_format    = 221
    Fpe_cgf_g16_format    = 222

    Fsoc_cgf_x            = 231
    Fsoc_cgf_y            = 232
    Fsoc_cgf_z            = 233
    Fsoc_cgf_atm_x        = 234
    Fsoc_cgf_atm_y        = 235
    Fsoc_cgf_atm_z        = 236
    Fsoc_cgf_x_g16_format = 237
    Fsoc_cgf_y_g16_format = 238
    Fsoc_cgf_z_g16_format = 239

    Fmocube               = 241

!   ----------
!   File names 
!   ----------

    Max_n_file = 1000
    ALLOCATE(Fname(1:Max_n_file))
    Fname = REPEAT(' ',30)

!   General information
    Fname(Finp)          = 'INPUT'
    Fname(Fcontrol_elec) = 'CONTROL_ELEC'
    Fname(Fcontrol_vib)  = 'CONTROL_VIB'
  
!   Information for nuclei  
    Fname(Fatmnum)       = 'ATMNUM'
    Fname(Fnuccharg)     = 'NUCCHARG'
    Fname(Fatmwt)        = 'ATMWT'
    Fname(Fcoord)        = 'COORD'
    Fname(Fnuccharg_soc) = 'ZEFF_SOC'

!   Basis functions  
    Fname(Fshtyp)        = 'SHELLTYP'
    Fname(Fnpgf_sh)      = 'NUMPGFSHELL'
    Fname(Fsh2atm)       = 'SHELL2ATM'
  
!   Molecular orbitals  
    Fname(Fao)           = 'ATMORB'
    Fname(Fcoef_soa_ref) = 'COEFAREF'
    Fname(Fcoef_sob_ref) = 'COEFBREF'
    Fname(Fcoef_soa_sys) = 'COEFASYS'
    Fname(Fcoef_sob_sys) = 'COEFBSYS'
  
!   Molecular vibrations  
    Fname(Ffreq)         = 'FREQ'
    Fname(Fvibmode)      = 'VIBMODE'
  
!   Electronic states   
    Fname(Fexstate)     = 'EXSTATE'
    Fname(Fxy)          = 'XY'
    Fname(Fcicoef)      = 'CICOEF'
    Fname(Fex_energy)   = 'EX_ENERGY'

    Fname(Fmocube)      = 'MOCUBE.sh'

!   Output files   
    Fname(Fhcore) = 'HCORE'

    Fname(Favcc)  = 'AVCC'
    Fname(Fvcc)   = 'VCC'
    !Fname(Fnac)   = 'NAC_Q'
    !Fname(Fpm)    = 'PM'
    !Fname(Fkic)   = 'KIC'

    Fname(Fdipole_au) = 'DM'
    Fname(Fsoc_au)   = 'SOC'

    Fname(Fdmat_soc_cgf_x)    = 'DMAT_SOC_CGF_X'
    Fname(Fdmat_soc_cgf_y)    = 'DMAT_SOC_CGF_Y'
    Fname(Fdmat_soc_cgf_z)    = 'DMAT_SOC_CGF_Z'

    Fname(Fke_cgf) = 'KE_CGF'
    Fname(Fke_cgf_g16_format) = 'KE_CGF_G16_FORMAT'
    
    Fname(Fpe_cgf_atm) = 'PE_CGF_ATM'
    Fname(Fpe_cgf_g16_format) = 'PE_CGF_G16_FORMAT'

    Fname(Felfld_cgf_atm) = 'ELFLD_CGF_ATM'
    Fname(Felfld_cgf_g16_format) = 'ELFLD_CGF_G16_FORMAT'

    Fname(Fsoc_cgf_x) = 'SOC_CGF_X'
    Fname(Fsoc_cgf_y) = 'SOC_CGF_Y'
    Fname(Fsoc_cgf_z) = 'SOC_CGF_Z'
    Fname(Fsoc_cgf_x_g16_format) = 'SOC_CGF_X_G16_FORMAT'
    Fname(Fsoc_cgf_y_g16_format) = 'SOC_CGF_Y_G16_FORMAT'
    Fname(Fsoc_cgf_z_g16_format) = 'SOC_CGF_Z_G16_FORMAT'
    Fname(Fsoc_cgf_atm_x) = 'SOC_CGF_ATM_X'
    Fname(Fsoc_cgf_atm_y) = 'SOC_CGF_ATM_Y'
    Fname(Fsoc_cgf_atm_z) = 'SOC_CGF_ATM_Z'
    
    Fname(Fopr_dmat) = 'OPR_DMAT'
    Fname(Fdmat_cgf) = 'DMAT_CGF'

  END SUBROUTINE assign_globals 
END MODULE global_constants







