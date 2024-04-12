! This module is part of dpp and defines global constants.
!
! dpp is a data preprocessing utility for d77.
!
! dpp is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE global_constants
  IMPLICIT NONE

  CHARACTER(LEN=20), SAVE :: Program_ver 
  CHARACTER(LEN=100), SAVE :: Text_blank 
  INTEGER, SAVE :: Max_ne, Max_n_atm, Max_n_mode, Max_n_hessian
  INTEGER, SAVE :: Max_n_cgf, Max_n_pgf, Max_n_so
  INTEGER, SAVE :: Max_n_state 
  INTEGER, SAVE :: Max_n_elec_config, Max_n_tuple

! ----------------------
! Mathematical constants
! ----------------------

  DOUBLE COMPLEX, SAVE :: Zeroc
  DOUBLE COMPLEX, SAVE :: Imag_unit
  DOUBLE PRECISION, SAVE :: Pi

! -------------------
! Universal constants
! -------------------

  DOUBLE PRECISION, SAVE :: Planck, Dirac, C_m, C_cm, Eps0

! -------------------------
! Electromagnetic constants
! -------------------------

  DOUBLE PRECISION, SAVE :: Echarg, Joule_to_ev, MuB

! ------------
! Atomic units
! ------------

  DOUBLE PRECISION, SAVE :: Bohr, Me, Hartree, Time_au
  DOUBLE PRECISION, SAVE :: Ang_to_bohr
  DOUBLE PRECISION, SAVE :: Hartree_to_ev, Ev_to_hartree
  DOUBLE PRECISION, SAVE :: Debye

! --------------------------
! Physico-chemical constants
! --------------------------

  DOUBLE PRECISION, SAVE :: Amc, Avogadro, Kb

  INTEGER, SAVE :: Threshold_tuple_default
  DOUBLE PRECISION, SAVE :: Threshold_ci_coef_default

! ------------
! Unit numbers 
! ------------

  INTEGER, SAVE :: Ffile_typ, &
                  &FR_control_elec, FW_control_elec, &
                  &FR_control_vib, FW_control_vib, &
                  &FR_atmnum, FW_atmnum, &
                  &FR_nuccharg, FW_nuccharg, &
                  &FR_atmwt, FW_atmwt, FW_atmwt_central, &
                  &FR_coord, FW_coord, FW_coord_central, &    
                  &FR_nuccharg_soc, FW_nuccharg_soc, &
                  &FR_shtyp, FW_shtyp, &
                  &FR_n_pgf_sh, FW_n_pgf_sh, &
                  &FR_sh2atm, FW_sh2atm, &
                  &FR_cntexp, FW_cntexp, &
                  &FR_cntcoef0, FW_cntcoef0, &
                  &FR_cntcoefp, FW_cntcoefp, &
                  &FW_basisset, FW_ao, &     
                  &FR_coefa_ref, FW_coefa_ref, &
                  &FR_coefb_ref, FW_coefb_ref, &
                  &FR_coefa_sys, FW_coefa_sys, &
                  &FR_coefb_sys, FW_coefb_sys, &
                  &FR_freq, FW_freq, &   
                  &FR_vibmode, FW_vibmode, FW_vibmode_modelsys_calc, &
                  &FR_hessian, FW_hessian, &
                  &FR_micopt, FW_micopt, &
                  &FR_cicoef_g16_cis, &
                  &FR_g16_cisd_log, &
                  &FR_g16_td_log_x, FR_g16_td_log_y, &
                  &FW_xy, FW_xy_count, &
                  &FR_ref_singlet_state, &
                  &FR_num_ex_states, &
                  &FR_orbnum_alpha, &
                  &FR_orbnum_beta, &
                  &FW_orbnum, &
                  &FR_threshold_values_ci, &
                  &FW_cicoef_as_read, &
                  &FW_cicoef, &
                  &FR_ex_energy_as_read, &
                  &FW_ex_energy, &
                  &FW_vib_fchk, &
                  &FW_atmwt_modelsys, FW_atmnum_modelsys, FW_coord_modelsys, &    
                  &FW_check_orthonorm_vibmode, &
                  &FR_duschinskymat, FW_duschinskymat, & 
                  &FR_amat, FW_amat, &
                  &FR_cmat, FW_cmat, & 
                  &FR_emat, FW_emat, & 
                  &FR_shiftvec, FW_shiftvec, &
                  &FR_bvec, FW_bvec, & 
                  &FR_dvec, FW_dvec


! ----------
! File names 
! ----------

  INTEGER, SAVE :: Max_n_file 
  CHARACTER(LEN=30), ALLOCATABLE, SAVE :: Fname(:)

CONTAINS

  SUBROUTINE assign_globals
    Program_ver = 'prepost_20190318'
    Text_blank = REPEAT(' ',100)
    Max_n_atm = 1000; Max_n_mode = 3000; Max_n_hessian = 9000000
    Max_ne = 500
    Max_n_cgf = 2000
    Max_n_so = 2*Max_n_cgf
    Max_n_pgf = 20
    Max_n_state = 50 
    Max_n_elec_config = 500000 
    Max_n_tuple = 8 
    Threshold_tuple_default = 8
    Threshold_ci_coef_default = 0.0D0


    Zeroc = (0.0D0, 0.0D0)
    Imag_unit = (0.0D0, 1.0D0)
    Pi = 4.0D0*DATAN(1.0D0)

    Planck = 6.62607015D-34  ! [J s]
    Dirac  = 1.054571817D-34 ! [J s]
    C_m    = 2.99792458D+8   ! [m s-1]
    C_cm   = 2.99792458D+10  ! [cm s-1]
    Eps0 = 8.8541878128D-12 ! [F m-1]

    Echarg = 1.6021766208D-19 ! [m]: 1 [e] = Echarg [C]
    Joule_to_ev = 1.0/Echarg ! [eV]: 1 [J] = Joule_to_ev [eV]
    MuB = 9.2740100783D-24 ! [J T-1]

    Me = 9.1093837015D-31 ! 1 [me] = Me [kg]
    Bohr = 0.529177210903D-10 ! 1 [bohr] = Bohr [m]
    Ang_to_bohr = 1.0D0/0.529177210903D0
  
    Hartree = 4.3597447222071E-18 ! 1 [hartree] = Hartree [J]
    Hartree_to_ev = Hartree * Joule_to_ev ! 1 [hartree] = Hartree_to_ev [eV]
    Ev_to_hartree = 1/Hartree_to_ev ! 1 [eV] = Ev_to_hartre [hartree]
 
    Time_au = Dirac/Hartree ! 1 [atomic unit of time] = Time_au [s]
    Debye = 2.5417464522531574

    Amc = 1.66053906660D-27 ! [kg]
    Avogadro = 6.02214076D+23 ! [mol-1]
    Kb = 1.380649D-23 ! [J K-1]

!   ------------
!   Unit numbers 
!   ------------

    Ffile_typ       = 100
    FR_control_elec = 101
    FW_control_elec = 102
    FR_control_vib  = 103
    FW_control_vib  = 104

    FR_atmnum   = 11 
    FW_atmnum   = 12
    FR_nuccharg = 13 
    FW_nuccharg = 14 
    FR_atmwt    = 15 
    FW_atmwt    = 16
    FR_coord    = 17 
    FW_coord    = 18
    FR_nuccharg_soc = 19 
    FW_nuccharg_soc = 20 

    FR_shtyp    = 21
    FW_shtyp    = 22
    FR_n_pgf_sh = 23
    FW_n_pgf_sh = 24
    FR_sh2atm   = 25
    FW_sh2atm   = 26
    FR_cntexp   = 27
    FW_cntexp   = 28
    FR_cntcoef0 = 29
    FW_cntcoef0 = 30
    FR_cntcoefp = 31
    FW_cntcoefp = 32
    FW_basisset = 39 
    FW_ao       = 40

    FR_coefa_ref = 41 
    FW_coefa_ref = 42
    FR_coefb_ref = 43 
    FW_coefb_ref = 44
    FR_coefa_sys = 45 
    FW_coefa_sys = 46
    FR_coefb_sys = 47 
    FW_coefb_sys = 48

    FR_freq     = 61
    FW_freq     = 62
    FR_vibmode  = 71 
    FW_vibmode  = 72
    FR_hessian  = 73 
    FW_hessian  = 74
    FR_micopt   = 75 
    FW_micopt   = 76
                                  
!   ras2sf
    FR_cicoef_g16_cis = 82
    FR_orbnum_alpha = 184
    FR_orbnum_beta = 185
    FW_orbnum = 186
    FR_threshold_values_ci = 191

    FR_ref_singlet_state = 302
    FR_num_ex_states = 303
    FR_ex_energy_as_read = 304
    FW_ex_energy = 305

    FR_g16_cisd_log = 192
    FR_g16_td_log_x = 193
    FR_g16_td_log_y = 194
    FW_xy = 195
    FW_xy_count = 196

    !FW_x_plus_y  = 92
    !FW_x_minus_y = 93
                      
    FW_cicoef_as_read = 201
    FW_cicoef = 202
    FW_vib_fchk = 203

    FW_atmnum_modelsys = 211
    FW_atmwt_modelsys  = 212
    FW_coord_modelsys  = 213
    FW_vibmode_modelsys_calc = 214

    FR_duschinskymat = 221
    FW_duschinskymat = 222
    FR_amat = 223
    FW_amat = 224
    FR_cmat = 225
    FW_cmat = 226
    FR_emat = 227
    FW_emat = 228
    FR_shiftvec = 229
    FW_shiftvec = 230
    FR_bvec = 231
    FW_bvec = 232
    FR_dvec = 233
    FW_dvec = 234



!   ----------
!   File names 
!   ----------

    Max_n_file = 1000
    ALLOCATE(Fname(1:Max_n_file))
    Fname = REPEAT(' ',30)

!   Control files
    Fname(Ffile_typ)       = 'FILE_TYP'
    Fname(FR_control_elec) = 'TEMP_CONTROL_ELEC'
    Fname(FW_control_elec) = 'CONTROL_ELEC'
    Fname(FR_control_vib)  = 'TEMP_CONTROL_VIB'
    Fname(FW_control_vib)  = 'CONTROL_VIB'
  
!   Information for nucleus  
    Fname(FR_atmnum)   = 'TEMP_ATMNUM'
    Fname(FW_atmnum)   = 'ATMNUM'
    Fname(FR_nuccharg) = 'TEMP_NUCCHARG'
    Fname(FW_nuccharg) = 'NUCCHARG'
    Fname(FR_atmwt)    = 'TEMP_ATMWT'
    Fname(FW_atmwt)    = 'ATMWT'
    Fname(FR_coord)    = 'TEMP_COORD'
    Fname(FW_coord)    = 'COORD'
    Fname(FR_nuccharg_soc) = 'TEMP_ZEFF_SOC'
    Fname(FW_nuccharg_soc) = 'ZEFF_SOC'

    Fname(FW_atmnum_modelsys) = 'ATMNUM_MODELSYS'
    Fname(FW_atmwt_modelsys)  = 'ATMWT_MODELSYS'
    Fname(FW_coord_modelsys)  = 'COORD_MODELSYS'
    Fname(FW_vibmode_modelsys_calc)  = 'VIBMODE_MODELSYS'
    Fname(FW_vib_fchk)  = 'VIBMODE_MODELSYS.fchk'

!   Basis functions  
    Fname(FR_shtyp)    = 'TEMP_SHELLTYP'
    Fname(FW_shtyp)    = 'SHELLTYP'
    Fname(FR_n_pgf_sh) = 'TEMP_NUMPGFSHELL'
    Fname(FW_n_pgf_sh) = 'NUMPGFSHELL'
    Fname(FR_sh2atm)   = 'TEMP_SHELL2ATM'
    Fname(FW_sh2atm)   = 'SHELL2ATM'
    Fname(FR_cntexp)   = 'TEMP_CNTEXP'
    Fname(FW_cntexp)   = 'CNTEXP'
    Fname(FR_cntcoef0) = 'TEMP_CNTCOEF0'
    Fname(FW_cntcoef0) = 'CNTCOEF0'
    Fname(FR_cntcoefp) = 'TEMP_CNTCOEFP'
    Fname(FW_cntcoefp) = 'CNTCOEFP'
    Fname(FW_basisset) = 'BASISSET'
    Fname(FW_ao)       = 'ATMORB'
  
!   Molecular orbitals  
    Fname(FR_coefa_ref) = 'TEMP_COEFAREF'
    Fname(FW_coefa_ref) = 'COEFAREF'
    Fname(FR_coefb_ref) = 'TEMP_COEFBREF'
    Fname(FW_coefb_ref) = 'COEFBREF'
    Fname(FR_coefa_sys) = 'TEMP_COEFASYS'
    Fname(FW_coefa_sys) = 'COEFASYS'
    Fname(FR_coefb_sys) = 'TEMP_COEFBSYS'
    Fname(FW_coefb_sys) = 'COEFBSYS'
  
!   Molecular vibrations  
    Fname(FR_freq)      = 'TEMP_FREQ'
    Fname(FW_freq)      = 'FREQ'
    Fname(FR_vibmode)   = 'TEMP_VIBMODE'
    Fname(FW_vibmode)   = 'VIBMODE'
    Fname(FR_hessian)   = 'TEMP_HESSIAN'
    Fname(FW_hessian)   = 'HESSIAN'
    Fname(FR_micopt)    = 'TEMP_MICOPT'
    Fname(FW_micopt)    = 'MICOPT'
  
!   Electronic states   
    Fname(FR_cicoef_g16_cis)  = 'G16_CIS'
    Fname(FR_ref_singlet_state) = 'REF_SINGLET_STATE'
    Fname(FR_num_ex_states) = 'TOTAL_NUM_EX_STATES'
    Fname(FR_orbnum_alpha) = 'ORBNUM_ALPHA'
    Fname(FR_orbnum_beta) = 'ORBNUM_BETA'
    Fname(FW_orbnum) = 'ORBNUM'
    Fname(FR_threshold_values_ci) = 'THRESHOLD_VALUES_CI'
    Fname(FW_cicoef_as_read) = 'CICOEF_AS_READ'
    Fname(FW_cicoef) = 'CICOEF'
    Fname(FR_ex_energy_as_read) = 'EX_ENERGY_AS_READ'
    Fname(FW_ex_energy) = 'EX_ENERGY'

!   g16 excited-states
    Fname(FR_g16_cisd_log) = 'G16_CISD'
    Fname(FR_g16_td_log_x) = 'G16_TD_X'
    Fname(FR_g16_td_log_y) = 'G16_TD_Y'
    Fname(FW_xy) = 'XY'
    Fname(FW_xy_count) = 'XY_COUNT'

!   int_vwf
    Fname(FR_duschinskymat) = 'TEMP_DUSCHINSKYMAT'
    Fname(FW_duschinskymat) = 'DUSCHINSKYMAT'
    Fname(FR_amat) = 'TEMP_AMAT'
    Fname(FW_amat) = 'AMAT'
    Fname(FR_cmat) = 'TEMP_CMAT'
    Fname(FW_cmat) = 'CMAT'
    Fname(FR_emat) = 'TEMP_EMAT'
    Fname(FW_emat) = 'EMAT'
    Fname(FR_shiftvec) = 'TEMP_SHIFTVEC'
    Fname(FW_shiftvec) = 'SHIFTVEC'
    Fname(FR_bvec) = 'TEMP_BVEC'
    Fname(FW_bvec) = 'BVEC'
    Fname(FR_dvec) = 'TEMP_DVEC'
    Fname(FW_dvec) = 'DVEC'



  END SUBROUTINE assign_globals 
END MODULE global_constants







