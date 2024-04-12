! Data preprocessing module of dpp
!
! dpp is a data preprocessing utility for d77.
!
! dpp is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE fchk2inputs
  USE global_constants
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE elec_fchk(qc_program)

    CHARACTER(LEN=10), INTENT(IN) :: qc_program
    INTEGER :: n_atm, n_coord, &
              &charg, mult, ne, nea, neb, &
              &n_basis_0, n_basis, n_sh, sp, n_pgf, &
              &betamo_ref, ncoefa_ref, ncoefb_ref, &
              &i_sh, &
              &shtyp(1:Max_n_cgf), &
              &n_pgf_sh(1:Max_n_cgf), &
              &sh2atm(1:Max_n_cgf), &
              &iao, nao, i_pgf_ao
    DOUBLE PRECISION :: cntexp(1:Max_n_pgf, 1:Max_n_cgf), &
                       &cntcoef0(1:Max_n_pgf, 1:Max_n_cgf), &
                       &cntcoefp(1:Max_n_pgf, 1:Max_n_cgf), &
                       &cntcoef(1:Max_n_pgf, 1:Max_n_cgf)
    CHARACTER(LEN=5) :: orbtyp(1:Max_n_cgf)
    INTEGER :: shtyp_ao(1:Max_n_cgf), n_pgf_ao(1:Max_n_cgf), &
              &ao2atm(1:Max_n_cgf)
    INTEGER :: fid

    shtyp=-10
    shtyp_ao=-10
    n_pgf_sh=0; n_pgf_ao=0; ao2atm=0
    cntcoef0=0.0D0; cntcoefp=0.0D0; cntcoef=0.0D0; cntexp=0.0D0
    nao = 0
  
    n_atm = 0; n_coord = 0
    charg = 0; mult = 0
    ne = 0; nea = 0; neb = 0
    n_basis_0 = 0; n_basis = 0; n_sh = 0; n_pgf = 0; sh2atm = 0
    ncoefa_ref = 0; ncoefb_ref = 0;
    sp = 1; betamo_ref = 1

    fid = FR_control_elec  
    OPEN(fid, FILE = Fname(fid))
      READ(fid,*) 
      READ(fid,*) 
      READ(fid,*) 
      READ(fid,*) 
      READ(fid,*) n_atm
      READ(fid,*) n_coord
      READ(fid,*) charg
      READ(fid,*) mult
      READ(fid,*) ne
      READ(fid,*) nea
      READ(fid,*) neb
      READ(fid,*) n_basis_0
      READ(fid,*) n_basis
      READ(fid,*) n_sh
      READ(fid,*) sp
      READ(fid,*) n_pgf
      READ(fid,*) betamo_ref
      SELECT CASE(betamo_ref)
        CASE(0)
          READ(fid,*) ncoefa_ref
          READ(fid,*) ncoefb_ref
        CASE DEFAULT ! No beta MOs
          READ(fid,*) ncoefa_ref
      END SELECT
    CLOSE(fid)

    fid = FW_control_elec
    OPEN(fid, FILE = Fname(fid))
      WRITE(fid, 1000) n_atm,&
      '! Number of atoms in reference system'
      WRITE(fid, 1000) charg,&
      '! Charge of reference system'
      WRITE(fid, 1000) mult,&
      '! Multiplicity of reference system'
      WRITE(fid, 1000) ne,&
      '! Number of electrons in reference system'
      WRITE(fid, 1000) n_basis_0,&
      '! Number of basis functions'
      WRITE(fid, 1000) n_basis,&
      '! Number of independent basis functions'
    CLOSE(fid)
    1000  FORMAT(I5,2X,A)
  
    CALL scan_data(FR_atmnum, Fname(FR_atmnum),&
                   FW_atmnum, Fname(FW_atmnum),&
                   n_atm, 6, 'INT ')
    CALL scan_data(FR_nuccharg, Fname(FR_nuccharg),&
                   FW_nuccharg, Fname(FW_nuccharg),&
                   n_atm, 5, 'REAL')
    CALL scan_data(FR_coord, Fname(FR_coord),&
                   FW_coord, Fname(FW_coord),&
                   n_coord, 5, 'REAL')
    IF(qc_program .EQ. 'g16') THEN 
      CALL scan_data(FR_nuccharg_soc, Fname(FR_nuccharg_soc),&
                     FW_nuccharg_soc, Fname(FW_nuccharg_soc),&
                     n_atm, 5, 'REAL')
    ELSE
    ENDIF

    CALL scan_data(FR_shtyp, Fname(FR_shtyp),&
                   FW_shtyp, Fname(FW_shtyp),&
                   n_sh, 6, 'INT ')

    CALL scan_data(FR_n_pgf_sh, Fname(FR_n_pgf_sh),&
                   FW_n_pgf_sh, Fname(FW_n_pgf_sh),&
                   n_sh, 6, 'INT ')

    CALL scan_data(FR_sh2atm, Fname(FR_sh2atm),&
                   FW_sh2atm, Fname(FW_sh2atm),&
                   n_sh, 6, 'INT ')

    CALL scan_data(FR_cntexp, Fname(FR_cntexp),&
                   FW_cntexp, Fname(FW_cntexp),&
                   n_pgf, 5, 'REAL')
    CALL scan_data(FR_cntcoef0, Fname(FR_cntcoef0),&
                   FW_cntcoef0, Fname(FW_cntcoef0),&
                   n_pgf, 5, 'REAL')
    SELECT CASE(sp)
      CASE(0)
        CALL scan_data(FR_cntcoefp, Fname(FR_cntcoefp),&
                       FW_cntcoefp, Fname(FW_cntcoefp),&
                       n_pgf, 5, 'REAL')
      CASE DEFAULT
    END SELECT
  
    CALL scan_data(FR_coefa_ref, Fname(FR_coefa_ref),&
                   FW_coefa_ref, Fname(FW_coefa_ref),&
                   ncoefa_ref, 5, 'REAL')
    SELECT CASE(betamo_ref)
      CASE(0)
        CALL scan_data(FR_coefb_ref, Fname(FR_coefb_ref),&
                       FW_coefb_ref, Fname(FW_coefb_ref),&
                       ncoefb_ref, 5, 'REAL')
      CASE DEFAULT
        !WRITE(*,*) 'No beta MOs in reference system'
    END SELECT
  
    OPEN(100, FILE = Fname(FW_shtyp))
    OPEN(110, FILE = Fname(FW_n_pgf_sh))
    OPEN(120, FILE = Fname(FW_sh2atm))
      DO i_sh = 1, n_sh
        READ(100,*) shtyp(i_sh)
        READ(110,*) n_pgf_sh(i_sh)
        READ(120,*) sh2atm(i_sh)
      ENDDO
    CLOSE(120)
    CLOSE(110)
    CLOSE(100)
  
!   Write atomic orbitals
    iao = 0
    SELECT CASE(sp)
      CASE(1) ! sp type = No
        OPEN(201, FILE = Fname(FW_cntexp))
        OPEN(202, FILE = Fname(FW_cntcoef0))
          DO i_sh = 1, n_sh
            SELECT CASE(shtyp(i_sh))
              CASE(0)
                iao = iao + 1
                orbtyp(iao) = 's'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  READ(201,*) cntexp(i_pgf_ao, iao)
                  READ(202,*) cntcoef0(i_pgf_ao, iao)
                ENDDO
              CASE(1)
                iao = iao + 1
                orbtyp(iao) = 'px'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  READ(201,*) cntexp(i_pgf_ao, iao)
                  READ(202,*) cntcoef0(i_pgf_ao, iao)
                ENDDO

                iao = iao + 1
                orbtyp(iao) = 'py'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-1)
                cntcoef0(:,iao) = cntcoef0(:,iao-1)

                iao = iao + 1
                orbtyp(iao) = 'pz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-2)
                cntcoef0(:,iao) = cntcoef0(:,iao-2)
              CASE(2)
                !6D
                iao = iao + 1
                orbtyp(iao) = 'dxx'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  READ(201,*) cntexp(i_pgf_ao, iao)
                  READ(202,*) cntcoef0(i_pgf_ao, iao)
                ENDDO
       
                iao = iao + 1
                orbtyp(iao) = 'dyy'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-1)
                cntcoef0(:,iao) = cntcoef0(:,iao-1)
       
                iao = iao + 1
                orbtyp(iao) = 'dzz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-2)
                cntcoef0(:,iao) = cntcoef0(:,iao-2)
  
                iao = iao + 1
                orbtyp(iao) = 'dxy'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-3)
                cntcoef0(:,iao) = cntcoef0(:,iao-3)
         
                iao = iao + 1
                orbtyp(iao) = 'dxz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-4)
                cntcoef0(:,iao) = cntcoef0(:,iao-4)
  
                iao = iao + 1
                orbtyp(iao) = 'dyz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-5)
                cntcoef0(:,iao) = cntcoef0(:,iao-5)
              CASE(3)
                !6D
                iao = iao + 1
                orbtyp(iao) = 'fxxx'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  READ(201,*) cntexp(i_pgf_ao, iao)
                  READ(202,*) cntcoef0(i_pgf_ao, iao)
                ENDDO

                iao = iao + 1
                orbtyp(iao) = 'fyyy'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-1)
                cntcoef0(:,iao) = cntcoef0(:,iao-1)

                iao = iao + 1
                orbtyp(iao) = 'fzzz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-2)
                cntcoef0(:,iao) = cntcoef0(:,iao-2)

                iao = iao + 1
                orbtyp(iao) = 'fxyy'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-3)
                cntcoef0(:,iao) = cntcoef0(:,iao-3)
                iao = iao + 1

                orbtyp(iao) = 'fxxy'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-4)
                cntcoef0(:,iao) = cntcoef0(:,iao-4)
                
                iao = iao + 1
                orbtyp(iao) = 'fxxz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-5)
                cntcoef0(:,iao) = cntcoef0(:,iao-5)

                iao = iao + 1
                orbtyp(iao) = 'fxzz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-6)
                cntcoef0(:,iao) = cntcoef0(:,iao-6)

                iao = iao + 1
                orbtyp(iao) = 'fyzz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-7)
                cntcoef0(:,iao) = cntcoef0(:,iao-7)

                iao = iao + 1
                orbtyp(iao) = 'fyyz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-8)
                cntcoef0(:,iao) = cntcoef0(:,iao-8)

                iao = iao + 1
                orbtyp(iao) = 'fxyz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-9)
                cntcoef0(:,iao) = cntcoef0(:,iao-9)

              CASE DEFAULT
                WRITE(*,*) 'Incorrect shell type detected'
                STOP
            END SELECT
          ENDDO
        CLOSE(202)
        CLOSE(201)
      CASE(0) ! sp type = yes
        OPEN(201, FILE = Fname(FW_cntexp))
        OPEN(202, FILE = Fname(FW_cntcoef0))
        OPEN(203, FILE = Fname(FW_cntcoefp))
          DO i_sh = 1, n_sh
            SELECT CASE(shtyp(i_sh))
              CASE(0)
                iao = iao + 1
                orbtyp(iao) = 's'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  READ(201,*) cntexp(i_pgf_ao, iao)
                  READ(202,*) cntcoef0(i_pgf_ao, iao)
                  READ(203,*) cntcoefp(i_pgf_ao, iao)
                ENDDO
              CASE(-1)
                iao = iao + 1
                orbtyp(iao) = 's'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  READ(201,*) cntexp(i_pgf_ao, iao)
                  READ(202,*) cntcoef0(i_pgf_ao, iao)
                  READ(203,*) cntcoefp(i_pgf_ao, iao)
                ENDDO
       
                iao = iao + 1
                orbtyp(iao) = 'px'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-1)
                cntcoef0(:,iao) = cntcoef0(:,iao-1)
                cntcoefp(:,iao) = cntcoefp(:,iao-1)
       
                iao = iao + 1
                orbtyp(iao) = 'py'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-2)
                cntcoef0(:,iao) = cntcoef0(:,iao-2)
                cntcoefp(:,iao) = cntcoefp(:,iao-2)
  
                iao = iao + 1
                orbtyp(iao) = 'pz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-3)
                cntcoef0(:,iao) = cntcoef0(:,iao-3)
                cntcoefp(:,iao) = cntcoefp(:,iao-3)
              CASE(1)
                iao = iao + 1
                orbtyp(iao) = 'px'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  READ(201,*) cntexp(i_pgf_ao, iao)
                  READ(202,*) cntcoef0(i_pgf_ao, iao)
                  READ(203,*) cntcoefp(i_pgf_ao, iao) ! 2021/1/18
                ENDDO

                iao = iao + 1
                orbtyp(iao) = 'py'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-1)
                cntcoef0(:,iao) = cntcoef0(:,iao-1)

                iao = iao + 1
                orbtyp(iao) = 'pz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-2)
                cntcoef0(:,iao) = cntcoef0(:,iao-2)
              CASE(2)
                !6D
                iao = iao + 1
                orbtyp(iao) = 'dxx'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  READ(201,*) cntexp(i_pgf_ao, iao)
                  READ(202,*) cntcoef0(i_pgf_ao, iao)
                  READ(203,*) cntcoefp(i_pgf_ao, iao)
                ENDDO
  
                iao = iao + 1
                orbtyp(iao) = 'dyy'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-1)
                cntcoef0(:,iao) = cntcoef0(:,iao-1)
                cntcoefp(:,iao) = cntcoefp(:,iao-1)
           
                iao = iao + 1
                orbtyp(iao) = 'dzz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-2)
                cntcoef0(:,iao) = cntcoef0(:,iao-2)
                cntcoefp(:,iao) = cntcoefp(:,iao-2)
           
                iao = iao + 1
                orbtyp(iao) = 'dxy'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-3)
                cntcoef0(:,iao) = cntcoef0(:,iao-3)
                cntcoefp(:,iao) = cntcoefp(:,iao-3)
  
                iao = iao + 1
                orbtyp(iao) = 'dxz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-4)
                cntcoef0(:,iao) = cntcoef0(:,iao-4)
                cntcoefp(:,iao) = cntcoefp(:,iao-4)
             
                iao = iao + 1
                orbtyp(iao) = 'dyz'
                shtyp_ao(iao) = shtyp(i_sh)
                n_pgf_ao(iao) = n_pgf_sh(i_sh)
                ao2atm(iao) = sh2atm(i_sh)
                cntexp(:,iao) = cntexp(:,iao-5)
                cntcoef0(:,iao) = cntcoef0(:,iao-5)
                cntcoefp(:,iao) = cntcoefp(:,iao-5)
              CASE DEFAULT
                WRITE(*,*) 'Incorrect shell type detected'
                WRITE(*,*) i_sh, shtyp_ao(iao)
                STOP
            END SELECT
            !WRITE(*,*) i_sh, shtyp_ao(iao)
          ENDDO
        CLOSE(203)
        CLOSE(202)
        CLOSE(201)
      CASE DEFAULT
        STOP
    END SELECT
    nao = iao
    IF(nao/=n_basis) THEN
      !WRITE(*,*) 'Error: incorrect basis set'
      WRITE(*,*) 'The number of AOs is not equal to'
      WRITE(*,*) 'the number of basis functions in FCHK.'
      WRITE(*,*) 'The number of basis functions in FCHK = ', n_basis
      WRITE(*,*) 'The number of AOs = ', nao
      !STOP
    ENDIF

    OPEN(FW_ao, FILE = Fname(FW_ao))
      DO iao = 1, nao
        SELECT CASE(shtyp_ao(iao))
          CASE(0) ! S orbital
            WRITE(FW_ao,9001) iao, '! iao'
            WRITE(FW_ao,9001) ao2atm(iao), '! i_atm'
            WRITE(FW_ao,9002) orbtyp(iao), n_pgf_ao(iao)
            DO i_pgf_ao = 1, n_pgf_ao(iao)
              WRITE(FW_ao,9102) cntexp(i_pgf_ao, iao),&
                               cntcoef0(i_pgf_ao, iao)
            ENDDO
          CASE(1)
            SELECT CASE(orbtyp(iao))
              CASE('px', 'py', 'pz') ! P orbital
                WRITE(FW_ao,9001) iao, '! iao'
                WRITE(FW_ao,9001) ao2atm(iao), '! i_atm'
                WRITE(FW_ao,9002) orbtyp(iao), n_pgf_ao(iao)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  WRITE(FW_ao,9102) cntexp(i_pgf_ao, iao),&
                                   cntcoef0(i_pgf_ao, iao)
                ENDDO
              CASE DEFAULT
                WRITE(*,*) 'Error'
                STOP
            END SELECT
          CASE(-1)
            SELECT CASE(orbtyp(iao))
              CASE('s') ! S orbital
                WRITE(FW_ao,9001) iao, '! iao'
                WRITE(FW_ao,9001) ao2atm(iao), '! i_atm'
                WRITE(FW_ao,9002) orbtyp(iao), n_pgf_ao(iao)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  WRITE(FW_ao,9102) cntexp(i_pgf_ao, iao),&
                                   cntcoef0(i_pgf_ao, iao)
                ENDDO
              CASE('px', 'py', 'pz') ! P orbital
                WRITE(FW_ao,9001) iao, '! iao'
                WRITE(FW_ao,9001) ao2atm(iao), '! i_atm'
                WRITE(FW_ao,9002) orbtyp(iao), n_pgf_ao(iao)
                DO i_pgf_ao = 1, n_pgf_ao(iao)
                  WRITE(FW_ao,9102) cntexp(i_pgf_ao, iao),&
                                   cntcoefp(i_pgf_ao, iao)
                                   !cntcoef0(i_pgf_ao, iao)
                ENDDO
              CASE DEFAULT
                WRITE(*,*) 'Error'
                STOP
            END SELECT
          CASE(2) ! D orbital
            WRITE(FW_ao,9001) iao, '! iao'
            WRITE(FW_ao,9001) ao2atm(iao), '! i_atm'
            WRITE(FW_ao,9002) orbtyp(iao), n_pgf_ao(iao)
            DO i_pgf_ao = 1, n_pgf_ao(iao)
              WRITE(FW_ao,9102) cntexp(i_pgf_ao, iao),&
                               cntcoef0(i_pgf_ao, iao)
            ENDDO
          CASE(3) ! F orbital
            WRITE(FW_ao,9001) iao, '! iao'
            WRITE(FW_ao,9001) ao2atm(iao), '! i_atm'
            WRITE(FW_ao,9002) orbtyp(iao), n_pgf_ao(iao)
            DO i_pgf_ao = 1, n_pgf_ao(iao)
              WRITE(FW_ao,9102) cntexp(i_pgf_ao, iao),&
                               cntcoef0(i_pgf_ao, iao)
            ENDDO
          CASE DEFAULT
            WRITE(*,*) 'Error'
            STOP
        END SELECT
      ENDDO
    CLOSE(FW_ao)
    9001 FORMAT(I0,2X,A)
    9002 FORMAT(A, I0)
    9102 FORMAT(2E18.10)
  
  END SUBROUTINE elec_fchk



  SUBROUTINE vib_fchk

    INTEGER n_atm
    INTEGER n_mode, n_coord_mode
  
    !n_atm = 0
  
    OPEN(FR_control_vib, FILE = Fname(FR_control_vib))
      READ(FR_control_vib,*) 
      READ(FR_control_vib,*) n_atm
      READ(FR_control_vib,*) n_mode
      READ(FR_control_vib,*) n_coord_mode
    CLOSE(FR_control_vib)
  
    OPEN(FW_control_vib, FILE = Fname(FW_control_vib))
      !WRITE(FW_control_vib,*) n_atm
      !WRITE(FW_control_vib,*) n_atm ! Added at 2020/4/16
      WRITE(FW_control_vib,*) n_mode
    CLOSE(FW_control_vib)
 
!   Vibrational modes
    CALL scan_data(FR_atmwt, Fname(FR_atmwt),&
                   FW_atmwt, Fname(FW_atmwt),&
                   n_atm, 5, 'REAL')
    CALL scan_data(FR_freq, Fname(FR_freq),&
                   FW_freq, Fname(FW_freq),&
                   n_mode, 5, 'REAL')
    CALL scan_data(FR_vibmode, Fname(FR_vibmode),&
                   FW_vibmode, Fname(FW_vibmode),&
                   n_coord_mode, 5, 'REAL')
  END SUBROUTINE vib_fchk


      SUBROUTINE read_atmwt(n_atm, atmnum, unit_mass, atmwt)

        INTEGER, INTENT(IN) :: n_atm
        INTEGER, INTENT(IN) :: atmnum(1:n_atm) 
        CHARACTER(5), INTENT(IN) :: unit_mass
        DOUBLE PRECISION, INTENT(OUT) :: atmwt(1:n_atm)

        INTEGER i_atm
        DOUBLE PRECISION atmwt_amu(1:n_atm)

        DO i_atm = 1, n_atm
          SELECT CASE(atmnum(i_atm))
            CASE(1);  atmwt_amu(i_atm) = 1.007825D0
            CASE(2);  atmwt_amu(i_atm) = 4.00260D0
            CASE(3);  atmwt_amu(i_atm) = 7.016003D0
            CASE(4);  atmwt_amu(i_atm) = 9.012182D0
            CASE(5);  atmwt_amu(i_atm) = 11.009305D0
            CASE(6);  atmwt_amu(i_atm) = 12.000000D0
            CASE(7);  atmwt_amu(i_atm) = 14.003074D0
            CASE(8);  atmwt_amu(i_atm) = 15.994915D0
            CASE(9);  atmwt_amu(i_atm) = 18.9984032D0
            CASE(10); atmwt_amu(i_atm) = 19.992435D0
            CASE(11); atmwt_amu(i_atm) = 22.989767D0
            CASE(12); atmwt_amu(i_atm) = 23.985042D0
            CASE(13); atmwt_amu(i_atm) = 26.981540D0
            CASE(14); atmwt_amu(i_atm) = 27.976927D0
            CASE(15); atmwt_amu(i_atm) = 30.973762D0
            CASE(16); atmwt_amu(i_atm) = 31.972070D0
            CASE(17); atmwt_amu(i_atm) = 34.968852D0
            CASE(18); atmwt_amu(i_atm) = 39.962384D0
            CASE(19); atmwt_amu(i_atm) = 38.963707D0
            CASE(20); atmwt_amu(i_atm) = 39.962591D0
            CASE(21); atmwt_amu(i_atm) = 44.955910D0
            CASE(22); atmwt_amu(i_atm) = 47.947947D0
            CASE(23); atmwt_amu(i_atm) = 50.943962D0
            CASE(24); atmwt_amu(i_atm) = 51.940509D0
            CASE(25); atmwt_amu(i_atm) = 54.938047D0
            CASE(26); atmwt_amu(i_atm) = 55.934939D0
            CASE(27); atmwt_amu(i_atm) = 58.933198D0
            CASE(28); atmwt_amu(i_atm) = 57.935346D0
            CASE(29); atmwt_amu(i_atm) = 62.939598D0
            CASE(30); atmwt_amu(i_atm) = 63.929145D0
            CASE(31); atmwt_amu(i_atm) = 68.925580D0
            CASE(32); atmwt_amu(i_atm) = 73.921177D0
            CASE(33); atmwt_amu(i_atm) = 74.921594D0
            CASE(34); atmwt_amu(i_atm) = 79.916520D0
            CASE(35); atmwt_amu(i_atm) = 78.918336D0
            CASE(36); atmwt_amu(i_atm) = 83.911507D0
            CASE DEFAULT
              WRITE(*,*) 'Invalid atomic number'
              WRITE(*,*) 'Atomic number = ', atmnum(i_atm)
              STOP
          END SELECT
        ENDDO
        
        SELECT CASE(unit_mass)
          CASE('amu')
            atmwt = atmwt_amu ! [u]
          CASE('au')
            ! atmwt [u]
            ! atmwt*Amc [kg]
            ! atmwt*Amc/Me [a.u.]
            atmwt = atmwt_amu*Amc/Me ! [a.u.]
          CASE DEFAULT
            WRITE(6,*) 'Invalid unit of mass'
            STOP
        END SELECT

      END SUBROUTINE read_atmwt


      SUBROUTINE scan_data(fid_in,  fname_in,&
                           fid_out, fname_out,&
                           n_elm, n_col, real_or_int)

        INTEGER, INTENT(IN) :: fid_in, fid_out
        CHARACTER(30), INTENT(IN) :: fname_in, fname_out
        INTEGER, INTENT(IN) :: n_elm, n_col
        CHARACTER(4), INTENT(IN) :: real_or_int
        INTEGER n_line, n_col_final
        INTEGER i_line, i_col
        DOUBLE PRECISION dummy_real(1:n_col)
        INTEGER dummy_int(1:n_col)

        n_line = 0; n_col_final = 0
        dummy_real = 0.0D0
        dummy_int = 0

        IF(n_elm <= 0) THEN
          WRITE(*,*) 'Error in scan_inp'
          WRITE(*,*) 'The number of elements is less than 0.'
          WRITE(*,*) 'n_elm = ', n_elm
          STOP
        ELSE  
        ENDIF

        OPEN(fid_in,  FILE = fname_in)
        OPEN(fid_out, FILE = fname_out)
          SELECT CASE(real_or_int)
            CASE('REAL')
              IF(MOD(n_elm, n_col) == 0) THEN
                n_line = n_elm/n_col
                DO i_line = 1, n_line
                  READ(fid_in,*) dummy_real(1:n_col)
                  DO i_col = 1, n_col
                    WRITE(fid_out,*) dummy_real(i_col)
                  ENDDO
                ENDDO
              ELSE
                IF(n_elm < n_col) THEN
                  READ(fid_in,*) dummy_real(1:n_elm)
                  DO i_col = 1, n_elm
                    WRITE(fid_out,*) dummy_real(i_col)
                  ENDDO
                ELSE
                  n_col_final = MOD(n_elm, n_col)
                  n_line = (n_elm - MOD(n_elm, n_col))/n_col + 1
                  DO i_line = 1, n_line - 1
                    READ(fid_in,*) dummy_real(1:n_col)
                    DO i_col = 1, n_col
                      WRITE(fid_out,*) dummy_real(i_col)
                    ENDDO
                  ENDDO
                  READ(fid_in,*) dummy_real(1:n_col_final)
                  DO i_col = 1, n_col_final
                    WRITE(fid_out,*) dummy_real(i_col)
                  ENDDO
                ENDIF
              ENDIF
            CASE('INT ')
              IF(MOD(n_elm, n_col) == 0) THEN
                n_line = n_elm/n_col
                DO i_line = 1, n_line
                  READ(fid_in,*) dummy_int(1:n_col)
                  DO i_col = 1, n_col
                    WRITE(fid_out,*) dummy_int(i_col)
                  ENDDO
                ENDDO
              ELSE
                IF(n_elm < n_col) THEN
                  READ(fid_in,*) dummy_int(1:n_elm)
                  DO i_col = 1, n_elm
                    WRITE(fid_out,*) dummy_int(i_col)
                  ENDDO
                ELSE
                  n_col_final = MOD(n_elm, n_col)
                  n_line = (n_elm - MOD(n_elm, n_col))/n_col + 1
                  DO i_line = 1, n_line - 1
                    READ(fid_in,*) dummy_int(1:n_col)
                    DO i_col = 1, n_col
                      WRITE(fid_out,*) dummy_int(i_col)
                    ENDDO
                  ENDDO
                  READ(fid_in,*) dummy_int(1:n_col_final)
                  DO i_col = 1, n_col_final
                    WRITE(fid_out,*) dummy_int(i_col)
                  ENDDO
                ENDIF
              ENDIF
            CASE DEFAULT
              STOP
          END SELECT
        CLOSE(fid_out)
        CLOSE(fid_in)
      END SUBROUTINE scan_data

      SUBROUTINE scan_data2(fid_in,  fname_in,&
                            fid_out, fname_out,&
                            n_elm, n_col, real_or_int)

        INTEGER, INTENT(IN) :: fid_in, fid_out
        CHARACTER(30), INTENT(IN) :: fname_in, fname_out
        INTEGER, INTENT(IN) :: n_elm, n_col
        CHARACTER(4), INTENT(IN) :: real_or_int
        INTEGER n_block, n_col_final
        INTEGER i_block, i_col
        INTEGER i_elm
        DOUBLE PRECISION dummy_real(1:n_elm, 1:n_col)
        INTEGER dummy_int(1:n_col)

        n_block = 0; n_col_final = 0
        dummy_real = 0.0D0
        dummy_int = 0

        IF(n_elm <= 0) THEN
          WRITE(*,*) 'Error in scan_inp'
          WRITE(*,*) 'The number of elements is less than 0.'
          WRITE(*,*) 'n_elm = ', n_elm
          STOP
        ELSE
        ENDIF

        OPEN(fid_in,  FILE = fname_in)
        OPEN(fid_out, FILE = fname_out)
          SELECT CASE(real_or_int)
            CASE('REAL')
              IF(MOD(n_elm, n_col) == 0) THEN
                n_block = n_elm/n_col
                DO i_block = 1, n_block
                  DO i_elm = 1, n_elm
                    READ(fid_in,*) (dummy_real(i_elm, i_col), i_col = 1, n_col)
                  ENDDO
                  DO i_col = 1, n_col
                    DO i_elm = 1, n_elm
                      WRITE(fid_out,*) dummy_real(i_elm, i_col)
                    ENDDO
                  ENDDO
                ENDDO
              ELSE
                IF(n_elm < n_col) THEN
                  DO i_elm = 1, n_elm
                    READ(fid_in,*) (dummy_real(i_elm, i_col), i_col = 1, n_elm)
                  ENDDO
                  DO i_col = 1, n_elm
                    DO i_elm = 1, n_elm
                      WRITE(fid_out,*) dummy_real(i_elm, i_col)
                    ENDDO
                  ENDDO
                ELSE
                  n_col_final = MOD(n_elm, n_col)
                  n_block = (n_elm - MOD(n_elm, n_col))/n_col + 1
                  DO i_block = 1, n_block - 1
                    DO i_elm = 1, n_elm
                      READ(fid_in,*) (dummy_real(i_elm, i_col), i_col = 1, n_col)
                    ENDDO
                    DO i_col = 1, n_col
                      DO i_elm = 1, n_elm
                        WRITE(fid_out,*) dummy_real(i_elm, i_col)
                      ENDDO
                    ENDDO
                  ENDDO
                  DO i_elm = 1, n_elm
                    READ(fid_in,*) (dummy_real(i_elm, i_col), i_col = 1, n_col_final)
                  ENDDO
                  DO i_col = 1, n_col_final
                    DO i_elm = 1, n_elm
                      WRITE(fid_out,*) dummy_real(i_elm, i_col)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
            CASE('INT ')
            CASE DEFAULT
              STOP
          END SELECT
        CLOSE(fid_out)
        CLOSE(fid_in)

      END SUBROUTINE scan_data2



      SUBROUTINE gen_vib_fchk(&
                 n_atm, charg, mult, ne, nea, neb, n_basis,&
                 xyznuc_dummy, atmnum, nuccharg,&
                 n_mode, n_hessian_cart,&
                 hessian_cart_dummy, freq, atmwt_amu,&
                 n_mode_coord, vibmode_cart_dummy)

        INTEGER, INTENT(IN) :: n_atm, charg, mult, ne, nea, neb, n_basis
        DOUBLE PRECISION, INTENT(IN) :: xyznuc_dummy(1:3*n_atm)
        INTEGER, INTENT(IN) :: atmnum(1:n_atm)
        DOUBLE PRECISION, INTENT(IN) :: nuccharg(1:n_atm)
        INTEGER, INTENT(IN) :: n_mode, n_hessian_cart
        DOUBLE PRECISION, INTENT(IN) :: hessian_cart_dummy&
                                        (1:n_hessian_cart) 
        DOUBLE PRECISION, INTENT(IN) :: atmwt_amu(1:n_atm)
        DOUBLE PRECISION, INTENT(IN) :: freq(1:n_mode)
        INTEGER, INTENT(IN) :: n_mode_coord
        DOUBLE PRECISION, INTENT(IN) :: vibmode_cart_dummy&
                                        (1:n_mode_coord)

        ! Local variables
        INTEGER n_elm_Vib_E2
        DOUBLE PRECISION, ALLOCATABLE :: dummy_real(:) 

        n_elm_Vib_E2 = 0

        OPEN(FW_vib_fchk, FILE = Fname(FW_vib_fchk))

          ! Write comment lines
          WRITE(FW_vib_fchk,'(1X, A)')&
          'fchk file generated using d77 prepost program    '
          WRITE(FW_vib_fchk,'(1X, A)')&
          'Freq      RHF                                    '

          ! Write number of atoms
          WRITE(FW_vib_fchk,9001)&
          'Number of atoms                            I     ',&
          n_atm       

          ! Write charge
          WRITE(FW_vib_fchk,9001)&
          'Charge                                     I     ',&
          charg       

          ! Write multiplicity
          WRITE(FW_vib_fchk,9001)&
          'Multiplicity                               I     ',&
          mult       

          ! Write number of electrons
          WRITE(FW_vib_fchk,9001)&
          'Number of electrons                        I     ',&
          ne

          ! Write number of alpha electrons
          WRITE(FW_vib_fchk,9001)&
          'Number of alpha electrons                  I     ',&
          nea       

          ! Write number of beta electrons
          WRITE(FW_vib_fchk,9001)&
          'Number of beta electrons                   I     ',&
          neb

          ! Write number of basis functions
          WRITE(FW_vib_fchk,9001)&
          'Number of basis functions                  I     ',&
          n_basis

          ! Write number of basis functions
          WRITE(FW_vib_fchk,9001)&
          'Number of independent functions            I     ',&
          n_basis

          ! Write number of point charges
          WRITE(FW_vib_fchk,9001)&
          'Number of point charges in /Mol/           I     ',&
          n_atm

          ! Write number of translation vectors
          WRITE(FW_vib_fchk,9001)&
          'Number of translation vectors              I     ',&
          0

          ! Write atomic numbers 
          WRITE(FW_vib_fchk,9001)&
          'Atomic numbers                             I   N=',&
          n_atm
          CALL write_data_int(FW_vib_fchk, n_atm, atmnum, 6)

          ! Write nuclear charges
          WRITE(FW_vib_fchk,9001)&
          'Nuclear charges                            R   N=',&
          n_atm
          CALL write_data_real(FW_vib_fchk, n_atm, nuccharg, 5)


          ! Write nuclear cooedinates
          WRITE(FW_vib_fchk,9001)&
          'Current cartesian coordinates              R   N=',&
          3*n_atm       
          CALL write_data_real(FW_vib_fchk, 3*n_atm,&
                               xyznuc_dummy(1:3*n_atm), 5)

          ! Write Cartesian gradient   
          WRITE(FW_vib_fchk,9001)&
          'Cartesian Gradient                         R   N=',& 
          n_mode
          ALLOCATE(dummy_real(1:n_mode))
          dummy_real = 0.0D0
          dummy_real(1:n_mode) = 0.0D0
          CALL write_data_real(FW_vib_fchk, n_mode, dummy_real, 5)
          DEALLOCATE(dummy_real)

          ! Write Cartesian force constant
          WRITE(FW_vib_fchk,9001)&
          'Cartesian Force Constants                  R   N=',& 
          n_hessian_cart
          CALL write_data_real(FW_vib_fchk, n_hessian_cart,&
                               hessian_Cart_dummy,5)

          ! Write number of vibrational modes
          WRITE(FW_vib_fchk,9001)&
          'Number of Normal Modes                     I     ',&
          n_mode
          WRITE(FW_vib_fchk,9001)&
          'Vib-NDFDPol                                I     ',&
          0
          WRITE(FW_vib_fchk,9001)&
          'Vib-NumROA                                 I     ',&
          0
          WRITE(FW_vib_fchk,9001)&
          'Vib-NDim                                   I     ',&
          6
          WRITE(FW_vib_fchk,9001)&
          'Vib-NDim0                                  I     ',&
          6

          ! Write Vib-AtMass
          WRITE(FW_vib_fchk,9001)&
          'Vib-AtMass                                 R   N=',&
          n_atm
          CALL write_data_real(FW_vib_fchk, n_atm, atmwt_amu, 5)
          ! Write freq
          n_elm_Vib_E2 = n_mode * 14
          WRITE(FW_vib_fchk,9001)&
          'Vib-E2                                     R   N=',&
          n_elm_Vib_E2
          ALLOCATE(dummy_real(1:n_elm_Vib_E2))
          dummy_real = 0.0D0
          dummy_real(1:n_mode) = freq
          CALL write_data_real(FW_vib_fchk, n_elm_Vib_E2, dummy_real,&
                          5)
          DEALLOCATE(dummy_real)


          WRITE(FW_vib_fchk,9001)&
          'Vib-Modes                                  R   N=',&
          n_mode_coord
          CALL write_data_real&
               (FW_vib_fchk, n_mode_coord, vibmode_cart_dummy, 5)

        CLOSE(FW_vib_fchk)
  9001 FORMAT(A, I12)
      END SUBROUTINE gen_vib_fchk

  SUBROUTINE gen_fchk_rot_vib(&
             n_atm, xyznuc_dummy, atmnum, nuccharg,&
             n_mode, n_hessian_cart, freq, atmwt_amu,&
             n_mode_coord, vibmode_cart_dummy)

    INTEGER, INTENT(IN) :: n_atm
    INTEGER :: charg, mult, ne, nea, neb, n_basis
    DOUBLE PRECISION, INTENT(IN) :: xyznuc_dummy(1:3*n_atm)
    INTEGER, INTENT(IN) :: atmnum(1:n_atm)
    DOUBLE PRECISION, INTENT(IN) :: nuccharg(1:n_atm)
    INTEGER, INTENT(IN) :: n_mode
    INTEGER :: n_hessian_cart
    DOUBLE PRECISION :: hessian_cart_dummy&
                        (1:n_hessian_cart)
    DOUBLE PRECISION, INTENT(IN) :: atmwt_amu(1:n_atm)
    DOUBLE PRECISION, INTENT(IN) :: freq(1:n_mode)
    INTEGER, INTENT(IN) :: n_mode_coord
    DOUBLE PRECISION, INTENT(IN) :: vibmode_cart_dummy&
                                    (1:n_mode_coord)

    ! Local variables
    INTEGER n_elm_Vib_E2
    DOUBLE PRECISION, ALLOCATABLE :: dummy_real(:)
    INTEGER fid

    n_elm_Vib_E2 = 0
    charg = 0; mult = 1; ne = 20; nea = 10; neb = 20; n_basis = 100

    fid = FW_vib_fchk
    OPEN(fid, FILE = Fname(fid))

      ! Write comment lines
      WRITE(fid,'(1X, A)')&
      'fchk file generated using d77 prepost program    '
      WRITE(fid,'(1X, A)')&
      'Freq      RHF                                    '

      ! Write number of atoms
      WRITE(fid,9001)&
      'Number of atoms                            I     ',&
      n_atm

      ! Write charge
      WRITE(fid,9001)&
      'Charge                                     I     ',&
      charg

      ! Write multiplicity
      WRITE(fid,9001)&
      'Multiplicity                               I     ',&
      mult

      ! Write number of electrons
      WRITE(fid,9001)&
      'Number of electrons                        I     ',&
      ne

      ! Write number of alpha electrons
      WRITE(fid,9001)&
      'Number of alpha electrons                  I     ',&
      nea

      ! Write number of beta electrons
      WRITE(fid,9001)&
      'Number of beta electrons                   I     ',&
      neb

      ! Write number of basis functions
      WRITE(fid,9001)&
      'Number of basis functions                  I     ',&
      n_basis

      ! Write number of basis functions
      WRITE(fid,9001)&
      'Number of independent functions            I     ',&
      n_basis

      ! Write number of point charges
      WRITE(fid,9001)&
      'Number of point charges in /Mol/           I     ',&
      n_atm

      ! Write number of translation vectors
      WRITE(fid,9001)&
      'Number of translation vectors              I     ',&
      0

      ! Write atomic numbers
      WRITE(fid,9001)&
      'Atomic numbers                             I   N=',&
      n_atm
      CALL write_data_int(fid, n_atm, atmnum, 6)

      ! Write nuclear charges
      WRITE(fid,9001)&
      'Nuclear charges                            R   N=',&
      n_atm
      CALL write_data_real(fid, n_atm, nuccharg, 5)


      ! Write nuclear cooedinates
      WRITE(fid,9001)&
      'Current cartesian coordinates              R   N=',&
      3*n_atm
      CALL write_data_real(fid, 3*n_atm,&
                           xyznuc_dummy(1:3*n_atm), 5)

      ! Write Cartesian gradient
      WRITE(fid,9001)&
      'Cartesian Gradient                         R   N=',&
      3*n_atm
      !n_mode
      ALLOCATE(dummy_real(1:n_mode))
      dummy_real = 0.0D0
      dummy_real(1:n_mode) = 0.0D0
      CALL write_data_real(fid, n_mode, dummy_real, 5)
      DEALLOCATE(dummy_real)

      ! Write Cartesian force constant
      hessian_cart_dummy=0.0D0
      WRITE(fid,9001)&
      'Cartesian Force Constants                  R   N=',&
      n_hessian_cart
      CALL write_data_real(fid, n_hessian_cart,&
                           hessian_Cart_dummy,5)

      ! Write number of vibrational modes
      WRITE(fid,9001)&
      'Number of Normal Modes                     I     ',&
      n_mode
      WRITE(fid,9001)&
      'Vib-NDFDPol                                I     ',&
      0
      WRITE(fid,9001)&
      'Vib-NumROA                                 I     ',&
      0
      WRITE(fid,9001)&
      'Vib-NDim                                   I     ',&
      6
      WRITE(fid,9001)&
      'Vib-NDim0                                  I     ',&
      6

      ! Write Vib-AtMass
      WRITE(fid,9001)&
      'Vib-AtMass                                 R   N=',&
      n_atm
      CALL write_data_real(fid, n_atm, atmwt_amu, 5)
      ! Write freq
      n_elm_Vib_E2 = n_mode * 14
      WRITE(fid,9001)&
      'Vib-E2                                     R   N=',&
      n_elm_Vib_E2
      ALLOCATE(dummy_real(1:n_elm_Vib_E2))
      dummy_real = 0.0D0
      dummy_real(1:n_mode) = freq
      CALL write_data_real(fid, n_elm_Vib_E2, dummy_real,&
                      5)
      DEALLOCATE(dummy_real)


      WRITE(fid,9001)&
      'Vib-Modes                                  R   N=',&
      n_mode_coord
      CALL write_data_real&
           (fid, n_mode_coord, vibmode_cart_dummy, 5)

    CLOSE(fid)

    9001 FORMAT(A, I12)
  END SUBROUTINE gen_fchk_rot_vib



      SUBROUTINE write_data_real(fid_out, n_elm, data_in,&
                                 n_col)

        INTEGER, INTENT(IN) :: fid_out
        INTEGER, INTENT(IN) :: n_elm, n_col
        DOUBLE PRECISION, INTENT(IN) :: data_in(1:n_elm)
        INTEGER n_line, n_col_final
        INTEGER i_line!, i_col
        INTEGER counter 

        n_line = 0; n_col_final = 0
        counter = 0
        
        IF(MOD(n_elm, n_col) == 0) THEN
          n_line = n_elm/n_col
          DO i_line = 1, n_line
            counter = (i_line-1) * n_col
            SELECT CASE(n_col)
              CASE(5)
                WRITE(fid_out,9015)& 
                data_in(counter + 1 : counter + n_col)
              CASE(6)
                WRITE(fid_out,9016)& 
                data_in(counter + 1 : counter + n_col)
              CASE DEFAULT
                STOP
            END SELECT
          ENDDO
        ELSE
          IF(n_elm < n_col) THEN
            SELECT CASE(n_elm)
              CASE(1)
                WRITE(fid_out,9011) data_in(1:n_elm)
              CASE(2)
                WRITE(fid_out,9012) data_in(1:n_elm)
              CASE(3)
                WRITE(fid_out,9013) data_in(1:n_elm)
              CASE(4)
                WRITE(fid_out,9014) data_in(1:n_elm)
              CASE(5)
                WRITE(fid_out,9015) data_in(1:n_elm)
              CASE(6)
                WRITE(fid_out,9016) data_in(1:n_elm)
              CASE DEFAULT
                STOP
            END SELECT
          ELSE
            n_col_final = MOD(n_elm, n_col)
            n_line = (n_elm - MOD(n_elm, n_col))/n_col + 1
            DO i_line = 1, n_line - 1
              counter = (i_line-1) * n_col
              SELECT CASE(n_col)
                CASE(5)
                  WRITE(fid_out,9015)& 
                  data_in(counter + 1 : counter + n_col)
                CASE(6)
                  WRITE(fid_out,9016)& 
                  data_in(counter + 1 : counter + n_col)
                CASE DEFAULT
                  STOP
              END SELECT
            ENDDO
            counter = (n_line-1) * n_col
            SELECT CASE(n_col_final)
              CASE(1)
                WRITE(fid_out,9011)& 
                data_in(counter + 1 : counter + n_col_final)
              CASE(2)
                WRITE(fid_out,9012)& 
                data_in(counter + 1 : counter + n_col_final)
              CASE(3)
                WRITE(fid_out,9013)&
                data_in(counter + 1 : counter + n_col_final)
              CASE(4)
                WRITE(fid_out,9014)&
                data_in(counter + 1 : counter + n_col_final)
              CASE(5)
                WRITE(fid_out,9015)&
                data_in(counter + 1 : counter + n_col_final)
              CASE(6)
                WRITE(fid_out,9016)&
                data_in(counter + 1 : counter + n_col_final)
              CASE DEFAULT
                STOP
            END SELECT
          ENDIF
        ENDIF

9011    FORMAT(E16.8)
9012    FORMAT(2E16.8)
9013    FORMAT(3E16.8)
9014    FORMAT(4E16.8)
9015    FORMAT(5E16.8)
9016    FORMAT(6E16.8)
      END SUBROUTINE write_data_real

      SUBROUTINE write_data_int(fid_out, n_elm, data_in,&
                                n_col)

        INTEGER, INTENT(IN) :: fid_out
        INTEGER, INTENT(IN) :: n_elm, n_col
        INTEGER, INTENT(IN) :: data_in(1:n_elm)
        INTEGER n_line, n_col_final
        INTEGER i_line!, i_col
        INTEGER counter 

        n_line = 0; n_col_final = 0
        counter = 0
        
        IF(MOD(n_elm, n_col) == 0) THEN
          n_line = n_elm/n_col
          DO i_line = 1, n_line
            counter = (i_line-1) * n_col
            SELECT CASE(n_col)
              CASE(5)
                WRITE(fid_out,9015)& 
                data_in(counter + 1 : counter + n_col)
              CASE(6)
                WRITE(fid_out,9016)& 
                data_in(counter + 1 : counter + n_col)
              CASE DEFAULT
                STOP
            END SELECT
          ENDDO
        ELSE
          IF(n_elm < n_col) THEN
            SELECT CASE(n_elm)
              CASE(1)
                WRITE(fid_out,9011) data_in(1:n_elm)
              CASE(2)
                WRITE(fid_out,9012) data_in(1:n_elm)
              CASE(3)
                WRITE(fid_out,9013) data_in(1:n_elm)
              CASE(4)
                WRITE(fid_out,9014) data_in(1:n_elm)
              CASE(5)
                WRITE(fid_out,9015) data_in(1:n_elm)
              CASE(6)
                WRITE(fid_out,9016) data_in(1:n_elm)
              CASE DEFAULT
                STOP
            END SELECT
          ELSE
            n_col_final = MOD(n_elm, n_col)
            n_line = (n_elm - MOD(n_elm, n_col))/n_col + 1
            DO i_line = 1, n_line - 1
              counter = (i_line-1) * n_col
              SELECT CASE(n_col)
                CASE(5)
                  WRITE(fid_out,9015)& 
                  data_in(counter + 1 : counter + n_col)
                CASE(6)
                  WRITE(fid_out,9016)& 
                  data_in(counter + 1 : counter + n_col)
                CASE DEFAULT
                  STOP
              END SELECT
            ENDDO
            counter = (n_line-1) * n_col
            SELECT CASE(n_col_final)
              CASE(1)
                WRITE(fid_out,9011) data_in(counter + 1 : counter + n_col_final)
              CASE(2)
                WRITE(fid_out,9012) data_in(counter + 1 : counter + n_col_final)
              CASE(3)
                WRITE(fid_out,9013) data_in(counter + 1 : counter + n_col_final)
              CASE(4)
                WRITE(fid_out,9014) data_in(counter + 1 : counter + n_col_final)
              CASE(5)
                WRITE(fid_out,9015) data_in(counter + 1 : counter + n_col_final)
              CASE(6)
                WRITE(fid_out,9016) data_in(counter + 1 : counter + n_col_final)
              CASE DEFAULT
                STOP
            END SELECT
          ENDIF
        ENDIF

9011    FORMAT(I12)
9012    FORMAT(2I12)
9013    FORMAT(3I12)
9014    FORMAT(4I12)
9015    FORMAT(5I12)
9016    FORMAT(6I12)
      END SUBROUTINE write_data_int

END MODULE fchk2inputs

