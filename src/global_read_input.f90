! This module is part of d77 and reads input files.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE global_read_input
  USE ifmod, ONLY: write_messages
  USE global_constants
  IMPLICIT NONE

! Gneral information of INPUT file
  CHARACTER(LEN=100), ALLOCATABLE, SAVE :: Text_input(:)
  INTEGER, SAVE :: N_line_elec_state, N_line_end_elec_state, &
                  &N_line_grid, N_line_end_grid, &
                  &N_line_mode_calc, N_line_end_mode_calc

! $CONTROL keywords  
  CHARACTER(LEN=100), SAVE :: Runtyp
  CHARACTER(LEN=100), SAVE :: Property
  CHARACTER(LEN=100), SAVE :: Method
  CHARACTER(LEN=100), SAVE :: Effcharg
  CHARACTER(LEN=100), SAVE :: Effcharg_soc
  DOUBLE PRECISION, SAVE   :: Temperature
  DOUBLE PRECISION, SAVE   :: Scale_fac_freq
  DOUBLE PRECISION, SAVE   :: Threshold_e_ab, &
                             &Threshold_distance, Threshold_s_cgf, Threshold_ci_ci, &
                             &Threshold_contribution
  CHARACTER(LEN=100), SAVE :: Check_orthonorm_s_cgf, Check_orthonorm_vibmode, &
                             &Active_space_only, Read_ci_coef, &
                             &Save_soc_cgf, &
                             &Save_avcc, &
                             &Debug
  CHARACTER(LEN=100), SAVE :: Fchk_mocube

! $ELEC_STATE keywords  
  INTEGER, SAVE :: Bra, Ket

  CHARACTER(LEN=100), SAVE :: Save_dmat_soc, Asoc
  CHARACTER(LEN=100), SAVE :: Gen_opr_dmat, Save_opr_dmat
  CHARACTER(LEN=100), SAVE :: Gen_dmat_cgf, Save_dmat_cgf
  CHARACTER(LEN=100), SAVE :: Gen_dmat_soc_cgf, Save_dmat_soc_cgf
  CHARACTER(LEN=100), SAVE :: Gen_soc_cgf, Gen_soc_cgf_atm
  CHARACTER(LEN=100), SAVE :: Gen_ke_cgf, Gen_pe_cgf_atm, Gen_elfld_cgf_atm
  CHARACTER(LEN=100), SAVE :: Save_g16_format

! $SPIN_ORBIT keywords  
  INTEGER, SAVE :: Ms

CONTAINS

  SUBROUTINE read_input_control

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'read_input_control'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

!   Local variables
    INTEGER :: i_line, n_line
    INTEGER :: N_line_control, N_line_end_control
    INTEGER :: N_line_spin_orbit, N_line_end_spin_orbit
    CHARACTER(LEN=100) :: text_temp
    INTEGER :: i_symbol
    CHARACTER(LEN=100) :: text_left
    CHARACTER(LEN=100) :: text_right, text_right_temp
    CHARACTER(LEN=4) :: real_int_char
    DOUBLE PRECISION :: value_real
    LOGICAL :: value_logc
    INTEGER :: value_int
    INTEGER :: i

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)

!   Global variables
    N_line_elec_state = 0; N_line_end_elec_state = 0
    N_line_grid       = 0; N_line_end_grid       = 0
    N_line_mode_calc  = 0; N_line_end_mode_calc  = 0

    Read_ci_coef       = 'Yes'
    Save_g16_format    = 'No'

    Gen_ke_cgf         = 'Calc' 
    Gen_pe_cgf_atm     = 'Calc'
    Gen_elfld_cgf_atm  = 'Calc'

    Fchk_mocube = REPEAT(' ',100)


!   -----------
!   Count lines  
!   -----------

    n_line = 0
    OPEN(Finp, FILE = Fname(Finp))
      DO
        READ(Finp, '()', END=100)
        n_line = n_line + 1
      ENDDO
    100 CLOSE(Finp)

!   -------------------------
!   Making copy of INPUT file
!   -------------------------

    ALLOCATE(Text_input(1:n_line)); Text_input = REPEAT(' ', 100)
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '===== Copy of your input file ====='
    WRITE(6,'(1X)')
    OPEN(Finp, FILE = Fname(Finp))
      DO i_line = 1, n_line
        READ(Finp, '(A)') Text_input(i_line)
        WRITE(6,'(1X, A)') TRIM(ADJUSTL(Text_input(i_line)))
      ENDDO
    CLOSE(Finp)
    WRITE(6,'(1X)')
    WRITE(6,'(1X, A)') '==================================='
    WRITE(6,'(1X)') 
    WRITE(6,'(1X, A)') 'Reading $CONTROL options in your input file'

!   Checking input options
    N_line_control    = 0; N_line_end_control    = 0
    N_line_elec_state = 0; N_line_end_elec_state = 0
    N_line_grid       = 0; N_line_end_grid       = 0
    N_line_mode_calc  = 0; N_line_end_mode_calc  = 0

    DO i_line = 1, n_line
      text_temp  = REPEAT(' ', 100)
      IF(INDEX(Text_input(i_line),'$')/=0) THEN
        DO i = 2, LEN_TRIM(Text_input(i_line))
          text_temp(i-1:i-1) = Text_input(i_line)(i:i)
        ENDDO
!       Capitalize the first letter of text_left and make the rest lowercase.
        CALL cap_text(text_temp)
        Text_input(i_line) = '$'//TRIM(text_temp)
      ELSE
      ENDIF
    ENDDO

    DO i_line = 1, n_line
      text_temp = Text_input(i_line)
      IF(INDEX(text_temp,'$Control')/=0) THEN
        N_line_control = i_line
      ELSEIF(INDEX(text_temp,'$End_control')/=0) THEN
        N_line_end_control = i_line
      ELSEIF(INDEX(text_temp,'$Elec_state')/=0) THEN
        N_line_elec_state = i_line
      ELSEIF(INDEX(text_temp,'$End_elec_state')/=0) THEN
        N_line_end_elec_state = i_line
      ELSEIF(INDEX(text_temp,'$Grid')/=0) THEN
        N_line_grid = i_line
      ELSEIF(INDEX(text_temp,'$End_grid')/=0) THEN
        N_line_end_grid = i_line
      ELSEIF(INDEX(text_temp,'$Vibmode')/=0) THEN
        N_line_mode_calc = i_line
      ELSEIF(INDEX(text_temp,'$End_vibmode')/=0) THEN
        N_line_end_mode_calc = i_line
      ELSEIF(INDEX(text_temp,'$Spin_orbit')/=0) THEN
        N_line_spin_orbit = i_line
      ELSEIF(INDEX(text_temp,'$End_spin_orbit')/=0) THEN
        N_line_end_spin_orbit = i_line
      ELSE
      ENDIF
    ENDDO

!   Default values
    Runtyp   = 'Int_pgf'
    Property = 'Dipole'
    Method   = 'Td'
    Debug    = 'No'

    Temperature    = 0.0D0
    Scale_fac_freq = 1.0D0
    
    Effcharg     = 'Nuccharg'
    Effcharg_soc = 'Read'

    Threshold_e_ab         = Threshold_e_ab_default
    Threshold_s_cgf        = Threshold_s_cgf_default
    Threshold_ci_ci        = Threshold_ci_ci_default
    Threshold_contribution = Threshold_contribution_default
    Threshold_distance     = Threshold_distance_default

    Check_orthonorm_s_cgf   = 'No'
    Check_orthonorm_vibmode = 'No'

    Gen_opr_dmat     = 'Calc'; Save_opr_dmat      = 'No'
    Gen_dmat_cgf     = 'Calc'; Save_dmat_cgf      = 'No'

    Active_space_only = 'No'
    Read_ci_coef      = 'Yes'
    
    Gen_elfld_cgf_atm = 'Read'
    Save_avcc         = 'No'

    IF(N_line_control /= 0 .AND. &
      &N_line_control + 1 < N_line_end_control) THEN
      DO i_line = N_line_control + 1, N_line_end_control - 1
        text_temp = REPEAT(' ', 100)
        text_temp = Text_input(i_line)
        i_symbol = INDEX(text_temp, '=')
        text_left = REPEAT(' ', 100)
        CALL get_text_left_symbol&
             (text_temp, i_symbol, text_left)
        text_right = REPEAT(' ', 100)
        CALL get_text_right_symbol&
             (text_temp, i_symbol, text_right)


!       Capitalize the first letter of text_left and make the rest lowercase.
        CALL cap_text(text_left)
!       Capitalize the first letter of text_right and make the rest lowercase.
        CALL cap_text(text_right)

        real_int_char = REPEAT(' ', 4)
        CALL get_typ(text_left, real_int_char)

        value_real = 0.0D0
        value_logc = .FALSE.
        SELECT CASE(real_int_char)
          CASE('REAL'); CALL get_real(text_right, value_real)
          CASE('LOGC'); CALL get_logc(text_right, value_logc)
          CASE DEFAULT
        END SELECT

        SELECT CASE(text_left)
          CASE('Runtyp');   Runtyp   = text_right
          CASE('Property'); Property = text_right
          CASE('Method');   Method   = text_right
          CASE('Debug');    Debug    = text_right

          CASE('Temperature');    Temperature    = value_real
          CASE('Scale_fac_freq'); Scale_fac_freq = value_real
          
          CASE('Effcharg');     Effcharg     = text_right
          CASE('Effcharg_soc'); Effcharg_soc = text_right

          CASE('Threshold_e_ab');         Threshold_e_ab         = value_real
          CASE('Threshold_s_cgf');        Threshold_s_cgf        = value_real
          CASE('Threshold_ci_ci');        Threshold_ci_ci        = value_real
          CASE('Threshold_contribution'); Threshold_contribution = value_real
          CASE('Threshold_distance');     Threshold_distance     = value_real

          CASE('Check_orthonorm_s_cgf');   Check_orthonorm_s_cgf   = text_right
          CASE('Check_orthonorm_vibmode'); Check_orthonorm_vibmode = text_right

          CASE('Save_avcc');  Save_avcc  = text_right

          CASE('Gen_opr_dmat');      Gen_opr_dmat      = text_right
          CASE('Save_opr_dmat');     Save_opr_dmat     = text_right
          CASE('Gen_dmat_cgf');      Gen_dmat_cgf      = text_right
          CASE('Save_dmat_cgf');     Save_dmat_cgf     = text_right

          CASE('Active_space_only'); Active_space_only = text_right
          CASE('Read_ci_coef');      Read_ci_coef      = text_right
          
          CASE('Save_soc_cgf'); Save_soc_cgf = text_right
          
          CASE('Gen_elfld_cgf_atm'); Gen_elfld_cgf_atm = text_right
          
          CASE('Fchk_mocube'); Fchk_mocube = text_right
          
          CASE DEFAULT
            WRITE(6,'(1X, A)') 'Invalid Control keyward: ', text_left
            CALL write_messages(-9999, Text_blank, type_program, name_program)
        END SELECT
      ENDDO
    ELSE
      WRITE(6,'(1X, A)') 'Missing Control option'
      CALL write_messages(-9999, Text_blank, type_program, name_program)
    ENDIF

    SELECT CASE(Property)
      CASE('Cicoef_only', 'Dvne_only', 'Rho', 'Ke', 'Pe', 'Dipole', 'Vc', 'Soc')
      CASE DEFAULT
        WRITE(6,'(1X, A, A)') 'Invalid property: ', Property
        CALL write_messages(-9999, Text_blank, type_program, name_program)
    END SELECT

    SELECT CASE(Runtyp)
      CASE('Calc_int_cgf', 'Int_pgf', 'Density')
      CASE DEFAULT
        WRITE(6,'(1X, A, A)') 'Invalid runtyp: ', Runtyp
        CALL write_messages(-9999, Text_blank, type_program, name_program)
    END SELECT  

    SELECT CASE(Runtyp)
      CASE('Int_pgf', 'Density')
        SELECT CASE(Method)
          CASE('Td', 'Cis')
            Read_ci_coef = 'Yes'
          CASE DEFAULT
            WRITE(6,'(1X, A, A)') 'Invalid method: ', method
            CALL write_messages(-9999, Text_blank, type_program, name_program)
        END SELECT  
      CASE DEFAULT
    END SELECT  

    WRITE(6,'(1X, A)') 'Done'

    SELECT CASE(Runtyp)
      CASE('Int_pgf', 'Density')
        SELECT CASE(Property)
          CASE('Electron_density', 'Rho', 'Dipole', 'Vc', 'Soc')
!           Reading Bra and Ket
            IF(N_line_elec_state /= 0 .AND.&
               N_line_elec_state + 1 < N_line_end_elec_state) THEN
              DO i_line = N_line_elec_state + 1, N_line_end_elec_state - 1
                text_temp = Text_input(i_line)
                i_symbol = INDEX(text_temp, '=')
                CALL get_text_left_symbol&
                     (text_temp, i_symbol, text_left)

!               Capitalize the first letter of text_left and make the rest lowercase.
                CALL cap_text(text_left)
  
                CALL get_typ(text_left, real_int_char)
                SELECT CASE(text_left)
                  CASE('Bra') !  text_left = bra
                    i_symbol = INDEX(text_temp, '<')
                    CALL get_text_right_symbol&
                         (text_temp, i_symbol, text_right_temp)
                    i_symbol = INDEX(text_right_temp, '|')
                    CALL get_text_left_symbol&
                         (text_right_temp, i_symbol, text_right)
                    value_int = 0
                    CALL get_int(text_right, value_int)
                    Bra = value_int
                  CASE('Ket') ! text_left = ket
                    i_symbol = INDEX(text_temp, '|')
                    CALL get_text_right_symbol&
                         (text_temp, i_symbol, text_right_temp)
                    i_symbol = INDEX(text_right_temp, '>')
                    CALL get_text_left_symbol&
                         (text_right_temp, i_symbol, text_right)
                    value_int = 0
                    CALL get_int(text_right, value_int)
                    Ket = value_int
                  CASE DEFAULT ! text_left /= bra, ket
                    WRITE(6,'(1X, A)') 'Incorrect ELEC_STATE keyword'
                    CALL write_messages(-9999, Text_blank, type_program, name_program)
                END SELECT
              ENDDO ! i_line = 1, n_line
            ELSE ! Default $ELEC_STATE options
              WRITE(6,'(1X, A)') 'Missing $ELEC_STATE option'
              CALL write_messages(-9999, Text_blank, type_program, name_program)
            ENDIF

            IF(Bra > Max_n_state) THEN
              WRITE(6,'(1X, A)') 'Too large Bra value'
              WRITE(6,'(1X, A, I0)') 'Bra value is less than or eqaul to ', Max_n_state
              CALL write_messages(-9999, Text_blank, type_program, name_program)
            ELSEIF(Ket > Max_n_state) THEN
              WRITE(6,'(1X, A)') 'Too large Ket value'
              WRITE(6,'(1X, A, I0)') 'Ket value is less than or eqaul to ', Max_n_state
              CALL write_messages(-9999, Text_blank, type_program, name_program)
            ELSE
            ENDIF
          CASE DEFAULT
        END SELECT
      CASE DEFAULT
    END SELECT

!   Reading $SPIN_ORBIT options
    SELECT CASE(Property)
      CASE('Soc')
        SELECT CASE(Runtyp)
          CASE('Int_pgf', 'Density')
            IF(N_line_spin_orbit /= 0 .AND.&
               N_line_spin_orbit + 1 < N_line_end_spin_orbit) THEN
              DO i_line = N_line_spin_orbit + 1, N_line_end_spin_orbit - 1
                text_temp = Text_input(i_line)
                WRITE(*,*) text_temp
                i_symbol = INDEX(text_temp, '=')
                CALL get_text_left_symbol&
                     (text_temp, i_symbol, text_left)
                CALL get_text_right_symbol&
                     (text_temp, i_symbol, text_right)
          
!               Capitalize the first letter of text_left and make the rest lowercase.
                CALL cap_text(text_left)
!               Capitalize the first letter of text_right and make the rest lowercase.
                CALL cap_text(text_right)
            
                real_int_char = REPEAT(' ', 4)
                CALL get_typ(text_left, real_int_char)
                Ms   = 0
                Asoc = 'No'
                Gen_soc_cgf      = 'Read'; Save_soc_cgf       = 'No'
                Gen_dmat_soc_cgf = 'Calc'; Save_dmat_soc_cgf  = 'No'
                value_int = 0
                SELECT CASE(real_int_char)
                  CASE('INT'); CALL get_int(text_right, value_int)
                  CASE DEFAULT
                END SELECT
                SELECT CASE(text_left)
                  CASE('Ms'); Ms = value_int
                  CASE('Gen_soc_cgf'); Gen_soc_cgf = text_right
                  CASE('Save_soc_cgf'); Save_soc_cgf = text_right
                  CASE('Gen_dmat_soc_cgf'); Gen_dmat_soc_cgf = text_right
                  CASE('Save_dmat_soc_cgf'); Save_dmat_soc_cgf = text_right
                  
                  CASE('Save_dmat_soc'); Save_dmat_soc = text_right
                  CASE DEFAULT
                    WRITE(6,'(1X, A)') 'Incorrect SPIN_ORBIT keyword: ', text_left
                    CALL write_messages(-9999, Text_blank, type_program, name_program)
                END SELECT
              ENDDO ! i_line = 1, n_line
            ELSE ! Missing SPIN_ORBIT option
              WRITE(6,'(1X, A)') 'Missing $SPIN_ORBIT option'
              CALL write_messages(-9999, Text_blank, type_program, name_program)
            ENDIF ! N_line_spin_orbit = 0
          CASE DEFAULT    
        END SELECT
      CASE DEFAULT
    END SELECT

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  
  END SUBROUTINE read_input_control


  SUBROUTINE get_real(value_text, value_real)  

!   Input variable
    CHARACTER(LEN=100), INTENT(IN) :: value_text ! value_text = -1.0D-5, -300.0

!   Output variable
    DOUBLE PRECISION, INTENT(OUT) :: value_real

!   Local variables
    INTEGER iE_or_D
    INTEGER sign_exp, sign_coef
    INTEGER nexp
    DOUBLE PRECISION coef
    INTEGER i, j
    INTEGER:: num09, ipunct

    value_real = 0.0D0

    iE_or_D = 0
    sign_exp = +1
    nexp=0   ! exponent
    sign_coef = +1

!   Reading exponent
    iE_or_D = INDEX(value_text, 'E') + INDEX(value_text,'e')&
            + INDEX(value_text, 'D') + INDEX(value_text,'d')
    IF(iE_or_D == 0) GOTO 4000 ! value_text = -300.0

    DO i = 1, LEN_TRIM(value_text) - iE_or_D ! For -1.6D-31, i = 1, 3
      j = LEN_TRIM(value_text) + 1 - i   ! j = 8, 6
      IF (value_text(j:j) == '+') EXIT   
      IF (value_text(j:j) == '-') THEN   
        sign_exp = -1             
      EXIT    
    ENDIF
    CALL char2num(value_text(j:j), num09)  ! convert to integer 0-9
    nexp = nexp + num09*(10**(i-1))   ! exponent number (integer)
    ENDDO
    nexp = sign_exp*nexp

 4000   CONTINUE

    !   Read coefficient
    ! -1.0D-5, -1.0, -1.D-5, -1., -1D-5, -1
    IF(iE_or_D == 0) THEN ! -1.0, -1., -1
      iE_or_D = INDEX(value_text,' ') ! iE_or_D = 5, 4, 3
    ENDIF
    ipunct = INDEX(value_text,'.')
    IF(ipunct .EQ. 0) THEN ! -1D-5, -1
      ipunct = iE_or_D ! ipunct = 3, 3
      GOTO 5000
    END IF
    IF(ipunct == LEN_TRIM(value_text)) GOTO 5000 ! -1.
    IF(ipunct + 1 == iE_or_D) GOTO 5000 ! -1.D-5
    !IF(ipunct == iE_or_D) GOTO 5000   ! -1.

    coef = 0.0D0
    DO i = 1, iE_or_D - ipunct - 1
      j = ipunct + i
      CALL char2num(value_text(j:j), num09)  
      coef=coef+DBLE(num09)*(DBLE(10)**DBLE(-i))
    ENDDO
    IF(ipunct == 1) GOTO 8000 
 5000   CONTINUE
    DO i = 1, ipunct - 1
    j = ipunct - i
    IF (value_text(j:j)=='+') CYCLE
    IF (value_text(j:j)=='-') THEN
        sign_coef = -1
        CYCLE
    END IF
    CALL char2num(value_text(j:j), num09)
    coef=coef+DBLE(num09*(10**(i-1)))
    ENDDO ! K=1,ipunct-1
 8000   CONTINUE
    value_real = DBLE(sign_coef)*coef*DBLE(10)**DBLE(nexp)
    CONTINUE
  END SUBROUTINE get_real

  SUBROUTINE get_int(value_text, value_int)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'get_int'
    CHARACTER(LEN=100), PARAmetER :: type_program = 'SUBROUTINE'

!   Input variable
    CHARACTER(LEN=100), INTENT(IN) :: value_text ! value_text = -1.0D-5, -300.0

!   Output variable
    INTEGER, INTENT(OUT) :: value_int

!   Local variables
    INTEGER iE_or_D
    INTEGER nexp, coef
    INTEGER sign_coef
    INTEGER i, j
    INTEGER num09

    value_int = 0

    iE_or_D = 0
    nexp = 0  ! exponent
    sign_coef = +1

!   Reading exponent
    DO i=1, LEN_TRIM(value_text)
      IF(value_text(i:i)=='.') THEN
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ENDIF
    ENDDO 
    iE_or_D = INDEX(value_text, 'E') + INDEX(value_text,'e')&
            + INDEX(value_text, 'D') + INDEX(value_text,'d')
    IF (iE_or_D == 0) THEN
      iE_or_D = INDEX(value_text,' ')  
        GOTO 4000
    ENDIF

    DO i = iE_or_D + 1, LEN_TRIM(value_text)
      IF(value_text(i:i) == '-') THEN
        CALL write_messages(-9999, Text_blank, type_program, name_program)
      ENDIF
    ENDDO
    DO i = 1, LEN_TRIM(value_text) - iE_or_D
      j = LEN_TRIM(value_text) + 1 - i
      IF (value_text(j:j)=='+') CYCLE  
        CALL char2num(value_text(j:j), num09)
        nexp=nexp+num09*(10**(i-1))   
    ENDDO
 4000   CONTINUE

    coef = 0    !  
    sign_coef = +1   
    DO i = 1, iE_or_D - 1  
      j = iE_or_D - i       
      IF(value_text(j:j)=='+')  CYCLE 
      IF(value_text(j:j)=='-')  THEN  
        sign_coef = - 1           
        CYCLE           
      ENDIF
      CALL char2num(value_text(j:j), num09)        
      coef = coef + num09*(10**(i-1))   
    ENDDO 
    coef = sign_coef * coef  
    value_int = coef * (10**nexp)

  END SUBROUTINE get_int


  SUBROUTINE get_logc(value_text, value_logc)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'get_logc'
    CHARACTER(LEN=100), PARAmetER :: type_program = 'SUBROUTINE'

!   Input variable
    CHARACTER(LEN=100), INTENT(IN) :: value_text ! value_text = TRUE

!   Output variable
    LOGICAL, INTENT(OUT) :: value_logc

    SELECT CASE(TRIM(ADJUSTL(value_text)))
      CASE('TRUE', 'T', 'true', 't')
        value_logc = .TRUE.
      CASE('FALSE', 'F', 'false', 'f')
        value_logc = .FALSE.
      CASE DEFAULT
        WRITE(6,'(A)') 'Error termination'
        WRITE(6,'(A)') 'Incorrect input file'
        CALL write_messages(-9999, Text_blank, type_program, name_program)
    END SELECT

  END SUBROUTINE get_logc


  SUBROUTINE char2num(char09, num09)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'char2num'
    CHARACTER(LEN=100), PARAmetER :: type_program = 'SUBROUTINE'

!   Input variable
    CHARACTER(1), INTENT(IN) :: char09

!   Output variable
    INTEGER, INTENT(OUT) :: num09

    num09 = 10 ! temporary values
    SELECT CASE(char09)
      CASE('0'); num09 = 0
      CASE('1'); num09 = 1
      CASE('2'); num09 = 2
      CASE('3'); num09 = 3
      CASE('4'); num09 = 4
      CASE('5'); num09 = 5
      CASE('6'); num09 = 6
      CASE('7'); num09 = 7
      CASE('8'); num09 = 8
      CASE('9'); num09 = 9
      CASE DEFAULT
        !CALL write_error(0)
        WRITE(6,*) char09, ' must be an integer'
        CALL write_messages(-9999, Text_blank, type_program, name_program)
    END SELECT
  
  END SUBROUTINE char2num

  SUBROUTINE get_text_left_symbol&
             (text, i_symbol, text_left)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'get_text_left_symbol'
    CHARACTER(LEN=100), PARAmetER :: type_program = 'SUBROUTINE'

!   Input variables
    CHARACTER(LEN=100), INTENT(IN) :: text
    INTEGER, INTENT(IN)            :: i_symbol

!   Output variable
    CHARACTER(LEN=100), INTENT(OUT) :: text_left

!   Local variables
    INTEGER i

    text_left = REPEAT(' ',100)

    ! text = '   temperature    =    298.0  ...'
    DO i = 1, i_symbol - 1
      text_left(i:i) = text(i:i)  ! text_left = '   temperature ...'
    ENDDO
    text_left = ADJUSTL(text_left)   ! text_left = 'temperature ...'
    IF(LEN_TRIM(text_left) == 0) THEN
      WRITE(6,'(1X)') 
      WRITE(6,'(1X, A)') '===== Error detected in your input file ====='
      WRITE(6,'(1X)') 
      WRITE(6,'(1X, A)') 'Missing input variable'
      CALL write_messages(-9999, Text_blank, type_program, name_program)
    ENDIF

  END SUBROUTINE get_text_left_symbol

  SUBROUTINE get_text_right_symbol&
             (text, i_symbol, text_right)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'get_text_right_symbol'
    CHARACTER(LEN=100), PARAmetER :: type_program = 'SUBROUTINE'

!   Input variables
    CHARACTER(LEN=100), INTENT(IN) :: text
    INTEGER, INTENT(IN)            :: i_symbol

!   Output variable
    CHARACTER(LEN=100), INTENT(OUT) :: text_right

!   Local variable
    INTEGER i

!   Initializing variables
    text_right = REPEAT(' ',100)

    ! text = '   temperature    =    298.0  ...'
    DO i = i_symbol + 1, 100
      text_right(i:i) = text(i:i) ! text_right = '... 298.0 ...'
    ENDDO
    text_right = ADJUSTL(text_right) ! text_right = '298.0       ...'
    IF(LEN_TRIM(text_right) == 0) THEN
      WRITE(6,'(1X, A)') '===== Error detected in your input file ====='
      WRITE(6,'(1X)') 
      WRITE(6,'(1X, A)') 'Missing input value'
      CALL write_messages(-9999, Text_blank, type_program, name_program)
    ELSE
    ENDIF

  END SUBROUTINE get_text_right_symbol


  SUBROUTINE get_typ(namevar, real_int_char)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'get_typ'
    CHARACTER(LEN=100), PARAmetER :: type_program = 'SUBROUTINE'

!   Input variable
    CHARACTER(LEN=100), INTENT(IN) :: namevar

!   Output variable
    CHARACTER(LEN=4), INTENT(OUT) :: real_int_char

    real_int_char = REPEAT(' ', 4) ! initialize real_int_char

    SELECT CASE(namevar)
      CASE('Runtyp');         real_int_char = 'CHAR'
      CASE('Property');       real_int_char = 'CHAR'
      CASE('Method');         real_int_char = 'CHAR'
      CASE('Effcharg');       real_int_char = 'CHAR'
      CASE('Effcharg_soc');   real_int_char = 'CHAR'
      CASE('Temperature');    real_int_char = 'REAL'
      CASE('Scale_fac_freq'); real_int_char = 'REAL'

      CASE('Threshold_s_cgf');        real_int_char = 'REAL'
      CASE('Threshold_e_ab');         real_int_char = 'REAL'
      CASE('Threshold_ci_ci');        real_int_char = 'REAL'
      CASE('Threshold_contribution'); real_int_char = 'REAL'
      CASE('Threshold_distance');     real_int_char = 'REAL'

      CASE('Check_orthonorm_s_cgf'); real_int_char = 'CHAR'
      CASE('Check_orthonorm_vibmode'); real_int_char = 'CHAR'
      
      
      CASE('Gen_ke_cgf'); real_int_char = 'CHAR' 
      CASE('Gen_pe_cgf_atm'); real_int_char = 'CHAR' 
      CASE('Gen_elfld_cgf_atm'); real_int_char = 'CHAR' 
      
      CASE('Save_soc_cgf'); real_int_char = 'CHAR'
      CASE('Save_avcc'); real_int_char = 'CHAR'

      CASE('Active_space_only'); real_int_char = 'CHAR'
      CASE('Read_ci_coef'); real_int_char = 'CHAR'
      CASE('Debug'); real_int_char = 'CHAR'
      CASE('Bra'); real_int_char = 'CHAR'
      CASE('Ket'); real_int_char = 'CHAR'
      CASE('Gen_dmat_cgf'); real_int_char = 'CHAR'
      CASE('Save_dmat_cgf'); real_int_char = 'CHAR'
      CASE('Gen_opr_dmat'); real_int_char = 'CHAR'
      CASE('Save_opr_dmat'); real_int_char = 'CHAR'
      CASE('Ms'); real_int_char = 'INT'
      CASE('Xmin'); real_int_char = 'REAL'
      CASE('Xmax'); real_int_char = 'REAL'
      CASE('Ymin'); real_int_char = 'REAL'
      CASE('Ymax'); real_int_char = 'REAL'
      CASE('Zmin'); real_int_char = 'REAL'
      CASE('Zmax'); real_int_char = 'REAL'
      CASE('Dx'); real_int_char = 'REAL'
      CASE('Dy'); real_int_char = 'REAL'
      CASE('Dz'); real_int_char = 'REAL'
      CASE('Fchk_mocube'); real_int_char = 'CHAR'
      CASE DEFAULT
        WRITE(6,'(1X, A)') 'Invalid CONTROL keyward: ', namevar
        CALL write_messages(-9999, Text_blank, type_program, name_program)
    END SELECT

  END SUBROUTINE get_typ


  SUBROUTINE cap_text(text)

!   Input and output variable
    CHARACTER(LEN=100), INTENT(INOUT) :: text

!   Local variables
    CHARACTER(LEN=100) :: text_temp
    INTEGER            :: i

!   Capitalize the first letter of text and make the rest lowercase.
    text_temp = REPEAT(' ', 100)
    text_temp = text
    IF(text_temp(1:1) >= 'a' .AND. text_temp(1:1) <= 'z') THEN
      text_temp(1:1) = CHAR(ICHAR(text_temp(1:1))-32)
    ENDIF
    DO i = 2, LEN_TRIM(text_temp)
      IF(text_temp(i:i) >= 'A' .AND. text_temp(i:i) <= 'Z') THEN
        text_temp(i:i) = CHAR(ICHAR(text_temp(i:i))+32)
      ENDIF
    ENDDO
    text = text_temp

  END SUBROUTINE cap_text


END MODULE global_read_input

