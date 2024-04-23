! This subroutine is part of d77 and computes
! one-particle reduced density matrix (gamma).

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

SUBROUTINE calc_opr_dmat&
          &(left_state, right_state, n_ci_ci, opr_dmat, trace_opr_dmat)

  USE ifmod, ONLY: write_messages
  USE global_constants
  USE global_read_input
  USE global_read_data
  USE func_int_slater
  IMPLICIT NONE

! ------------------------
! Declaration of variables
! ------------------------

! Program name
  CHARACTER(LEN=100), PARAMETER :: name_program = 'calc_opr_dmat'
  CHARACTER(LEN=100), PARAMETER :: type_program = 'SUBROUTINE'

! Input variables
  INTEGER, INTENT(IN) :: left_state, right_state

! Output variables
  INTEGER, INTENT(OUT)          :: n_ci_ci
  DOUBLE PRECISION, INTENT(OUT) :: opr_dmat(:,:)
  DOUBLE PRECISION, INTENT(OUT) :: trace_opr_dmat 

! Local variables
  INTEGER :: i, j 
  INTEGER :: left_elec_config, right_elec_config

  INTEGER :: tuple_temp, tuple_left, tuple_right
  INTEGER :: a_ras_temp(1:Max_n_tuple), r_ras_temp(1:Max_n_tuple)
  INTEGER :: a(1:Max_n_tuple), r(1:Max_n_tuple), &
            &b(1:Max_n_tuple), s(1:Max_n_tuple)
  INTEGER :: min_a, max_r, min_b, max_s
  INTEGER :: min_so_temp, max_so_temp
  INTEGER :: min_so(1:Max_n_ci_ci), max_so(1:Max_n_ci_ci)
 
  INTEGER :: p, q
  INTEGER :: counter, i_ci_ci
  INTEGER :: tuple_left_ci_ci(1:Max_n_ci_ci), &
            &tuple_right_ci_ci(1:Max_n_ci_ci)
  INTEGER :: a_ci_ci(1:Max_n_tuple, 1:Max_n_ci_ci), &
            &r_ci_ci(1:Max_n_tuple, 1:Max_n_ci_ci), &
            &b_ci_ci(1:Max_n_tuple, 1:Max_n_ci_ci), &
            &s_ci_ci(1:Max_n_tuple, 1:Max_n_ci_ci)

  DOUBLE PRECISION :: ci_coef_temp
  DOUBLE PRECISION :: ci_coef_left, ci_coef_right
  DOUBLE PRECISION :: ci_ci_temp 
  DOUBLE PRECISION :: ci_ci(1:Max_n_ci_ci)

  IF(Debug == 'Yes') &
 &CALL write_messages(2, Text_blank, type_program, name_program)

! ----------------------
! Initializing variables
! ----------------------

! Output variables
  n_ci_ci = 0
  opr_dmat = 0.0D0
  trace_opr_dmat = 0.0D0

! Local variables
  left_elec_config = 0; right_elec_config = 0

  tuple_temp = 0; tuple_left = 0; tuple_right = 0
  a_ras_temp = 0; r_ras_temp = 0 
  a = 0; r = 0; b = 0; s = 0

  ci_coef_temp = 0.0D0
  ci_coef_left = 0.0D0; ci_coef_right = 0.0D0
  min_a = 0; max_r = 0; min_b = 0; max_s = 0
  min_so_temp = 0; max_so_temp = 0
  min_so = 0; max_so = 0

  ci_ci_temp = 0.0D0 
  ci_ci = 0.0D0

  tuple_left_ci_ci  = 0
  tuple_right_ci_ci = 0 
  a_ci_ci = 0
  r_ci_ci = 0
  b_ci_ci = 0
  s_ci_ci = 0

! --------------------
! Calculating opr_dmat
! --------------------

! Calculating for the case in which left_state <= right_state 
  SELECT CASE(method)
    CASE('Cis')
      IF(left_state == 0 .AND. right_state == 0) THEN
        IF(Active_space_only == 'Yes') THEN
        ELSE
          n_ci_ci = 1 
          DO p = 1, ne!, 2!, 2!1, n_cgf
            DO q = 1, ne!, 2!, 2!1, n_cgf
              opr_dmat(q, p) = func_int_slater_0_pq_0(ne, p, q)
            ENDDO ! q
          ENDDO ! p
        ENDIF 
      ELSEIF(left_state == 0 .AND. right_state /= 0) THEN
        counter = 0 
        a = 0; r = 0
        min_a = 0; max_r = N_so + 1
        tuple_left = 0
        ci_coef_left = 1.0D0
        DO right_elec_config = 1, n_elec_config(right_state)
          b = 0; s = 0
          min_b = 0; max_s = N_so + 1
          min_so_temp = 0; max_so_temp = N_so + 1
          tuple_right = tuple(right_elec_config, right_state)
          ci_coef_right = ci_coef(right_elec_config, right_state)
          SELECT CASE(tuple_right)
            CASE(0)
            CASE DEFAULT
              DO j = 1, tuple_right
                b(j) = a_ras(j, right_elec_config, right_state)
                s(j) = r_ras(j, right_elec_config, right_state)
              ENDDO
              min_b = MINVAL(b(1:tuple_right))
              max_s = MAXVAL(s(1:tuple_right))
          END SELECT
  
          ci_ci_temp = ci_coef_left * ci_coef_right
          IF(tuple_left == 0 .AND. tuple_right == 0) THEN
            min_so_temp = 1 
            max_so_temp = ne
          ELSEIF(tuple_left == 0 .AND. tuple_right /= 0) THEN
            IF(Active_space_only == 'Yes') THEN
              min_so_temp = min_b 
            ELSE
              min_so_temp = 1 
            ENDIF
            max_so_temp = max_s
          ELSE
          ENDIF
          IF(min_so_temp == 0 .OR. max_so_temp > N_so) THEN
          ELSE
          ENDIF
  
          IF(ABS(ci_ci_temp) > threshold_ci_ci) THEN
            IF(Active_space_only == 'Yes') THEN
              IF(tuple_left == 0 .AND. tuple_right == 0) CYCLE
            ELSE
            ENDIF
            counter = counter + 1
            ci_ci(counter) = ci_ci_temp
            tuple_left_ci_ci(counter)  = tuple_left
            tuple_right_ci_ci(counter) = tuple_right
            a_ci_ci(:, counter) = a
            r_ci_ci(:, counter) = r
            b_ci_ci(:, counter) = b
            s_ci_ci(:, counter) = s
            min_so(counter) = min_so_temp
            max_so(counter) = max_so_temp
          ELSE
            CYCLE
          ENDIF
        ENDDO ! right_elec_config = 1, n_elec_config(right_state)
        CALL FLUSH(6) 
        IF(counter > Max_n_ci_ci) THEN
          WRITE(6,'(1X, A)') 'Too many electronic configurations'
          WRITE(6,'(1X, A, I0)') 'The number of &
         &n_ci_ci &
         &must be less than or qual to ', Max_n_ci_ci
          WRITE(6,'(1X, A, I0)') 'n_ci_ci = ', counter
        ELSEIF(counter <= Max_n_ci_ci .AND.&
               counter >= 1) THEN
          n_ci_ci = counter
  
          DO i_ci_ci = 1, n_ci_ci
            !DO p = 1, max_so(i_ci_ci)!, 2!, 2!1, n_cgf
            !  DO q = 1, max_so(i_ci_ci)!, 2!, 2!1, n_cgf
            DO p = min_so(i_ci_ci), max_so(i_ci_ci)!, 2!, 2!1, n_cgf
              DO q = min_so(i_ci_ci), max_so(i_ci_ci)!, 2!, 2!1, n_cgf
                opr_dmat(q, p) &
               &= opr_dmat(q, p) &
               &+ ci_ci(i_ci_ci) &
               &* func_int_slater_m_pq_n&
                 &(ne, tuple_left_ci_ci(i_ci_ci), &
                  &tuple_right_ci_ci(i_ci_ci), &
                  &a_ci_ci(:, i_ci_ci), r_ci_ci(:, i_ci_ci), &
                  &p, q, &
                  &b_ci_ci(:, i_ci_ci), s_ci_ci(:, i_ci_ci))
              ENDDO
            ENDDO
          ENDDO
        ELSE
        ENDIF
      ELSEIF(left_state >= 1 .AND. right_state >= 1) THEN
        counter = 0 
        DO left_elec_config = 1, n_elec_config(left_state)
          a = 0; r = 0
          min_a = 0; max_r = N_so + 1
          tuple_left=tuple(left_elec_config, left_state)
          ci_coef_left=ci_coef(left_elec_config, left_state)
          SELECT CASE(tuple_left)
            CASE(0)
            CASE DEFAULT
              DO i = 1, tuple_left
                a(i) = a_ras(i, left_elec_config, left_state)
                r(i) = r_ras(i, left_elec_config, left_state)
              ENDDO
              min_a = MINVAL(a(1:tuple_left))
              max_r = MAXVAL(r(1:tuple_left))
          END SELECT
          DO right_elec_config = 1, n_elec_config(right_state)
            b = 0; s = 0
            min_b = 0; max_s = N_so + 1
            min_so_temp = 0; max_so_temp = N_so + 1
            tuple_right = tuple(right_elec_config, right_state)
            ci_coef_right = ci_coef(right_elec_config, right_state)
            SELECT CASE(tuple_right)
              CASE(0)
              CASE DEFAULT
                DO j = 1, tuple_right
                  b(j) = a_ras(j, right_elec_config, right_state)
                  s(j) = r_ras(j, right_elec_config, right_state)
                ENDDO
                min_b = MINVAL(b(1:tuple_right))
                max_s = MAXVAL(s(1:tuple_right))
            END SELECT
  
            ci_ci_temp = ci_coef_left * ci_coef_right
            IF(tuple_left == 0 .AND. tuple_right == 0) THEN
              min_so_temp = 1 
              max_so_temp = ne
            ELSEIF(tuple_left == 0 .AND. tuple_right /= 0) THEN
              IF(Active_space_only == 'Yes') THEN
                min_so_temp = min_b 
              ELSE
                min_so_temp = 1 
              ENDIF
              max_so_temp = max_s
            ELSEIF(tuple_left /= 0 .AND. tuple_right == 0) THEN
              IF(Active_space_only == 'Yes') THEN
                min_so_temp = min_a 
              ELSE
                min_so_temp = 1 
              ENDIF
              max_so_temp = max_r
            ELSEIF(tuple_left /= 0 .AND. tuple_right /= 0) THEN
              IF(Active_space_only == 'Yes') THEN
                min_so_temp = MIN(min_a, min_b) 
              ELSE
                min_so_temp = 1 
              ENDIF
              max_so_temp = MAX(max_r, max_s)
            ELSE
            ENDIF
            IF(min_so_temp == 0 .OR. min_so_temp > N_so) THEN
            ELSE
            ENDIF
  
            IF(ABS(ci_ci_temp) > threshold_ci_ci) THEN
              IF(Active_space_only == 'Yes') THEN
                IF(tuple_left == 0 .AND. tuple_right == 0) CYCLE
              ELSE
              ENDIF
              counter = counter + 1
              ci_ci(counter) = ci_ci_temp
              tuple_left_ci_ci(counter)  = tuple_left
              tuple_right_ci_ci(counter) = tuple_right
              a_ci_ci(:, counter) = a
              r_ci_ci(:, counter) = r
              b_ci_ci(:, counter) = b
              s_ci_ci(:, counter) = s
              min_so(counter) = min_so_temp
              max_so(counter) = max_so_temp
            ELSE
              CYCLE
            ENDIF
          ENDDO
        ENDDO
  
        IF(counter > Max_n_ci_ci) THEN
        ELSEIF(counter <= Max_n_ci_ci .AND.&
               counter >= 1) THEN
          n_ci_ci = counter
  
          DO i_ci_ci = 1, n_ci_ci
            DO p = min_so(i_ci_ci), max_so(i_ci_ci)!, 2!, 2!1, n_cgf
              DO q = min_so(i_ci_ci), max_so(i_ci_ci)!, 2!, 2!1, n_cgf
                opr_dmat(q, p) &
               &= opr_dmat(q, p) &
               &+ ci_ci(i_ci_ci) &
               &* func_int_slater_m_pq_n&
                &(ne, tuple_left_ci_ci(i_ci_ci), &
                 &tuple_right_ci_ci(i_ci_ci), &
                 &a_ci_ci(:, i_ci_ci), r_ci_ci(:, i_ci_ci), &
                 &p, q, &
                 &b_ci_ci(:, i_ci_ci), s_ci_ci(:, i_ci_ci))
              ENDDO
            ENDDO
          ENDDO
        ELSE
        ENDIF
      ENDIF

    CASE('Td')
      IF(left_state == 0 .AND. right_state == 0) THEN
        IF(Active_space_only == 'Yes') THEN
        ELSE
          n_ci_ci = 1
          DO p = 1, ne!, 2!, 2!1, n_cgf
            DO q = 1, ne!, 2!, 2!1, n_cgf
              opr_dmat(q, p) = func_int_slater_0_pq_0(ne, p, q)
            ENDDO ! q
          ENDDO ! p
        ENDIF
      ELSEIF(left_state == 0 .AND. right_state /= 0) THEN
        counter = 0
        a = 0; r = 0
        min_a = 0; max_r = N_so + 1
        ci_coef_left = 1.0D0
        DO right_elec_config = 1, n_elec_config(right_state)
          b = 0; s = 0
          min_b = 0; max_s = N_so + 1
          min_so_temp = 0; max_so_temp = N_so + 1
          !ci_coef_right = x_minus_y(right_elec_config, right_state)
          ci_coef_right = x_plus_y(right_elec_config, right_state)
          b(1) = a_ras(1, right_elec_config, right_state)
          s(1) = r_ras(1, right_elec_config, right_state)
          min_b = b(1)
          max_s = s(1)

          ci_ci_temp = ci_coef_left * ci_coef_right
          IF(Active_space_only == 'Yes') THEN
            min_so_temp = min_b
            max_so_temp = max_s
          ELSE
            min_so_temp = 1
            max_so_temp = N_so
          ENDIF
          !max_so_temp = max_s

          IF(ABS(ci_ci_temp) > threshold_ci_ci) THEN
            counter = counter + 1
            ci_ci(counter) = ci_ci_temp
            tuple_left_ci_ci(counter)  = 0
            tuple_right_ci_ci(counter) = 1
            a_ci_ci(:, counter) = a
            r_ci_ci(:, counter) = r
            b_ci_ci(:, counter) = b
            s_ci_ci(:, counter) = s
            min_so(counter) = min_so_temp
            max_so(counter) = max_so_temp
          ELSE
            CYCLE
          ENDIF
        ENDDO ! right_elec_config = 1, n_elec_config(right_state)
        CALL FLUSH(6)
        IF(counter > Max_n_ci_ci) THEN
          WRITE(6,'(1X, A)') 'Too many electronic configurations'
          WRITE(6,'(1X, A, I0)') 'The number of &
         &n_ci_ci &
         &must be less than or qual to ', Max_n_ci_ci
          WRITE(6,'(1X, A, I0)') 'n_ci_ci = ', counter
        ELSEIF(1 <= counter .AND. counter <= Max_n_ci_ci) THEN
          n_ci_ci = counter

          DO i_ci_ci = 1, n_ci_ci
            !DO p = 1, max_so(i_ci_ci)!, 2!, 2!1, n_cgf
            !  DO q = 1, max_so(i_ci_ci)!, 2!, 2!1, n_cgf
            DO p = min_so(i_ci_ci), max_so(i_ci_ci)!, 2!, 2!1, n_cgf
              DO q = min_so(i_ci_ci), max_so(i_ci_ci)!, 2!, 2!1, n_cgf
                opr_dmat(q, p) &
               &= opr_dmat(q, p) &
               &+ ci_ci(i_ci_ci) &
               &* func_int_slater_m_pq_n&
                 &(ne, tuple_left_ci_ci(i_ci_ci), &
                  &tuple_right_ci_ci(i_ci_ci), &
                  &a_ci_ci(:, i_ci_ci), r_ci_ci(:, i_ci_ci), &
                  &p, q, &
                  &b_ci_ci(:, i_ci_ci), s_ci_ci(:, i_ci_ci))
              ENDDO
            ENDDO
          ENDDO
        ELSE
        ENDIF
      ELSEIF(left_state >= 1 .AND. right_state >= 1) THEN
        counter = 0
        DO left_elec_config = 1, n_elec_config(left_state)
          a = 0; r = 0
          min_a = 0; max_r = N_so + 1
          ci_coef_left=X_minus_y(left_elec_config, left_state)
          !ci_coef_left=X_plus_y(left_elec_config, left_state)
          a(1) = a_ras(1, left_elec_config, left_state)
          r(1) = r_ras(1, left_elec_config, left_state)
          min_a = a(1)
          max_r = r(1)
          DO right_elec_config = 1, n_elec_config(right_state)
            b = 0; s = 0
            min_b = 0; max_s = N_so + 1
            min_so_temp = 0; max_so_temp = N_so + 1
            ci_coef_right = X_plus_y(right_elec_config, right_state)
            b(1) = a_ras(1, right_elec_config, right_state)
            s(1) = r_ras(1, right_elec_config, right_state)
            min_b = b(1)
            max_s = s(1)

            ci_ci_temp = ci_coef_left * ci_coef_right
            IF(Active_space_only == 'Yes') THEN
              min_so_temp = MIN(min_a, min_b)
              max_so_temp = MAX(max_r, max_s)
            ELSE
              min_so_temp = 1
              max_so_temp = N_so
            ENDIF
            !max_so_temp = MAX(max_r, max_s)

            IF(ABS(ci_ci_temp) > threshold_ci_ci) THEN
              counter = counter + 1
              ci_ci(counter) = ci_ci_temp
              tuple_left_ci_ci(counter)  = 1
              tuple_right_ci_ci(counter) = 1
              a_ci_ci(:, counter) = a
              r_ci_ci(:, counter) = r
              b_ci_ci(:, counter) = b
              s_ci_ci(:, counter) = s
              min_so(counter) = min_so_temp
              max_so(counter) = max_so_temp
            ELSE
              CYCLE
            ENDIF
          ENDDO
        ENDDO

        IF(counter > Max_n_ci_ci) THEN
          WRITE(6,'(1X, A)') 'Too many electronic configurations'
          WRITE(6,'(1X, A, I0)') 'The number of &
         &n_ci_ci &
         &must be less than or qual to ', Max_n_ci_ci
          WRITE(6,'(1X, A, I0)') 'n_ci_ci = ', counter
        ELSEIF(1 <= counter .AND. counter <= Max_n_ci_ci) THEN
          n_ci_ci = counter

          DO i_ci_ci = 1, n_ci_ci
            DO p = min_so(i_ci_ci), max_so(i_ci_ci)!, 2!, 2!1, n_cgf
              DO q = min_so(i_ci_ci), max_so(i_ci_ci)!, 2!, 2!1, n_cgf
                opr_dmat(q, p) &
               &= opr_dmat(q, p) &
               &+ ci_ci(i_ci_ci) &
               &* func_int_slater_m_pq_n&
                 &(ne, tuple_left_ci_ci(i_ci_ci), &
                  &tuple_right_ci_ci(i_ci_ci), &
                  &a_ci_ci(:, i_ci_ci), r_ci_ci(:, i_ci_ci), &
                  &p, q, &
                  &b_ci_ci(:, i_ci_ci), s_ci_ci(:, i_ci_ci))
              ENDDO
            ENDDO
          ENDDO
        ELSE
        ENDIF
      ENDIF

    CASE DEFAULT
      CALL write_messages(-5, Text_blank, type_program, name_program)
  END SELECT    

  DO p = 1, N_so!, 2!, 2!1, n_cgf
    trace_opr_dmat = trace_opr_dmat + opr_dmat(p, p)
  ENDDO ! p  

  IF(Save_opr_dmat == 'Yes') THEN
    OPEN(Fopr_dmat, FILE=Fname(Fopr_dmat))
      DO p = 1, N_so!, 2!, 2!1, n_cgf
        DO q = 1, N_so
          WRITE(Fopr_dmat, *) p, q, opr_dmat(p, q)
        ENDDO ! q
      ENDDO ! p
    CLOSE(Fopr_dmat)
  ELSE
  ENDIF

  IF(Debug == 'Yes') &
 &CALL write_messages(3, Text_blank, type_program, name_program)

END SUBROUTINE calc_opr_dmat
