! This module is part of d77 and computes
! integrals between Slater determinants.
!
! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU Lesser General Public License v3.0
! as published by the Free Software Foundation.
! http://www.gnu.org/licenses
!
! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE func_int_slater
  USE ifmod, ONLY: write_messages
  USE global_constants
  USE global_read_data
  IMPLICIT NONE
  CONTAINS

  FUNCTION func_kronecker_delta(index_1, index_2) &
   &RESULT(kronecker_delta)
    ! Input variables
    INTEGER, INTENT(IN) :: index_1, index_2
    ! Output variable
    DOUBLE PRECISION :: kronecker_delta
    ! Initializing output variables
    kronecker_delta = 0.0D0
    IF(index_1 == index_2) kronecker_delta = 1.0D0
  END FUNCTION func_kronecker_delta


  !------------
  ! int_slater
  !------------

  FUNCTION func_int_slater_m_pq_n&
          &(ne, tuple_left, tuple_right, a, r, p, q, b, s) &      
   &RESULT(f)
    ! Constants
    CHARACTER(100) name_program 
    CHARACTER(100) type_program 
    ! Input variables
    INTEGER, INTENT(IN) :: ne
    INTEGER, INTENT(IN) :: tuple_left, tuple_right
    INTEGER, INTENT(IN) :: a(1:Max_n_tuple), r(1:Max_n_tuple)
    INTEGER, INTENT(IN) :: p, q
    INTEGER, INTENT(IN) :: b(1:Max_n_tuple), s(1:Max_n_tuple)
    ! Output variable
    DOUBLE PRECISION :: f

    ! Constants
    name_program = 'func_int_slater_m_pq_n'
    type_program = 'FUNCTION'
    ! Initializing output variables
    f = 0.0D0

    SELECT CASE(tuple_left)
      CASE(0)
        SELECT CASE(tuple_right)
          CASE(0)
            f = func_int_slater_0_pq_0&
               &(ne, p, q)
          CASE(1)
            f = func_int_slater_0_pq_1&
               &(ne, p, q, b(1), s(1))
          CASE DEFAULT
        END SELECT
      CASE(1)
        SELECT CASE(tuple_right)
          CASE(0)
            f = func_int_slater_1_pq_0&
               &(ne, a(1), r(1), p, q)
          CASE(1)
            f = func_int_slater_1_pq_1&
               &(ne, a(1), r(1), p, q, b(1), s(1))
          CASE DEFAULT
        END SELECT
      CASE DEFAULT
        !!CALL error_termination(type_program, name_program)
    END SELECT
  END FUNCTION func_int_slater_m_pq_n


  !----------------------------------
  ! < Phi_n | p'q | Phi_n > (n = 0-4)
  !----------------------------------

  ! < Phi_0 | p'q | Phi_0 >
  FUNCTION func_int_slater_0_pq_0(N, p, q) &
   &RESULT(f)
    ! Input variables
    INTEGER, INTENT(IN) :: N, p, q
    ! Output variable
    DOUBLE PRECISION :: f
    ! Initializing output variables
    f = 0.0D0
    IF(p <= N .AND. q <= N) f = func_kronecker_delta(p, q)
  END FUNCTION func_int_slater_0_pq_0


  ! < Phi_1 | p'q | Phi_1 >
  FUNCTION func_int_slater_1_pq_1(N, a, r, p, q, b, s) &
   &RESULT(f)
    ! Constants
    CHARACTER(72) name_program 
    CHARACTER(72) type_program 
    ! Input variables
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: a, r
    INTEGER, INTENT(IN) :: p, q
    INTEGER, INTENT(IN) :: b, s
    ! Output variable
    DOUBLE PRECISION :: f
    ! Local variables
    DOUBLE PRECISION r_p, r_s

    ! Constants
    name_program = 'func_int_slater_1_pq_1'
    type_program = 'FUNCTION'
    ! Initializing output variables
    f = 0.0D0
    ! Initializing local variables
    r_p = 0.0D0; r_s = 0.0D0

    !IF(a >= N+1 .OR. b >= N+1) THEN
      !!CALL write_error(31)
    !  WRITE(6,9001) 'a/b is greater than or euqal to ', N+1
      !!CALL error_termination(type_program, name_program)
    !ELSEIF(r <= N .OR. s <= N) THEN
      !!CALL write_error(31)
    !  WRITE(6,9001) 'r/s is less than or euqal to ', N
      !!CALL error_termination(type_program, name_program)
    !ELSE
    !ENDIF
    !9001 FORMAT(1X, A)

    r_p = func_kronecker_delta(r, p)
    r_s = func_kronecker_delta(r, s)
    f = r_p * func_int_slater_0_pq_1(N, a, q, b, s) &
     &+ r_s * func_int_slater_0_cpqd_0(N, a, p, q, b)
    !f = Identity(r, p) * Identity(s, q) * Identity(a, b) &
    ! &+ Identity(r, s) &
    ! &* (Identity(a, b) * Identity(p, q) - Identity(a, q) * Identity(p, b))
  END FUNCTION func_int_slater_1_pq_1




  !------------------------------------
  ! < Phi_n | p'q | Phi_n-1 > (n = 1-4)
  !------------------------------------

  ! < Phi_1 | p'q | Phi_0 >
  FUNCTION func_int_slater_1_pq_0(N, a, r, p, q) &
   &RESULT(f)
    ! Constants
    CHARACTER(72) name_program 
    CHARACTER(72) type_program 
    ! Input variables
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: a, r
    INTEGER, INTENT(IN) :: p, q
    ! Output variable
    DOUBLE PRECISION :: f
    ! Local variables
    DOUBLE PRECISION r_p

    ! Constants
    name_program = 'func_int_slater_1_pq_0'
    type_program = 'FUNCTION'
    ! Initializing output variables
    f = 0.0D0
    ! Initializing local variables
    r_p = 0.0D0

    !IF(a >= N+1) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 'a is greater than or euqal to ', N+1
      !CALL error_termination(type_program, name_program)
    !ELSEIF(r <= N) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 'r is less than or euqal to ', N
      !CALL error_termination(type_program, name_program)
    !ELSE
    !ENDIF
    !9001 FORMAT(1X, A)

    r_p = func_kronecker_delta(r, p)
    f = r_p * func_int_slater_0_pq_0(N, a, q)
  END FUNCTION func_int_slater_1_pq_0



  !------------------------------------
  ! < Phi_n-1 | p'q | Phi_n > (n = 1-4)
  !------------------------------------

  ! < Phi_0 | p'q | Phi_1 >
  FUNCTION func_int_slater_0_pq_1(N, p, q, b, s) &
   &RESULT(f)
    ! Constants
    CHARACTER(72) name_program 
    CHARACTER(72) type_program 
    ! Input variables
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: p, q
    INTEGER, INTENT(IN) :: b, s
    ! Output variable
    DOUBLE PRECISION :: f
    ! Local variables
    DOUBLE PRECISION s_q

    ! Constants
    name_program = 'func_int_slater_0_pq_1'
    type_program = 'FUNCTION'
    ! Initializing output variables
    f = 0.0D0
    ! Initializing local variables
    s_q = 0.0D0

    !IF(b >= N+1) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 'b is greater than or euqal to ', N+1
      !CALL error_termination(type_program, name_program)
    !ELSEIF(s <= N) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 's is less than or euqal to ', N
      !CALL error_termination(type_program, name_program)
    !ELSE
    !ENDIF
    !9001 FORMAT(1X, A)

    s_q = func_kronecker_delta(s, q)
    f = s_q * func_int_slater_0_pq_0(N, p, b)
  END FUNCTION func_int_slater_0_pq_1



  !------------------------------------
  ! < Phi_n | cp'qd | Phi_n > (n = 0-3)
  !------------------------------------

  ! < Phi_0 | c'p'qd | Phi_0 >
  FUNCTION func_int_slater_0_cpqd_0(N, c, p, q, d) &
   &RESULT(f)
    ! Constants
    CHARACTER(72) name_program 
    CHARACTER(72) type_program 
    ! Input variables
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: c, p, q, d
    ! Output variable
    DOUBLE PRECISION :: f
    ! Local variables
    DOUBLE PRECISION c_d, c_q

    ! Constants
    name_program = 'func_int_slater_0_cpqd_0'
    type_program = 'FUNCTION'
    ! Initializing output variables
    f = 0.0D0
    ! Initializing local variables
    c_d = 0.0D0; c_q = 0.0D0

    !IF(c >= N+1 .OR. d >= N+1) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 'c/d is greater than or euqal to ', N+1
      !CALL error_termination(type_program, name_program)
    !ELSE
    !ENDIF
    !9001 FORMAT(1X, A)

    c_d = func_kronecker_delta(c, d)
    c_q = func_kronecker_delta(c, q)
    f = c_d * func_int_slater_0_pq_0(N, p, q) &
     &- c_q * func_int_slater_0_pq_0(N, p, d)
  END FUNCTION func_int_slater_0_cpqd_0


  ! < Phi_1 | c'p'qd | Phi_1 >
  FUNCTION func_int_slater_1_cpqd_1(N, a, r, c, p, q, d, b, s)&
   &RESULT(f)
    ! Constants
    CHARACTER(72) name_program
    CHARACTER(72) type_program
    ! Input variables
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: a, r
    INTEGER, INTENT(IN) :: c, p, q, d
    INTEGER, INTENT(IN) :: b, s
    ! Output variable
    DOUBLE PRECISION :: f
    ! Local variables
    DOUBLE PRECISION c_d, c_q, c_b

    ! Constants
    name_program = 'func_int_slater_1_cpqd_1'
    type_program = 'FUNCTION'
    ! Initializing output variables
    f = 0.0D0
    ! Initializing local variables
    c_d = 0.0D0; c_q = 0.0D0; c_b = 0.0D0

    !IF(a >= N+1 .OR. b >= N+1 .OR. c >= N+1 .OR. d >= N+1) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 'a/b/c/d is greater than or euqal to ', N+1
      !CALL error_termination(type_program, name_program)
    !ELSEIF(r <= N .OR. s <= N) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 'r/s is less than or euqal to ', N
      !CALL error_termination(type_program, name_program)
    !ELSE
    !ENDIF
    !9001 FORMAT(1X, A)

    c_d = func_kronecker_delta(c, d)
    c_q = func_kronecker_delta(c, q)
    c_b = func_kronecker_delta(c, b)
    f = c_d * func_int_slater_1_pq_1&
            &(N, a, r, p, q, b, s) &
     &- c_q * func_int_slater_1_pq_1&
            &(N, a, r, p, d, b, s) &
     &- c_b * func_int_slater_1_pq_1&
            &(N, a, r, p, q, d, s)
  END FUNCTION func_int_slater_1_cpqd_1




  !--------------------------------------
  ! < Phi_n-1 | cp'qd | Phi_n > (n = 1-3)
  !--------------------------------------

  ! < Phi_0 | c'p'qd | Phi_1 >
  FUNCTION func_int_slater_0_cpqd_1(N, c, p, q, d, b, s) &
   &RESULT(f)
    ! Constants
    CHARACTER(72) name_program
    CHARACTER(72) type_program
    ! Input variables
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: c, p, q, d
    INTEGER, INTENT(IN) :: b, s
    ! Output variable
    DOUBLE PRECISION :: f
    ! Local variables
    DOUBLE PRECISION c_d, c_q, c_b

    ! Constants
    name_program = 'func_int_slater_1_cpqd_1'
    type_program = 'FUNCTION'
    ! Initializing output variables
    f = 0.0D0
    ! Initializing local variables
    c_d = 0.0D0; c_q = 0.0D0; c_b = 0.0D0
   
    !IF(b >= N+1 .OR. c >= N+1 .OR. d >= N+1) THEN
    !CALL write_error(31)
    !WRITE(6,9001) 'b/c/d is greater than or euqal to ', N+1
    !CALL error_termination(type_program, name_program)
    !ELSEIF(s <= N) THEN
    !CALL write_error(31)
    !WRITE(6,9001) 's is less than or euqal to ', N
    !CALL error_termination(type_program, name_program)
    !ELSE
    !ENDIF
    !9001 FORMAT(1X, A)
   
    c_d = func_kronecker_delta(c, d)
    c_q = func_kronecker_delta(c, q)
    c_b = func_kronecker_delta(c, b)
    f = c_d * func_int_slater_0_pq_1&
             &(N, p, q, b, s) &
     &- c_q * func_int_slater_0_pq_1&
             &(N, p, d, b, s) &
     &- c_b * func_int_slater_0_pq_1&
             &(N, p, q, d, s)
  END FUNCTION func_int_slater_0_cpqd_1


  !--------------------------------------
  ! < Phi_n | cp'qd | Phi_n-1 > (n = 1-3)
  !--------------------------------------

  ! < Phi_1 | c'p'qd | Phi_0 >
  FUNCTION func_int_slater_1_cpqd_0(N, a, r, c, p, q, d) &
   &RESULT(f)
    CHARACTER(72) name_program 
    CHARACTER(72) type_program 
    !--------------------------
    ! Input variables
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: a, r
    INTEGER, INTENT(IN) :: c, p, q, d
    ! Output variable
    DOUBLE PRECISION :: f
    ! Local variables
    DOUBLE PRECISION c_d, c_q

    ! Constants
    name_program = 'func_int_slater_1_cpqd_0'
    type_program = 'FUNCTION'
    ! Initializing output variables
    f = 0.0D0
    ! Initializing local variables
    c_d = 0.0D0; c_q = 0.0D0

    !IF(a >= N+1 .OR. c >= N+1 .OR. d >= N+1) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 'a/c/d is greater than or euqal to ', N+1
      !CALL error_termination(type_program, name_program)
    !ELSEIF(r <= N) THEN
      !CALL write_error(31)
    !  WRITE(6,9001) 'r is less than or euqal to ', N
      !CALL error_termination(type_program, name_program)
    !ELSE
    !ENDIF
    !9001 FORMAT(1X, A)

    c_d = func_kronecker_delta(c, d)
    c_q = func_kronecker_delta(c, q)
    f = c_d * func_int_slater_1_pq_0(N, a, r, p, q) &
     &- c_q * func_int_slater_1_pq_0(N, a, r, p, d)
  END FUNCTION func_int_slater_1_cpqd_0

END MODULE func_int_slater
