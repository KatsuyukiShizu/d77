! This subroutine is part of d77 and computes density.   

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE func_rho
  USE ifmod, ONLY: write_messages
  USE global_read_input
  USE global_read_data
  IMPLICIT NONE
  CONTAINS


  FUNCTION func_phi_phi_cgf(xyz, num_cgf) &
   &RESULT(f)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'func_phi_phi_cgf'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'FUNCTION'

!   Input variable
    DOUBLE PRECISION, INTENT(IN) :: xyz(:)
    INTEGER, INTENT(IN)          :: num_cgf

!   Output variable
    DOUBLE PRECISION :: f(1:num_cgf, 1:num_cgf)

!   Local variables
    INTEGER          :: i_cgf, j_cgf, i_pgf, j_pgf
    DOUBLE PRECISION :: dxyz(1:3), dxyz_i(1:3), dxyz_j(1:3)
    DOUBLE PRECISION :: dxyz2, dxyz_i2, dxyz_j2
    DOUBLE PRECISION :: distance
    DOUBLE PRECISION :: alpha, beta, alpha_plus_beta 
    DOUBLE PRECISION :: e_ab
    DOUBLE PRECISION :: xyz_ab(1:3)
    DOUBLE PRECISION :: xyz_ab2 ! Distance between CGF centers
    DOUBLE PRECISION :: xyz_p_ab(1:3)
    
    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)
    
    f = 0.0D0
    dxyz = 0.0D0; dxyz_i = 0.0D0; dxyz_j = 0.0D0
    dxyz2 = 0.0D0; dxyz_i2 = 0.0D0; dxyz_j2 = 0.0D0
    distance = 0.0D0

    DO i_cgf = 1, num_cgf - 1
      dxyz_i = xyz - Xyzao(1:3, i_cgf)
      dxyz_i2 = DOT_PRODUCT(dxyz_i, dxyz_i)
      DO j_cgf = i_cgf + 1, num_cgf
        dxyz_j = xyz - Xyzao(1:3, j_cgf)
        dxyz_j2 = DOT_PRODUCT(dxyz_j, dxyz_j)
        distance = DMIN1(DSQRT(dxyz_i2), DSQRT(dxyz_j2))
        IF(distance > Threshold_distance) CYCLE

        xyz_ab = xyzao(1:3, j_cgf) - xyzao(1:3, i_cgf)
        xyz_ab2 = DOT_PRODUCT(xyz_ab, xyz_ab)
        DO i_pgf = 1, N_pgf(i_cgf)
          alpha = cntexp(i_pgf, i_cgf)
          DO j_pgf = 1, N_pgf(j_cgf)
            beta = cntexp(j_pgf, j_cgf)
            alpha_plus_beta = alpha + beta
            e_ab = DEXP(-xyz_ab2*alpha*beta/alpha_plus_beta)

            IF(e_ab > Threshold_e_ab) THEN
              xyz_p_ab &
             &= (alpha * xyzao(1:3, i_cgf) + beta * xyzao(1:3, j_cgf)) &
             &/ alpha_plus_beta
              dxyz = xyz - xyz_p_ab
              dxyz2 = DOT_PRODUCT(dxyz, dxyz)
              f(i_cgf, j_cgf) &
             &= f(i_cgf, j_cgf) &
             &+ cntcoef(i_pgf, i_cgf) &
             &* cntcoef(j_pgf, j_cgf) &
             &* cntnormfac(i_pgf, i_cgf) &
             &* cntnormfac(j_pgf, j_cgf) &
             &* e_ab &
             &* DEXP(-alpha_plus_beta*dxyz2)
            ELSE
            ENDIF

          ENDDO ! j_pgf
        ENDDO ! i_pgf

        f(i_cgf, j_cgf) &
        &= f(i_cgf, j_cgf) &
        &* (dxyz_i(1)**DBLE(Lmn(1, i_cgf))) &
        &* (dxyz_i(2)**DBLE(Lmn(2, i_cgf))) &
        &* (dxyz_i(3)**DBLE(Lmn(3, i_cgf))) &
        &* (dxyz_j(1)**DBLE(Lmn(1, j_cgf))) &
        &* (dxyz_j(2)**DBLE(Lmn(2, j_cgf))) &
        &* (dxyz_j(3)**DBLE(Lmn(3, j_cgf)))
        f(j_cgf, i_cgf) =  f(i_cgf, j_cgf)

      ENDDO ! j_cgf
    ENDDO ! i_cgf
    
    DO i_cgf = 1, num_cgf
      dxyz_i = xyz - Xyzao(1:3, i_cgf)
      dxyz_i2 = DOT_PRODUCT(dxyz_i, dxyz_i)
      dxyz2 = dxyz_i2
      DO i_pgf = 1, N_pgf(i_cgf)
        alpha = cntexp(i_pgf, i_cgf)
        DO j_pgf = 1, N_pgf(i_cgf)
          beta = cntexp(j_pgf, i_cgf)
          alpha_plus_beta = alpha + beta

          f(i_cgf, i_cgf) &
         &= f(i_cgf, i_cgf) &
         &+ cntcoef(i_pgf, i_cgf) &
         &* cntcoef(j_pgf, i_cgf) &
         &* cntnormfac(i_pgf, i_cgf) &
         &* cntnormfac(j_pgf, i_cgf) &
         &* DEXP(-alpha_plus_beta*dxyz2)
        ENDDO
      ENDDO

      f(i_cgf, i_cgf) &
     &= f(i_cgf, i_cgf) &
     &* dxyz_i(1)**(2*Lmn(1, i_cgf)) &
     &* dxyz_i(2)**(2*Lmn(2, i_cgf)) &
     &* dxyz_i(3)**(2*Lmn(3, i_cgf))
    ENDDO
    
    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  
  END FUNCTION func_phi_phi_cgf


  FUNCTION func_rho_0(dmat_cgf, phi_phi_cgf) &
   &RESULT(rho)
    USE global_read_data
    IMPLICIT NONE
    
!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'func_rho_0'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'FUNCTION'
    
!   Input variables
    DOUBLE PRECISION, INTENT(IN) :: dmat_cgf(:,:), phi_phi_cgf(:,:)
    
!   Output variable
    DOUBLE PRECISION :: rho
    
!   Local variables
    INTEGER :: num_cgf
    INTEGER :: i_cgf, j_cgf
    
    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)
    
    num_cgf = 0
    num_cgf = SIZE(dmat_cgf,1)
    rho = 0.0D0
    DO i_cgf = 1, num_cgf - 1
      DO j_cgf = i_cgf + 1, num_cgf
        rho = rho&
        + (dmat_cgf(i_cgf,j_cgf) + dmat_cgf(j_cgf,i_cgf))&
        * phi_phi_cgf(i_cgf,j_cgf)
      ENDDO
    ENDDO
    DO i_cgf = 1, num_cgf
      rho = rho + dmat_cgf(i_cgf, i_cgf) * phi_phi_cgf(i_cgf, i_cgf)
    ENDDO
    
    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  
  END FUNCTION func_rho_0


! --------------
! Kinetic energy
! --------------

  FUNCTION func_grad_phi_grad_phi(xyz, num_cgf) RESULT(f)

!   ------------------------
!   Declaration of variables
!   ------------------------

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'func_grad_phi_grad_phi'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'FUNCTION'

!   Input variable
    DOUBLE PRECISION, INTENT(IN) :: xyz(:)
    INTEGER, INTENT(IN)          :: num_cgf

!   Output variable
    DOUBLE PRECISION :: f(1:num_cgf, 1:num_cgf, 1:3)

!   Local variables
    INTEGER :: i_cgf, j_cgf, i_pgf, j_pgf
    DOUBLE PRECISION :: dxyz(1:3), dxyz_i(1:3), dxyz_j(1:3)
    DOUBLE PRECISION :: dxyz2, dxyz_i2, dxyz_j2
    DOUBLE PRECISION :: distance
    DOUBLE PRECISION :: alpha, beta, alpha_plus_beta 
    DOUBLE PRECISION :: xyz_p_ab(1:3)
    DOUBLE PRECISION :: e_ab
    DOUBLE PRECISION :: xyz_ab(1:3)
    DOUBLE PRECISION :: xyz_ab2 ! Distance between CGF centers
    DOUBLE PRECISION :: const

    INTEGER :: l_a, m_a, n_a, &
              &l_b, m_b, n_b
    DOUBLE PRECISION :: alpha_a_2, alpha_b_2, alpha_ab_4
    INTEGER :: lmn_a_temp(1:3), lmn_b_temp(1:3)
    DOUBLE PRECISION :: temp

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)
    
    f = 0.0D0
    dxyz = 0.0D0; dxyz_i = 0.0D0; dxyz_j = 0.0D0
    dxyz2 = 0.0D0; dxyz_i2 = 0.0D0; dxyz_j2 = 0.0D0
    distance = 0.0D0
    const = 0.0D0

    DO i_cgf = 1, num_cgf - 1
      l_a = Lmn(1, i_cgf)
      m_a = Lmn(2, i_cgf)
      n_a = Lmn(3, i_cgf)
      dxyz_i = xyz - Xyzao(1:3, i_cgf)
      dxyz_i2 = DOT_PRODUCT(dxyz_i, dxyz_i)
      DO j_cgf = i_cgf + 1, num_cgf
        l_b = Lmn(1, j_cgf)
        m_b = Lmn(2, j_cgf)
        n_b = Lmn(3, j_cgf)
        dxyz_j = xyz - Xyzao(1:3, j_cgf)
        dxyz_j2 = DOT_PRODUCT(dxyz_j, dxyz_j)
        distance = DMIN1(DSQRT(dxyz_i2), DSQRT(dxyz_j2))
        IF(distance > Threshold_distance) CYCLE

        xyz_ab = xyzao(1:3, j_cgf) - xyzao(1:3, i_cgf)
        xyz_ab2 = DOT_PRODUCT(xyz_ab, xyz_ab)
        DO i_pgf = 1, N_pgf(i_cgf)
          alpha = cntexp(i_pgf, i_cgf)
          DO j_pgf = 1, N_pgf(j_cgf)
            beta = cntexp(j_pgf, j_cgf)
            alpha_plus_beta = alpha + beta
            e_ab = DEXP(-xyz_ab2*alpha*beta/alpha_plus_beta)

            IF(e_ab > Threshold_e_ab) THEN
              xyz_p_ab &
             &= (alpha * xyzao(1:3, i_cgf) + beta * xyzao(1:3, j_cgf)) &
             &/ alpha_plus_beta

              dxyz = xyz - xyz_p_ab
              dxyz2 = DOT_PRODUCT(dxyz, dxyz)
              alpha_a_2 = 2.0D0*alpha
              alpha_b_2 = 2.0D0*beta
              alpha_ab_4 = alpha_a_2 * alpha_b_2
              const = cntcoef(i_pgf, i_cgf) &
                   &* cntcoef(j_pgf, j_cgf) &
                   &* cntnormfac(i_pgf, i_cgf) &
                   &* cntnormfac(j_pgf, j_cgf) &
                   &* e_ab &
                   &* DEXP(-alpha_plus_beta*dxyz2)

!             -----------------------------------------
!             Calculating x component of kinetic energy
!             -----------------------------------------

              temp = 0.0D0
              lmn_a_temp(1:3) = Lmn(1:3, i_cgf) 
              lmn_b_temp(1:3) = Lmn(1:3, j_cgf) 

              lmn_a_temp(1) = l_a - 1
              lmn_b_temp(1) = l_b - 1
              IF(l_a >= 1 .AND. l_b >= 1) &
             &temp = temp + DBLE(l_a*l_b) &
             &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 1)

              lmn_a_temp(1) = l_a - 1
              lmn_b_temp(1) = l_b + 1
              IF(l_a >= 1) &
             &temp = temp - alpha_b_2 * DBLE(l_a) &
                  &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 1)
        
              lmn_a_temp(1) = l_a + 1
              lmn_b_temp(1) = l_b - 1
              IF(l_b >= 1) &
             &temp = temp - alpha_a_2 * DBLE(l_b) &
                  &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 1) 
        
              lmn_a_temp(1) = l_a + 1
              lmn_b_temp(1) = l_b + 1
              temp = temp + alpha_ab_4 &
                  &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 1) 
        
              temp = temp * const
              f(i_cgf, j_cgf, 1) = f(i_cgf, j_cgf, 1) + temp

!             -----------------------------------------
!             Calculating y component of Kinetic energy
!             -----------------------------------------
        
              temp = 0.0D0
              lmn_a_temp(1:3) = Lmn(1:3, i_cgf)
              lmn_b_temp(1:3) = Lmn(1:3, j_cgf)
        
              lmn_a_temp(2) = m_a - 1
              lmn_b_temp(2) = m_b - 1
              IF(m_a >= 1 .AND. m_b >= 1) &
             &temp = temp + DBLE(m_a*m_b) &
                  &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 2)
        
              lmn_a_temp(2) = m_a - 1
              lmn_b_temp(2) = m_b + 1
              IF(m_a >= 1) &
             &temp = temp - alpha_b_2 * DBLE(m_a) &
                  &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 2)
        
              lmn_a_temp(2) = m_a + 1
              lmn_b_temp(2) = m_b - 1
              IF(m_b >= 1) &
             &temp = temp - alpha_a_2 * DBLE(m_b) &
                  &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 2)
        
              lmn_a_temp(2) = m_a + 1
              lmn_b_temp(2) = m_b + 1
              temp = temp + alpha_ab_4 &
                  &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 2)
        
              temp = temp * const
              f(i_cgf, j_cgf, 2) = f(i_cgf, j_cgf, 2) + temp
        
!             -----------------------------------------
!             Calculating z component of Kinetic energy
!             -----------------------------------------

              temp = 0.0D0 
              lmn_a_temp(1:3) = Lmn(1:3, i_cgf)
              lmn_b_temp(1:3) = Lmn(1:3, j_cgf)
         
              lmn_a_temp(3) = n_a - 1
              lmn_b_temp(3) = n_b - 1
              IF(n_a >= 1 .AND. n_b >= 1) &
             &temp = temp + DBLE(n_a*n_b) &
             &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 3)
         
              lmn_a_temp(3) = n_a - 1
              lmn_b_temp(3) = n_b + 1
              IF(n_a >= 1) &
             &temp = temp - alpha_b_2 * DBLE(n_a) &
             &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 3)
         
              lmn_a_temp(3) = n_a + 1
              lmn_b_temp(3) = n_b - 1
              IF(n_b >= 1) &
             &temp = temp - alpha_a_2 * DBLE(n_b) &
             &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 3)
         
              lmn_a_temp(3) = n_a + 1
              lmn_b_temp(3) = n_b + 1
              temp = temp + alpha_ab_4 &
             &* func_xyz_ke(dxyz_i, dxyz_j, lmn_a_temp, lmn_b_temp, 3)
         
              temp = temp * const
              f(i_cgf, j_cgf, 3) = f(i_cgf, j_cgf, 3) + temp

              !f(i_cgf, j_cgf, 1:3) = f(i_cgf, j_cgf, 1:3) * const
            ELSE
            ENDIF
          ENDDO ! j_pgf
        ENDDO ! i_pgf
        
        f(i_cgf, j_cgf, 1) &
       &= f(i_cgf, j_cgf, 1) & 
       &* (dxyz_i(2)**DBLE(m_a)) &
       &* (dxyz_i(3)**DBLE(n_a)) &
       &* (dxyz_j(2)**DBLE(m_b)) &
       &* (dxyz_j(3)**DBLE(n_b)) 
     
        f(i_cgf, j_cgf, 2) &
       &= f(i_cgf, j_cgf, 2) &
       &* (dxyz_i(3)**DBLE(n_a)) &
       &* (dxyz_i(1)**DBLE(l_a)) &
       &* (dxyz_j(3)**DBLE(n_b)) &
       &* (dxyz_j(1)**DBLE(l_b))
     
        f(i_cgf, j_cgf, 3) &
       &= f(i_cgf, j_cgf, 3) &
       &* (dxyz_i(1)**DBLE(l_a)) &
       &* (dxyz_i(2)**DBLE(m_a)) &
       &* (dxyz_j(1)**DBLE(l_b)) &
       &* (dxyz_j(2)**DBLE(m_b))
  
        f(j_cgf, i_cgf, 1:3) = f(i_cgf, j_cgf, 1:3)

      ENDDO ! j_cgf
    ENDDO ! i_cgf

    DO i_cgf = 1, Num_cgf
      dxyz_i = xyz - Xyzao(1:3, i_cgf)
      dxyz_i2 = DOT_PRODUCT(dxyz_i, dxyz_i)
      distance = DSQRT(dxyz_i2)
      dxyz2 = dxyz_i2
      IF(distance > Threshold_distance) CYCLE

      DO i_pgf = 1, N_pgf(i_cgf)
        alpha = cntexp(i_pgf, i_cgf)
        DO j_pgf = 1, N_pgf(i_cgf)
          beta = cntexp(j_pgf, i_cgf)
          alpha_plus_beta = alpha + beta

          alpha_a_2 = 2.0D0*alpha
          alpha_b_2 = 2.0D0*beta
          alpha_ab_4 = alpha_a_2 * alpha_b_2

          const = cntcoef(i_pgf, i_cgf) &
               &* cntcoef(j_pgf, i_cgf) &
               &* cntnormfac(i_pgf, i_cgf) &
               &* cntnormfac(j_pgf, i_cgf) &
               &* DEXP(-alpha_plus_beta*dxyz2)

!         -----------------------------------------
!         Calculating x component of kinetic energy
!         -----------------------------------------

          temp = 0.0D0
          lmn_a_temp(1:3) = Lmn(1:3, i_cgf)
          lmn_b_temp(1:3) = Lmn(1:3, i_cgf)
      
          lmn_a_temp(1) = l_a - 1
          lmn_b_temp(1) = l_b - 1
          IF(l_a >= 1 .AND. l_b >= 1) &
         &temp = temp + DBLE(l_a*l_b) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 1)
      
          lmn_a_temp(1) = l_a - 1
          lmn_b_temp(1) = l_b + 1
          IF(l_a >= 1) &
         &temp = temp - alpha_b_2 * DBLE(l_a) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 1)
      
          lmn_a_temp(1) = l_a + 1
          lmn_b_temp(1) = l_b - 1
          IF(l_b >= 1) &
         &temp = temp - alpha_a_2 * DBLE(l_b) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 1)
      
          lmn_a_temp(1) = l_a + 1
          lmn_b_temp(1) = l_b + 1
          temp = temp + alpha_ab_4 &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 1)
      
          temp = temp * const
          f(i_cgf, i_cgf, 1) = f(i_cgf, i_cgf, 1) + temp

!         -----------------------------------------
!         Calculating y component of kinetic energy
!         -----------------------------------------
      
          temp = 0.0D0
          lmn_a_temp(1:3) = Lmn(1:3, i_cgf)
          lmn_b_temp(1:3) = Lmn(1:3, i_cgf)
      
          lmn_a_temp(2) = m_a - 1
          lmn_b_temp(2) = m_b - 1
          IF(m_a >= 1 .AND. m_b >= 1) &
         &temp = temp + DBLE(m_a*m_b) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 2)
      
          lmn_a_temp(2) = m_a - 1
          lmn_b_temp(2) = m_b + 1
          IF(m_a >= 1) &
         &temp = temp - alpha_b_2 * DBLE(m_a) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 2)
      
          lmn_a_temp(2) = m_a + 1
          lmn_b_temp(2) = m_b - 1
          IF(m_b >= 1) &
         &temp = temp - alpha_a_2 * DBLE(m_b) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 2)
      
          lmn_a_temp(2) = m_a + 1
          lmn_b_temp(2) = m_b + 1
          temp = temp + alpha_ab_4 &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 2)
      
          temp = temp * const
          f(i_cgf, i_cgf, 2) = f(i_cgf, i_cgf, 2) + temp
      
      
!         -----------------------------------------
!         Calculating z component of kinetic energy
!         -----------------------------------------
      
          temp = 0.0D0
          lmn_a_temp(1:3) = Lmn(1:3, i_cgf)
          lmn_b_temp(1:3) = Lmn(1:3, i_cgf)
      
          lmn_a_temp(3) = n_a - 1
          lmn_b_temp(3) = n_b - 1
          IF(n_a >= 1 .AND. n_b >= 1) &
         &temp = temp + DBLE(n_a*n_b) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 3)
      
          lmn_a_temp(3) = n_a - 1
          lmn_b_temp(3) = n_b + 1
          IF(n_a >= 1) &
         &temp = temp - alpha_b_2 * DBLE(n_a) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 3)
      
          lmn_a_temp(3) = n_a + 1
          lmn_b_temp(3) = n_b - 1
          IF(n_b >= 1) &
         &temp = temp - alpha_a_2 * DBLE(n_b) &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 3)
      
          lmn_a_temp(3) = n_a + 1
          lmn_b_temp(3) = n_b + 1
          temp = temp + alpha_ab_4 &
         &* func_xyz_ke(dxyz_i, dxyz_i, lmn_a_temp, lmn_b_temp, 3)
      
      
          temp = temp * const
          f(i_cgf, i_cgf, 3) = f(i_cgf, i_cgf, 3) + temp

        ENDDO ! j_pgf
      ENDDO ! i_pgf

      f(i_cgf, i_cgf, 1) &
     &= f(i_cgf, i_cgf, 1) &
     &* (dxyz_i(2)**DBLE(2*m_a)) &
     &* (dxyz_i(3)**DBLE(2*n_a)) 

      f(i_cgf, i_cgf, 2) &
     &= f(i_cgf, i_cgf, 2) &
     &* (dxyz_i(3)**DBLE(2*n_a)) &
     &* (dxyz_i(1)**DBLE(2*l_a)) 

      f(i_cgf, i_cgf, 3) &
     &= f(i_cgf, i_cgf, 3) &
     &* (dxyz_i(1)**DBLE(2*l_a)) &
     &* (dxyz_i(2)**DBLE(2*m_a)) 
    ENDDO ! i_cgf

    f = 0.50D0*f

    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  
  END FUNCTION func_grad_phi_grad_phi

  FUNCTION func_ked(num_cgf, dmat_cgf, grad_phi_grad_phi_cgf) &
   &RESULT(ked)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'func_ked'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'FUNCTION'

!   Input variables
    INTEGER, INTENT(IN)          :: num_cgf
    DOUBLE PRECISION, INTENT(IN) :: dmat_cgf(1:num_cgf, 1:num_cgf), &
                                   &grad_phi_grad_phi_cgf(1:num_cgf,1:num_cgf,1:3)

!   Output variable
    DOUBLE PRECISION :: ked(1:3)

!   DO loop variables
    INTEGER i_cgf, j_cgf, i_xyz

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)
    
    ked = 0.0D0
    DO i_cgf = 1, num_cgf - 1
      DO j_cgf = i_cgf + 1, num_cgf
        DO i_xyz = 1, 3
          ked(i_xyz) = ked(i_xyz) &
         &+ dmat_cgf(i_cgf,j_cgf) * grad_phi_grad_phi_cgf(i_cgf,j_cgf,i_xyz) &
         &+ dmat_cgf(j_cgf,i_cgf) * grad_phi_grad_phi_cgf(j_cgf,i_cgf,i_xyz)
        ENDDO
      ENDDO
    ENDDO
    DO i_cgf = 1, num_cgf
      DO i_xyz = 1, 3
        ked(i_xyz) = ked(i_xyz) &
       &+ dmat_cgf(i_cgf,i_cgf) * grad_phi_grad_phi_cgf(i_cgf,i_cgf,i_xyz)
      ENDDO
    ENDDO
    
    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  
  END FUNCTION func_ked


  FUNCTION func_xyz_ke(dxyz_a, dxyz_b, lmn_a, lmn_b, i_xyz) &
  RESULT(f)

!   Program name
    CHARACTER(LEN=100), PARAMETER :: name_program = 'func_xyz_ke'
    CHARACTER(LEN=100), PARAMETER :: type_program = 'FUNCTION'

!   Input variables
    DOUBLE PRECISION, INTENT(IN) :: dxyz_a(1:3), dxyz_b(1:3)
    INTEGER, INTENT(IN) :: lmn_a(1:3), lmn_b(1:3)
    INTEGER, INTENT(IN) :: i_xyz

!   Output variable
    DOUBLE PRECISION :: f

    IF(Debug == 'Yes') &
   &CALL write_messages(2, Text_blank, type_program, name_program)
    
    SELECT CASE(i_xyz)
      CASE(1)
        f = dxyz_a(1)**DBLE(lmn_a(1)) &
         &* dxyz_b(1)**DBLE(lmn_b(1)) 
      CASE(2)
        f = dxyz_a(2)**DBLE(lmn_a(2)) &
         &* dxyz_b(2)**DBLE(lmn_b(2)) 
      CASE(3)
        f = dxyz_a(3)**DBLE(lmn_a(3)) &
         &* dxyz_b(3)**DBLE(lmn_b(3)) 
      CASE DEFAULT
        STOP
    END SELECT
    
    IF(Debug == 'Yes') &
   &CALL write_messages(3, Text_blank, type_program, name_program)
  
  END FUNCTION func_xyz_ke

END MODULE func_rho
