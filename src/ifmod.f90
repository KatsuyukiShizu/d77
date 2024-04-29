! This interface module is part of d77.

! d77 is free software and can be redistributed and/or modified
! under the terms of the GNU General Public License v3.0
! as published by the Free Software Foundation.
! https://www.gnu.org/licenses/gpl-3.0.html

! For bug reports, e-mail to shizu@scl.kyoto-u.ac.jp

MODULE ifmod
  INTERFACE

    SUBROUTINE calc_int_cgf
    END SUBROUTINE calc_int_cgf
    
    SUBROUTINE rho
    END SUBROUTINE rho
    
    SUBROUTINE dipole_int_pgf
    END SUBROUTINE dipole_int_pgf

    SUBROUTINE dipole_density
    END SUBROUTINE dipole_density
    
    SUBROUTINE vc_int_pgf
    END SUBROUTINE vc_int_pgf
    
    SUBROUTINE vc_density
    END SUBROUTINE vc_density
    
    SUBROUTINE soc_int_pgf
    END SUBROUTINE soc_int_pgf

    SUBROUTINE calc_def&
          &(xyz_a, xyz_b, xyz_p, &
           &lmn_a, lmn_b, cntexp_p, &
           &d, e, f)
      DOUBLE PRECISION, INTENT(IN)  :: xyz_a(1:3), xyz_b(1:3), xyz_p(1:3)
      INTEGER, INTENT(IN)           :: lmn_a(1:3), lmn_b(1:3)
      DOUBLE PRECISION, INTENT(IN)  :: cntexp_p
      DOUBLE PRECISION, INTENT(OUT) :: d(0:lmn_a(1)+lmn_b(1), 0:lmn_a(1), 0:lmn_b(1))
      DOUBLE PRECISION, INTENT(OUT) :: e(0:lmn_a(2)+lmn_b(2), 0:lmn_a(2), 0:lmn_b(2))
      DOUBLE PRECISION, INTENT(OUT) :: f(0:lmn_a(3)+lmn_b(3), 0:lmn_a(3), 0:lmn_b(3))
    END SUBROUTINE calc_def

    SUBROUTINE calc_rlmn_0&
              &(l, m, n, max_j, &
               &a, b, c, t, &
               &alpha, &
               &rlmn)
      DOUBLE PRECISION, INTENT(IN)  :: a, b, c, t
      DOUBLE PRECISION, INTENT(IN)  :: alpha
      INTEGER, INTENT(IN)           :: l, m, n, max_j
      DOUBLE PRECISION, INTENT(OUT) :: rlmn
    END SUBROUTINE calc_rlmn_0

!   ----------------
!   Density matrices
!   ----------------

    SUBROUTINE calc_opr_dmat(left_state, right_state, n_ci_ci, opr_dmat, trace_opr_dmat)
      INTEGER, INTENT(IN)           :: left_state, right_state
      INTEGER, INTENT(OUT)          :: n_ci_ci
      DOUBLE PRECISION, INTENT(OUT) :: opr_dmat(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: trace_opr_dmat
    END SUBROUTINE calc_opr_dmat

    SUBROUTINE calc_opr_dmat_td(left_state, right_state, n_ci_ci, opr_dmat, trace_opr_dmat)
      INTEGER, INTENT(IN)           :: left_state, right_state
      INTEGER, INTENT(OUT)          :: n_ci_ci
      DOUBLE PRECISION, INTENT(OUT) :: opr_dmat(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: trace_opr_dmat
    END SUBROUTINE calc_opr_dmat_td

    SUBROUTINE calc_dmat_cgf(opr_dmat, dmat_cgf)
      DOUBLE PRECISION, INTENT(IN)  :: opr_dmat(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: dmat_cgf(:,:)
    END SUBROUTINE calc_dmat_cgf

    SUBROUTINE calc_dmat_soc_cgf(opr_dmat, dmat_soc_cgf)
      DOUBLE PRECISION, INTENT(IN)  :: opr_dmat(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: dmat_soc_cgf(:,:,:)
    END SUBROUTINE calc_dmat_soc_cgf

!   ----------
!   Cube files
!   ----------

    SUBROUTINE cube_rho(dmat_cgf, int_rho)
      DOUBLE PRECISION, INTENT(IN)  :: dmat_cgf(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: int_rho
    END SUBROUTINE cube_rho

    SUBROUTINE cube_dipole_density(projection, dmat_cgf, int_rho, int_dd, int_dd_p)
      CHARACTER(LEN=4), INTENT(IN)    :: projection
      DOUBLE PRECISION, INTENT(IN)    :: dmat_cgf(:,:)
      DOUBLE PRECISION, INTENT(INOUT) :: int_rho
      DOUBLE PRECISION, INTENT(INOUT) :: int_dd(1:3), int_dd_p
    END SUBROUTINE cube_dipole_density

    SUBROUTINE cube_vc_density(dmat_cgf, int_rho, int_vcd)
      DOUBLE PRECISION, INTENT(IN)  :: dmat_cgf(:,:)
      DOUBLE PRECISION, INTENT(OUT) :: int_rho
      DOUBLE PRECISION, INTENT(OUT) :: int_vcd(:)
    END SUBROUTINE cube_vc_density

!   ----------

    SUBROUTINE write_messages(message_number, text, type_program, name_program)
      INTEGER, INTENT(IN)            :: message_number
      CHARACTER(LEN=100), INTENT(IN) :: text, type_program, name_program
    END SUBROUTINE write_messages

    SUBROUTINE write_g16_log_format(fid, data)
      INTEGER, INTENT(IN)          :: fid
      DOUBLE PRECISION, INTENT(IN) :: data(:,:)
    END SUBROUTINE write_g16_log_format

    SUBROUTINE write_mat(fid, data)
      INTEGER, INTENT(IN)          :: fid
      DOUBLE PRECISION, INTENT(IN) :: data(:,:)
    END SUBROUTINE write_mat

  END INTERFACE
END MODULE ifmod
