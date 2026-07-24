!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief A 'vertical W2' mass matrix used in the damping layer term.
!>
!> @details The kernel modifies the two rows of the velocity mass matrix
!> corresponding to the degrees of freedom of the vertical component of velocity
!> to account for Rayleigh damping in the absorbing layer region for use in runs
!> with non-flat bottom boundary.
!>
module semi_lump_kernel_mod

  use argument_mod,              only: arg_type,          &
                                       GH_OPERATOR,       &
                                       GH_SCALAR,         &
                                       GH_READ, GH_WRITE, &
                                       GH_REAL,           &
                                       GH_LOGICAL,        &
                                       CELL_COLUMN
  use constants_mod,             only: i_def, r_def, l_def
  use fs_continuity_mod,         only: W2
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------

  type, public, extends(kernel_type) :: semi_lump_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                       &
         arg_type(GH_OPERATOR, GH_REAL,    GH_WRITE, W2, W2), &
         arg_type(GH_OPERATOR, GH_REAL,    GH_READ,  W2, W2), &
         arg_type(GH_SCALAR,   GH_LOGICAL, GH_READ)           &
         /)
        integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: semi_lump_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: semi_lump_code

contains

  !> @brief Computes the reduced mass matrix for the damping layer term in the momentum equation.
  !!
  !! @param[in] cell     Identifying number of cell
  !! @param[in] nlayers  Number of layers
  !! @param[in] ncell_3d ncell*ndf
  !! @param[in,out] mm   Local stencil or mass matrix
  !! @param[in] chi1     1st coordinate field in Wchi
  !! @param[in] chi2     2nd coordinate field in Wchi
  !! @param[in] chi3     3rd coordinate field in Wchi
  !! @param[in] panel_id Field giving the ID for mesh panels
  !! @param[in] dl_base_height
  !!                     Base height of damping layer
  !! @param[in] dl_strength
  !!                     Strength of damping layer
  !! @param[in] domain_height
  !!                     The model domain height
  !! @param[in] radius   The planet radius
  !! @param[in] element_order_h The model finite element order in the horizontal
  !!                            direction
  !! @param[in] element_order_v The model finite element order in the vertical
  !!                            direction
  !! @param[in] dt       The model timestep length
  !! @param[in] ndf_w2   Degrees of freedom per cell
  !! @param[in] basis_w2 Vector basis functions evaluated at quadrature points.
  !! @param[in] ndf_chi  Degrees of freedom per cell for chi field
  !! @param[in] undf_chi Unique degrees of freedom for chi field
  !! @param[in] map_chi  Dofmap for the cell at the base of the column, for the
  !!                     space on which the chi field lives
  !! @param[in] basis_chi Vector basis functions evaluated at quadrature points
  !! @param[in] diff_basis_chi Vector differential basis functions evaluated at
  !!                           quadrature points
  !! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
  !! @param[in] undf_pid Number of unique degrees of freedom for panel_id
  !! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
  !! @param[in] nqp_h    Number of horizontal quadrature points
  !! @param[in] nqp_v    Number of vertical quadrature points
  !! @param[in] wqp_h    Horizontal quadrature weights
  !! @param[in] wqp_v    Vertical quadrature weights
  subroutine semi_lump_code(cell, nlayers,          & 
                            ncell_3d1, mm_lump,     &
                            ncell_3d2, mm,          &
                            hlump,                  &
                            ndf_w2)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: cell
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_3d1, ncell_3d2
    integer(kind=i_def), intent(in)    :: ndf_w2

    real(kind=r_def),    intent(inout) :: mm_lump(ncell_3d1,ndf_w2,ndf_w2)
    real(kind=r_def),    intent(in)    :: mm(ncell_3d2,ndf_w2,ndf_w2)

    logical(kind=l_def), intent(in)    :: hlump

    ! Internal variables
    integer(kind=i_def) :: df1, df2, k, ik, i1, i2

    ! Copy matrix
    do df2 = 1, ndf_w2
      do df1 = 1, ndf_w2
        do k = 1, nlayers
          ik = k + (cell-1)*nlayers
          mm_lump(ik,df1,df2) = mm(ik,df1,df2)
        end do
      end do
    end do

    ! Lump aligned terms
    if ( hlump ) then
      do k = 1, nlayers
        ik = k + (cell-1)*nlayers
        mm_lump(ik,1,1) = mm_lump(ik,1,1) + mm_lump(ik,1,3)
        mm_lump(ik,1,3) = 0.0_r_def
        mm_lump(ik,2,2) = mm_lump(ik,2,2) + mm_lump(ik,2,4)
        mm_lump(ik,2,4) = 0.0_r_def
        mm_lump(ik,3,3) = mm_lump(ik,3,3) + mm_lump(ik,3,1)
        mm_lump(ik,3,1) = 0.0_r_def
        mm_lump(ik,4,4) = mm_lump(ik,4,4) + mm_lump(ik,4,2)
        mm_lump(ik,4,2) = 0.0_r_def
      end do
   else
      if ( ndf_w2 == 6 ) then
        i1 = 5
        i2 = 6
      else
        ! Wtheta
        i1 = 1
        i2 = 2
      end if
      do k = 1, nlayers
        ik = k + (cell-1)*nlayers
        mm_lump(ik,i1,i1) = mm_lump(ik,i1,i1) + mm_lump(ik,i1,i2)
        mm_lump(ik,i1,i2) = 0.0_r_def
        mm_lump(ik,i2,i2) = mm_lump(ik,i2,i2) + mm_lump(ik,i2,i1)
        mm_lump(ik,i2,i1) = 0.0_r_def
      end do

   end if
    
  end subroutine semi_lump_code

end module semi_lump_kernel_mod
