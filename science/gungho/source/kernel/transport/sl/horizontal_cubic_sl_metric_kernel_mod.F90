!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates the advective increments in x and y at time n+1 using
!!          cubic semi-Lagrangian transport.
!> @details This kernel using cubic interpolation to solve the one-dimensional
!!          advection equation in both x and y, giving advective increments
!!          in both directions. This is the second part of the COSMIC splitting,
!!          so the x increment works on the field previously advected in the
!!          y-direction(and vice versa).
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest order.

module horizontal_cubic_sl_metric_kernel_mod

  use argument_mod,          only: arg_type,                     &
                                   GH_FIELD, GH_REAL,            &
                                   CELL_COLUMN, GH_WRITE,        &
                                   GH_READ, GH_SCALAR,           &
                                   STENCIL, CROSS2D, GH_INTEGER, &
                                   ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,         only: r_tran, i_def, l_def
  use fs_continuity_mod,     only: W2H
  use kernel_mod,            only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_cubic_sl_metric_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  & ! increment
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,   &
                                                            STENCIL(CROSS2D)), & ! z
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2H)                         & ! wind
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_cubic_sl_metric_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: horizontal_cubic_sl_metric_code

contains

  !> @brief Compute advective transport in x and y directions using 1D
  !!        Semi-Lagrangian schemes, with a cubic reconstruction. This is the
  !!        "outer" step of a COSMIC splitting scheme.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] increment_x       Advective increment in x direction
  !> @param[in,out] increment_y       Advective increment in y direction
  !> @param[in]     field_x           Field from x direction
  !> @param[in]     stencil_sizes_x   Sizes of the branches of the cross stencil
  !> @param[in]     stencil_max_x     Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_x     Dofmap for the field_x stencil
  !> @param[in]     field_y           Field from y direction
  !> @param[in]     stencil_sizes_y   Sizes of the branches of the cross stencil
  !> @param[in]     stencil_max_y     Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_y     Dofmap for the field_y stencil
  !> @param[in]     dep_pts           Departure points
  !> @param[in]     monotone          Horizontal monotone option for cubic SL
  !> @param[in]     ndf_wf            Num of DoFs for field per cell
  !> @param[in]     undf_wf           Num of DoFs for this partition for field
  !> @param[in]     map_wf            Map for Wf
  !> @param[in]     ndf_w2h           Num of DoFs for W2H per cell
  !> @param[in]     undf_w2h          Num of DoFs for this partition for W2H
  !> @param[in]     map_w2h           Map for W2H
  subroutine horizontal_cubic_sl_metric_code( nlayers,         &
                                              increment,       &
                                              z,               &
                                              stencil_sizes,   &
                                              stencil_max,     &
                                              stencil_map,     &
                                              wind,            &
                                              ndf_wf,          &
                                              undf_wf,         &
                                              map_wf,          &
                                              ndf_w2h,         &
                                              undf_w2h,        &
                                              map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wf
    integer(kind=i_def), intent(in) :: ndf_wf
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_max
    integer(kind=i_def), intent(in) :: stencil_sizes(4)

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: stencil_map(ndf_wf,stencil_max,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: increment(undf_wf)
    real(kind=r_tran),   intent(in)    :: z(undf_wf)
    real(kind=r_tran),   intent(in)    :: wind(undf_w2h)

    ! Local scalars
    integer(kind=i_def) :: k, kp, km

    real(kind=r_tran)   :: dzdx, dzdy
    real(kind=r_tran)   :: z_l, z_r, up, um, vp, vm

    ! Interpolation coefficients
    real(kind=r_tran), parameter :: b0 = -1.0/6.0
    real(kind=r_tran), parameter :: b1 = 5.0/6.0
    real(kind=r_tran), parameter :: b2 = 2.0/6.0
    real(kind=r_tran), parameter :: c1 = 2.0/6.0
    real(kind=r_tran), parameter :: c2 = 5.0/6.0
    real(kind=r_tran), parameter :: c3 = -1.0/6.0
   
    do k = 0, nlayers
      km = max(0, k-1)
      kp = min(nlayers-1, k)

      um =  0.5_r_tran*(wind(map_w2h(1)+km) + wind(map_w2h(1)+kp))
      up =  0.5_r_tran*(wind(map_w2h(3)+km) + wind(map_w2h(3)+kp))
      vm = -0.5_r_tran*(wind(map_w2h(2)+km) + wind(map_w2h(2)+kp))
      vp = -0.5_r_tran*(wind(map_w2h(4)+km) + wind(map_w2h(4)+kp))

    ! dzdx
    ! Compute upwind Z on the left and right sides of the cell
    if ( um > 0.0_r_tran ) then
      z_l = b0*z(stencil_map(1,3,1)+k) + b1*z(stencil_map(1,2,1)+k) + b2*z(stencil_map(1,1,1)+k)
    else
      z_l = c1*z(stencil_map(1,2,1)+k) + c2*z(stencil_map(1,1,1)+k) + c3*z(stencil_map(1,2,3)+k)
    end if
    if ( up > 0.0_r_tran ) then
      z_r = b0*z(stencil_map(1,2,1)+k) + b1*z(stencil_map(1,1,1)+k) + b2*z(stencil_map(1,2,3)+k)
    else
      z_r = c1*z(stencil_map(1,1,3)+k) + c2*z(stencil_map(1,2,3)+k) + c3*z(stencil_map(1,3,3)+k)
    end if

    dzdx = (z_r- z_l)

    ! dzdy
    ! Compute upwind Z on the left and right sides of the cell
    if ( vm > 0.0_r_tran ) then
      z_l = b0*z(stencil_map(1,3,2)+k) + b1*z(stencil_map(1,2,2)+k) + b2*z(stencil_map(1,1,2)+k)
    else
      z_l = c1*z(stencil_map(1,2,2)+k) + c2*z(stencil_map(1,1,2)+k) + c3*z(stencil_map(1,2,4)+k)
    end if
    if ( vp > 0.0_r_tran ) then
      z_r = b0*z(stencil_map(1,2,2)+k) + b1*z(stencil_map(1,1,2)+k) + b2*z(stencil_map(1,2,4)+k)
    else
      z_r = c1*z(stencil_map(1,1,4)+k) + c2*z(stencil_map(1,2,4)+k) + c3*z(stencil_map(1,3,4)+k)
    end if

    dzdy = (z_r- z_l)

    increment(map_wf(1)+k) = 0.5_r_tran*(um+up)*dzdx + 0.5_r_tran*(vm+vp)*dzdy
  end do       

  end subroutine horizontal_cubic_sl_metric_code

end module horizontal_cubic_sl_metric_kernel_mod
