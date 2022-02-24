!-----------------------------------------------------------------------
! Defines Precision and Constants
!-----------------------------------------------------------------------

module Constants

    implicit none

    ! Precision
    integer, parameter           :: rp = selected_real_kind(15)

    ! Constants
    real(kind=rp), parameter     :: pi = acos(-1.0_rp)
    real(kind=rp), parameter     :: R = 287.05_rp
    real(kind=rp), parameter     :: g = 9.81_rp
    real(kind=rp), parameter     :: Cp = 1004.0_rp
    real(kind=rp), parameter     :: Lv = 2.26e+06      ! Latent heat of Vapourization   (J/kg)

    ! Parameters
    integer, parameter           :: nlev = 100
    integer, parameter           :: MaxIter = 6
    real(kind=rp),parameter      :: tol = 5.0_rp
    ! integer, parameter           :: MaxIter = 10
    ! real(kind=rp),parameter      :: tol = 15.0_rp

    type :: properties
        real(kind=rp) :: p          ! Pressure                                    (hPa)
        real(kind=rp) :: T          ! Temperature                                 (K)
        real(kind=rp) :: q          ! Specific Humidity                           (gm/gm)
        real(kind=rp) :: z          ! Height                                      (m)
    end type properties

    type :: sounding
        real(kind=rp), dimension(nlev) :: p        ! Pressure                     (hPa)
        real(kind=rp), dimension(nlev) :: T        ! Temperature                  (K)
        real(kind=rp), dimension(nlev) :: Td       ! Dew point Temperature        (K)
        real(kind=rp), dimension(nlev) :: q        ! Specific Humidity            (gm/gm)
        real(kind=rp), dimension(nlev) :: qs       ! Saturation Specific Humidity (gm/gm)
        real(kind=rp), dimension(nlev) :: z        ! Height                       (m)
        !real(kind=rp), dimension(nlev) :: omega    ! dp/dt                       (Pa/s)
    end type sounding

end module Constants
