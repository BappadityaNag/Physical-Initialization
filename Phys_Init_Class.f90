!-----------------------------------------------------------------------
! Module : Physical Initialization Class
! Author : Bappaditya Nag
! Date   : 08/18/2017
!
! Contains
! 1. Disp_of_SuperSat
! 2. Moist_Adiabat
! 3. Tetens
! 4. Read_Var
! 5. Write_Var
!-----------------------------------------------------------------------

module Phys_Init_Class
implicit none
contains

!-----------------------------------------------------------------------
! Computes Disposition of Super-saturation
!
! Usage : rain_mod = Disp_of_SuperSat(real rain_mod, real rain_GPM, 
!                    real qv, real qcl, real pres, real T, real znu,
!                    real z, real omg)
!-----------------------------------------------------------------------
     
     subroutine Disp_of_SuperSat (rain_mod,rain_GPM,qv,qcl,pres,T,znu,z,omg)
        use Constants
        implicit none

        real(kind=rp), dimension(1:,1:,1:)            :: rain_mod, rain_GPM
        real(kind=rp), dimension(1:,1:,1:)            :: z, qv, qcl, pres, T, omg
        real(kind=rp), dimension(1:ubound(T,3))       :: znu, p, zm, Tm, qm
        integer                                       :: nlevs, nlats, nlons
        integer                                       :: ilat, ilon, ilev
        real(kind=rp)                                 :: dT, alpha
        
        nlons = ubound(T,1); nlats = ubound(T,2); nlevs = ubound(T,3)
        dT = 0.01_rp

        do ilat = 1, nlats
            do ilon = 1, nlons
                do ilev = 1, nlevs
                    p(ilev) = pres(ilon,ilat,nlevs) + znu(ilev) * (pres(ilon,ilat,1)-pres(ilon,ilat,nlevs))
                enddo
                call Moist_Adiabat (z(ilon,ilat,1),T(ilon,ilat,1),qv(ilon,ilat,1),pres(ilon,ilat,1),dT, &
                                    pres(ilon,ilat,:),nlevs,zm,Tm,qm)

                if (rain_mod(ilon,ilat,1) .gt. 1.0e-03_rp .and. rain_GPM(ilon,ilat,1) .le. 1.0e-03_rp) &
                    qcl(ilon,ilat,:) = 0.2_rp * qm(:)
                if (rain_mod(ilon,ilat,1) .le. 1.0e-03_rp .and. rain_GPM(ilon,ilat,1) .gt. 1.0e-03_rp) &
                    qcl(ilon,ilat,:) = 1.8_rp * qm(:)
                if (abs(rain_mod(ilon,ilat,1)-rain_GPM(ilon,ilat,1)) .gt. 1.0e-04_rp .and. &
                      rain_GPM(ilon,ilat,1) .gt. 1.0e-03_rp .and. rain_mod(ilon,ilat,1) .gt. 1.0e-02_rp) then
                        alpha = ( rain_GPM(ilon,ilat,1) - rain_mod(ilon,ilat,1) ) / rain_mod(ilon,ilat,1)
                        qcl(ilon,ilat,:) = ( 1.0_rp + alpha ) * qcl(ilon,ilat,:)
                endif
                
                rain_mod(ilon,ilat,1) = 0.0_rp
                do ilev = 2, nlevs
                    if (omg(ilon,ilat,ilev) .le. 0.0_rp) then
                        rain_mod(ilon,ilat,1) = rain_mod(ilon,ilat,1) &
                        + ( omg(ilon,ilat,ilev)+omg(ilon,ilat,ilev-1) )/2.0_rp * (qm(ilev-1)-qm(ilev))
                    endif
                enddo
                rain_mod(ilon,ilat,1) = -0.5*60*60/(g*100) * rain_mod(ilon,ilat,1)
                ! rain_mod(ilon,ilat,1) = -3.0*60*60/(g*100) * rain_mod(ilon,ilat,1)
                do ilev = 2, nlevs-2
!                     if (omg(ilon,ilat,ilev) .le. 0.0_rp .and. qcl(ilon,ilat,ilev) .gt. qm(ilev)) then
!                         rain_mod(ilon,ilat,1) = rain_mod(ilon,ilat,1) + qcl(ilon,ilat,ilev) - qm(ilev)
!                     endif
                    if (qcl(ilon,ilat,ilev) .gt. qm(ilev)) &
                    rain_mod(ilon,ilat,1) = rain_mod(ilon,ilat,1) + qcl(ilon,ilat,ilev) - qm(ilev)
                enddo
                if (rain_GPM(ilon,ilat,1) .le. 1.0e-03_rp) &
                    rain_mod(ilon,ilat,1) = 0.0_rp
            enddo
            write(*,*) "Computing for ilat: ", ilat
        enddo
    end subroutine Disp_of_SuperSat
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Computes a Moist Adiabat from Data (z,T,q) at level p0 via the 
! gz + CpT + Lq Method assuming it to be invariant along the constructed
! sounding.
! Note : The Subroutine is written to start at p0 and work upward
! through decreasing Pressure. One may work down through increasing 
! Pressure by setting dp negative. Also note that if this is done, the 
! order of the "if (ET-EM)" should be reversed by making it "(EM-ET)".
!
! Usage : call Moist_Adiabat(z0,T0,q0,p0,dT,dp,nlevs,z,T,q)
!
! Input Data
! z0    : Height (m) of p0 Surface
! T0    : Temperature (K) at p0 Surface
! q0    : Specific Humidity (gm/gm) at p0
! p0    : Pressure (hPa) of Reference Surface
! dT    : Temperature Increment (K) for use in guessing Temperature at p
!        Surface
! dp    : Pressure change between desired Pressure Surfaces
! nlevs : Number of desired Pressure Surfaces (including p0)
!
! Output Data
! z     : Height (m) at each Pressure Surface
! T     : Temperature (C) at each Pressure Surface
! q     : Absolute Humidity at each Pressure Surface
!-----------------------------------------------------------------------

     subroutine Moist_Adiabat(z0,T0,qv0,p0,dT,p,nlevs,z,T,qm)
        use Constants
        implicit none

        integer, intent(in)             :: nlevs
        real(kind=rp)                   :: p0, z0, T0, qv0
        real(kind=rp)                   :: pBar, zBar, TBar, qmBar
        real(kind=rp), dimension(nlevs) :: z, T, qm, p
        real(kind=rp)                   :: dT, dz, dp, es, ET, EM
        real(kind=rp), parameter        :: HL = 2.49e6
        integer                         :: L, k

        ! Start from given level
        z(1) = z0; T(1) = T0; qm(1) = qv0
        EM   = g*z(1)+Cp*T(1)+HL*qm(1)     ! Compute Reference State Energy
        do L = 2, nlevs                    ! Begin Main Processing Loop
            dp = p(L-1)-p(L)
            do k = 1, 10000                ! Begin Temp, Humidity, Height Loop for this Surface
                T(L) = T(L-1)-(k*dT)       ! Guess Temperature
                ! Compute Saturation Vapour Pressure and Specific Humidity for this Temperature
                es   = Tetens(T(L))
                qm(L) = 0.622_rp * es/(p(L)-0.378_rp*es)
                ! Compute averaged Quantities for use in Height (Thickness) Calculations
                qmBar = (qm(L)+qm(L-1))/2.0_rp
                TBar = (T(L)+T(L-1))/2.0_rp
                pBar = (p(L)+p(L-1))/2.0_rp
                ! Compute Thickness and Height of Pressure Surface
                dz   = (R*TBar*(1.0_rp+0.61_rp*qmBar)*dp)/(g*pBar)
                z(L) = z(L-1)+dz
                ET   = g*z(L)+Cp*T(L)+HL*qm(L)                          ! Compute Static Energy
                ! If Static Energy is not the same as that at p0, guess a new Temperature
                ! If it is same, go to next Pressure Surface
                if ((ET-EM) .le. 0.01_rp) exit
            enddo
        enddo
    end subroutine Moist_Adiabat
    
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Moist Adiabat (z, T, qv) for a fixed interval of Pressure
!-----------------------------------------------------------------------
     
     subroutine Moist_Adiabat_fixed_dp(z0,T0,qv0,p0,dT,dp,nlevs,z,T,qv)
        use Constants
        implicit none

        integer, intent(in)             :: nlevs
        real(kind=rp)                   :: p0, z0, T0, qv0
        real(kind=rp)                   :: pBar, zBar, TBar, qvBar
        real(kind=rp), dimension(nlevs) :: z, T, qv
        real(kind=rp)                   :: p, dT, dz, dp, es, ET, EM
        real(kind=rp), parameter        :: HL = 2.49e6
        integer                         :: L, k

        ! Start from given level
        z(1) = z0; T(1) = T0; qv(1) = qv0
        EM   = g*z(1)+Cp*T(1)+HL*qv(1)      ! Compute Reference State Energy
        do L = 2, nlevs                    ! Begin Main Processing Loop
            p = p0-((L-1)*dp)              ! Increment Pressure
            do k = 1, 1000                 ! Begin Temp, Humidity, Height Loop for this Surface
                T(L) = T(L-1)-(k*dT)       ! Guess Temperature
                ! Compute Saturation Vapour Pressure and Specific Humidity for this Temperature
                es   = Tetens(T(L))
                qv(L) = 0.622_rp * es/(p-0.378_rp*es)
                ! Compute averaged Quantities for use in Height (Thickness) Calculations
                qvBar = (qv(L)+qv(L-1))/2.0_rp
                TBar = (T(L)+T(L-1))/2.0_rp
                pBar = (p+(p0-(L-2)*dp))/2.0_rp
                ! Compute Thickness and Height of Pressure Surface
                dz   = (R*TBar*(1.0_rp+0.61_rp*qvBar)*dp)/(g*pBar)
                z(L) = z(L-1)+dz
                ET   = g*z(L)+Cp*T(L)+HL*qv(L)                           ! Compute Static Energy
                ! If Static Energy is not the same as that at p0, guess a new Temperature
                ! If it is same, go to next Pressure Surface
                if ((ET-EM) .le. 0.01_rp) exit
            enddo
        enddo
    end subroutine Moist_Adiabat_fixed_dp

!-----------------------------------------------------------------------
! Temperature (T) -> Saturation Evaporation Ratio (es)
! Usage : var = Read_Var (real var, string VAR_NAME, string FILE_NAME, 
!                         real lons, real lats, real levs, int nlons, int nlats, int nlevs)
!-----------------------------------------------------------------------

    function Tetens(T) result(es)
        use Constants
        implicit none

        real(kind=rp) :: a, b                                           ! Constants in Tetens 
        real(kind=rp) :: T                                              ! Temperature                     (K)
        real(kind=rp) :: es                                             ! Saturation Vapour Pressure      (hPa)

        if ( T .gt. 263.0_rp ) then
             a = 17.26;  b = 35.86
        else
             a = 21.87;  b = 7.66
        endif
        es = 6.112_rp * exp ( a*(T-273.15_rp)/(T-b) )
    end function Tetens

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Read a variable from a netCDF file.
! Usage : var = Read_Var (real var, string VAR_NAME, string FILE_NAME, 
!                         real lons, real lats, real levs, int nlons, int nlats, int nlevs)
!-----------------------------------------------------------------------

    subroutine Read_Var (var, VAR_NAME, FILE_NAME, lons, lats, levs, nlons, nlats, nlevs)
        use Constants
        use netcdf
        implicit none

        character (len = *) :: FILE_NAME
        integer             :: ncid

        ! Reading 4D data, lat-lon grid, with 1 timesteps of data.
        integer, parameter :: ndims = 4, nrecs = 1
        integer            :: nlons, nlats, nlevs
        character (len = *), parameter :: LEV_NAME = "lev"
        character (len = *), parameter :: LAT_NAME = "lat"
        character (len = *), parameter :: LON_NAME = "lon"
        character (len = *), parameter :: REC_NAME = "time"
        integer :: lev_dimid, lon_dimid, lat_dimid, rec_dimid           ! File Dimension

        ! The start and count arrays will tell the netCDF library where to
        ! read our data.
        integer :: start(ndims), count(ndims)

        ! In addition to the latitude and longitude dimensions, we will also
        ! create latitude and longitude variables which will hold the actual
        ! latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real(kind=rp) :: lons(nlons), lats(nlats), levs(nlevs)
        integer       :: lon_varid, lat_varid, lev_varid

        ! "Variables."
        character (len = *) :: VAR_NAME
        integer             :: varid
        integer             :: dimids(ndims)

        ! We recommend that each variable carry a "units" attribute.
        character (len = *), parameter :: units     = "units"
        character (len = *), parameter :: LON_UNITS = "degrees_east"
        character (len = *), parameter :: LAT_UNITS = "degrees_north"
        character (len = *), parameter :: LEV_UNITS = ""

        ! Program variables to hold the data we will read in. We will only
        ! need enough space to hold one timestep of data; one record.
        real(kind=rp) :: var(nlons, nlats, nlevs)

        integer :: lev, lat, lon, rec, i                                ! Loop indices

        ! Open the file.
        call check (nf90_open(FILE_NAME, nf90_nowrite, ncid))

        ! Get the varids of the latitude and longitude coordinate variables.
        call check (nf90_inq_varid(ncid, LAT_NAME, lat_varid))
        call check (nf90_inq_varid(ncid, LON_NAME, lon_varid))
        if (nlevs .gt. 1) call check (nf90_inq_varid(ncid, LEV_NAME, lev_varid))

        ! Read the latitude and longitude data.
        call check (nf90_get_var(ncid, lon_varid, lons))
        call check (nf90_get_var(ncid, lat_varid, lats))
        if (nlevs .gt. 1) call check (nf90_get_var(ncid, lev_varid, levs))

        ! Get the varids of the netCDF variables.
        call check (nf90_inq_varid(ncid, VAR_NAME, varid))

        ! Read 1 record of nlevs*nlats*nlons values, starting at the beginning 
        ! of the record (the (1, 1, 1, rec) element in the netCDF file).
        count = (/ nlons, nlats, nlevs, 1 /)
        start = (/ 1, 1, 1, 1 /)

        ! Read the data from the file
        do rec = 1, nrecs
            start(4) = rec
            call check (nf90_get_var(ncid, varid, var, start = start, &
                                     count = count))
        end do

        ! Close the file. This frees up any internal netCDF resources
        ! associated with the file.
        call check (nf90_close(ncid))
    end subroutine Read_Var
    
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Write a variable to a netCDF file.
! Usage : Write_Var (real var, string VAR_NAME, string VAR_UNITS, string FILE_NAME, 
!                    real lons, real lats, real levs, int nlons, int nlats, int nlevs)
!-----------------------------------------------------------------------

    subroutine Write_Var (var, VAR_NAME, VAR_UNITS, FILE_NAME, lons, lats, levs, nlons, nlats, nlevs)
        use Constants
        use netcdf
        implicit none

        character (len = *) :: FILE_NAME
        integer             :: ncid

        ! Reading 4D data, lat-lon grid, with 1 timesteps of data.
        integer, parameter  :: ndims = 4, nrecs = 1

        ! We create latitude and longitude variables which will hold the 
        ! actual latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real(kind=rp), dimension(1:)   :: lons, lats, levs
        integer                        :: lon_varid, lat_varid, lev_varid
        integer                        :: nlons, nlats, nlevs
        character (len = *), parameter :: LEV_NAME = "level"
        character (len = *), parameter :: LAT_NAME = "lat"
        character (len = *), parameter :: LON_NAME = "lon"
        character (len = *), parameter :: REC_NAME = "time"
        integer :: lev_dimid, lon_dimid, lat_dimid, rec_dimid           ! File Dimension

        ! The start and count arrays will tell the netCDF library where to
        ! read our data.
        integer :: start(ndims), count(ndims)

        ! "Variables."
        character (len = *) :: VAR_NAME
        integer             :: varid
        integer             :: dimids(ndims)

        ! We recommend that each variable carry a "units" attribute.
        character (len = *), parameter :: UNITS     = "units"
        character (len = *)            :: VAR_UNITS
        character (len = *), parameter :: LAT_UNITS = "degrees_north"
        character (len = *), parameter :: LON_UNITS = "degrees_east"
        character (len = *), parameter :: LEV_UNITS = ""

        ! Program variables to hold the data we will read in. We will only
        ! need enough space to hold one timestep of data; one record.
        real(kind=rp)                  :: var(1:,1:,1:)

        integer :: lev, lat, lon, rec, i                                ! Loop indices

        ! Create the file. The nf90_clobber parameter tells netCDF to
        ! overwrite this file, if it already exists.
        call check (nf90_create(FILE_NAME, nf90_clobber, ncid))

        call check (nf90_def_dim(ncid, LEV_NAME, nlevs, lev_dimid))
        call check (nf90_def_dim(ncid, LAT_NAME, nlats, lat_dimid))
        call check (nf90_def_dim(ncid, LON_NAME, nlons, lon_dimid))
        call check (nf90_def_dim(ncid, REC_NAME, nf90_unlimited, rec_dimid))
    
        ! Define the coordinate variables. We will only define coordinate
        ! variables for lat and lon.  Ordinarily we would need to provide
        ! an array of dimension IDs for each variable's dimensions, but
        ! since coordinate variables only have one dimension, we can
        ! simply provide the address of that dimension ID (lat_dimid) and
        ! similarly for (lon_dimid).
        call check (nf90_def_var(ncid, LEV_NAME, NF90_REAL, lev_dimid, lev_varid))
        call check (nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid))
        call check (nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid))

        ! Assign units attributes to coordinate variables.
        call check (nf90_put_att(ncid, lev_varid, UNITS, LEV_UNITS))
        call check (nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS))
        call check (nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS))

        ! The dimids array is used to pass the dimids of the dimensions of
        ! the netCDF variables. Both of the netCDF variables we are creating
        ! share the same four dimensions. In Fortran, the unlimited
        ! dimension must come last on the list of dimids.
        dimids = (/ lon_dimid, lat_dimid, lev_dimid, rec_dimid /)

        ! Define the netCDF variables for the data.
        call check (nf90_def_var(ncid, VAR_NAME, NF90_REAL, dimids, varid))

        ! Assign units attributes to the netCDF variables.
        call check (nf90_put_att(ncid, varid, UNITS, VAR_UNITS))

        ! End define mode. This tells netCDF we are done defining metadata.
        call check (nf90_enddef(ncid))

        ! Read the latitude and longitude data.
        call check (nf90_put_var(ncid, lev_varid, levs))
        call check (nf90_put_var(ncid, lat_varid, lats))
        call check (nf90_put_var(ncid, lon_varid, lons))

        ! Read 1 record of nlevs*nlats*nlons values, starting at the beginning 
        ! of the record (the (1, 1, 1, rec) element in the netCDF file).
        count = (/ nlons, nlats, nlevs, 1 /)
        start = (/ 1, 1, 1, 1 /)

        ! Read the surface pressure and temperature data from the file, one
        ! record at a time.
        do rec = 1, nrecs
            start(4) = rec
            call check( nf90_put_var(ncid, varid, var, start = start, &
                                     count = count) )
        end do

        ! Close the file. This frees up any internal netCDF resources
        ! associated with the file.
        call check (nf90_close(ncid))
    end subroutine Write_Var

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Always check the return code of every netCDF function call. In
! this example program, wrapping netCDF calls with "call check()"
! makes sure that any return which is not equal to nf90_noerr (0)
! will print a netCDF error message and exit.
!-----------------------------------------------------------------------
    subroutine check(status)
        use netcdf
        integer, intent(in) :: status
    
        if (status /= nf90_noerr) then 
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check
!-----------------------------------------------------------------------

end module Phys_Init_Class
