!-----------------------------------------------------------------------
! Program : Physical Initialization
! Author  : Bappaditya Nag
! Date    : 08/18/2017
!
! Test Program for Physical Initialization
!-----------------------------------------------------------------------

program Phys_Init_Main
    use Constants
    use netcdf
    ! use mpi
    use Phys_Init_Class
    implicit none
    
!-----------------------------------------------------------------------
!     integer, parameter                          :: nlevs = 29, nlats = 449, nlons = 406
!     integer, parameter                          :: nlevs = 29, nlats = 108, nlons = 105
    integer, parameter                          :: nlevs = 29, nlats = 224, nlons = 192
!     integer, parameter                          :: nlevs = 29, nlats = 108, nlons = 209
!     integer, parameter                          :: nlevs = 29, nlats = 454, nlons = 404
    real(kind=rp)                               :: lons(nlons), lats(nlats), levs(nlevs)
    real(kind=rp), dimension(nlons,nlats,1)     :: rain_Obvs, rain_WRF, rain_mod
    real(kind=rp), dimension(nlons,nlats,nlevs) :: qv, qcl, pres, tc, z, omg
    real(kind=rp), dimension(nlevs)             :: znu, p
    integer                                     :: ilat, ilon, ilev
    integer                                     :: iter
    real(kind=rp), dimension(nlons,nlats)       :: error
    real(kind=rp)                               :: maxerr, dT
    character(100)                              :: filename, qfilename, infile2D, infile3D
    character (len=5)                           :: charI
    integer                                     :: idx_lat, idx_lon
    real(kind=rp), dimension(nlevs)             :: zm, Tm, qm

    real(kind=rp)                               :: start, finish
    integer(kind=4)                             :: procid, ierr, numprocs
!-----------------------------------------------------------------------

    infile2D = "wrfout2D_2017-08-28-18:00:00.nc"
    infile3D = "wrfout3D_2017-08-28-18:00:00.nc"
    write(*,*) "Reading netCDF ..."
    call Read_Var (rain_Obvs, "RAIN"  , "./pre_Phys_Init/GPM_2017-08-28-00:0:00.nc" , lons, lats, (/1.0_rp/), nlons, nlats, 1)
!     call Read_Var (rain_Obvs, "RAIN"  , "./pre_Phys_Init/TRMM_2013-06-16-12:00:00.nc" , lons, lats, (/1.0_rp/), nlons, nlats, 1)
    call Read_Var (rain_WRF , "RAIN"  , "./pre_Phys_Init/"//infile2D                  , lons, lats, (/1.0_rp/), nlons, nlats, 1)
    call Read_Var (znu      , "ZNU"   , "./pre_Phys_Init/"//infile3D                  , lons, lats, levs, 1,     1,     nlevs)
    call Read_Var (z        , "Z"     , "./pre_Phys_Init/"//infile3D                  , lons, lats, levs, nlons, nlats, nlevs)
    call Read_Var (qv       , "QVAPOR", "./pre_Phys_Init/"//infile3D                  , lons, lats, levs, nlons, nlats, nlevs)
    call Read_Var (qcl      , "QCLOUD", "./pre_Phys_Init/"//infile3D                  , lons, lats, levs, nlons, nlats, nlevs)
    call Read_Var (pres     , "PRES"  , "./pre_Phys_Init/"//infile3D                  , lons, lats, levs, nlons, nlats, nlevs)
    call Read_Var (tc       , "TC"    , "./pre_Phys_Init/"//infile3D                  , lons, lats, levs, nlons, nlats, nlevs)
    call Read_Var (omg      , "OMG"   , "./pre_Phys_Init/"//infile3D                  , lons, lats, levs, nlons, nlats, nlevs)

    open (unit=10, file='timing.txt')
    open (unit=20, file='error.txt')
    
    tc = tc + 273.15_rp
    pres = pres/100.0_rp
    ! qcl = qcl*1000_rp
    qcl = qcl*1000_rp + 0.1_rp
    
!-----------------------------------------------------------------------
! Single Processor

    maxerr = 20.0_rp
    rain_mod = rain_WRF
    iter = 0
    do while (iter .le. MaxIter .and. maxerr .gt. tol)
        call cpu_time(start)
        write(*,*) "Iterations : ", iter
        maxerr = 0.0_rp
        call Disp_of_SuperSat (rain_mod,rain_Obvs,qv,qcl,pres,tc,znu,z,omg)
        do ilat = 10, nlats-10
            do ilon = 10, nlons-10
                error(ilon,ilat) = abs(rain_mod(ilon,ilat,1) - rain_Obvs(ilon,ilat,1))
                maxerr = max(error(ilon,ilat), maxerr)
                if (maxerr .eq. error(ilon,ilat)) write(20,'(A40,4F20.10)') "Lat, Lon, Rain_Mod, rain_Obvs .......", &
                lats(ilat), lons(ilon), rain_mod(ilon,ilat,1), rain_Obvs(ilon,ilat,1)
            enddo
        enddo
        write(*,*) "MaxError : ", maxerr
        iter = iter + 1
        write(charI,"(I5)") iter
        filename = 'wrfmod2D_2017-08-28-18:00:00_'//trim(adjustl(charI))//'.nc'
        qfilename = 'wrfmod3D_2017-08-28-18:00:00_'//trim(adjustl(charI))//'.nc'
        call Write_Var (rain_mod, "RAIN_MOD","mm",filename, lons, lats, (/1.0_rp/), nlons, nlats, 1)
        ! call Write_Var (qcl/1000.0_rp, "QCLOUD_MOD","kg/kg", qfilename, lons, lats, levs, nlons, nlats, nlevs)
        call Write_Var (qcl/(1000.0_rp+qcl), "SPHUM_MOD"," ", qfilename, lons, lats, levs, nlons, nlats, nlevs)
        call cpu_time(finish)
        write(10,'(I3,F10.4)') iter, finish-start
        write(20,'(I3)') iter
    enddo

!!-----------------------------------------------------------------------
!! Parallel Processing
!!-----------------------------------------------------------------------

!        call MPI_Init ( ierr )
!        call MPI_Comm_size ( MPI_COMM_WORLD, numprocs, ierr )
!        call MPI_Comm_rank ( MPI_COMM_WORLD, procid, ierr )
!        ! MPI_Bcast(ilat, COUNT, MPI_real, ROOT, COMM, IERROR)

!       !$OMP PARALLEL &
!       !$OMP PRIVATE ( ilon, ilat, ilev )
!       !$OMP DO


!    maxerr = 20.0_rp
!    rain_mod = rain_WRF
!    iter = 0
!    do while (iter .le. MaxIter .and. maxerr .gt. tol)
!        call cpu_time(start)
!        write(*,*) "Iterations : ", iter
!        maxerr = 0.0_rp
!        call Disp_of_SuperSat (rain_mod,rain_Obvs,qv,qcl,pres,tc,znu,z,omg)

!                        if (ilat .eq. 1 .and. ilon .eq. 3) &
!                write(*,*) "Computing for procid, ilat, ilon: ", procid, ilat, ilon

!        do ilat = 10, nlats-10
!            do ilon = 10, nlons-10
!                error(ilon,ilat) = abs(rain_mod(ilon,ilat,1) - rain_Obvs(ilon,ilat,1))
!                maxerr = max(error(ilon,ilat), maxerr)
!                if (maxerr .eq. error(ilon,ilat)) write(20,'(A40,4F20.10)') "Lat, Lon, Rain_Mod, rain_Obvs .......", &
!                lats(ilat), lons(ilon), rain_mod(ilon,ilat,1), rain_Obvs(ilon,ilat,1)
!            enddo
!        enddo
!        write(*,*) "MaxError : ", maxerr
!        iter = iter + 1
!        write(charI,"(I5)") iter
!        filename = 'wrfmod2D_2017-08-28-00:30:00_'//trim(adjustl(charI))//'.nc'
!        call Write_Var (rain_mod, "RAIN_MOD","mm",filename, lons, lats, levs, nlons, nlats, 1)
!        call cpu_time(finish)
!        write(10,'(I3,F10.4)') iter, finish-start
!        write(20,'(I3)') iter
!    enddo

!        !$omp end do
!        !$omp end parallel

!-----------------------------------------------------------------------

! Test for Computation of Moist Adiabat
!
!    dT = 0.01
!    do ilat = 1, nlats
!        do ilon = 1, nlons
!            if (rain_WRF(ilon,ilat,1) .gt. 0.1_rp .and. rain_Obvs(ilon,ilat,1) .gt. 0.1_rp) then
!                if (abs(rain_WRF(ilon,ilat,1)-rain_Obvs(ilon,ilat,1)) .gt. 1.0_rp) then
!                    idx_lat = ilat
!                    idx_lon = ilon
!                endif
!            endif
!        enddo
!    enddo
!    call Moist_Adiabat (z(idx_lon,idx_lat,1),tc(idx_lon,idx_lat,1),qv(idx_lon,idx_lat,1), &
!                        pres(idx_lon,idx_lat,1),dT,pres(idx_lon,idx_lat,:),nlevs,zm,Tm,qm)
!    ! call Moist_Adiabat_fixed_dp(z(idx_lon,idx_lat,1),tc(idx_lon,idx_lat,1),qv(idx_lon,idx_lat,1), &
!    !                     pres(idx_lon,idx_lat,1),dT,20.0_rp,nlevs,zm,Tm,qm)
!    open (unit=10, file='moist_adiabat.txt')
!    do ilev = 1, nlevs
!        write(10,'(5F20.10)') zm(ilev), Tm(ilev), qcl(idx_lon,idx_lat,ilev), qm(ilev), pres(idx_lon,idx_lat,ilev)
!    enddo

!-----------------------------------------------------------------------

end program Phys_Init_Main
