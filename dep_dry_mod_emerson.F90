!>\file dep_dry_mod.F90
!! This file is for the dry depostion driver.
!-------------REVISION HISTORY---------------!
! XX/XX/XXXX : original implementation (Ravan Ahmadov)
! 08/17/2023 : modified to follow Emerson et al., (2020) (Jordan Schnell)
! 08/17/2023 : gravitational settling folowing the coarse pm settling driver (Jordan Schnell)

module dep_dry_mod_emerson

  use mpas_kind_types
  use dep_data_mod     
  use mpas_smoke_config, only : n_dbg_lines
  use mpas_smoke_init
  use mpas_timer, only : mpas_timer_start, mpas_timer_stop

  implicit none

  private

  public :: dry_dep_driver_emerson, particle_settling_wrapper

contains
    subroutine dry_dep_driver_emerson(rmol,ustar,znt,num_chem,ddvel,         &
               chem,delz,snowh,t_phy,p_phy,rho_phy,ivgtyp,gravity,dt,        &
               drydep_flux,tend_chem_settle,dbg_opt,settling_opt,vg,         &
               ids,ide, jds,jde, kds,kde,                                    &
               ims,ime, jms,jme, kms,kme,                                    &
               its,ite, jts,jte, kts,kte, curr_secs                          )
!
! compute dry deposition velocity for aerosol particles
! Based on Emerson et al. (2020), PNAS,
! www.pnas.org/cgi/doi/10.1073/pnas.2014761117
! Code adapted from Hee-Ryu and Min, (2022):
! https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2021MS002792
!----------------------------------------------------------------------
  IMPLICIT NONE
 
       INTEGER,  INTENT(IN   ) ::  num_chem,                               &
                                   ids,ide, jds,jde, kds,kde,              &
                                   ims,ime, jms,jme, kms,kme,              &
                                   its,ite, jts,jte, kts,kte

       REAL(RKIND) :: curr_secs

       REAL(RKIND), DIMENSION( ims:ime , jms:jme )        ,            &
           INTENT(IN) :: ustar, rmol, znt, snowh
       REAL(RKIND), DIMENSION( ims:ime , kms:kme , jms:jme ),          &
           INTENT(IN   ) :: t_phy, rho_phy, p_phy, delz                    
       INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN) ::  ivgtyp          
       REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ),  &
                                             INTENT(IN ) :: chem
       REAL(RKIND), INTENT(IN) :: gravity,dt
       LOGICAL, INTENT(IN) :: dbg_opt
       INTEGER, INTENT(IN) :: settling_opt
 !
 ! Output arrays 
       REAL(RKIND), DIMENSION( ims:ime, jms:jme, 1:num_chem ), INTENT(INOUT)       :: drydep_flux
       REAL(RKIND), DIMENSION( ims:ime, jms:jme, 1:num_chem ), INTENT(OUT)         :: ddvel
       REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem), INTENT(OUT) :: tend_chem_settle, vg
 ! Local
       real(RKIND), DIMENSION( its:ite, kts:kte, jts:jte)     :: airkinvisc        ! Air kinetic viscosity [cm2/s]
       real(RKIND), DIMENSION( its:ite, kts:kte, jts:jte)     :: freepath          ! Air molecular freepath [cm]
       real(RKIND), DIMENSION( its:ite, jts:jte)              :: aer_res           ! Aerodynmaic resistance
       real(RKIND), DIMENSION( its:ite, jts:jte)              :: A, eps0           ! land surface params [-]
       real(RKIND) :: Cc                      ! Cunningham/slip correction factor [-]
       real(RKIND) :: DDp, Eb                 ! Brownian diffusion []
       real(RKIND) :: Eim                     ! Impaction []
       real(RKIND) :: Ein                     ! Interception []
       real(RKIND) :: Sc                      ! Schmit number []
       real(RKIND) :: St                      ! Stokes number []
       real(RKIND) :: dp                      ! aerosol diameter [cm]
       real(RKIND) :: aerodens                ! aerosol density [g/cm3] 
       real(RKIND) :: Rs                      ! Surface resistance
       real(RKIND) :: growth_fac,vsettl,dtmax,conver,converi,dzmin
       real(RKIND) :: rmol_local
       integer :: i, j, k, ntdt, nv
       integer :: icall=0
!> -- Gas constant
       real(RKIND), parameter :: RSI = 8.314510_RKIND
       real(RKIND) :: four_ninths, g100

       logical, parameter :: do_timing = .true.

       integer, parameter :: max_iter_settle = 5

       real(RKIND) :: ust

       tend_chem_settle(:,:,:,:) = 0._RKIND
       ddvel(:,:,:)              = 0._RKIND
       vg(:,:,:,:)               = 0._RKIND
 !      growth_fac        = 1._RKIND
       conver            = 1.e-9_RKIND
       converi           = 1.e9_RKIND
 !      four_ninths       = 4._RKIND / 9._RKIND
       g100             = gravity * 1.e2_RKIND
       ntdt=INT(dt)
 
       if (mod(int(curr_secs),1800) .eq. 0) then
           icall = 0
       endif
       
! Set up 2D vars
       ! MODIS type lu, large roughness lengths (e.g., urban or forest)
       ! -----------------------------------------------------------------------
       ! *** TO DO -- set A and eps0 for all land surface types *** !!!
       ! -----------------------------------------------------------------------
       ! Set if snow greater than 1 cm
       do j = jts, jte
       do i = its, ite
          if ( ivgtyp(i,j) .eq. 13 .or. ivgtyp(i,j) .le. 5 ) then ! Forest
             A(i,j) = A_for
             eps0(i,j) = eps0_for
          else if ( ivgtyp(i,j) .eq. 17 ) then ! water
             A(i,j) = A_wat
             eps0(i,j) = eps0_wat
          else ! otherwise
             A(i,j) = A_grs
             eps0(i,j) = eps0_grs
          end if
          rmol_local = rmol(i,j)
          call depvel( rmol_local, dep_ref_hgt, znt(i,j), ustar(i,j), aer_res(i,j) )
       enddo
       enddo
!
! Set up 3D vars
!       if  (do_timing) call mpas_timer_start('dep_prep')
       do j = jts, jte
       do k = kts, kte
       do i = its, ite
          airkinvisc(i,k,j)  = ( 1.8325e-4_RKIND * ( 416.16_RKIND / ( t_phy(i,k,j) + 120._RKIND) ) *   &
                       ( ( t_phy(i,k,j) / 296.16_RKIND )**1.5_RKIND ) ) / (rho_phy(i,k,j) * 1.e-3_RKIND) ! Convert density to mol/cm^3
          ! Air molecular freepath (cm)  
          freepath(i,k,j)    = 7.39758e-4_RKIND * airkinvisc(i,k,j) / sqrt( t_phy(i,k,j) )
!          delz_flip(i,k,j)   = delz(i,kte-kts+k,j)
       enddo
       enddo
       enddo
!      
! 3D + chem --> vg
!       if  (do_timing) call mpas_timer_start('vg_and_ddvel_calc')
       do nv = 1, num_chem
          if (aero_diam(nv) .lt. 0) cycle  ! At some point we'll do something different for gasses
          ! Convert diameter to cm and aerodens to g/cm3
          dp       = aero_diam(nv) * 100._RKIND
          aerodens = aero_dens(nv) * 1.e-3_RKIND
          do j = jts, jte
          do k = kts, kte
          do i = its, ite
             ! Cunningham correction factor
             Cc = 1._RKIND + 2._RKIND * freepath(i,k,j) / dp * ( 1.257_RKIND + 0.4_RKIND*exp( -0.55_RKIND * dp / freepath(i,k,j) ) )
             ! Gravitational Settling
             vg(i,k,j,nv) = aerodens * dp * dp * g100 * Cc / &       ! Convert gravity to cm/s^2
                    ( 18._RKIND * airkinvisc(i,k,j) * (rho_phy(i,k,j)*1.e-3_RKIND) ) ! Convert density to mol/cm^
          enddo
          enddo
          enddo
       enddo

! 2D + chem (surface dep)
       k=kts 
       do nv = 1, num_chem
          if (aero_diam(nv) .lt. 0) cycle  ! At some point we'll do something different for gasses
          ! Convert diameter to cm and aerodens to g/cm3
          dp       = aero_diam(nv) * 100._RKIND
          aerodens = aero_dens(nv) * 1.e-3_RKIND
          do j = jts, jte
          do i = its, ite
                ust        = ustar(i,j)*1.e2_RKIND
                ! Cunningham correction factor
                Cc = 1._RKIND + 2._RKIND * freepath(i,k,j) / dp * ( 1.257_RKIND + 0.4_RKIND*exp( -0.55_RKIND * dp / freepath(i,k,j) ) )
                ! Brownian Diffusion
                DDp = ( boltzmann * t_phy(i,k,j) ) * Cc / (3._RKIND * pi * airkinvisc(i,k,j) * (rho_phy(i,k,j)*1.e-3_RKIND)  * dp) ! Convert density to mol/cm^3
                ! Schmit number
                Sc = airkinvisc(i,k,j) / DDp
                ! Brownian Diffusion
                Eb = Cb * Sc**(-0.666666667_RKIND)
                ! Stokes number
                St = ust*ust * vg(i,k,j,nv) / airkinvisc(i,k,j) / g100 ! Convert ustar to cm/s, gravity to cm/s^2
                ! Impaction 
                Eim = Cim * ( St / ( alpha + St ) )**1.7_RKIND
                ! Interception
                Ein = Cin * ( dp / A(i,j) )**vv
                ! Surface resistance
                Rs = 1._RKIND / ( ust * ( Eb + Eim + Ein) * eps0(i,j) ) ! Convert ustar to cm/s
                ! Compute final ddvel = aer_res + RS, set max at max_dep_vel in dep_data_mod.F[ m/s]
                ! The /100. term converts from cm/s to m/s, required for MYNN.
                if ( settling_opt .gt. 0 ) then
                      ddvel(i,j,nv) = max(min( ( vg(i,k,j,nv) + 1._RKIND / (aer_res(i,j)+Rs) )*1.e-2_RKIND, max_dep_vel),0._RKIND)
                else
                      ddvel(i,j,nv) = max(min( ( 1._RKIND / (aer_res(i,j)+Rs) )*1.e-2_RKIND, max_dep_vel),0._RKIND)
                endif
                if ( dbg_opt .and. (icall .le. n_dbg_lines) ) then
                   icall = icall + 1
                endif
                drydep_flux(i,j,nv) = drydep_flux(i,j,nv) + chem(i,k,j,nv)*rho_phy(i,k,j)*ddvel(i,j,nv)*dt*10._RKIND
          enddo
          enddo
       enddo
        
end subroutine dry_dep_driver_emerson
!
!--------------------------------------------------------------------------------
!
subroutine depvel( rmol, zr, z0, ustar, aer_res )
!--------------------------------------------------
!     THIS FUNCTION HAS BEEN DESIGNED TO EVALUATE AN UPPER LIMIT
!     FOR THE POLLUTANT DEPOSITION VELOCITY AS A FUNCTION OF THE
!     SURFACE ROUGHNESS AND METEOROLOGICAL CONDITIONS.
!     PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
!         Modified by Darrell A. Winner  (Feb. 1991)
!                  by Winfried Seidl     (Aug. 1997)
!.....PROGRAM VARIABLES...
!     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
!     ZR       - REFERENCE HEIGHT
!     Z0       - SURFACE ROUGHNESS HEIGHT
!     USTAR    - FRICTION VELOCITY U*
!     AER_RES  - AERODYNAMIC RESISTANCE
!.....REFERENCES...
!     MCRAE, G.J. ET AL. (1983) MATHEMATICAL MODELING OF PHOTOCHEMICAL
!       AIR POLLUTION, ENVIRONMENTAL QUALITY LABORATORY REPORT 18,
!       CALIFORNIA INSTITUTE OF TECHNOLOGY, PASADENA, CALIFORNIA.
!.....RESTRICTIONS...
!     1. THE MODEL EDDY DIFFUSIVITIES ARE BASED ON MONIN-OBUKHOV
!        SIMILARITY THEORY AND SO ARE ONLY APPLICABLE IN THE
!        SURFACE LAYER, A HEIGHT OF O(30M).
!     2. ALL INPUT UNITS MUST BE CONSISTENT
!     3. THE PHI FUNCTIONS USED TO CALCULATE THE FRICTION
!        VELOCITY U* AND THE POLLUTANT INTEGRALS ARE BASED
!        ON THE WORK OF BUSINGER ET AL.(1971).
!     4. THE MOMENTUM AND POLLUTANT DIFFUSIVITIES ARE NOT
!        THE SAME FOR THE CASES L<0 AND L>0.
!--------------------------------------------------
! .. Scalar Arguments ..
!--------------------------------------------------
        REAL(RKIND), intent(in)    :: ustar, z0, zr
        REAL(RKIND), intent(out)   :: aer_res
        REAL(RKIND), intent(inout) :: rmol
!--------------------------------------------------
! .. Local Scalars ..
!--------------------------------------------------
        INTEGER :: l
        REAL(RKIND)    :: ao, ar, polint, vk
!--------------------------------------------------
! .. Intrinsic Functions ..
!--------------------------------------------------
        INTRINSIC alog
!--------------------------------------------------
!     Set the von Karman constant
!--------------------------------------------------
        vk = karman
!--------------------------------------------------
!     DETERMINE THE STABILITY BASED ON THE CONDITIONS
!             1/L < 0 UNSTABLE
!             1/L = 0 NEUTRAL
!             1/L > 0 STABLE
!--------------------------------------------------
        if(abs(rmol) < 1.E-6 ) rmol = 0._RKIND
        IF ( rmol < 0._RKIND ) THEN
          ar = ((1._RKIND-9._RKIND*zr*rmol)**(0.25)+0.001_RKIND)**2.
          ao = ((1._RKIND-9._RKIND*z0*rmol)**(0.25)+0.001_RKIND)**2.
          polint = 0.74_RKIND*(alog((ar-1._RKIND)/(ar+1._RKIND))-alog((ao-1._RKIND)/(ao+1._RKIND)))
        ELSEIF ( rmol==0._RKIND ) THEN
          polint = 0.74_RKIND*alog(zr/z0)
        ELSE
          polint = 0.74_RKIND*alog(zr/z0) + 4.7_RKIND*rmol*(zr-z0)
        END IF
        !vgpart = ustar*vk/polint
        aer_res = polint/(karman*max(ustar,1.0e-4_RKIND)) * 1.e-2_RKIND ! convert to s/cm
end subroutine depvel
!
!--------------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------------
!
subroutine particle_settling_wrapper(tend_chem_settle,chem,rho_phy,delz_flip,vg,  &
                             dt, kts,kte,its,ite,jts,jte,num_chem,ims,ime, jms,jme, kms,kme       )
     IMPLICIT NONE
     
     INTEGER, INTENT(IN ) :: kts, kte,its,ite,jts,jte,num_chem,ims,ime, jms,jme, kms,kme 
     REAL(RKIND), INTENT(IN) :: dt
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT (IN)  :: rho_phy
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT (IN)  :: delz_flip
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(IN) :: chem
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(IN) :: vg
     REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT) :: tend_chem_settle
     
     REAL(RKIND) :: dt_settl, growth_fac, four_ninths, dzmin, vsettl, dtmax
     INTEGER     :: ntdt, ndt_settl
!
!--- Local------
     integer, parameter :: max_iter_settle = 10
     real(RKIND), parameter :: one_over_dyn_visc = 1.e5_RKIND ! 5.5248E5_RKIND
     INTEGER :: k,n,l2,i,j,nv
     REAL(RKIND) :: temp_tc, transfer_to_below_level, vd_wrk1
  
     growth_fac        = 1._RKIND
     four_ninths       = 4._RKIND / 9._RKIND
     ntdt = INT(dt)

     do nv = 1,num_chem
     if (aero_diam(nv) .lt. 0) cycle  ! At some point we'll do something different for gasses
     ! -- NOTE, diameters and densities are NOT converted to cm and g/cm3 like in Emerson
     vsettl = four_ninths * gravity * aero_dens(nv) * ( growth_fac * ( 0.5_RKIND * aero_diam(nv) ))**2.0_RKIND * one_over_dyn_visc

     do j = jts,jte
     do i = its,ite
     dzmin = delz_flip(i,kts,j)
     do k = kts,kte
        dzmin = min(dzmin,delz_flip(i,k,j))
     enddo
     ! -- Get necessary info for settling
     ! -- Determine the maximum time-step satisying the CFL condition:
     ! -- dt_settl calculations (from original coarsepm_settling)
     ! 1.5E-5 = dyn_visc --> dust_data_mod.F90

     dtmax = dzmin / vsettl
 
     ndt_settl = MAX( 1, INT( ntdt /dtmax) )
     ! Limit maximum number of iterations
     IF (ndt_settl > max_iter_settle) ndt_settl = max_iter_settle
     dt_settl = REAL(ntdt,kind=RKIND) /REAL(ndt_settl,kind=RKIND)

     do n = 1,ndt_settl
     transfer_to_below_level = 0._RKIND

     do k = kte,kts,-1

        l2 = kte - k + 1
     
        temp_tc = chem(i,k,j,nv)
     
        vd_wrk1 = dt_settl * vg(i,k,j,nv)*1.e-2_RKIND / delz_flip(i,l2,j)               ! convert vg to m/s
     
        temp_tc = temp_tc * (1._RKIND - vd_wrk1) + transfer_to_below_level
        if ( k .gt. kts ) then
        transfer_to_below_level =(temp_tc*vd_wrk1)*((delz_flip(i,l2,j) &
                        *rho_phy(i,k,j))/(delz_flip(i,l2+1,j)*rho_phy(i,k-1,j)))          ! [ug/kg]
        endif
        tend_chem_settle(i,k,j,nv) = tend_chem_settle(i,k,j,nv) + (temp_tc - chem(i,k,j,nv))

     enddo ! k
     enddo ! n
     enddo ! i
     enddo ! j
     enddo ! nv



end subroutine particle_settling_wrapper
end module dep_dry_mod_emerson
