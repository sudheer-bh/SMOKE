!>\file  pollen_mod.F90
!! This file contains the MPAS pollen emission scheme.

module pollen_mod
!
!  This module developed by Jordan Schnell (NOAA GSL)
!  Allison Steiner, Tamanna Subba, and Yingxiao Zhang
!  For serious questions contact jordan.schnell@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init
  use dep_data_mod 

  implicit none

  private

  public :: pollen_driver

contains

  subroutine pollen_driver       (                               &
       num_chem, chem,                                           &
       dt, u10, v10, rho, dz8w, t, z_at_w, ktop2d,               &
       xland, rainc, rainnc, relhum,                             &
       swdown, total_flashrate, cldfrac,                         &
       num_pols_per_polp, pollen_emis_scale_factor,              &
       tree_pollen_emis_scale_factor,                            &
       grass_pollen_emis_scale_factor,                           &
       weed_pollen_emis_scale_factor,                            &
       e_bio_in, e_bio_out, kbio,                                &
       num_e_bio_in, num_e_bio_out,                              &
       index_e_bio_in_polp_tree,                                 &
       index_e_bio_in_polp_grass,                                &
       index_e_bio_in_polp_weed,                                 &
       index_e_bio_out_polp_tree,                                &
       index_e_bio_out_polp_grass,                               &
       index_e_bio_out_polp_weed,                                &
       ids,ide, jds,jde, kds,kde,                                &
       ims,ime, jms,jme, kms,kme,                                &
       its,ite, jts,jte, kts,kte                                 )
  
    IMPLICIT NONE

    INTEGER,      INTENT(IN   ) ::                             &
         ids,ide, jds,jde, kds,kde,                            &
         ims,ime, jms,jme, kms,kme,                            &
         its,ite, jts,jte, kts,kte   ! Dimensions
    INTEGER, INTENT(IN) :: num_chem  ! Size of chemsitry/scalar array
    INTEGER, INTENT(IN) :: num_e_bio_in, num_e_bio_out 
    INTEGER, INTENT(IN) :: kbio
    INTEGER, INTENT(IN) :: &
       index_e_bio_in_polp_tree, &
       index_e_bio_in_polp_grass, &
       index_e_bio_in_polp_weed, &
       index_e_bio_out_polp_tree, &
       index_e_bio_out_polp_grass, &
       index_e_bio_out_polp_weed
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ), &
           INTENT(INOUT) :: chem     ! Chemistry array
    REAL(RKIND), INTENT(IN) :: dt    ! Timestep (s)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme),                & 
           INTENT(IN) :: v10, u10    ! 10 m wind speed (m/s)
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ),       &
           INTENT(IN) :: rho         ! Dry air density (kg/m3) 
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ),       &
           INTENT(IN) :: dz8w        ! Layer Thickness (m)
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ),       &
           INTENT(IN) :: t           ! Temperature (K)
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ),       &
           INTENT(IN) :: z_at_w      ! Layer Height (m)
    INTEGER, DIMENSION( ims:ime , jms:jme ),                   &
           INTENT(IN) :: ktop2d      ! Index of cloud top
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ),               &
           INTENT(IN) :: xland       ! land/water
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ),                         &
           INTENT(IN) :: rainc       ! Convective precip this timestep (m)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ),                         &
           INTENT(IN) :: rainnc      ! Resolved precip this timestep  (m)
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ),       &
           INTENT(IN) :: relhum      ! Relative humidty (frac, 0-1)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ),               &
           INTENT(IN) :: swdown      ! Shortwave radiation at the surface (W/m2)
    REAL(RKIND), DIMENSION( ims:ime , jms:jme ),              &
           INTENT(IN) :: total_flashrate ! lightning flash rate (s-1)
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ),      &
           INTENT(IN) :: cldfrac
    REAL(RKIND), DIMENSION( ims:ime , 1:kbio, jms:jme, 1:num_e_bio_in ),              & 
           INTENT(IN) :: e_bio_in ! Emission potential for pollen species (grains/m2) 
    REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_e_bio_out ),               &
           INTENT(INOUT) :: e_bio_out                                  
! namelist
    REAL(RKIND), INTENT(IN) :: pollen_emis_scale_factor, num_pols_per_polp,    &
                               tree_pollen_emis_scale_factor,                  &
                               grass_pollen_emis_scale_factor,                 &
                               weed_pollen_emis_scale_factor

! local variables
    INTEGER     :: i, j, k, l, l_oc,ibin                                     
    REAL(RKIND) :: factaa, emis
    REAL(RKIND) :: fr, fh, fw, fa
    REAL(RKIND) :: ppemfact_numb_tree, ppemfact_numb_grass, ppemfact_numb_weed
    REAL(RKIND) :: ppemfact_mass_tree, ppemfact_mass_grass, ppemfact_mass_weed
    REAL(RKIND) :: rain, wind
    REAL(RKIND) :: flashrate_for_rupture
    INTEGER     :: kfreeze 
    REAL(RKIND) :: depth, ratio, cgfrac
    REAL(RKIND), DIMENSION( its:ite , jts:jte ) :: ic_flashrate, cg_flashrate
    INTEGER     :: ierr
    REAL(RKIND) :: diam_polp_tree, diam_polp_grass, diam_polp_weed
    REAL(RKIND) :: rho_polp_tree, rho_polp_grass, rho_polp_weed
    REAL(RKIND) :: day_to_sec, piover6
    REAL(RKIND) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9, fac10, fac11, fac12, fac13

! Define some constants (TODO, move to data mod or namelist?)
    REAL(RKIND), PARAMETER :: pr_low    = 0._RKIND     ! [mm/day] lowest rainfall
    REAL(RKIND), PARAMETER :: pr_high   = 0.5_RKIND    ! [mm/day] highest rainfall                                                     
    REAL(RKIND), PARAMETER :: rh_low    = 50._RKIND    ! [%] lowest relative humidity
    REAL(RKIND), PARAMETER :: rh_high   = 70._RKIND    ! [%] highest relative humidity
    REAL(RKIND), PARAMETER :: rh_rupt   = 80._RKIND    ! [%} Rupture humidity
    REAL(RKIND), PARAMETER :: converi   = 1.e-9_RKIND  ! Conversion factor
    REAL(RKIND), PARAMETER :: sph       = 3600._RKIND  ! Conversion from hours to seconds
    REAL(RKIND), PARAMETER :: hpd       = 24._RKIND    ! Conversion from days to hours
    REAL(RKIND), PARAMETER :: rho_polp  = 1200._RKIND  ! Pollen (coarse) density (kg/m3)
    REAL(RKIND), PARAMETER :: diam_polp = 25._RKIND    ! Diameter (coarse) pollen (um)
    REAL(RKIND), PARAMETER :: rho_pols  = 1425._RKIND  ! Pollen (coarse)density (kg/m3)
    REAL(RKIND), PARAMETER :: diam_pols = 0.15_RKIND   ! Diameter(coarse) pollen (um)
    REAL(RKIND), PARAMETER :: pols_to_polp_frac_rh = 0.7_RKIND  ! Fraction of polp to pols for humidity rupture
    REAL(RKIND), PARAMETER :: pols_to_polp_frac_lt = 0.1_RKIND  ! Fraction of polp to pols for lightning rupture
! Lightning Parameters
    REAL(RKIND), PARAMETER :: dH_min = 5.5_RKIND
    REAL(RKIND), PARAMETER :: dH_max = 14._RKIND
    REAL(RKIND), PARAMETER :: coef_A = 0.021_RKIND
    REAL(RKIND), PARAMETER :: coef_B = -0.648_RKIND
    REAL(RKIND), PARAMETER :: coef_C = 7.493_RKIND
    REAL(RKIND), PARAMETER :: coef_D = -36.54_RKIND
    REAL(RKIND), PARAMETER :: coef_E = 63.09_RKIND
    REAL(RKIND), PARAMETER :: cldtop_adjustment = 0._RKIND

    LOGICAL, PARAMETER :: do_pollen_lightning_rupture = .false.
    LOGICAL, PARAMETER :: do_pollen_rh_rupture = .true.

    day_to_sec = 1._RKIND / hpd / sph
    piover6    = pi / 6._RKIND

  ! Set the densities / diameters based on what is available
  if ( p_polp_tree > 0 ) then
     diam_polp_tree = 1.E6_RKIND * aero_diam(p_polp_tree)
     rho_polp_tree  = aero_dens(p_polp_tree)
  else
     diam_polp_tree = diam_polp
     rho_polp_tree  = rho_polp
  endif 
!
  if ( p_polp_grass > 0 ) then
     diam_polp_grass = 1.E6_RKIND * aero_diam(p_polp_grass)
     rho_polp_grass  = aero_dens(p_polp_grass)
  else
     diam_polp_grass = diam_polp
     rho_polp_grass  = rho_polp
  endif 
!
  if ( p_polp_weed > 0 ) then
     diam_polp_weed = 1.E6_RKIND * aero_diam(p_polp_weed)
     rho_polp_weed  = aero_dens(p_polp_weed)
  else
     diam_polp_weed = diam_polp
     rho_polp_weed  = rho_polp
  endif

  fac1 = piover6 * rho_polp_tree  * diam_polp_tree**3._RKIND * converi
  fac2 = piover6 * rho_polp_grass * diam_polp_grass**3._RKIND * converi
  fac3 = piover6 * rho_polp_weed  * diam_polp_weed**3._RKIND * converi

  ! Loop over the grid cells to calculate and add the emissions
    do j = jts, jte
    do i = its, ite
   
      ! No emissions over water or at night
    if ( xland(i,j) .lt. 1.5_RKIND .and. swdown(i,j) .gt. 0._RKIND ) then ! land

     ! Compute the rainfall at the current time step (and convert to mm)
       rain  = (rainc(i,j) + rainnc(i,j)) * 1.e3_RKIND
     ! Determine the rainfall factor
       if ( rain < 0.5_RKIND ) then
          fr = 1._RKIND
       else if ( rain > 1.0_RKIND) then
          fr = 0._RKIND
       else
          fr = ( 1._RKIND - rain ) / ( 1._RKIND - 0.5_RKIND )
       end if

     ! Compute the wind and wind factor
       wind = ( sqrt( u10(i,j)**2._RKIND + v10(i,j)**2._RKIND ) )
       fw = 1.5_RKIND * ( 1._RKIND -  exp( -1._RKIND * wind * 0.2_RKIND ))

     ! Compute the relative humidity factor
       if ( (relhum(i,kts,j) * 100._RKIND) < rh_low ) then
          fh = 1._RKIND
       else if ( (relhum(i,kts,j) * 100._RKIND) > rh_high ) then
          fh = 0.1_RKIND
       else
          fh = (rh_high - (relhum(i,kts,j) * 100._RKIND )) / (rh_high - rh_low)
       end if
 
     ! Combine the factors
     ! Emissions are described / day, convert to /sec 
       fa = fh * fw * fr * day_to_sec
  
     ! Compute the number emissions
       ppemfact_numb_tree  = e_bio_in(i,1,j,index_e_bio_in_polp_tree)  * fa !, 0._RKIND)
       ppemfact_numb_grass = e_bio_in(i,1,j,index_e_bio_in_polp_grass) * fa !, 0._RKIND)
       ppemfact_numb_weed  = e_bio_in(i,1,j,index_e_bio_in_polp_weed)  * fa !, 0._RKIND)
   
     ! Convert number emissions to mass emissions 
     ! fac[1-3] = # * pi/6 * dens * diam^3 * 1.e-9
       ppemfact_mass_tree  = ppemfact_numb_tree  * fac1
       ppemfact_mass_grass = ppemfact_numb_grass * fac2
       ppemfact_mass_weed  = ppemfact_numb_weed  * fac3 
       ! Calculate the conversion factor
       factaa = dt / ( dz8w(i,kts,j) * rho(i,kts,j) )

     ! Compute the mass emissions, update the diagnostic and chemistry arrays
       emis = tree_pollen_emis_scale_factor * pollen_emis_scale_factor * factaa * ppemfact_mass_tree
       e_bio_out(i,kts,j,index_e_bio_out_polp_tree)  = e_bio_out(i,kts,j,index_e_bio_out_polp_tree) + emis
       if (p_polp_tree .gt. 0)   chem(i,kts,j,p_polp_tree)  = chem(i,kts,j,p_polp_tree) + emis

       emis = grass_pollen_emis_scale_factor * pollen_emis_scale_factor * factaa * ppemfact_mass_grass
       e_bio_out(i,kts,j,index_e_bio_out_polp_grass)  = e_bio_out(i,kts,j,index_e_bio_out_polp_grass) + emis
       if (p_polp_grass .gt. 0)  chem(i,kts,j,p_polp_grass)  = chem(i,kts,j,p_polp_grass) + emis

       emis = weed_pollen_emis_scale_factor * pollen_emis_scale_factor * factaa * ppemfact_mass_weed
       e_bio_out(i,kts,j,index_e_bio_out_polp_weed)  = e_bio_out(i,kts,j,index_e_bio_out_polp_weed) + emis
       if (p_polp_weed .gt. 0)  chem(i,kts,j,p_polp_weed)  = chem(i,kts,j,p_polp_weed) + emis

       if (p_polp_all .gt. 0) then
          emis = pollen_emis_scale_factor * factaa * &
                (tree_pollen_emis_scale_factor * ppemfact_mass_tree + &
                 grass_pollen_emis_scale_factor * ppemfact_mass_grass + &
                 weed_pollen_emis_scale_factor * ppemfact_mass_weed) 
          chem(i,kts,j,p_polp_all)   = chem(i,kts,j,p_polp_all)  + emis

     endif ! if land

   enddo
   enddo

  fac4 = pols_to_polp_frac_rh * rho_pols/rho_polp       * &
         ((diam_pols)**3._RKIND)/((diam_polp)**3._RKIND)       * num_pols_per_polp
  fac5 = pols_to_polp_frac_rh * rho_pols/rho_polp_tree  * &
         ((diam_pols)**3._RKIND)/((diam_polp_tree)**3._RKIND)  * num_pols_per_polp
  fac6 = pols_to_polp_frac_rh * rho_pols/rho_polp_grass * &
         ((diam_pols)**3._RKIND)/((diam_polp_grass)**3._RKIND) * num_pols_per_polp
  fac7 = pols_to_polp_frac_rh * rho_pols/rho_polp_weed * & 
         ((diam_pols)**3._RKIND)/((diam_polp_weed)**3._RKIND)  * num_pols_per_polp
  fac8 = (1._RKIND - pols_to_polp_frac_rh)
  
  if ( do_pollen_rh_rupture ) then
  ! Loop over the grid cells to simulate pollen rupture due to humidity, polp --> pols
    do j = jts, jte
    do k = kts, kte
    do i = its, ite
     ! Rupture only occurs for RH > 80%
       if ( (relhum(i,k,j) * 100._RKIND) < rh_rupt ) cycle
     ! Convert polp to pols due to humidty rupture
       if ( p_pols_all .gt. 0 ) then
          if ( p_polp_all .gt. 0) then
             chem(i,k,j,p_pols_all)  = chem(i,k,j,p_pols_all) + chem(i,k,j,p_polp_all)    * fac4
          else

   
             chem(i,k,j,p_pols_all)  = chem(i,k,j,p_pols_all) + ((chem(i,k,j,p_polp_tree) *fac5)  +   &
                                                                 (chem(i,k,j,p_polp_grass)*fac6) +   &
                                                                 (chem(i,k,j,p_polp_weed) *fac7)      )
          endif
       else
          
          chem(i,k,j,p_pols_tree)  = chem(i,k,j,p_pols_tree)  + chem(i,k,j,p_polp_tree)  *fac5
          chem(i,k,j,p_pols_grass) = chem(i,k,j,p_pols_grass) + chem(i,k,j,p_polp_grass) *fac6 
          chem(i,k,j,p_pols_weed)  = chem(i,k,j,p_pols_weed)  + chem(i,k,j,p_polp_weed)  *fac7
       endif
       ! Remove the converted amount from polp
       if ( p_polp_all .gt. 0 ) then
          chem(i,k,j,p_polp_all)   = chem(i,k,j,p_polp_all)   * fac8
       else
          chem(i,k,j,p_polp_tree)  = chem(i,k,j,p_polp_tree)  * fac8
          chem(i,k,j,p_polp_grass) = chem(i,k,j,p_polp_grass) * fac8
          chem(i,k,j,p_polp_weed)  = chem(i,k,j,p_polp_weed)  * fac8
       endif

    enddo
    enddo
    enddo
   endif ! do_pollen_rh_rupture

   if (do_pollen_lightning_rupture) then
! Compute the lightning flash rates, first intitialize
    do j=jts,jte
    do i=its,ite
       ic_flashrate(i,j) = 0._RKIND
       cg_flashrate(i,j) = 0._RKIND
    enddo
    enddo
! ----------
    do j=jts,jte
    do i=its,ite
    if ( total_flashrate(i,j) .gt. 0.) then
      ! Look for freezing level
        kfreeze = ktop2d(i,j)
        do while ( t(i,kfreeze,j) .lt. 273.15 .and. kfreeze .gt. 1 )
            kfreeze = kfreeze - 1
        enddo
      ! Calculate the depth between ktop and kfreeze (km)
        depth = ( z_at_w(i,ktop2d(i,j),j) - z_at_w(i,kfreeze,j) ) * 1.E-3_RKIND + cldtop_adjustment
        if (depth .le. 0.) continue
        depth = max( dH_min, min( dH_max, depth ))

        ratio = (((coef_A*depth+coef_B )*depth+coef_C)*depth+coef_D)*depth+coef_E
        cgfrac = 1._RKIND / (ratio+1._RKIND)

        cg_flashrate(i,j) = total_flashrate(i,j) * cgfrac
        ic_flashrate(i,j) = total_flashrate(i,j) - cg_flashrate(i,j)
    endif
    enddo
    enddo

    fac9  =  pols_to_polp_frac_lt * num_pols_per_polp * rho_pols/rho_polp       * &
             ((diam_pols)**3._RKIND)/((diam_polp)**3._RKIND) * flashrate_for_rupture
    fac10 =  pols_to_polp_frac_lt * num_pols_per_polp * rho_pols/rho_polp_tree  * &
             ((diam_pols)**3._RKIND)/((diam_polp_tree)**3._RKIND) * flashrate_for_rupture
    fac11 =  pols_to_polp_frac_lt * num_pols_per_polp * rho_pols/rho_polp_grass * &
             ((diam_pols)**3._RKIND)/((diam_polp_grass)**3._RKIND) * flashrate_for_rupture
    fac12 =  pols_to_polp_frac_lt * num_pols_per_polp * rho_pols/rho_polp_weed  * &
             ((diam_pols)**3._RKIND)/((diam_polp_weed)**3._RKIND) * flashrate_for_rupture
    fac13 =  1._RKIND - (pols_to_polp_frac_lt * flashrate_for_rupture)

  ! Loop over the grid cells to simulate pollen rupture due to lightning, polp --> pols
    do j = jte, jte
    do k = kte, kts, -1
    do i = its, ite
       if (cldfrac(i,k,j)>0.3) then
          flashrate_for_rupture = ic_flashrate(i,j)
       else
          flashrate_for_rupture = cg_flashrate(i,j)
       endif
     ! Convert polp->pols due to lightning
       if ( p_pols_all .gt. 0 ) then
         if (p_polp_all .gt. 0 ) then
          chem(i,k,j,p_pols_all)  = chem(i,k,j,p_pols_all) +  chem(i,k,j,p_polp_all)*fac9
         else
          chem(i,k,j,p_pols_all)  = chem(i,k,j,p_pols_all) + ((chem(i,k,j,p_polp_tree) *fac10)  +   &
                                                                 (chem(i,k,j,p_polp_grass)*fac11) +   &
                                                                 (chem(i,k,j,p_polp_weed) *fac12)      )
         endif
       else
          chem(i,k,j,p_pols_tree)  = chem(i,k,j,p_pols_tree)  + chem(i,k,j,p_polp_tree)  *fac10
          chem(i,k,j,p_pols_grass) = chem(i,k,j,p_pols_grass) + chem(i,k,j,p_polp_grass) *fac11
          chem(i,k,j,p_pols_weed)  = chem(i,k,j,p_pols_weed)  + chem(i,k,j,p_polp_weed)  *fac12
       endif 
     ! Remove the converted amount from polp
       if ( p_polp_all .gt. 0 ) then
          chem(i,k,j,p_polp_all)   = chem(i,k,j,p_polp_all)   * fac13
       else
          chem(i,k,j,p_polp_tree)  = chem(i,k,j,p_polp_tree)  * fac13
          chem(i,k,j,p_polp_grass) = chem(i,k,j,p_polp_grass) * fac13
          chem(i,k,j,p_polp_weed)  = chem(i,k,j,p_polp_weed)  * fac13
       endif
    enddo
    enddo
    enddo
  endif ! if do_pollen_lightning_rupture
   end subroutine pollen_driver

  
end module pollen_mod
