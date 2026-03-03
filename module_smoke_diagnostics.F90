!>\file  module_smoke_diagnostics.F90
!! This file contains the MPAS-Aerosols/RRFS diagnostics module

module module_smoke_diagnostics
!
!  Jordan Schnell (NOAA GSL)
!  For serious questions contact jordan.schnell@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init
  use rad_data_mod

  implicit none

  private

  public :: mpas_aod_diag, mpas_visibility_diag, mpas_int_sm_diag

contains
  subroutine mpas_aod_diag(chem,aod3d,rho_phy,dz8w,num_chem,        &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte         )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte,         &
                                  num_chem

  REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: rho_phy, dz8w
  REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(IN) :: chem
  REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: aod3d

  real(RKIND) :: ext

  integer:: i,k,j,nv
  
  aod3d(:,:,:) = 0._RKIND
  do nv = 1, num_chem
  if ( nv .eq. p_ch4 ) cycle
  ext = sc_eff(nv) + ab_eff(nv)
  do j = jts, jte
  do k = kts, kte
  do i = its, ite
     aod3d(i,k,j)= aod3d(i,k,j) +  1.e-6_RKIND * ext * chem(i,k,j,nv)*rho_phy(i,k,j)*dz8w(i,k,j)
  enddo
  enddo
  enddo
  enddo

  end subroutine mpas_aod_diag

  subroutine mpas_int_sm_diag(chem,int_sm,centroid_h,sm_top,sm_bot,  &
                                  rho_phy,dz8w,num_chem,z_at_w,     &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte         )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte,         &
                                  num_chem

  REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: rho_phy, dz8w, z_at_w
  REAL(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(IN) :: chem
  REAL(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: int_sm, centroid_h, sm_top, sm_bot

  real(RKIND), parameter :: smoke_threshold = 1.0e-9

  integer:: i,k,j,nv

  sm_top(:,:)     = 0.0_RKIND
  sm_bot(:,:)     = -1.0_RKIND
  int_sm(:,:)     = 0._RKIND
  centroid_h(:,:) = 0._RKIND

  real(RKIND) :: layer_mass

  do nv = 1, num_chem
    if ( nv == p_smoke_fine .or. nv == p_smoke_coarse ) then
      do j = jts, jte
        do i = its, ite
          do k = kts, kte

            if (chem(i,k,j,nv) > smoke_threshold) then
             if (sm_bot(i,j) < 0.0_RKIND) then
                sm_bot(i,j) = z_at_w(i,k,j) ! Smoke layer bottom MSL
             endif
             sm_top(i,j) = z_at_w(i,k,j) ! Smoke layer top MSL
            endif

            ! mass increment for the layer
            layer_mass = 1.e-6_RKIND * chem(i,k,j,nv) * rho_phy(i,k,j) * dz8w(i,k,j)
            
            ! column integrated smoke
            int_sm(i,j) = int_sm(i,j) + layer_mass
            
            ! Accumulate mass-weighted height
            centroid_h(i,j) = centroid_h(i,j) + (layer_mass * z_at_w(i,k,j))
          enddo
        enddo
      enddo
    endif
  enddo

  ! Sanity check
  where (sm_bot < 0.0_RKIND) sm_bot = 0.0_RKIND

  do j = jts, jte
    do i = its, ite
      if (int_sm(i,j) > 1.e-12_RKIND) then 
         centroid_h(i,j) = centroid_h(i,j) / int_sm(i,j)
      else ! zero for clear air
         centroid_h(i,j) = 0.0_RKIND
      endif
    enddo
  enddo

  end subroutine mpas_intsm_diag

  subroutine mpas_visibility_diag(qcloud,qrain,qice,qsnow,qgrpl,    &
                                  blcldw,blcldi,                    &
                                  rho_phy,wind10m,wind,             &
                                  rh2m,rh,qv,t2m,t,                 &
                                  coszen,extcoef55,vis,             &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte         )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte

   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN)   :: qcloud,qrain,qice,qsnow,qgrpl
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN)   :: blcldi,blcldw
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN)   :: rho_phy,wind,rh,qv,t,extcoef55
   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN)     :: wind10m,rh2m,t2m,coszen
   REAL(RKIND),DIMENSION(ims:ime,jms:jme), INTENT(OUT)   :: vis
                                                                               
  ! local
   INTEGER :: i,j,k,d
   REAL(RKIND), PARAMETER :: visfactor = 3.912_RKIND
   REAL(RKIND) :: bc, br, bi, bs, hydro_extcoeff, extcoeff, vis_haze
   REAL(RKIND) :: tvd, qrh, prob_ext_coeff_gt_p29, haze_ext_coeff
   REAL(RKIND) :: vis_hydlith, alpha_haze,rhmax
   REAL(RKIND) :: bg,qcloud2,blcldw2,qrain2,qice2,blcldi2,qsnow2,qgrpl2,extcoeff552,vis_night,zen_fac
   REAL(RKIND) :: tiny_number
   REAL(RKIND), DIMENSION(its:ite,jts:jte) :: vis_alpha

   tiny_number = 1.e-12_RKIND

   do j = jts,jte
   do i = its,ite
      !Initialize
      qcloud2 = 0._RKIND
      blcldw2 = 0._RKIND
      qrain2  = 0._RKIND
      qice2   = 0._RKIND
      blcldi2 = 0._RKIND
      qsnow2  = 0._RKIND
      qgrpl2  = 0._RKIND
      extcoeff552 = 0._RKIND
      ! Follwowing UPP: CALVIS_GSD.f, take max of hydrometeors in lowest 3 levels   
      ! - in UPP, only bottom rho_phy is used, shouldn't we use rho_phy from that level (as below)?
      k = kts
!      do k = 1,3
         qcloud2    = qcloud(i,k,j)*rho_phy(i,k,j)*1000._RKIND !max(qcloud2,qcloud(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
         blcldw2    = blcldw(i,k,j)*rho_phy(i,k,j)*1000._RKIND !max(blcldw2,blcldw(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
         qrain2     = qrain(i,k,j)*rho_phy(i,k,j)*1000._RKIND !max(qrain2,qrain(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
         qice2      = qice(i,k,j)*rho_phy(i,k,j)*1000._RKIND ! max(qice2,qice(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
         blcldi2    = blcldi(i,k,j)*rho_phy(i,k,j)*1000._RKIND !max(blcldi2,blcldi(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
         qsnow2     = qsnow(i,k,j)*rho_phy(i,k,j)*1000._RKIND !max(qsnow2,qsnow(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
         qgrpl2     = qgrpl(i,k,j)*rho_phy(i,k,j)*1000._RKIND !max(qgrpl2,qgrpl(i,k,j)*rho_phy(i,k,j)*1000._RKIND)
         extcoeff552= extcoef55(i,k,j)*1.e-3_RKIND !max(extcoeff552,extcoeff55(i,k,j)) ! JLS - EXT55 is in units = 1/km, covert to 1/m
!      enddo

      bc = 144.7_RKIND * (qcloud2+blcldw2)  ** 0.88_RKIND
      br = 2.24_RKIND  * qrain2  ** 0.75_RKIND
      bi = 327.8_RKIND * (qice2+blcldi2)  ** 1.00_RKIND
      bs = 10.36_RKIND * qsnow2  ** 0.7776_RKIND
      bg = 8.0_RKIND   * qgrpl2  ** 0.7500_RKIND

      hydro_extcoeff=(bc+br+bi+bs+bg)*1.e-3_RKIND ! m^-1

! TODO, JLS, bring in QV 2M
      vis_haze=999999._RKIND
      IF (qv(i,1,j) .gt. 0._RKIND) THEN
        vis_haze=1500._RKIND*(105._RKIND-rh2m(i,j))*(5._RKIND/min(1000._RKIND*qv(i,1,j),5._RKIND))
      ENDIF

      ! Follow UPP/CALVIS_GSD, First compute max RH of lowest 2 layers
      rhmax = max( rh(i,1,j), rh(i,2,j) )
      qrh = max(0._RKIND,min(0.8_RKIND,(rhmax*0.01_RKIND-0.15_RKIND)))
      vis_haze = 90000._RKIND * exp(-2.5_RKIND*qrh) 

      ! Calculate a Weibull function "alpha" term.  This can be
      ! used later with visibility (which acts as the "beta" term
      ! in the Weibull function) to create a probability distribution
      ! for visibility. Alpha can be thought of as the "level of
      ! certainty" that beta (model visibility) is correct. Fog is
      ! notoriously difficult to model. In the below algorithm,
      ! the alpha value (certainty) decreases as PWAT, mixing ratio,
      ! or winds decrease (possibly foggy conditions), but increases
      ! if RH decreases (more certainly not foggy).  If PWAT is lower
      ! then there is a higher chance of radiation fog because there
      ! is less insulating cloud above.
      ! -------------------------------------------------------------

!      alpha_haze=3.6
!      IF (q2m(i,j) .gt. 0.) THEN
!        alpha_haze=0.1 + pwater(i,j)/25.     + wind125m(i,j)/3. + &
!                        (100.-rh2m(i,j))/10. + 1./(1000.*q2m(i,j))
!        alpha_haze=min(alpha_haze,3.6)
!      ENDIF

      ! Calculate visibility from hydro/lithometeors
      ! Maximum visibility -> 999999 meters
      ! --------------------------------------------
      extcoeff=hydro_extcoeff+extcoeff552
      IF (extcoeff .gt. tiny_number ) THEN
        vis_hydlith=min(visfactor/extcoeff, 999999._RKIND)
      ELSE
        vis_hydlith=999999._RKIND
      ENDIF

      ! -- Dec 2003 - Roy Rasmussen (NCAR) expression for night vs. day vis
      !   1.609 factor is number of km in mile.
      vis_night = 1.69 * ((vis_hydlith/1.609)**0.86) * 1.609
      zen_fac = min(0.1,max(coszen(i,j),0.)) * 10._RKIND
      vis_hydlith = zen_fac * vis_hydlith + (1.-zen_fac)*vis_night

      ! Calculate total visibility
      ! Take minimum visibility from hydro/lithometeors and haze
      ! Set alpha to be alpha_haze if haze dominates, or 3.6
      ! (a Guassian distribution) when hydro/lithometeors dominate
      ! ----------------------------------------------------------
      IF (vis_hydlith < vis_haze) THEN
         vis(i,j)=vis_hydlith
         !vis_alpha(i,j)=3.6
      ELSE
         vis(i,j)=vis_haze
         !vis_alpha(i,j)=alpha_haze
      ENDIF

   enddo ! i
   enddo ! j
  
  end subroutine mpas_visibility_diag

end module module_smoke_diagnostics
