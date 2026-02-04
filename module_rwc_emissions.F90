!>\file  module_wildfire_smoke_emissions.F90
!! This file contains the MPAS-Aerosols/RRFS wildfire emission module

module module_rwc_emissions
!
!  This module developed by Johana Romero-Alvarez and Jordan Schnell (NOAA GSL)
!  For serious questions contact johana.romero-alvarez@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init

  implicit none

  private

  public :: mpas_smoke_rwc_emis_driver

contains


  subroutine mpas_smoke_rwc_emis_driver(ktau,dt,gmt,julday,krwc,                     &
                           xlat,xlong,xland,chem,num_chem,dz8w,t_phy,rho_phy,        &
                           rwc_emis_scale_factor,                                    &
                           RWC_denominator,RWC_annual_sum,           &
                           RWC_annual_sum_smoke_fine, RWC_annual_sum_smoke_coarse,   &
                           RWC_annual_sum_unspc_fine, RWC_annual_sum_unspc_coarse,   &
                           e_ant_out,  num_e_ant_out,         &
                           index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
                           index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse,   &
                           index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
                           index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse, &
                           ids,ide, jds,jde, kds,kde,                                &
                           ims,ime, jms,jme, kms,kme,                                &
                           its,ite, jts,jte, kts,kte                                 )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ktau,julday, num_chem, krwc,           &
                                  ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte,         &
                                  num_e_ant_out,       &
           index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
           index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse,   &
           index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
           index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse

   REAL(RKIND), INTENT(IN    ) :: dt,gmt,rwc_emis_scale_factor

   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: xlat,xlong,xland
   REAL(RKIND),DIMENSION(ims:ime,1:krwc,jms:jme),INTENT(IN) :: RWC_annual_sum_smoke_fine,RWC_annual_sum_smoke_coarse, &
                                                        RWC_annual_sum_unspc_fine,RWC_annual_sum_unspc_coarse, &
                                                        RWC_annual_sum
   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: RWC_denominator
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN) :: dz8w,rho_phy,t_phy
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_e_ant_out),INTENT(INOUT) :: e_ant_out
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT)     :: chem
                                                                               
  ! local
   INTEGER :: i,j,k,n,kemit
   INTEGER, PARAMETER :: offset_from_kts = 1
   REAL(RKIND) :: conv_aer, conv_gas, emis, t_phy_f, frac
   REAL(RKIND), DIMENSION(its:ite,jts:jte) :: rwc_t_thresh_grid

   REAL(RKIND), PARAMETER :: rwc_t_thresh = 283.15_RKIND ! [ 50 F]
   REAL(RKIND), PARAMETER :: spd_r = 1./86400._RKIND
   REAL(RKIND), PARAMETER :: emis_max = 1.0_RKIND

! TODO, read in TBL or define otherwise  
   rwc_t_thresh_grid(:,:) = rwc_t_thresh

! For now, just the surface
   k = kts
   kemit = kts + offset_from_kts

   do j = jts, jte
   do i = its, ite
     ! Is it cold enough to emit wood burning emissions? And are we over land?
      if ( t_phy(i,k,j) .lt. rwc_t_thresh_grid(i,j) .and. (( xland(i,j)-1.5) .lt. 0.) .and. RWC_denominator(i,j) .gt. 0.) then
        ! Conversion factor for aerosol emissions (ug/m2/s) --> ug/kg
         conv_aer = dt / (rho_phy(i,k,j) *  dz8w(i,k,j))
        ! Conversion factor for gas phase emissions (mol/m2/s) --> ppm/ppm
         conv_gas = 60._RKIND * 1.E6_RKIND * 4.828E-4_RKIND * dt / ( rho_phy(i,k,j) * dz8w(i,k,j) )
        ! Convert temperature to Fahrenheit
         t_phy_f = 1.8_RKIND * (t_phy(i,k,j)-273.15_RKIND) + 32._RKIND
        ! Calculate the fraction of total emisisons based on the linear equation, convert from /day to /sec
         frac = (42.12_RKIND - 0.79_RKIND*t_phy_f) / RWC_denominator(i,j) * spd_r
         if ( p_smoke_fine .gt. 0 .and. index_e_ant_out_smoke_fine .gt. 0 ) then
            emis = min(rwc_emis_scale_factor * conv_aer * frac * RWC_annual_sum_smoke_fine(i,1,j),emis_max)
            chem(i,kemit,j,p_smoke_fine) = chem(i,kemit,j,p_smoke_fine) + emis
            e_ant_out(i,kemit,j,index_e_ant_out_smoke_fine) = e_ant_out(i,kemit,j,index_e_ant_out_smoke_fine) + emis
         endif
         if ( (p_smoke_coarse .gt. 0 .or. p_unspc_coarse .gt. 0 ) .and. index_e_ant_out_smoke_coarse .gt. 0 ) then
            emis = min(rwc_emis_scale_factor * conv_aer * frac * RWC_annual_sum_smoke_coarse(i,1,j),emis_max)
            if ( p_smoke_coarse .gt. 0 ) then
                    chem(i,kemit,j,p_smoke_coarse) = chem(i,kemit,j,p_smoke_coarse) + emis
            elseif ( p_unspc_coarse .gt. 0 ) then
                    chem(i,kemit,j,p_unspc_coarse) = chem(i,kemit,j,p_unspc_coarse) + emis
            endif
            e_ant_out(i,kemit,j,index_e_ant_out_smoke_coarse) =  e_ant_out(i,kemit,j,index_e_ant_out_smoke_coarse) + emis
         endif
         if ( p_unspc_fine .gt. 0 .and. index_e_ant_out_unspc_fine .gt. 0 ) then
              emis = min(rwc_emis_scale_factor * conv_aer * frac * RWC_annual_sum_unspc_fine(i,1,j),emis_max)
              chem(i,kemit,j,p_unspc_fine) = chem(i,kemit,j,p_unspc_fine) + emis
              e_ant_out(i,kemit,j,index_e_ant_out_unspc_fine) = e_ant_out(i,kemit,j,index_e_ant_out_unspc_fine) + emis
         endif
         if ( p_unspc_coarse .gt. 0 .and. index_e_ant_out_unspc_coarse .gt. 0 ) then
              emis = min(rwc_emis_scale_factor * conv_aer * frac * RWC_annual_sum_unspc_coarse(i,1,j),emis_max)
              chem(i,kemit,j,p_unspc_coarse) = chem(i,kemit,j,p_unspc_coarse) + emis
              e_ant_out(i,kemit,j,index_e_ant_out_unspc_coarse) = e_ant_out(i,kemit,j,index_e_ant_out_unspc_coarse) + emis
         endif
      endif ! Temp check
   enddo
   enddo        

  end subroutine mpas_smoke_rwc_emis_driver

end module module_rwc_emissions
