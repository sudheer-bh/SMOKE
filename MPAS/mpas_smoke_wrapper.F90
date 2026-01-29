!>\file mpas_smoke_wrapper.F90
!! This file is MPAS Smoke wrapper
!! Haiqin.Li@noaa.gov 09/2024

module mpas_smoke_wrapper

   use mpas_kind_types
   use mpas_pool_routines
   use mpas_constants,        only: cp
   use mpas_timer, only : mpas_timer_start, mpas_timer_stop
   use mpas_smoke_config
   use mpas_smoke_init 
   use module_plumerise,      only : ebu_driver
   use module_fire_emissions
   use module_add_emiss_burn, only : add_emis_burn
   use dep_dry_simple_mod,    only : dry_dep_driver_simple
   use dep_dry_mod_emerson,   only : dry_dep_driver_emerson, particle_settling_wrapper
   use dep_data_mod,          only : aero_dry_dep_init, aero_wet_dep_init 
   use rad_data_mod,          only : aero_rad_init
   use module_wetdep_ls,      only : wetdep_ls
   use dust_fengsha_mod,      only : gocart_dust_fengsha_driver
   use pollen_mod,            only : pollen_driver
   use module_anthro_emissions
   use module_rwc_emissions
   use module_tactic_sna
   use ssalt_mod
   use module_smoke_diagnostics

   implicit none

   private

   public :: mpas_smoke_driver

contains

    subroutine mpas_smoke_driver(                                                            &
           num_chem              , chemistry_start             , chem           ,            &
           config_extra_chemical_tracers,                                                    &
           kanthro    , kbio, kfire, kvol, krwc,                                             &
           config_ultrafine, config_coarse,                                                  &
           index_smoke_ultrafine , index_smoke_fine            , index_smoke_coarse,         &
           index_dust_ultrafine  , index_dust_fine             , index_dust_coarse,          &
           index_ssalt_fine      , index_ssalt_coarse          ,                             &
           index_polp_tree       , index_polp_grass            , index_polp_weed,            &
           index_polp_all        , index_pols_all              ,                             &
           index_pols_tree       , index_pols_grass            , index_pols_weed,            &
           index_unspc_ultrafine , index_unspc_fine            , index_unspc_coarse,         &
           index_no3_a_fine      , index_so4_a_fine            , index_nh4_a_fine,           &
           index_nh3             , index_so2                   , index_ch4,                  &
           index_co              , index_nox                   , index_bact_fine,            &
           index_e_bb_in_smoke_ultrafine, index_e_bb_in_smoke_fine, index_e_bb_in_smoke_coarse, &
           index_e_bb_in_co, index_e_bb_in_nh3, index_e_bb_in_ch4,                           &
           index_e_bb_in_nox, index_e_bb_in_so2, &
           index_e_ant_in_unspc_ultrafine, index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,  &
           index_e_ant_in_smoke_ultrafine,index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse,   &
           index_e_ant_in_no3_a_fine, index_e_ant_in_so4_a_fine, &
           index_e_ant_in_nh4_a_fine, index_e_ant_in_so2, &
           index_e_ant_in_nh3, index_e_ant_in_ch4,                                           &
           index_e_ant_in_co, index_e_ant_in_nox,                                            &
           index_e_bio_in_polp_tree, index_e_bio_in_polp_grass, index_e_bio_in_polp_weed,    &
           index_e_vol_in_vash_fine, index_e_vol_in_vash_coarse,                             &
           index_e_bb_out_smoke_ultrafine, index_e_bb_out_smoke_fine, index_e_bb_out_smoke_coarse, & 
           index_e_bb_out_nox, index_e_bb_out_ch4,index_e_bb_out_co,                         &
           index_e_bb_out_so2,index_e_bb_out_nh3,       &
           index_e_ant_out_unspc_ultrafine, index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse,   &
           index_e_ant_out_smoke_ultrafine, index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse,   &
           index_e_ant_out_no3_a_fine, index_e_ant_out_so4_a_fine, &
           index_e_ant_out_nh4_a_fine, index_e_ant_out_so2, &
           index_e_ant_out_nh3, index_e_ant_out_ch4, &
           index_e_ant_out_co, index_e_ant_out_nox, &
           index_e_bio_out_polp_tree, index_e_bio_out_polp_grass, index_e_bio_out_polp_weed, &
           index_e_vol_out_vash_fine,  index_e_vol_out_vash_coarse,                          &
           index_e_dust_out_dust_ultrafine, index_e_dust_out_dust_fine, index_e_dust_out_dust_coarse,  &
           index_e_ss_out_ssalt_fine, index_e_ss_out_ssalt_coarse,                           &
           frp_in                , frp_out,    fre_in, fre_out, hwp,  coef_bb_dc          ,  &
           totprcp_prev24        , hwp_avg     , frp_avg,    fre_avg, fire_end_hr,           &
           nblocks               , EFs_map,                                                  &
           eco_id, efs_smold, efs_flam, efs_rsmold,fmc_avg,                                  &
           hfx_bb                , qfx_bb         ,  frac_grid_burned    ,                   &
           min_bb_plume          , max_bb_plume,                                             &
           sandfrac_in           , clayfrac_in           , uthres_in            ,            &
           uthres_sg_in          , albedo_drag_in        , feff_in              ,            &
           sep_in                ,                                                           &
           e_ant_in, e_bb_in, e_bio_in, e_vol_in,                                            &
           e_ant_out, e_bb_out, e_bio_out, e_dust_out, e_ss_out, e_vol_out,                  &
           num_e_ant_in, num_e_bb_in, num_e_bio_in, num_e_vol_in,                            &
           num_e_ant_out, num_e_bb_out, num_e_bio_out,                                       &
           num_e_dust_out, num_e_ss_out, num_e_vol_out   , drydep_flux           ,           &
           ddvel                 , wetdep_resolved       , tend_chem_settle      ,           &
           do_mpas_smoke         , do_mpas_dust          , do_mpas_pollen        ,           &
           do_mpas_anthro        , do_mpas_ssalt         , do_mpas_volc          ,           &
           do_mpas_sna           , do_mpas_methane       , do_mpas_hab           ,           &
           do_mpas_rwc           ,                                                           &
           calc_bb_emis_online   , bb_beta               ,                                   &
           hwp_method            , hwp_alpha             , wetdep_ls_opt        ,            &
           wetdep_ls_alpha       , plumerise_opt         , plume_wind_eff       ,            &
           plume_alpha           , bb_emis_scale_factor, ebb_dcycle             ,            &
           drydep_opt            , pm_settling           , add_fire_heat_flux   ,            &
           add_fire_moist_flux   , plumerisefire_frq     , bb_qv_scale_factor   ,            &
           dust_alpha            , dust_gamma            , dust_drylimit_factor ,            &
           dust_moist_correction ,                                                           &
           num_pols_per_polp     , pollen_emis_scale_factor,                                 &
           tree_pollen_emis_scale_factor, grass_pollen_emis_scale_factor        ,            &
           weed_pollen_emis_scale_factor,                                                    &
           bb_input_prevh        , rwc_emis_scale_factor,                 &
           RWC_denominator       , RWC_annual_sum       ,                                    &
           RWC_annual_sum_smoke_fine, RWC_annual_sum_smoke_coarse,                           &
           RWC_annual_sum_unspc_fine, RWC_annual_sum_unspc_coarse,                           &
           nwfa                  , nifa                 ,  vis                  ,            &
           qc_vis, qr_vis, qi_vis, qs_vis, qg_vis, blcldw_vis, blcldi_vis,                   &
           hno3_bkgd             , coszen                , aod3d_smoke, aod3d   ,            &
           ktau                  , dt                    , dxcell               ,            &
           area                  ,                                                           &
           xland                 , u10                   , v10                  ,            &
           ust                   , xlat                  , xlong                ,            &
           tskin                 , pblh                  , t2m                  ,            &
           p8w                   , dz8w                  , z_at_w               ,            &
           p_phy                 , t_phy                 , u_phy                ,            &
           v_phy                 , qv                    , vvel                 ,            &
           pi_phy                , rho_phy               , kpbl                 ,            &
           nsoil                 , smois                 , tslb                 ,            &
           ivgtyp                , isltyp                , nlcat                ,            &
           swdown                , z0                    , snowh                ,            &
           julian                , rmol                  , raincv               ,            &
           rainncv               , dpt2m                 , znt                  ,            &
           mavail                , g                     , vegfra               ,            &
           landusef              , cldfrac               , ktop_deep            ,            &
           cp                    , rd                    , gmt                  ,            &
           ids       , ide       , jds       , jde       , kds       , kde      ,            &
           ims       , ime       , jms       , jme       , kms       , kme      ,            &
           its       , ite       , jts       , jte       , kts       , kte                   &
                                                                                             )
    implicit none

! intent arguments:
! array indexes
    integer,intent(in):: ids,ide,jds,jde,kds,kde,        &
                         ims,ime,jms,jme,kms,kme,        &
                         its,ite,jts,jte,kts,kte
! Timestep, day, constants
    real(RKIND),intent(in):: dt, julian, g, cp, rd, gmt
! Time step #
    integer,intent(in):: ktau
    integer,intent(in)::nblocks
! Dimensions and indexes
    integer,intent(in):: nsoil, nlcat, num_chem, chemistry_start
    integer,intent(in):: kanthro, kbio, kfire, kvol, krwc
    integer,intent(in):: num_e_ant_in,  num_e_bb_in,  num_e_bio_in,  num_e_vol_in
    integer,intent(in):: num_e_ant_out, num_e_bb_out, num_e_bio_out, num_e_dust_out, num_e_ss_out, num_e_vol_out
! 2D mesh arguments
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)             :: xlat, xlong, dxcell, area, xland   ! grid
! 2D Met input
    integer,intent(in), dimension(ims:ime, jms:jme)                :: isltyp, ivgtyp ! domainant soil, vegetation type
    integer,intent(in), dimension(ims:ime, jms:jme)                :: kpbl          ! k-index of PBLH
    integer,intent(in), dimension(ims:ime, jms:jme),optional       :: ktop_deep
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)            :: u10, v10      ! 10-m winds
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)            :: tskin, t2m, dpt2m            ! temperature
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)            :: pblh              ! PBL height [m]
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)            :: vegfra
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)            :: swdown, z0, snowh, znt
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)            :: coszen
    real(RKIND),intent(in), dimension(ims:ime, jms:jme)            :: raincv, rainncv, mavail                    
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme)         :: rmol, ust
! 2D Fire Input
    real(RKIND),intent(in), dimension(ims:ims, jms:jme), optional      :: totprcp_prev24, fire_end_hr,fmc_avg,     &
                                                                          efs_smold, efs_flam, efs_rsmold
    integer,intent(in), dimension(ims:ime,jms:jme),optional            :: eco_id
    real(RKIND),intent(in),dimension(ims:ime, jms:jme),optional        :: frp_in, fre_in      ! Fire input
! 2D + Time Fire Input
    real(RKIND),intent(in), dimension(ims:ims, jms:jme, nblocks),        &
                                                   optional      :: hwp_avg, fre_avg, frp_avg
! Residential Wood burning
    real(RKIND),intent(in), dimension(ims:ims, jms:jme),optional    :: RWC_denominator
    real(RKIND),intent(in), dimension(ims:ims, 1:krwc,jms:jme),optional :: RWC_annual_sum,                        &
                                                                           RWC_annual_sum_smoke_fine, RWC_annual_sum_smoke_coarse, &
                                                                           RWC_annual_sum_unspc_fine, RWC_annual_sum_unspc_coarse
! 3D Met input 
    real(RKIND),intent(in), dimension(ims:ime, kms:kme, jms:jme)   :: p8w,    dz8w,    z_at_w, cldfrac,   &
                                                                       p_phy,  t_phy,   u_phy,  v_phy,     &
                                                                       pi_phy, rho_phy, vvel  
    real(RKIND),intent(inout),dimension(ims:ime, kms:kme, jms:jme),optional  :: nifa, nwfa, hno3_bkgd 
! 3D (2D + first 3 levels) for visibility calculations
    real(RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional :: qc_vis, qr_vis, qi_vis, qs_vis, qg_vis, blcldw_vis, blcldi_vis
! 3D emission input
    real(RKIND),intent(in), dimension(ims:ime,1:kanthro,jms:jme,1:num_e_ant_in),optional :: e_ant_in
    real(RKIND),intent(in), dimension(ims:ime,1:kfire,jms:jme,1:num_e_bb_in),optional  :: e_bb_in
    real(RKIND),intent(in), dimension(ims:ime,1:kbio,jms:jme,1:num_e_bio_in),optional  :: e_bio_in
    real(RKIND),intent(in), dimension(ims:ime,1:kvol,jms:jme,1:num_e_vol_in),optional  :: e_vol_in
! JLS - TODO, if we update QV via moist flux, we will need to update the scalar in the driver
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme),optional           :: qv
    real(RKIND),intent(in), dimension(ims:ime,1:nsoil, jms:jme)   ,optional            :: smois, tslb
    real(RKIND),intent(in), dimension(ims:ime,1:nlcat, jms:jme)  ,optional             :: landusef
! Chemistry indexes into MPAS scalar array
    integer, intent(in),optional :: &
                           index_smoke_ultrafine, index_smoke_fine, index_smoke_coarse,                    &
                           index_dust_ultrafine, index_dust_fine,  index_dust_coarse,                     &
                           index_polp_tree,  index_polp_grass,   index_polp_weed,   &
                           index_pols_tree,  index_pols_grass,   index_pols_weed,   &
                           index_pols_all,   index_polp_all,                        &
                           index_unspc_ultrafine, index_unspc_fine, index_unspc_coarse,                    &
                           index_ssalt_fine, index_ssalt_coarse,                    &
                           index_so4_a_fine, index_no3_a_fine, index_nh4_a_fine,    &
                           index_so2, index_nh3, index_ch4,                         &
                           index_co,  index_nox, index_bact_fine    
    integer, intent(in),optional :: &
                           index_e_bb_in_smoke_ultrafine, index_e_bb_in_smoke_fine, index_e_bb_in_smoke_coarse, &
                           index_e_bb_in_co, index_e_bb_in_ch4, index_e_bb_in_nox, &
                           index_e_bb_in_so2, index_e_bb_in_nh3, &
                           index_e_ant_in_unspc_ultrafine, index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse, &
                           index_e_ant_in_smoke_ultrafine, index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse, &
                           index_e_ant_in_no3_a_fine, index_e_ant_in_so4_a_fine, &
                           index_e_ant_in_nh4_a_fine, index_e_ant_in_so2, &
                           index_e_ant_in_nh3, index_e_ant_in_ch4,    &
                           index_e_ant_in_co,  index_e_ant_in_nox,    &
                           index_e_bio_in_polp_tree, index_e_bio_in_polp_grass, index_e_bio_in_polp_weed, &
                           index_e_vol_in_vash_fine,  index_e_vol_in_vash_coarse
    integer, intent(in),optional :: &
                           index_e_bb_out_smoke_ultrafine, index_e_bb_out_smoke_fine, index_e_bb_out_smoke_coarse, &
                           index_e_bb_out_ch4, index_e_bb_out_co, index_e_bb_out_nox, index_e_bb_out_nh3, index_e_bb_out_so2, &
                           index_e_ant_out_unspc_ultrafine, index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
                           index_e_ant_out_smoke_ultrafine, index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse, &
                           index_e_ant_out_no3_a_fine, index_e_ant_out_so4_a_fine, &
                           index_e_ant_out_nh4_a_fine, index_e_ant_out_so2, &
                           index_e_ant_out_nh3, index_e_ant_out_ch4, &
                           index_e_ant_out_co, index_e_ant_out_nox, &
                           index_e_bio_out_polp_tree, index_e_bio_out_polp_grass, index_e_bio_out_polp_weed, &
                           index_e_vol_out_vash_fine,  index_e_vol_out_vash_coarse, &
                           index_e_dust_out_dust_ultrafine, index_e_dust_out_dust_fine, index_e_dust_out_dust_coarse, &
                           index_e_ss_out_ssalt_fine, index_e_ss_out_ssalt_coarse
! 2D dust input arrays 
    real(RKIND),intent(in), dimension(ims:ime, jms:jme),optional  :: sandfrac_in, clayfrac_in, uthres_in, &        ! dust (FENGSHA) input
                                                                     uthres_sg_in, albedo_drag_in, feff_in, sep_in ! dust (FENGSHA) input
! 2D input/output arrays
    real(RKIND),intent(inout),dimension(ims:ime, jms:jme),optional :: frp_out, fre_out, EFs_map
    real(RKIND),intent(inout),dimension(ims:ime, jms:jme),optional :: hwp, coef_bb_dc
    real(RKIND),intent(inout),dimension(ims:ime, jms:jme),optional :: hfx_bb, qfx_bb, frac_grid_burned
    integer,intent(inout),dimension(ims:ime,jms:jme),optional      :: min_bb_plume, max_bb_plume
! 2D + chem output arrays
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme, 1:num_chem),optional :: wetdep_resolved
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme, 1:num_chem),optional :: drydep_flux
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme, 1:num_chem),optional :: ddvel
! 2D output arrays
    real(RKIND),intent(inout), dimension(ims:ime, jms:jme),optional                   :: vis
! 3D output emissions
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme,1:num_e_ant_out),optional   :: e_ant_out
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme,1:num_e_bb_out),optional    :: e_bb_out
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme,1:num_e_bio_out),optional   :: e_bio_out
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme,1:num_e_dust_out),optional  :: e_dust_out
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme,1:num_e_ss_out),optional    :: e_ss_out
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme,1:num_e_vol_out),optional   :: e_vol_out
! 3D + chem output arrays
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem)                :: chem
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),optional       :: tend_chem_settle
    real(RKIND),intent(inout), dimension(ims:ime, kms:kme, jms:jme),optional                   :: aod3d_smoke, aod3d
!>-- Namelist options
     logical,intent(in)               :: do_mpas_smoke
     logical,intent(in)               :: do_mpas_dust
     logical,intent(in)               :: do_mpas_pollen
     logical,intent(in)               :: do_mpas_anthro
     logical,intent(in)               :: do_mpas_ssalt
     logical,intent(in)               :: do_mpas_volc
     logical,intent(in)               :: do_mpas_sna
     logical,intent(in)               :: do_mpas_methane
     logical,intent(in)               :: do_mpas_hab
     logical,intent(in)               :: do_mpas_rwc
     character(len=*),intent(in)      :: config_extra_chemical_tracers
     logical,intent(in)               :: config_ultrafine, config_coarse
     logical,intent(in)               :: calc_bb_emis_online
     integer,intent(in)               :: hwp_method
     real(RKIND),intent(in)           :: hwp_alpha
     integer,intent(in)               :: wetdep_ls_opt
     real(kind=RKIND),intent(in)      :: wetdep_ls_alpha
     integer,intent(in)               :: plumerise_opt
     integer,intent(in)               :: plume_wind_eff
     real(kind=RKIND),intent(in)      :: plume_alpha
     real(kind=RKIND),intent(in)      :: bb_emis_scale_factor, bb_qv_scale_factor
     real(kind=RKIND),intent(in)      :: rwc_emis_scale_factor
     integer,intent(in)               :: ebb_dcycle
     integer,intent(in)               :: drydep_opt
     integer,intent(in)               :: pm_settling
     logical,intent(in)               :: add_fire_heat_flux
     logical,intent(in)               :: add_fire_moist_flux
     integer,intent(in)               :: plumerisefire_frq
     integer,intent(in)               :: bb_beta
     real(RKIND),intent(in)           :: dust_alpha, dust_gamma
     real(RKIND),intent(in)           :: dust_drylimit_factor, dust_moist_correction
     integer,intent(in)               :: bb_input_prevh
     real(RKIND),intent(in),optional  :: pollen_emis_scale_factor, num_pols_per_polp 
     real(RKIND),intent(in),optional  :: tree_pollen_emis_scale_factor, &
                                         grass_pollen_emis_scale_factor, &
                                         weed_pollen_emis_scale_factor

!----------------------------------
!>-- Local Variables
!>-- 3D met
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: rri,     &
                     wind_phy,theta_phy,zmid,t8w,relhum
    real(RKIND), dimension(ims:ime,jms:jme) :: rh2m
!>-- indexes, time
    integer :: julday
!>- dust & chemistry variables
    real(RKIND), dimension(ims:ime, 1:nlcat, jms:jme) :: vegfrac
    ! JLS, temporary, need to read in like SMOKE_RRFS/MPAS
    real(RKIND), dimension(ims:ime, jms:jme) :: total_flashrate
!>- plume variables
    ! -- buffers
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme, num_e_bb_in) :: ebu
    real(RKIND), dimension(ims:ime, jms:jme)          :: flam_frac,                               &
                                                         fire_hist, peak_hr,                      &
                                                         hwp_day_avg,                             &
                                                         uspdavg2d, hpbl2d, wind10m
    real(RKIND), dimension(ims:ime, jms:jme)          :: lu_nofire, lu_qfire, lu_sfire
    integer,     dimension(ims:ime, jms:jme)          :: fire_type
    integer,     dimension(ims:ime, jms:jme)          :: kpbl_thetav
    logical                                           :: call_plume
    real(RKIND), parameter                            :: conv_frpi   = 1.e-06_RKIND  ! FRP conversion factor, MW to W
    real(RKIND), parameter                            :: conv_frei   = 1.e-06_RKIND  ! FRE conversion factor, MW-s to W-s
!>- Dry deposition - temporary - move to output 
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem) :: vg     ! gravitational settling velocity
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme)             :: thetav, dz8w_flip
!> -- other
    real(RKIND)    :: theta
    real(RKIND)    :: curr_secs
    integer        :: nbegin, nv
    integer        :: i, j, k, kp, n
    character(100) :: errmsg
    integer        :: errflg
    logical        :: do_plumerise
! FRP/plumerise related thresholds
    real(RKIND), parameter :: frp_min        = 1.e+7     ! Minimum FRP (Watts) to distribute smoke in PBL, 10MW
    real(RKIND), parameter :: frp_max        = 2.e+10    ! Maximum FRP over 3km Pixel, 20,000 MW
    real(RKIND), parameter :: fre_min        = -999.     ! Minimum FRE (Watt seconds) to distrubute smoke in PBL, TODO
    real(RKIND), parameter :: fre_max        = -999.     ! Maximumm FRE (Watt seconds) "                       ", TODO     
    real(RKIND), parameter :: zpbl_threshold = 2.e+3     ! Minimum PBL depth to have plume rise 
    real(RKIND), parameter :: uspd_threshold = 5.        ! Wind speed averaged across PBL depth to control smoke release levels 
    real(RKIND), parameter :: frp_wthreshold = 1.e+9     ! Minimum FRP (Watts) to have plume rise in windy conditions
    real(RKIND), parameter :: fre_wthreshold = 1.e+9     ! Minimum FRE (Watt-seconds) to have plume rise in windy conditions
    real(RKIND), parameter :: ebb_min        = 1.e-3     ! Minimum smoke emissions (ug/m2/s)

    logical, parameter :: do_timing = .false.

    errmsg = ''
    errflg = 0
 
    call mpas_log_write( ' Beginning Aerosol Driver')

  ! If not simulating smoke or pollen, get outta here...
    if ( (.not. do_mpas_smoke) .and. (.not. do_mpas_pollen) .and. &
         (.not. do_mpas_dust ) .and. (.not. do_mpas_anthro) .and. &
         (.not. do_mpas_ssalt) .and. (.not. do_mpas_sna)    .and. &
         (.not. do_mpas_methane) .and. (.not. do_mpas_hab)  .and. &
         (.not. do_mpas_rwc) )  return
! 
   if (ktau == 1) then
!   Reorder chemistry indices
      call set_scalar_indices(chemistry_start,                                     &
                    index_smoke_ultrafine, index_smoke_fine, index_smoke_coarse,   &
                    index_dust_ultrafine, index_dust_fine, index_dust_coarse,      &
                    index_polp_tree, index_polp_grass, index_polp_weed,            &
                    index_pols_tree, index_pols_grass, index_pols_weed,            &
                    index_pols_all,  index_polp_all,                               &
                    index_unspc_ultrafine, index_unspc_fine, index_unspc_coarse,   &
                    index_ssalt_fine, index_ssalt_coarse,                          &
                    index_no3_a_fine, index_so4_a_fine, index_nh4_a_fine,          &
                    index_so2, index_nh3, index_ch4, index_nox, index_co,          &
                    index_bact_fine                                                )
                                                                          
      call mpas_log_write( ' Initializing dry deposition parameterss ')
      call aero_dry_dep_init()
      call mpas_log_write( ' Initializing wet deposition parameterss ')
      call aero_wet_dep_init()
      call mpas_log_write( ' Initializing radiation feedback parameterss ')
      call aero_rad_init()
   endif
!
    uspdavg2d   = 0._RKIND
    hpbl2d      = 0._RKIND
    peak_hr     = 0._RKIND
    flam_frac   = 0._RKIND
    fire_type   = 0

    curr_secs = ktau * dt
    julday = int(julian)

    call mpas_log_write( ' Calling smoke prep')
    if  (do_timing) call mpas_timer_start('smoke_prep')
    !>- get ready for chemistry run
    call mpas_smoke_prep(                                                   &
        do_mpas_smoke,                                                      &
        ktau, nlcat,cp,ebb_dcycle,ebb_min,                                  &
        xland,xlat,xlong,ivgtyp,isltyp,landusef,                            & ! JLS TODO LANDUSEF /= VEGFRA?
        snowh,u10,v10,wind10m,t2m,dpt2m,mavail,hwp,hwp_day_avg,             &
        index_e_bb_in_smoke_fine,num_e_bb_in,kfire,e_bb_in,                 &
        t_phy,u_phy,v_phy,p_phy,pi_phy,z_at_w,                              &
        dz8w,dz8w_flip,                                                     &
        rho_phy,qv,relhum,rh2m,rri,                                         &
        total_flashrate,                                                    &
        wind_phy,theta_phy,zmid,kpbl_thetav,                                &
        peak_hr,coef_bb_dc,fire_hist,                                       &
        lu_nofire, lu_qfire, lu_sfire, fire_type,                           &
        ids,ide, jds,jde, kds,kde,                                          &
        ims,ime, jms,jme, kms,kme,                                          &
        its,ite, jts,jte, kts,kte)
   if  (do_timing) call mpas_timer_stop('smoke_prep')

        
   if ( do_mpas_smoke ) then
! Are we calculating emissions online?
      if ( calc_bb_emis_online ) then
         call mpas_log_write( ' Calling module wildfire smoke emissions ')
           ! Calculate the wildfire emission factors
             call compute_emission_factors( EFs_map, vegfrac, bb_beta,                       &
                                            eco_id, efs_smold, efs_flam, efs_rsmold,         &
                                            fmc_avg, mavail,                                 &
                                            nlcat, dt, gmt,                                  &
                                            ids, ide, jds, jde, kds, kde,                    &
                                            ims, ime, jms, jme, kms, kme,                    &
                                            its, ite, jts, jte, kts, kte)
         !TODO, JR replace by dynamic EFs 
         if ( ebb_dcycle .eq. 2 ) then
           call calculate_smoke_emissions(   dt, julday, nlcat, EFs_map, fre_avg,            &
                                             ebb_dcycle, area, nblocks, ktau,                &
                                             bb_input_prevh,ebu,num_e_bb_in,                 &
                                             ids, ide, jds, jde, kds, kde,                   &
                                             ims, ime, jms, jme, kms, kme,                   &
                                             its, ite, jts, jte, kts, kte)
         else
           call calculate_smoke_emissions(   dt, julday, nlcat, EFs_map, fre_in,             &
                                             ebb_dcycle, area, 1, ktau,                      &
                                             1,ebu,num_e_bb_in,                              &
                                             ids, ide, jds, jde, kds, kde,                   &
                                             ims, ime, jms, jme, kms, kme,                   &
                                             its, ite, jts, jte, kts, kte)

         endif
      else ! i
         if (ktau==1) then
           do nv=1,num_e_bb_in
           do j=jts,jte
           do i=its,ite
           
             ebu(i,kts,j,nv)= e_bb_in(i,kts,j,nv)
             do k=kts+1,kte
              ebu(i,k,j,nv)= 0._RKIND
             enddo
           enddo
           enddo
           enddo
         else
           do nv=1,num_e_bb_in
           do j=jts,jte
           do k=kts,kte
           do i=its,ite
           ! ebu is divided by coef_bb_dc since it is applied in the output
             ebu(i,k,j,nv) = e_bb_out(i,k,j,nv) / MAX(1.E-4_RKIND,coef_bb_dc(i,j))
           enddo
           enddo
           enddo
           enddo
         endif ! ktau = 1
      endif ! calc emis_online

  ! Compute the heat/moisture fluxes
    if ( add_fire_heat_flux ) then
     do j = jts,jte
     do i = its,ite
       if ( coef_bb_dc(i,j)*frp_in(i,j) .ge. 1.E7_RKIND ) then
          hfx_bb(i,j)           = min(max(0._RKIND,0.88_RKIND * coef_bb_dc(i,j)*frp_in(i,j) / &
                                  0.55_RKIND / area(i,j)) ,5000._RKIND) ! W m-2 [0 - 10,000]
          frac_grid_burned(i,j) = min(max(0._RKIND, 1.3_RKIND * 0.0006_RKIND * &
                                  coef_bb_dc(i,j)*frp_in(i,j)/area(i,j) ), &
                                  1._RKIND)
       else
          hfx_bb(i,j)           = 0._RKIND
          frac_grid_burned(i,j) = 0._RKIND
       endif
     enddo
     enddo
    endif
   ! JLS, input emissions or scale?
    if (add_fire_moist_flux) then
      do j = jts,jte
      do i = its,ite
        if ( coef_bb_dc(i,j)*frp_in(i,j) .ge. 1.E7_RKIND ) then
           qfx_bb(i,j)           = 0._RKIND
        else
           qfx_bb(i,j)           = 0._RKIND
        endif
      enddo
      enddo
    endif

    if ( ebb_dcycle .eq. 2 ) then
       call  diurnal_cycle (  dt,dz8w,rho_phy,pi,ebb_min,              &
                              chem,num_chem,julday,gmt,xlat,xlong,     &
                              fire_end_hr,peak_hr,curr_secs,coef_bb_dc,&
                              fire_hist,hwp,hwp_avg,hwp_day_avg,       &  !I think fire_hist replaced sc_factor
                              vegfrac, eco_id, nblocks,                &
                              lu_nofire, lu_qfire, lu_sfire,           &
                              swdown,ebb_dcycle,ebu,num_e_bb_in,       &
                              index_e_bb_in_smoke_fine,                &
                              fire_type, qv, add_fire_moist_flux,      &
                              bb_qv_scale_factor, hwp_alpha,           &
                              ids,ide, jds,jde, kds,kde,               &
                              ims,ime, jms,jme, kms,kme,               &
                              its,ite, jts,jte, kts,kte                )
    endif

    ! Apply the diurnal cycle coefficient to frp_out ()
    do j=jts,jte
    do i=its,ite
      if ( fire_type(i,j) .eq. 4 ) then ! only apply scaling factor to wildfires
         frp_out(i,j) = min(bb_emis_scale_factor*frp_in(i,j)*coef_bb_dc(i,j),frp_max)
      else
         frp_out(i,j) = min(frp_in(i,j)*coef_bb_dc(i,j),frp_max)
      endif
    enddo
    enddo

    ! Deterimine if this is a plumerise timestep
    do_plumerise = .false.
    if (plumerise_opt .gt. 0 ) do_plumerise = .true.
    ! plumerise frequency in minutes set up by the namelist input
    call_plume       = (do_plumerise .and. (plumerisefire_frq > 0))
    if (call_plume) call_plume = (mod(int(curr_secs), max(1, 60*plumerisefire_frq)) == 0) .or. (ktau == 2)

    ! compute wild-fire plumes
    if  (do_timing) call mpas_timer_start('ebu_driver')
    if (call_plume) then
      call mpas_log_write( ' Calling ebu_driver')
      call ebu_driver (                                               &
                 flam_frac,kfire,                                     &
                 theta_phy,qv,                                        &
                 rho_phy,vvel,u_phy,v_phy,pi_phy,wind_phy,            &
                 z_at_w,zmid,g,cp,rd,                                 &
                 frp_out, min_bb_plume, max_bb_plume,                 &
                 plume_wind_eff,                                      &
                 do_plumerise,                                        &
                 kpbl_thetav,kpbl,curr_secs,                          &
                 xlat, xlong, uspdavg2d, hpbl2d, plume_alpha,         &
                 frp_min, frp_wthreshold,                             &
                 zpbl_threshold, uspd_threshold,                      &
                 e_bb_in, ebu, num_e_bb_in,                           &
                 ids,ide, jds,jde, kds,kde,                           &
                 ims,ime, jms,jme, kms,kme,                           &
                 its,ite, jts,jte, kts,kte, errmsg, errflg            )
      if(errflg/=0) return
    end if
    if  (do_timing) call mpas_timer_stop('ebu_driver')

    ! -- add biomass burning emissions at every timestep
    if (addsmoke_flag == 1) then
    if  (do_timing) call mpas_timer_start('add_emiss_burn')
     call mpas_log_write( ' Calling add_emis_burn')
     call add_emis_burn(dt,dz8w,rho_phy,pi,ebb_min,                   &
                        chem,num_chem,                                &
                        julday,gmt,xlat,xlong,                        &
                        fire_end_hr, peak_hr,curr_secs,               &
                        coef_bb_dc,fire_hist,hwp,hwp_day_avg,         &
                        swdown,ebb_dcycle,ebu,num_e_bb_in,            &
                        fire_type,                                    &
                        qv, add_fire_moist_flux,                      &
                        bb_emis_scale_factor, aod3d_smoke,            &
                        index_e_bb_in_smoke_ultrafine,                &
                        index_e_bb_in_smoke_fine,                     &
                        index_e_bb_in_smoke_coarse,                   &
                        index_e_bb_in_co, index_e_bb_in_nh3,          &
                        index_e_bb_in_ch4,                            &
                        index_e_bb_in_nox, index_e_bb_in_so2,         &
                        ids,ide, jds,jde, kds,kde,                    &
                        ims,ime, jms,jme, kms,kme,                    &
                        its,ite, jts,jte, kts,kte                     )
    endif
    if  (do_timing) call mpas_timer_stop('add_emiss_burn')
  endif ! if do_mpas_smoke

    ! -- add pollen emissions at every timestep
    if ( do_mpas_pollen ) then
    if  (do_timing) call mpas_timer_start('pollen_driver')
    call mpas_log_write( ' Calling pollen driver')
    call  pollen_driver       (                                       &
       num_chem, chem,                                                &
       dt, u10, v10, rho_phy, dz8w, t_phy, z_at_w, ktop_deep,         &
       xland, raincv, rainncv, relhum, swdown, total_flashrate,       &
       cldfrac,                                                       &
       num_pols_per_polp,pollen_emis_scale_factor,                    &
       tree_pollen_emis_scale_factor,                                 &
       grass_pollen_emis_scale_factor,                                &
       weed_pollen_emis_scale_factor,                                 &
       e_bio_in, e_bio_out, kbio,                                     &
       num_e_bio_in, num_e_bio_out,                                   &
       index_e_bio_in_polp_tree,                                      &
       index_e_bio_in_polp_grass,                                     &
       index_e_bio_in_polp_weed,                                      &
       index_e_bio_out_polp_tree,                                     &
       index_e_bio_out_polp_grass,                                    &
       index_e_bio_out_polp_weed,                                     &
       ids,ide, jds,jde, kds,kde,                                     &
       ims,ime, jms,jme, kms,kme,                                     &
       its,ite, jts,jte, kts,kte                                      )
    if  (do_timing) call mpas_timer_stop('pollen_driver')
    endif

  ! -- add sea salt emissions
  if (do_mpas_ssalt .or. do_mpas_hab) then
    if  (do_timing) call mpas_timer_start('seasalt_driver')
     call gocart_seasalt_driver (                                     &
             dt,rri,t_phy,u_phy,v_phy,                                &
             num_chem,chem,rho_phy,dz8w,u10,v10,                      &
             ust,p8w,tskin,xland,xlat,xlong,area,g,                   &
             e_ss_out,num_e_ss_out,                                   &
             index_e_ss_out_ssalt_fine,                               &
             index_e_ss_out_ssalt_coarse,pi,                          &
             num_emis_seas,seas_opt,                                  &
             ids,ide, jds,jde, kds,kde,                               &
             ims,ime, jms,jme, kms,kme,                               &
             its,ite, jts,jte, kts,kte                                )
    if  (do_timing) call mpas_timer_stop('seasalt_driver')
    endif

    if ( do_mpas_dust ) then
    if  (do_timing) call mpas_timer_start('dust_driver')
    call mpas_log_write( ' Calling dust driver')
!    if ( dust_opt .eq. 5 ) then
    !-- compute dust (FENGSHA)
       call gocart_dust_fengsha_driver(dt,ktau,chem,rho_phy,               &
            smois,tslb,p8w,                                           &
            isltyp,snowh,xland,area,g,                                &
            ust,znt,                                                  &
            clayfrac_in,sandfrac_in,                                  &
            uthres_in, uthres_sg_in,                                  &
            albedo_drag_in, feff_in, sep_in,                          &
            e_dust_out, num_e_dust_out,                               &
            index_e_dust_out_dust_fine,                               &
            index_e_dust_out_dust_coarse,                             &
            num_emis_dust,num_chem,nsoil,                             &
            dust_alpha, dust_gamma, dust_drylimit_factor,             &
            dust_moist_correction,                                    &
            ids,ide, jds,jde, kds,kde,                                &
            ims,ime, jms,jme, kms,kme,                                &
            its,ite, jts,jte, kts,kte                                 )
!    end if
    if  (do_timing) call mpas_timer_stop('dust_driver')
    end if

    if ( do_mpas_anthro ) then
    if  (do_timing) call mpas_timer_start('anthro_driver')
       call mpas_log_write( ' Calling anthro emis driver')
       call mpas_smoke_anthro_emis_driver(dt,gmt,julday,kanthro,      &
            xlat,xlong, chem,num_chem,dz8w,t_phy,rho_phy,             &
            e_ant_in, e_ant_out, num_e_ant_in, num_e_ant_out,         &
            index_e_ant_in_unspc_ultrafine,                           &
            index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
            index_e_ant_in_no3_a_fine, index_e_ant_in_so4_a_fine,     &
            index_e_ant_in_nh4_a_fine,                                &
            index_e_ant_in_so2, index_e_ant_in_nh3,index_e_ant_in_ch4,&
            index_e_ant_in_nox, index_e_ant_in_co,                    &
            index_e_ant_out_unspc_ultrafine,                          &
            index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
            index_e_ant_out_no3_a_fine, index_e_ant_out_so4_a_fine,   &
            index_e_ant_out_nh4_a_fine,                               &
            index_e_ant_out_so2, index_e_ant_out_nh3,                 &
            index_e_ant_out_ch4, index_e_ant_out_nox,                 &
            index_e_ant_out_co,                                       &
            ids,ide, jds,jde, kds,kde,                                &
            ims,ime, jms,jme, kms,kme,                                &
            its,ite, jts,jte, kts,kte                                 )
    if  (do_timing) call mpas_timer_stop('anthro_driver')
    endif

    if ( do_mpas_rwc ) then
       call mpas_log_write( ' Calling online residential wood combustion driver')
       call mpas_smoke_rwc_emis_driver(ktau,dt,gmt,julday,krwc,            &
            xlat,xlong, xland, chem,num_chem,dz8w,t_phy,rho_phy,      &
            rwc_emis_scale_factor,                                    &
            RWC_denominator, RWC_annual_sum,         &
            RWC_annual_sum_smoke_fine, RWC_annual_sum_smoke_coarse,   &
            RWC_annual_sum_unspc_fine, RWC_annual_sum_unspc_coarse,   &
            e_ant_out, num_e_ant_out,         &
            index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
            index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse,   &
            index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
            index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse, &
            ids,ide, jds,jde, kds,kde,                                &
            ims,ime, jms,jme, kms,kme,                                &
            its,ite, jts,jte, kts,kte                                 )
    endif

    !>-- compute dry deposition, based on Emerson et al., (2020)
    if (drydep_opt == 1) then
    if  (do_timing) call mpas_timer_start('drydep_driver')
     ! Set up the arrays if this is the first time through
     call mpas_log_write( ' Calling dry_dep_driver_emerson')
     call dry_dep_driver_emerson(rmol,ust,znt,num_chem,ddvel,         &
        chem,dz8w,snowh,t_phy,p_phy,rho_phy,ivgtyp,g,dt,              &
        drydep_flux,tend_chem_settle,dbg_opt,                         &
        pm_settling,vg,                                               &
        ids,ide, jds,jde, kds,kde,                                    &
        ims,ime, jms,jme, kms,kme,                                    &
        its,ite, jts,jte, kts,kte, curr_secs                          )
    if  (do_timing) call mpas_timer_stop('drydep_driver')
    if  (do_timing) call mpas_timer_start('settling')
    if (pm_settling .gt. 0 ) then
     call particle_settling_wrapper(tend_chem_settle,chem,rho_phy,dz8w_flip,vg,   &
                 dt,kts,kte,its,ite,jts,jte,num_chem,ims,ime, jms,jme, kms,kme    )
    endif
    if  (do_timing) call mpas_timer_stop('settling')
    !>-- compute dry deposition based on simple parameterization (HRRR-Smoke)
    elseif (drydep_opt == 2) then
     call dry_dep_driver_simple(rmol,ust,ndvel,ddvel,                 &
        ids,ide, jds,jde, kds,kde,                                    &
        ims,ime, jms,jme, kms,kme,                                    &
        its,ite, jts,jte, kts,kte                                     )
    else
        call mpas_log_write( 'Not performing dry deposition')
    endif

 !>- large-scale wet deposition
    if (wetdep_ls_opt == 1) then
    if  (do_timing) call mpas_timer_start('wetdep_ls')
       call mpas_log_write( ' Calling wetdep_ls')
       call  wetdep_ls(dt,g,chem,rainncv,qv,                          &
                     rho_phy,num_chem,dz8w,vvel,p_phy,                &
                     wetdep_ls_alpha,                                 &
                     wetdep_resolved,                                 &
                     ids,ide, jds,jde, kds,kde,                       &
                     ims,ime, jms,jme, kms,kme,                       &
                     its,ite, jts,jte, kts,kte                        )
    if  (do_timing) call mpas_timer_stop('wetdep_ls')
    endif

    !>-- output of MPAS-Smoke
    call mpas_log_write( ' Calculating AOD ')
    call mpas_aod_diag(           chem,aod3d,rho_phy,dz8w,num_chem,        &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte         )
    call mpas_log_write( ' Calculating VIS ')
    call mpas_visibility_diag(    qc_vis,qr_vis,qi_vis,qs_vis,qg_vis,    &
                                  blcldw_vis,blcldi_vis,                 &
                                  rho_phy,wind10m,wind_phy,              &
                                  rh2m,relhum,qv, &
                                  t2m,t_phy, &
                                  coszen,aod3d,vis,                 &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte         )
    

    if (do_mpas_smoke) then
    ! UPP/MPASSIT expects FRP in MW
    do j=jts,jte
    do i=its,ite
       if ( fire_type(i,j) .eq. 4 ) then ! only apply scaling factor to wildfires
          frp_out(i,j) = 1.e-6_RKIND * min(bb_emis_scale_factor*frp_in(i,j)*coef_bb_dc(i,j),frp_max)
       else
          frp_out(i,j) = 1.e-6_RKIND * min(frp_in(i,j)*coef_bb_dc(i,j),frp_max)
       endif
    enddo
    enddo
    ! Assign the coef to the emissions on the way out
    do nv=1,num_e_bb_in
    do j=jts,jte
    do k=kts,kte
    do i=its,ite
          e_bb_out(i,k,j,nv)=ebu(i,k,j,nv) * coef_bb_dc(i,j)
    enddo
    enddo
    enddo
    enddo
    endif
    
 end subroutine mpas_smoke_driver

 subroutine mpas_smoke_prep(                                                &
        do_mpas_smoke,                                                      &
        ktau, nlcat,cp,ebb_dcycle,ebb_min,                                  &
        xland,xlat,xlong,ivgtyp,isltyp,vegfrac,                             &
        snowh,u10,v10,wind10m,t2m,dpt2m,wetness,hwp,hwp_day_avg,            &
        index_e_bb_in_smoke_fine,num_e_bb_in,kfire,e_bb_in,                 &
        t_phy,u_phy,v_phy,p_phy,pi_phy,z_at_w,                              &
        dz8w,dz8w_flip,                                                     &
        rho_phy,qv,relhum,rh2m,rri,                                         &
        total_flashrate,                                                    &
        wind_phy,theta_phy,zmid,kpbl_thetav,                                &
        peak_hr,coef_bb_dc,fire_hist,                                       &
        lu_nofire, lu_qfire, lu_sfire, fire_type,                           &
        ids,ide, jds,jde, kds,kde,                                          &
        ims,ime, jms,jme, kms,kme,                                          &
        its,ite, jts,jte, kts,kte                                           )

    !intent arguments:
    integer,intent(in):: ids,ide,jds,jde,kds,kde,                           &
                         ims,ime,jms,jme,kms,kme,                           &
                         its,ite,jts,jte,kts,kte

    integer,intent(in)     :: ktau, nlcat,ebb_dcycle,                       &
                              kfire, num_e_bb_in, index_e_bb_in_smoke_fine
    logical,intent(in)     :: do_mpas_smoke
    real(RKIND),intent(in) :: cp, ebb_min
    integer,intent(in), dimension(ims:ime, jms:jme) :: isltyp, ivgtyp
    real(RKIND),intent(in),   dimension(ims:ime, nlcat, jms:jme) :: vegfrac

    real(RKIND),intent(in),   dimension(ims:ime, jms:jme) :: xland, xlat, xlong,                   &
                                        snowh, u10, v10, t2m, dpt2m, wetness
    real(RKIND),intent(in),   dimension(ims:ime, kms:kme, jms:jme) :: qv, z_at_w,                  &
                                        p_phy, t_phy, u_phy, v_phy, pi_phy, rho_phy, dz8w
    real(RKIND),intent(out),  dimension(ims:ime, kms:kme, jms:jme) :: zmid, wind_phy,theta_phy,   &
                                       relhum, rri, dz8w_flip
    integer    ,intent(out),  dimension(ims:ime, jms:jme) :: fire_type, kpbl_thetav
    real(RKIND),intent(out),  dimension(ims:ime, jms:jme) :: lu_nofire, lu_qfire, lu_sfire
    real(RKIND),intent(out),  dimension(ims:ime, jms:jme) :: fire_hist, peak_hr, hwp_day_avg
    real(RKIND),intent(inout),dimension(ims:ime,jms:jme),optional :: hwp, coef_bb_dc
    real(RKIND),intent(out),  dimension(ims:ime, jms:jme) :: total_flashrate, wind10m, rh2m
    real(RKIND),intent(in),   dimension(ims:ime,1:kfire,jms:jme,1:num_e_bb_in),optional :: e_bb_in

    !local variables
    real(RKIND), parameter :: delta_theta4gust = 0.5
    real(RKIND) :: theta, wdgust, snoweq, term2, term3
!    real(RKIND), dimension(ims:ime, jms:jme) :: 
    real(RKIND), dimension(ims:ime, kms:kme, jms:jme) :: thetav

    integer :: i, j, k, k1, nv

    if ( do_mpas_smoke ) then
       if (ktau==1) then
         do j=jts,jte
         do i=its,ite
           fire_hist   (i,j) = 1._RKIND
           coef_bb_dc  (i,j) = 1._RKIND
           if (xlong(i,j)<230.) then
               peak_hr(i,j)= 0.0* 3600.     ! peak at 24 UTC, fires in Alaska
           elseif(xlong(i,j)<245.) then
               peak_hr(i,j)= 23.0* 3600.
           elseif (xlong(i,j)<260.) then
               peak_hr(i,j)= 22.0* 3600.    ! peak at 22 UTC, fires in the western US
           elseif (xlong(i,j)<275.) then
               peak_hr(i,j)= 21.0* 3600.
           elseif (xlong(i,j)<290.) then    ! peak at 20 UTC, fires in the eastern US
               peak_hr(i,j)= 20.0* 3600.
           else
               peak_hr(i,j)= 19.0* 3600.
           endif
         enddo
         enddo
       endif
    endif

    do j=jts,jte
    do k=kts,kte
    do i=its,ite
      zmid(i,k,j)=z_at_w(i,k,j)
    enddo
    enddo
    enddo

    do j=jts,jte
    do k=kts,kte
    do i=its,ite
      theta_phy(i,k,j) = t_phy(i,k,j)/pi_phy(i,k,j)*cp
      wind_phy(i,k,j) = sqrt(u_phy(i,k,j)**2. + v_phy(i,k,j)**2.)
      rri(i,k,j) = 1._RKIND/rho_phy(i,k,j)
    enddo
    enddo
    enddo

    ! Calculate relative humidity [0.1 -- 0.95]
    do j=jts,jte
    do k=kts,kte
    do i=its,ite
         relhum(i,k,j) = max(.1_RKIND,MIN( .95_RKIND, qv(i,k,j) / &
               (3.80_RKIND*exp(17.27_RKIND*(t_phy(i,k,j)-273._RKIND)/ &
               (t_phy(i,k,j)-36._RKIND))/(.01_RKIND*p_phy(i,k,j)))))
    enddo
    enddo
    enddo

    ! JLS, TODO, lightning flashrate
    do j=jts,jte
    do i=its,ite
       total_flashrate(i,j) = 0._RKIND
       rh2m(i,j) = relhum(i,1,j)
! TODO, JLS, check
!       rh2m(i,j) = EXP((17.625_RKIND*(dpt2m(i,j)+273.15_RKIND))/(243.04_RKIND+(dpt2m(i,j)+273.15_RKIND))) / &
!                   EXP((17.625_RKIND*(t2m(i,j)+273.15_RKIND))/(243.04_RKIND+(t2m(i,j)+273.15_RKIND)))
    enddo
    enddo
    

    !---- Calculate PBLH and K-PBL based on virtual potential temperature profile
    !---- First calculate THETAV
    do j = jts,jte
    do k = kts,kte
    do i = its,ite
       theta = t_phy(i,k,j) * (1.E5_RKIND/p_phy(i,k,j))**0.286_RKIND
       thetav(i,k,j) = theta * (1._RKIND + 0.61_RKIND * qv(i,k,j))
       dz8w_flip(i,k,j)   = dz8w(i,kte-k+kts,j)
    enddo
    enddo
    enddo
    !---- Now use the UPP code to deterimine the height and level
    do i = its, ite
    do j = jts, jte
       if ( thetav(i,kts+1,j) .lt. ( thetav(i,kts,j) + delta_theta4gust) ) then
          do k = kts+1, kte
             k1 = k
!--- give theta-v at the sfc a 0.5K boost in the PBLH definition
             if ( thetav(i,kts+k-1,j) .gt. ( thetav(i,kts,j) + delta_theta4gust) ) then
                exit
             endif
          enddo
          kpbl_thetav(i,j) = k1
       else
          kpbl_thetav(i,j) = kts + 1
       endif
   enddo
   enddo

    if (do_mpas_smoke) then
      !RAR: change this to the fractional LU type; fire_type: 0- no fires, 1- Ag
      ! or urban fires, 2- prescribed fires in wooded area, 3- wildfires
       if (ebb_dcycle==2) then
         do j=jts,jte
         do i=its,ite
           if (e_bb_in(i,1,j,index_e_bb_in_smoke_fine)<ebb_min) then
              fire_type(i,j) = 0
              lu_nofire(i,j) = 1.0
           else
             ! Permanent wetlands, snow/ice, water, barren tundra:
             lu_nofire(i,j)= vegfrac(i,11,j) + vegfrac(i,15,j) + vegfrac(i,17,j) + vegfrac(i,20,j)
             ! cropland, urban, cropland/natural mosaic, barren and sparsely
             ! vegetated and non-vegetation areas:
             lu_qfire(i,j) = lu_nofire(i,j) + vegfrac(i,12,j) + vegfrac(i,13,j) + vegfrac(i,14,j) + vegfrac(i,16,j)
             ! Savannas and grassland fires, these fires last longer than the Ag fires:
             lu_sfire(i,j) = lu_qfire(i,j) + vegfrac(i,8,j) + vegfrac(i,9,j) + vegfrac(i,10,j)
             if (lu_nofire(i,j)>0.95) then ! no fires
               fire_type(i,j) = 0
             else if (lu_qfire(i,j)>0.9) then   ! Ag. and urban fires
               fire_type(i,j) = 1
             else if (xlong(i,j)>260. .AND. xlat(i,j)>25. .AND. xlat(i,j)<41.) then
               fire_type(i,j) = 2    ! slash burn and wildfires in the east, eastern temperate forest ecosystem
             else if (lu_sfire(i,j)>0.8) then
               fire_type(i,j) = 3    ! savanna and grassland fires
             else
               fire_type(i,j) = 4    ! potential wildfires
             end if
           end if
         end do
         end do
       endif ! ebb_dycycle == 2
   

! TODO - JLS add HWP options 
       !>-- HWP: Pre-release of RRFSv1 method - using wind gust calculated via UPP Method
       do i=its, ite
       do j=jts, jte
         wind10m(i,j) = sqrt(u10(i,j)**2.+v10(i,j)**2.)
         wdgust  =max(wind10m(i,j),3._RKIND)
         snoweq  =max((25._RKIND - snowh(i,j))/25._RKIND,0._RKIND)
!         term2   = max(t2m(i,j)-dpt2m(i,j),15._RKIND)**1.03 ! TODO, floating point exception
         term2   = 15._RKIND
         term3   = (1._RKIND-wetness(i,j))**0.4
         hwp(i,j)= 0.177_RKIND * wdgust**0.97 * term2 * ((1._RKIND-wetness(i,j))**0.4) * snoweq
         hwp_day_avg(i,j)=hwp(i,j)
       enddo
       enddo
    endif ! do_mpas_smoke
   

  end subroutine mpas_smoke_prep
  

!> @}
  end module mpas_smoke_wrapper
