!>\file  mpas_smoke_init.F90
!! This file initializes and saves the chemistry pointers so they are available to lower routines
!
! Jordan.Schnell@noaa.gov  
! 01/2025
!
module mpas_smoke_init

  use mpas_kind_types
  
  implicit none

   integer, save ::     &
   p_smoke_fine   = -1, &
   p_smoke_coarse = -1, &
   p_dust_fine    = -1, &
   p_dust_coarse  = -1, &
   p_polp_all     = -1, &
   p_polp_tree    = -1, &
   p_polp_grass   = -1, & 
   p_polp_weed    = -1, &
   p_pols_tree    = -1, &
   p_pols_grass   = -1, & 
   p_pols_weed    = -1, &
   p_pols_all     = -1, &
   p_unspc_fine   = -1, &
   p_unspc_coarse = -1, &
   p_ssalt_fine   = -1, &
   p_ssalt_coarse = -1, &
   p_no3_a_fine   = -1, &
   p_so4_a_fine   = -1, &
   p_nh4_a_fine   = -1, &
   p_so2          = -1, &
   p_nh3          = -1, &
   p_ch4          = -1, &
   p_mold_fine    = -1, &
   p_bact_fine    = -1, &
   p_plastic_fine = -1, &
   p_dms          = -1

   contains

  subroutine set_scalar_indices( chemistry_start,                                     &
                                 index_smoke_fine, index_smoke_coarse,                &
                                 index_dust_fine, index_dust_coarse,                  &                  
                                 index_polp_tree, index_polp_grass, index_polp_weed,  &
                                 index_pols_tree, index_pols_grass, index_pols_weed,  &
                                 index_pols_all,  index_polp_all,                     &
                                 index_unspc_fine, index_unspc_coarse,                &
                                 index_ssalt_fine, index_ssalt_coarse,                &
                                 index_no3_a_fine, index_so4_a_fine, index_nh4_a_fine,&
                                 index_so2, index_nh3, index_ch4                      )

    implicit none

    integer, intent(in), optional :: chemistry_start,                           &
                           index_smoke_fine, index_smoke_coarse,                &
                           index_dust_fine, index_dust_coarse,                  &
                           index_polp_tree, index_polp_grass, index_polp_weed,  &
                           index_pols_tree, index_pols_grass, index_pols_weed,  &
                           index_pols_all,  index_polp_all,                     & 
                           index_unspc_fine, index_unspc_coarse,                &
                           index_ssalt_fine, index_ssalt_coarse,                &
                           index_no3_a_fine, index_so4_a_fine, index_nh4_a_fine,&
                           index_so2, index_nh3, index_ch4

        if(present(index_smoke_fine))   p_smoke_fine    = index_smoke_fine - chemistry_start + 1
        if(present(index_smoke_coarse)) p_smoke_coarse  = index_smoke_coarse - chemistry_start + 1
        if(present(index_dust_fine))    p_dust_fine     = index_dust_fine - chemistry_start + 1
        if(present(index_dust_coarse))  p_dust_coarse   = index_dust_coarse - chemistry_start + 1
        if(present(index_polp_all))     p_polp_all      = index_polp_all - chemistry_start + 1
        if(present(index_polp_tree))    p_polp_tree     = index_polp_tree - chemistry_start + 1
        if(present(index_polp_grass))   p_polp_grass    = index_polp_grass - chemistry_start + 1
        if(present(index_polp_weed))    p_polp_weed     = index_polp_weed - chemistry_start + 1
        if(present(index_pols_tree))    p_pols_tree     = index_pols_tree - chemistry_start + 1
        if(present(index_pols_grass))   p_pols_grass    = index_pols_grass - chemistry_start + 1
        if(present(index_pols_weed))    p_pols_weed     = index_pols_weed - chemistry_start + 1
        if(present(index_pols_all))     p_pols_all      = index_pols_all - chemistry_start + 1
        if(present(index_unspc_fine))   p_unspc_fine    = index_unspc_fine - chemistry_start + 1
        if(present(index_unspc_coarse)) p_unspc_coarse  = index_unspc_coarse - chemistry_start + 1
        if(present(index_ssalt_fine))   p_ssalt_fine    = index_ssalt_fine - chemistry_start + 1
        if(present(index_ssalt_coarse)) p_ssalt_coarse  = index_ssalt_coarse - chemistry_start + 1
        if(present(index_no3_a_fine))   p_no3_a_fine    = index_no3_a_fine - chemistry_start + 1
        if(present(index_so4_a_fine))   p_so4_a_fine    = index_so4_a_fine - chemistry_start + 1
        if(present(index_nh4_a_fine))   p_nh4_a_fine    = index_nh4_a_fine - chemistry_start + 1
        if(present(index_so2))          p_so2           = index_so2 - chemistry_start + 1
        if(present(index_nh3))          p_nh3           = index_nh3 - chemistry_start + 1
        if(present(index_ch4))          p_ch4           = index_ch4 - chemistry_start + 1
 
   end subroutine set_scalar_indices


end module mpas_smoke_init
