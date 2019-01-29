MODULE dtadyn
   !!======================================================================
   !!                       ***  MODULE  dtadyn  ***
   !! Off-line : interpolation of the physical fields
   !!======================================================================
   !! History :   OPA  ! 1992-01 (M. Imbard) Original code
   !!             8.0  ! 1998-04 (L.Bopp MA Foujols) slopes for isopyc. 
   !!              -   ! 1998-05 (L. Bopp) read output of coupled run
   !!             8.2  ! 2001-01 (M. Levy et M. Benjelloul) add netcdf FORMAT
   !!   NEMO      1.0  ! 2005-03 (O. Aumont and A. El Moussaoui) F90
   !!              -   ! 2005-12 (C. Ethe) Adapted for DEGINT
   !!             3.0  ! 2007-06 (C. Ethe) use of iom module
   !!             3.3  ! 2010-11 (C. Ethe) Full reorganization of the off-line: phasing with the on-line
   !!             3.4  ! 2011-05 (C. Ethe) Use of fldread
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dta_dyn_init : initialization, namelist read, and SAVEs control
   !!   dta_dyn      : Interpolation of the fields
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE c1d             ! 1D configuration: lk_c1d
   USE dom_oce         ! ocean domain: variables
   USE zdf_oce         ! ocean vertical physics: variables
   USE sbc_oce         ! surface module: variables
   USE trc_oce         ! share ocean/biogeo variables
   USE phycst          ! physical constants
   USE trabbl          ! active tracer: bottom boundary layer
   USE ldfslp          ! lateral diffusion: iso-neutral slopes
   USE ldfeiv          ! eddy induced velocity coef. 
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE zdfmxl          ! vertical physics: mixed layer depth
   USE eosbn2          ! equation of state - Brunt Vaisala frequency
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE zpshde          ! z-coord. with partial steps: horizontal derivatives
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! distributed memory computing library
   USE prtctl          ! print control
   USE fldread         ! read input fields 
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dta_dyn_init   ! called by opa.F90
   PUBLIC   dta_dyn        ! called by step.F90

   CHARACTER(len=100) ::   cn_dir     = './'    !: Root directory for location of ssr files
   LOGICAL            ::   ln_dynwzv  = .true.  !: vertical velocity read in a file (T) or computed from u/v (F)
   LOGICAL            ::   ln_dynbbl  = .true.  !: bbl coef read in a file (T) or computed (F)
   LOGICAL            ::   ln_degrad  = .false. !: degradation option enabled or not 

   INTEGER  , PARAMETER ::   jpfld = 19     ! maximum number of files to read
   INTEGER  , SAVE      ::   jf_tem         ! index of temperature
   INTEGER  , SAVE      ::   jf_sal         ! index of salinity
   INTEGER  , SAVE      ::   jf_uwd         ! index of u-wind
   INTEGER  , SAVE      ::   jf_vwd         ! index of v-wind
   INTEGER  , SAVE      ::   jf_wwd         ! index of w-wind
   INTEGER  , SAVE      ::   jf_avt         ! index of Kz
   INTEGER  , SAVE      ::   jf_mld         ! index of mixed layer deptht
   INTEGER  , SAVE      ::   jf_emp         ! index of water flux
   INTEGER  , SAVE      ::   jf_qsr         ! index of solar radiation
   INTEGER  , SAVE      ::   jf_wnd         ! index of wind speed
   INTEGER  , SAVE      ::   jf_ice         ! index of sea ice cover
   INTEGER  , SAVE      ::   jf_ubl         ! index of u-bbl coef
   INTEGER  , SAVE      ::   jf_vbl         ! index of v-bbl coef
   INTEGER  , SAVE      ::   jf_ahu         ! index of u-diffusivity coef
   INTEGER  , SAVE      ::   jf_ahv         ! index of v-diffusivity coef 
   INTEGER  , SAVE      ::   jf_ahw         ! index of w-diffusivity coef
   INTEGER  , SAVE      ::   jf_eiu         ! index of u-eiv
   INTEGER  , SAVE      ::   jf_eiv         ! index of v-eiv
   INTEGER  , SAVE      ::   jf_eiw         ! index of w-eiv

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_dyn  ! structure of input fields (file informations, fields read)
   !                                               ! 
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: wdta       ! vertical velocity at 2 time step
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:  ) :: wnow       ! vertical velocity at 2 time step
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: uslpdta    ! zonal isopycnal slopes
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: vslpdta    ! meridional isopycnal slopes
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: wslpidta   ! zonal diapycnal slopes
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: wslpjdta   ! meridional diapycnal slopes
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: uslpnow    ! zonal isopycnal slopes
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: vslpnow    ! meridional isopycnal slopes
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: wslpinow   ! zonal diapycnal slopes
   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: wslpjnow   ! meridional diapycnal slopes

   INTEGER :: nrecprev_tem , nrecprev_uwd

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OFF 3.3 , NEMO Consortium (2010)
   !! $Id: dtadyn.F90 3270 2012-01-20 11:44:56Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dta_dyn( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dta_dyn  ***
      !!
      !! ** Purpose :  Prepares dynamics and physics fields from a NEMO run
      !!               for an off-line simulation of passive tracers
      !!
      !! ** Method : calculates the position of data 
      !!             - computes slopes if needed
      !!             - interpolates data if needed
      !!----------------------------------------------------------------------
      !
      USE oce, ONLY:  zts    => tsa 
      USE oce, ONLY:  zuslp  => ua   , zvslp  => va
      USE oce, ONLY:  zwslpi => rotb , zwslpj => rotn
      USE oce, ONLY:  zu     => ub   , zv     => vb,  zw => hdivb
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      INTEGER  ::   isecsbc    ! number of seconds between Jan. 1st 00h of nit000 year and the middle of time step
      REAL(wp) ::   ztinta     ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb     ! ratio applied to before records when doing time interpolation
      INTEGER  ::   iswap_tem, iswap_uwd    ! 
      !!----------------------------------------------------------------------
      
      !
      IF( nn_timing == 1 )  CALL timing_start( 'dta_dyn')
      !
      isecsbc = nsec_year + nsec1jan000 
      !
      IF( kt == nit000 ) THEN
         nrecprev_tem = 0
         nrecprev_uwd = 0
         !
         CALL fld_read( kt, 1, sf_dyn )      !==   read data at kt time step   ==!
         !
         IF( lk_ldfslp .AND. .NOT.lk_c1d .AND. sf_dyn(jf_tem)%ln_tint ) THEN    ! Computes slopes (here avt is used as workspace)                       
            zts(:,:,:,jp_tem) = sf_dyn(jf_tem)%fdta(:,:,:,1) * tmask(:,:,:)   ! temperature
            zts(:,:,:,jp_sal) = sf_dyn(jf_sal)%fdta(:,:,:,1) * tmask(:,:,:)   ! salinity 
            avt(:,:,:)        = sf_dyn(jf_avt)%fdta(:,:,:,1) * tmask(:,:,:)   ! vertical diffusive coef.
            CALL dta_dyn_slp( kt, zts, zuslp, zvslp, zwslpi, zwslpj )
            uslpdta (:,:,:,1) = zuslp (:,:,:) 
            vslpdta (:,:,:,1) = zvslp (:,:,:) 
            wslpidta(:,:,:,1) = zwslpi(:,:,:) 
            wslpjdta(:,:,:,1) = zwslpj(:,:,:) 
         ENDIF
         IF( ln_dynwzv .AND. sf_dyn(jf_uwd)%ln_tint )  THEN    ! compute vertical velocity from u/v
            zu(:,:,:) = sf_dyn(jf_uwd)%fdta(:,:,:,1)
            zv(:,:,:) = sf_dyn(jf_vwd)%fdta(:,:,:,1)
            CALL dta_dyn_wzv( zu, zv, zw )
            wdta(:,:,:,1) = zw(:,:,:) * tmask(:,:,:)
         ENDIF
      ELSE
         nrecprev_tem = sf_dyn(jf_tem)%nrec_a(2)
         nrecprev_uwd = sf_dyn(jf_uwd)%nrec_a(2)
         !
         CALL fld_read( kt, 1, sf_dyn )      !==   read data at kt time step   ==!
         !
      ENDIF
      ! 
      IF( lk_ldfslp .AND. .NOT.lk_c1d ) THEN    ! Computes slopes (here avt is used as workspace)                       
         iswap_tem = 0
         IF(  kt /= nit000 .AND. ( sf_dyn(jf_tem)%nrec_a(2) - nrecprev_tem ) /= 0 )  iswap_tem = 1
         IF( ( isecsbc > sf_dyn(jf_tem)%nrec_b(2) .AND. iswap_tem == 1 ) .OR. kt == nit000 )  THEN    ! read/update the after data
            write(numout,*)
            write(numout,*) ' Compute new slopes at kt = ', kt
            IF( sf_dyn(jf_tem)%ln_tint ) THEN                 ! time interpolation of data
               IF( kt /= nit000 ) THEN
                  uslpdta (:,:,:,1) =  uslpdta (:,:,:,2)         ! swap the data
                  vslpdta (:,:,:,1) =  vslpdta (:,:,:,2)  
                  wslpidta(:,:,:,1) =  wslpidta(:,:,:,2) 
                  wslpjdta(:,:,:,1) =  wslpjdta(:,:,:,2) 
               ENDIF
               !
               zts(:,:,:,jp_tem) = sf_dyn(jf_tem)%fdta(:,:,:,2) * tmask(:,:,:)   ! temperature
               zts(:,:,:,jp_sal) = sf_dyn(jf_sal)%fdta(:,:,:,2) * tmask(:,:,:)   ! salinity 
               avt(:,:,:)        = sf_dyn(jf_avt)%fdta(:,:,:,2) * tmask(:,:,:)   ! vertical diffusive coef.
               CALL dta_dyn_slp( kt, zts, zuslp, zvslp, zwslpi, zwslpj )
               !
               uslpdta (:,:,:,2) = zuslp (:,:,:) 
               vslpdta (:,:,:,2) = zvslp (:,:,:) 
               wslpidta(:,:,:,2) = zwslpi(:,:,:) 
               wslpjdta(:,:,:,2) = zwslpj(:,:,:) 
            ELSE
               zts(:,:,:,jp_tem) = sf_dyn(jf_tem)%fnow(:,:,:) * tmask(:,:,:)
               zts(:,:,:,jp_sal) = sf_dyn(jf_sal)%fnow(:,:,:) * tmask(:,:,:)
               avt(:,:,:)        = sf_dyn(jf_avt)%fnow(:,:,:) * tmask(:,:,:)
               CALL dta_dyn_slp( kt, zts, zuslp, zvslp, zwslpi, zwslpj )
               uslpnow (:,:,:)   = zuslp (:,:,:) 
               vslpnow (:,:,:)   = zvslp (:,:,:) 
               wslpinow(:,:,:)   = zwslpi(:,:,:) 
               wslpjnow(:,:,:)   = zwslpj(:,:,:) 
            ENDIF
         ENDIF
         IF( sf_dyn(jf_tem)%ln_tint )  THEN
            ztinta =  REAL( isecsbc - sf_dyn(jf_tem)%nrec_b(2), wp )  &
               &    / REAL( sf_dyn(jf_tem)%nrec_a(2) - sf_dyn(jf_tem)%nrec_b(2), wp )
            ztintb =  1. - ztinta
            uslp (:,:,:) = ztintb * uslpdta (:,:,:,1)  + ztinta * uslpdta (:,:,:,2)  
            vslp (:,:,:) = ztintb * vslpdta (:,:,:,1)  + ztinta * vslpdta (:,:,:,2)  
            wslpi(:,:,:) = ztintb * wslpidta(:,:,:,1)  + ztinta * wslpidta(:,:,:,2)  
            wslpj(:,:,:) = ztintb * wslpjdta(:,:,:,1)  + ztinta * wslpjdta(:,:,:,2)  
         ELSE
            uslp (:,:,:) = uslpnow (:,:,:)
            vslp (:,:,:) = vslpnow (:,:,:)
            wslpi(:,:,:) = wslpinow(:,:,:)
            wslpj(:,:,:) = wslpjnow(:,:,:)
         ENDIF
      ENDIF
      !
      IF( ln_dynwzv )  THEN    ! compute vertical velocity from u/v
         iswap_uwd = 0
         IF(  kt /= nit000 .AND. ( sf_dyn(jf_uwd)%nrec_a(2) - nrecprev_uwd ) /= 0 )  iswap_uwd = 1
         IF( ( isecsbc > sf_dyn(jf_uwd)%nrec_b(2) .AND. iswap_uwd == 1 ) .OR. kt == nit000 )  THEN    ! read/update the after data
            write(numout,*)
            write(numout,*) ' Compute new vertical velocity at kt = ', kt
            write(numout,*)
            IF( sf_dyn(jf_uwd)%ln_tint ) THEN                 ! time interpolation of data
               IF( kt /= nit000 )  THEN
                  wdta(:,:,:,1) =  wdta(:,:,:,2)     ! swap the data for initialisation
               ENDIF
               zu(:,:,:) = sf_dyn(jf_uwd)%fdta(:,:,:,2)
               zv(:,:,:) = sf_dyn(jf_vwd)%fdta(:,:,:,2)
               CALL dta_dyn_wzv( zu, zv, zw )
               wdta(:,:,:,2) = zw(:,:,:) * tmask(:,:,:)
            ELSE
               zu(:,:,:) = sf_dyn(jf_uwd)%fnow(:,:,:) 
               zv(:,:,:) = sf_dyn(jf_vwd)%fnow(:,:,:)
               CALL dta_dyn_wzv( zu, zv, zw )
               wnow(:,:,:)  = zw(:,:,:) * tmask(:,:,:)
            ENDIF
         ENDIF
         IF( sf_dyn(jf_uwd)%ln_tint )  THEN
            ztinta =  REAL( isecsbc - sf_dyn(jf_uwd)%nrec_b(2), wp )  &
               &    / REAL( sf_dyn(jf_uwd)%nrec_a(2) - sf_dyn(jf_uwd)%nrec_b(2), wp )
            ztintb =  1. - ztinta
            wn(:,:,:) = ztintb * wdta(:,:,:,1)  + ztinta * wdta(:,:,:,2)  
         ELSE
            wn(:,:,:) = wnow(:,:,:)
         ENDIF
      ENDIF
      !
      tsn(:,:,:,jp_tem) = sf_dyn(jf_tem)%fnow(:,:,:) * tmask(:,:,:)    ! temperature
      tsn(:,:,:,jp_sal) = sf_dyn(jf_sal)%fnow(:,:,:) * tmask(:,:,:)    ! salinity
      !
      CALL eos    ( tsn, rhd, rhop )                                       ! In any case, we need rhop
      CALL zdf_mxl( kt )                                                   ! In any case, we need mxl 
      !
      avt(:,:,:)       = sf_dyn(jf_avt)%fnow(:,:,:) * tmask(:,:,:)    ! vertical diffusive coefficient 
      un (:,:,:)       = sf_dyn(jf_uwd)%fnow(:,:,:) * umask(:,:,:)    ! u-velocity
      vn (:,:,:)       = sf_dyn(jf_vwd)%fnow(:,:,:) * vmask(:,:,:)    ! v-velocity 
      IF( .NOT.ln_dynwzv ) &                                           ! w-velocity read in file 
         wn (:,:,:)    = sf_dyn(jf_wwd)%fnow(:,:,:) * tmask(:,:,:)    
      hmld(:,:)        = sf_dyn(jf_mld)%fnow(:,:,1) * tmask(:,:,1)    ! mixed layer depht
      wndm(:,:)        = sf_dyn(jf_wnd)%fnow(:,:,1) * tmask(:,:,1)    ! wind speed - needed for gas exchange
      emp (:,:)        = sf_dyn(jf_emp)%fnow(:,:,1) * tmask(:,:,1)    ! E-P
      emps(:,:)        = emp(:,:) 
      fr_i(:,:)        = sf_dyn(jf_ice)%fnow(:,:,1) * tmask(:,:,1)     ! Sea-ice fraction
      qsr (:,:)        = sf_dyn(jf_qsr)%fnow(:,:,1) * tmask(:,:,1)    ! solar radiation

      !                                                      ! bbl diffusive coef
#if defined key_trabbl && ! defined key_c1d
      IF( ln_dynbbl ) THEN                                        ! read in a file
         ahu_bbl(:,:)  = sf_dyn(jf_ubl)%fnow(:,:,1) * umask(:,:,1)
         ahv_bbl(:,:)  = sf_dyn(jf_vbl)%fnow(:,:,1) * vmask(:,:,1)
      ELSE                                                        ! Compute bbl coefficients if needed
         tsb(:,:,:,:) = tsn(:,:,:,:)
         CALL bbl( kt, nit000, 'TRC')
      END IF
#endif
#if ( ! defined key_degrad && defined key_traldf_c2d && defined key_traldf_eiv ) && ! defined key_c1d 
      aeiw(:,:)        = sf_dyn(jf_eiw)%fnow(:,:,1) * tmask(:,:,1)    ! w-eiv
      !                                                           ! Computes the horizontal values from the vertical value
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            aeiu(ji,jj) = .5 * ( aeiw(ji,jj) + aeiw(ji+1,jj  ) )  ! Average the diffusive coefficient at u- v- points
            aeiv(ji,jj) = .5 * ( aeiw(ji,jj) + aeiw(ji  ,jj+1) )  ! at u- v- points
         END DO
      END DO
      CALL lbc_lnk( aeiu, 'U', 1. )   ;   CALL lbc_lnk( aeiv, 'V', 1. )    ! lateral boundary condition
#endif
      
#if defined key_degrad && ! defined key_c1d 
      !                                          ! degrad option : diffusive and eiv coef are 3D
      ahtu(:,:,:) = sf_dyn(jf_ahu)%fnow(:,:,:) * umask(:,:,:)
      ahtv(:,:,:) = sf_dyn(jf_ahv)%fnow(:,:,:) * vmask(:,:,:)
      ahtw(:,:,:) = sf_dyn(jf_ahw)%fnow(:,:,:) * tmask(:,:,:)
#  if defined key_traldf_eiv 
      aeiu(:,:,:) = sf_dyn(jf_eiu)%fnow(:,:,:) * umask(:,:,:)
      aeiv(:,:,:) = sf_dyn(jf_eiv)%fnow(:,:,:) * vmask(:,:,:)
      aeiw(:,:,:) = sf_dyn(jf_eiw)%fnow(:,:,:) * tmask(:,:,:)
#  endif
#endif
      !
      !
      IF(ln_ctl) THEN                  ! print control
         CALL prt_ctl(tab3d_1=tsn(:,:,:,jp_tem), clinfo1=' tn      - : ', mask1=tmask, ovlap=1, kdim=jpk   )
         CALL prt_ctl(tab3d_1=tsn(:,:,:,jp_sal), clinfo1=' sn      - : ', mask1=tmask, ovlap=1, kdim=jpk   )
         CALL prt_ctl(tab3d_1=un               , clinfo1=' un      - : ', mask1=umask, ovlap=1, kdim=jpk   )
         CALL prt_ctl(tab3d_1=vn               , clinfo1=' vn      - : ', mask1=vmask, ovlap=1, kdim=jpk   )
         CALL prt_ctl(tab3d_1=wn               , clinfo1=' wn      - : ', mask1=tmask, ovlap=1, kdim=jpk   )
         CALL prt_ctl(tab3d_1=avt              , clinfo1=' kz      - : ', mask1=tmask, ovlap=1, kdim=jpk   )
         CALL prt_ctl(tab2d_1=fr_i             , clinfo1=' fr_i    - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab2d_1=hmld             , clinfo1=' hmld    - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab2d_1=emps             , clinfo1=' emps    - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab2d_1=wndm             , clinfo1=' wspd    - : ', mask1=tmask, ovlap=1 )
         CALL prt_ctl(tab2d_1=qsr              , clinfo1=' qsr     - : ', mask1=tmask, ovlap=1 )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'dta_dyn')
      !
   END SUBROUTINE dta_dyn


   SUBROUTINE dta_dyn_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dta_dyn_init  ***
      !!
      !! ** Purpose :   Initialisation of the dynamical data     
      !! ** Method  : - read the data namdta_dyn namelist
      !!
      !! ** Action  : - read parameters
      !!----------------------------------------------------------------------
      INTEGER  :: ierr, ierr0, ierr1, ierr2, ierr3   ! return error code
      INTEGER  :: ifpr                               ! dummy loop indice
      INTEGER  :: jfld                               ! dummy loop arguments
      INTEGER  :: inum, idv, idimv                   ! local integer
      !!
      CHARACTER(len=100)            ::  cn_dir   !   Root directory for location of core files
      TYPE(FLD_N), DIMENSION(jpfld) ::  slf_d    ! array of namelist informations on the fields to read
      TYPE(FLD_N) :: sn_tem, sn_sal, sn_mld, sn_emp, sn_ice, sn_qsr, sn_wnd  ! informations about the fields to be read
      TYPE(FLD_N) :: sn_uwd, sn_vwd, sn_wwd, sn_avt, sn_ubl, sn_vbl          !   "                                 "
      TYPE(FLD_N) :: sn_ahu, sn_ahv, sn_ahw, sn_eiu, sn_eiv, sn_eiw          !   "                                 "
      !
      NAMELIST/namdta_dyn/cn_dir, ln_dynwzv, ln_dynbbl, ln_degrad,    &
         &                sn_tem, sn_sal, sn_mld, sn_emp, sn_ice, sn_qsr, sn_wnd,  &
         &                sn_uwd, sn_vwd, sn_wwd, sn_avt, sn_ubl, sn_vbl,          &
         &                sn_ahu, sn_ahv, sn_ahw, sn_eiu, sn_eiv, sn_eiw       

      !!----------------------------------------------------------------------
      !                                   ! ============
      !                                   !   Namelist
      !                                   ! ============
      ! (NB: frequency positive => hours, negative => months)
      !                !   file      ! frequency !  variable  ! time intep !  clim  ! 'yearly' or ! weights  ! rotation   !
      !                !   name      !  (hours)  !   name     !   (T/F)    !  (T/F) !  'monthly'  ! filename ! pairs      !
      sn_tem  = FLD_N( 'dyna_grid_T' ,    120    , 'votemper' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_sal  = FLD_N( 'dyna_grid_T' ,    120    , 'vosaline' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_mld  = FLD_N( 'dyna_grid_T' ,    120    , 'somixght' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_emp  = FLD_N( 'dyna_grid_T' ,    120    , 'sowaflcd' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_ice  = FLD_N( 'dyna_grid_T' ,    120    , 'soicecov' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_qsr  = FLD_N( 'dyna_grid_T' ,    120    , 'soshfldo' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_wnd  = FLD_N( 'dyna_grid_T' ,    120    , 'sowindsp' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_uwd  = FLD_N( 'dyna_grid_U' ,    120    , 'vozocrtx' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_vwd  = FLD_N( 'dyna_grid_V' ,    120    , 'vomecrty' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_wwd  = FLD_N( 'dyna_grid_W' ,    120    , 'vovecrtz' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_avt  = FLD_N( 'dyna_grid_W' ,    120    , 'votkeavt' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_ubl  = FLD_N( 'dyna_grid_U' ,    120    , 'sobblcox' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_vbl  = FLD_N( 'dyna_grid_V' ,    120    , 'sobblcoy' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_ahu  = FLD_N( 'dyna_grid_U' ,    120    , 'vozoahtu' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_ahv  = FLD_N( 'dyna_grid_V' ,    120    , 'vomeahtv' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_ahw  = FLD_N( 'dyna_grid_W' ,    120    , 'voveahtz' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_eiu  = FLD_N( 'dyna_grid_U' ,    120    , 'vozoaeiu' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_eiv  = FLD_N( 'dyna_grid_V' ,    120    , 'vomeaeiv' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_eiw  = FLD_N( 'dyna_grid_W' ,    120    , 'voveaeiw' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      !
      REWIND( numnam )                          ! read in namlist namdta_dyn
      READ  ( numnam, namdta_dyn )
      !                                         ! store namelist information in an array
      !                                         ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dta_dyn : offline dynamics '
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namdta_dyn'
         WRITE(numout,*) '      vertical velocity read from file (T) or computed (F) ln_dynwzv  = ', ln_dynwzv
         WRITE(numout,*) '      bbl coef read from file (T) or computed (F)          ln_dynbbl  = ', ln_dynbbl
         WRITE(numout,*) '      degradation option enabled (T) or not (F)            ln_degrad  = ', ln_degrad
         WRITE(numout,*)
      ENDIF
      ! 
      IF( ln_degrad .AND. .NOT.lk_degrad ) THEN
         CALL ctl_warn( 'dta_dyn_init: degradation option requires key_degrad activated ; force ln_degrad to false' )
         ln_degrad = .FALSE.
      ENDIF
      IF( ln_dynbbl .AND. ( .NOT.lk_trabbl .OR. lk_c1d ) ) THEN
         CALL ctl_warn( 'dta_dyn_init: bbl option requires key_trabbl activated ; force ln_dynbbl to false' )
         ln_dynbbl = .FALSE.
      ENDIF

      jf_tem = 1   ;   jf_sal = 2   ;  jf_mld = 3   ;  jf_emp = 4   ;   jf_ice = 5   ;   jf_qsr = 6 
      jf_wnd = 7   ;   jf_uwd = 8   ;  jf_vwd = 9   ;  jf_wwd = 10  ;   jf_avt = 11  ;   jfld  = 11
      !
      slf_d(jf_tem) = sn_tem   ;   slf_d(jf_sal) = sn_sal   ;   slf_d(jf_mld) = sn_mld
      slf_d(jf_emp) = sn_emp   ;   slf_d(jf_ice) = sn_ice   ;   slf_d(jf_qsr) = sn_qsr
      slf_d(jf_wnd) = sn_wnd   ;   slf_d(jf_uwd) = sn_uwd   ;   slf_d(jf_vwd) = sn_vwd
      slf_d(jf_wwd) = sn_wwd   ;   slf_d(jf_avt) = sn_avt 
      !
      IF( .NOT.ln_degrad ) THEN     ! no degrad option
         IF( lk_traldf_eiv .AND. ln_dynbbl ) THEN        ! eiv & bbl
                 jf_ubl  = 12      ;         jf_vbl  = 13      ;         jf_eiw  = 14   ;   jfld = 14
           slf_d(jf_ubl) = sn_ubl  ;   slf_d(jf_vbl) = sn_vbl  ;   slf_d(jf_eiw) = sn_eiw
         ENDIF
         IF( .NOT.lk_traldf_eiv .AND. ln_dynbbl ) THEN   ! no eiv & bbl
                 jf_ubl  = 12      ;         jf_vbl  = 13      ;   jfld = 13
           slf_d(jf_ubl) = sn_ubl  ;   slf_d(jf_vbl) = sn_vbl
         ENDIF
         IF( lk_traldf_eiv .AND. .NOT.ln_dynbbl ) THEN   ! eiv & no bbl
           jf_eiw = 12   ;   jfld = 12   ;   slf_d(jf_eiw) = sn_eiw
         ENDIF
      ELSE
              jf_ahu  = 12      ;         jf_ahv  = 13      ;         jf_ahw  = 14   ;   jfld = 14
        slf_d(jf_ahu) = sn_ahu  ;   slf_d(jf_ahv) = sn_ahv  ;   slf_d(jf_ahw) = sn_ahw
        IF( lk_traldf_eiv .AND. ln_dynbbl ) THEN         ! eiv & bbl
                 jf_ubl  = 15      ;         jf_vbl  = 16      
           slf_d(jf_ubl) = sn_ubl  ;   slf_d(jf_vbl) = sn_vbl  
                 jf_eiu  = 17      ;         jf_eiv  = 18      ;          jf_eiw  = 19   ;   jfld = 19
           slf_d(jf_eiu) = sn_eiu  ;   slf_d(jf_eiv) = sn_eiv  ;    slf_d(jf_eiw) = sn_eiw
        ENDIF
        IF( .NOT.lk_traldf_eiv .AND. ln_dynbbl ) THEN    ! no eiv & bbl
                 jf_ubl  = 15      ;         jf_vbl  = 16      ;   jfld = 16
           slf_d(jf_ubl) = sn_ubl  ;   slf_d(jf_vbl) = sn_vbl
        ENDIF
        IF( lk_traldf_eiv .AND. .NOT.ln_dynbbl ) THEN    ! eiv & no bbl
                 jf_eiu  = 15      ;         jf_eiv  = 16      ;         jf_eiw  = 17   ;   jfld = 17
           slf_d(jf_eiu) = sn_eiu  ;   slf_d(jf_eiv) = sn_eiv  ;   slf_d(jf_eiw) = sn_eiw
        ENDIF
      ENDIF
  
      ALLOCATE( sf_dyn(jfld), STAT=ierr )         ! set sf structure
      IF( ierr > 0 ) THEN
         CALL ctl_stop( 'dta_dyn: unable to allocate sf structure' )   ;   RETURN
      ENDIF
      ! Open file for each variable to get his number of dimension
      DO ifpr = 1, jfld
         CALL iom_open( slf_d(ifpr)%clname, inum )
         idv   = iom_varid( inum , slf_d(ifpr)%clvar )  ! id of the variable sdjf%clvar
         idimv = iom_file ( inum )%ndims(idv)             ! number of dimension for variable sdjf%clvar
         IF( inum /= 0 )   CALL iom_close( inum )       ! close file if already open
         IF( idimv == 3 ) THEN    ! 2D variable
                                      ALLOCATE( sf_dyn(ifpr)%fnow(jpi,jpj,1)    , STAT=ierr0 )
            IF( slf_d(ifpr)%ln_tint ) ALLOCATE( sf_dyn(ifpr)%fdta(jpi,jpj,1,2)  , STAT=ierr1 )
         ELSE                     ! 3D variable
                                      ALLOCATE( sf_dyn(ifpr)%fnow(jpi,jpj,jpk)  , STAT=ierr0 )
            IF( slf_d(ifpr)%ln_tint ) ALLOCATE( sf_dyn(ifpr)%fdta(jpi,jpj,jpk,2), STAT=ierr1 )
         ENDIF
         IF( ierr0 + ierr1 > 0 ) THEN
            CALL ctl_stop( 'dta_dyn_init : unable to allocate sf_dyn array structure' )   ;   RETURN
         ENDIF
      END DO
      !                                         ! fill sf with slf_i and control print
      CALL fld_fill( sf_dyn, slf_d, cn_dir, 'dta_dyn_init', 'Data in file', 'namdta_dyn' )
      !
      IF( lk_ldfslp .AND. .NOT.lk_c1d ) THEN                  ! slopes 
         IF( sf_dyn(jf_tem)%ln_tint ) THEN      ! time interpolation
            ALLOCATE( uslpdta (jpi,jpj,jpk,2), vslpdta (jpi,jpj,jpk,2),    &
            &         wslpidta(jpi,jpj,jpk,2), wslpjdta(jpi,jpj,jpk,2), STAT=ierr2 )
         ELSE
            ALLOCATE( uslpnow (jpi,jpj,jpk)  , vslpnow (jpi,jpj,jpk)  ,    &
            &         wslpinow(jpi,jpj,jpk)  , wslpjnow(jpi,jpj,jpk)  , STAT=ierr2 )
         ENDIF 
         IF( ierr2 > 0 ) THEN
            CALL ctl_stop( 'dta_dyn_init : unable to allocate slope arrays' )   ;   RETURN
         ENDIF
      ENDIF
      IF( ln_dynwzv ) THEN                  ! slopes 
         IF( sf_dyn(jf_uwd)%ln_tint ) THEN      ! time interpolation
            ALLOCATE( wdta(jpi,jpj,jpk,2), STAT=ierr3 )
         ELSE
            ALLOCATE( wnow(jpi,jpj,jpk)  , STAT=ierr3 )
         ENDIF 
         IF( ierr3 > 0 ) THEN
            CALL ctl_stop( 'dta_dyn_init : unable to allocate wdta arrays' )   ;   RETURN
         ENDIF
      ENDIF
      !
      CALL dta_dyn( nit000 )
      !
   END SUBROUTINE dta_dyn_init

   SUBROUTINE dta_dyn_wzv( pu, pv, pw )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE wzv  ***
      !!
      !! ** Purpose :   Compute the now vertical velocity after the array swap
      !!
      !! ** Method  : - compute the now divergence given by :
      !!         * z-coordinate ONLY !!!!
      !!         hdiv = 1/(e1t*e2t) [ di(e2u  u) + dj(e1v  v) ]
      !!     - Using the incompressibility hypothesis, the vertical
      !!      velocity is computed by integrating the horizontal divergence
      !!      from the bottom to the surface.
      !!        The boundary conditions are w=0 at the bottom (no flux).
      !!----------------------------------------------------------------------
      USE oce, ONLY:  zhdiv => hdivn
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) :: pu, pv    !:  horizontal velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) :: pw        !:  vertical velocity
      !!
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zu, zu1, zv, zv1, zet
      !!----------------------------------------------------------------------
      !
      ! Computation of vertical velocity using horizontal divergence
      zhdiv(:,:,:) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zu  = pu(ji  ,jj  ,jk) * umask(ji  ,jj  ,jk) * e2u(ji  ,jj  ) * fse3u(ji  ,jj  ,jk)
               zu1 = pu(ji-1,jj  ,jk) * umask(ji-1,jj  ,jk) * e2u(ji-1,jj  ) * fse3u(ji-1,jj  ,jk)
               zv  = pv(ji  ,jj  ,jk) * vmask(ji  ,jj  ,jk) * e1v(ji  ,jj  ) * fse3v(ji  ,jj  ,jk)
               zv1 = pv(ji  ,jj-1,jk) * vmask(ji  ,jj-1,jk) * e1v(ji  ,jj-1) * fse3v(ji  ,jj-1,jk)
               zet = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
               zhdiv(ji,jj,jk) = ( zu - zu1 + zv - zv1 ) * zet 
            END DO
         END DO
      END DO
      CALL lbc_lnk( zhdiv, 'T', 1. )      ! Lateral boundary conditions on zhdiv
      !
      ! computation of vertical velocity from the bottom
      pw(:,:,jpk) = 0._wp
      DO jk = jpkm1, 1, -1
         pw(:,:,jk) = pw(:,:,jk+1) - fse3t(:,:,jk) * zhdiv(:,:,jk)
      END DO
      !
   END SUBROUTINE dta_dyn_wzv

   SUBROUTINE dta_dyn_slp( kt, pts, puslp, pvslp, pwslpi, pwslpj )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE dta_dyn_slp  ***
      !!
      !! ** Purpose : Computation of slope
      !!
      !!---------------------------------------------------------------------
      INTEGER ,                              INTENT(in ) :: kt       ! time step
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in ) :: pts      ! temperature/salinity
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(out) :: puslp    ! zonal isopycnal slopes
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(out) :: pvslp    ! meridional isopycnal slopes
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(out) :: pwslpi   ! zonal diapycnal slopes
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(out) :: pwslpj   ! meridional diapycnal slopes
      !!---------------------------------------------------------------------
#if defined key_ldfslp && ! defined key_c1d
      CALL eos( pts, rhd, rhop )   ! Time-filtered in situ density 
      CALL bn2( pts, rn2 )         ! before Brunt-Vaisala frequency
      IF( ln_zps )   &
         &  CALL zps_hde( kt, jpts, pts, gtsu, gtsv, rhd, gru, grv )  ! Partial steps: before Horizontal DErivative
         !                                                            ! of t, s, rd at the bottom ocean level
      CALL zdf_mxl( kt )            ! mixed layer depth
      CALL ldf_slp( kt, rhd, rn2 )  ! slopes
      puslp (:,:,:) = uslp (:,:,:) 
      pvslp (:,:,:) = vslp (:,:,:) 
      pwslpi(:,:,:) = wslpi(:,:,:) 
      pwslpj(:,:,:) = wslpj(:,:,:) 
#else
      puslp (:,:,:) = 0.            ! to avoid warning when compiling
      pvslp (:,:,:) = 0.
      pwslpi(:,:,:) = 0.
      pwslpj(:,:,:) = 0.
#endif
      !
   END SUBROUTINE dta_dyn_slp
   !!======================================================================
END MODULE dtadyn
