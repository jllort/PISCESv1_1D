MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_sink_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wscal    !: Calcite and BSi sinking speeds

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking, sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal, sinksil   !: CaCO3 and BSi sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer            !: Small BFe sinking fluxes
#if ! defined key_kriest
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer2           !: Big iron sinking fluxes
#endif

   INTEGER  :: iksed  = 10

#if  defined key_kriest
   REAL(wp) ::  xkr_sfact    = 250.     !: Sinking factor
   REAL(wp) ::  xkr_stick    = 0.2      !: Stickiness
   REAL(wp) ::  xkr_nnano    = 2.337    !: Nbr of cell in nano size class
   REAL(wp) ::  xkr_ndiat    = 3.718    !: Nbr of cell in diatoms size class
   REAL(wp) ::  xkr_nmeso    = 7.147    !: Nbr of cell in mesozoo  size class
   REAL(wp) ::  xkr_naggr    = 9.877    !: Nbr of cell in aggregates  size class

   REAL(wp) ::  xkr_frac 

   REAL(wp), PUBLIC ::  xkr_dnano       !: Size of particles in nano pool
   REAL(wp), PUBLIC ::  xkr_ddiat       !: Size of particles in diatoms pool
   REAL(wp), PUBLIC ::  xkr_dmeso       !: Size of particles in mesozoo pool
   REAL(wp), PUBLIC ::  xkr_daggr       !: Size of particles in aggregates pool
   REAL(wp), PUBLIC ::  xkr_wsbio_min   !: min vertical particle speed
   REAL(wp), PUBLIC ::  xkr_wsbio_max   !: max vertical particle speed

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   xnumm   !:  maximum number of particles in aggregates
#endif

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsink.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_kriest
   !!----------------------------------------------------------------------
   !!   'key_kriest'                                                    ???
   !!----------------------------------------------------------------------

   SUBROUTINE p4z_sink ( kt, jnt )
      !!---------------------------------------------------------------------
      !!                ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to
      !!              gravitational sinking - Kriest parameterization
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, jnt
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zagg1, zagg2, zagg3, zagg4, zagg5, zaggsi, zaggsh
      REAL(wp) :: zagg , zaggdoc, znumdoc
      REAL(wp) :: znum , zeps, zfm, zgm, zsm
      REAL(wp) :: zdiv , zdiv1, zdiv2, zdiv3, zdiv4, zdiv5
      REAL(wp) :: zval1, zval2, zval3, zval4
      REAL(wp) :: zrfact2
      INTEGER  :: ik1
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: znum3d 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')
      !
      CALL wrk_alloc( jpi, jpj, jpk, znum3d )
      !
      !     Initialisation of variables used to compute Sinking Speed
      !     ---------------------------------------------------------

      znum3d(:,:,:) = 0.e0
      zval1 = 1. + xkr_zeta
      zval2 = 1. + xkr_zeta + xkr_eta
      zval3 = 1. + xkr_eta

      !     Computation of the vertical sinking speed : Kriest et Evans, 2000
      !     -----------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) /= 0.e0 ) THEN
                  znum = trn(ji,jj,jk,jppoc) / ( trn(ji,jj,jk,jpnum) + rtrn ) / xkr_massp
                  ! -------------- To avoid sinking speed over 50 m/day -------
                  znum  = MIN( xnumm(jk), znum )
                  znum  = MAX( 1.1      , znum )
                  znum3d(ji,jj,jk) = znum
                  !------------------------------------------------------------
                  zeps  = ( zval1 * znum - 1. )/ ( znum - 1. )
                  zfm   = xkr_frac**( 1. - zeps )
                  zgm   = xkr_frac**( zval1 - zeps )
                  zdiv  = MAX( 1.e-4, ABS( zeps - zval2 ) ) * SIGN( 1., ( zeps - zval2 ) )
                  zdiv1 = zeps - zval3
                  wsbio3(ji,jj,jk) = xkr_wsbio_min * ( zeps - zval1 ) / zdiv    &
                     &             - xkr_wsbio_max *   zgm * xkr_eta  / zdiv
                  wsbio4(ji,jj,jk) = xkr_wsbio_min *   ( zeps-1. )    / zdiv1   &
                     &             - xkr_wsbio_max *   zfm * xkr_eta  / zdiv1
                  IF( znum == 1.1)   wsbio3(ji,jj,jk) = wsbio4(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO

      wscal(:,:,:) = MAX( wsbio3(:,:,:), 50._wp )

      !   INITIALIZE TO ZERO ALL THE SINKING ARRAYS
      !   -----------------------------------------

      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0

     !   Compute the sedimentation term using p4zsink2 for all the sinking particles
     !   -----------------------------------------------------

      CALL p4z_sink2( wsbio3, sinking , jppoc )
      CALL p4z_sink2( wsbio4, sinking2, jpnum )
      CALL p4z_sink2( wsbio3, sinkfer , jpsfe )
      CALL p4z_sink2( wscal , sinksil , jpdsi )
      CALL p4z_sink2( wscal , sinkcal , jpcal )

     !  Exchange between organic matter compartments due to coagulation/disaggregation
     !  ---------------------------------------------------

      zval1 = 1. + xkr_zeta
      zval2 = 1. + xkr_eta
      zval3 = 3. + xkr_eta
      zval4 = 4. + xkr_eta

      DO jk = 1,jpkm1
         DO jj = 1,jpj
            DO ji = 1,jpi
               IF( tmask(ji,jj,jk) /= 0.e0 ) THEN

                  znum = trn(ji,jj,jk,jppoc)/(trn(ji,jj,jk,jpnum)+rtrn) / xkr_massp
                  !-------------- To avoid sinking speed over 50 m/day -------
                  znum  = min(xnumm(jk),znum)
                  znum  = MAX( 1.1,znum)
                  !------------------------------------------------------------
                  zeps  = ( zval1 * znum - 1.) / ( znum - 1.)
                  zdiv  = MAX( 1.e-4, ABS( zeps - zval3) ) * SIGN( 1., zeps - zval3 )
                  zdiv1 = MAX( 1.e-4, ABS( zeps - 4.   ) ) * SIGN( 1., zeps - 4.    )
                  zdiv2 = zeps - 2.
                  zdiv3 = zeps - 3.
                  zdiv4 = zeps - zval2
                  zdiv5 = 2.* zeps - zval4
                  zfm   = xkr_frac**( 1.- zeps )
                  zsm   = xkr_frac**xkr_eta

                  !    Part I : Coagulation dependant on turbulence
                  !    ----------------------------------------------

                  zagg1 = ( 0.163 * trn(ji,jj,jk,jpnum)**2               &
                     &            * 2.*( (zfm-1.)*(zfm*xkr_mass_max**3-xkr_mass_min**3)    &
                     &            * (zeps-1)/zdiv1 + 3.*(zfm*xkr_mass_max-xkr_mass_min)    &
                     &            * (zfm*xkr_mass_max**2-xkr_mass_min**2)                  &
                     &            * (zeps-1.)**2/(zdiv2*zdiv3)) 
                  zagg2 =  2*0.163*trn(ji,jj,jk,jpnum)**2*zfm*                       &
                     &                   ((xkr_mass_max**3+3.*(xkr_mass_max**2          &
                     &                    *xkr_mass_min*(zeps-1.)/zdiv2                 &
                     &                    +xkr_mass_max*xkr_mass_min**2*(zeps-1.)/zdiv3)    &
                     &                    +xkr_mass_min**3*(zeps-1)/zdiv1)                  &
                     &                    -zfm*xkr_mass_max**3*(1.+3.*((zeps-1.)/           &
                     &                    (zeps-2.)+(zeps-1.)/zdiv3)+(zeps-1.)/zdiv1))    

                  zagg3 =  0.163*trn(ji,jj,jk,jpnum)**2*zfm**2*8. * xkr_mass_max**3  
                  
                 !    Aggregation of small into large particles
                 !    Part II : Differential settling
                 !    ----------------------------------------------

                  zagg4 =  2.*3.141*0.125*trn(ji,jj,jk,jpnum)**2*                       &
                     &                 xkr_wsbio_min*(zeps-1.)**2                         &
                     &                 *(xkr_mass_min**2*((1.-zsm*zfm)/(zdiv3*zdiv4)      &
                     &                 -(1.-zfm)/(zdiv*(zeps-1.)))-                       &
                     &                 ((zfm*zfm*xkr_mass_max**2*zsm-xkr_mass_min**2)     &
                     &                 *xkr_eta)/(zdiv*zdiv3*zdiv5) )   

                  zagg5 =   2.*3.141*0.125*trn(ji,jj,jk,jpnum)**2                         &
                     &                 *(zeps-1.)*zfm*xkr_wsbio_min                        &
                     &                 *(zsm*(xkr_mass_min**2-zfm*xkr_mass_max**2)         &
                     &                 /zdiv3-(xkr_mass_min**2-zfm*zsm*xkr_mass_max**2)    &
                     &                 /zdiv)  
                  zaggsi = ( zagg4 + zagg5 ) * xstep / 10.

                  zagg = 0.5 * xkr_stick * ( zaggsh + zaggsi )

                  !     Aggregation of DOC to small particles
                  !     --------------------------------------

                  zaggdoc = ( 0.4 * trn(ji,jj,jk,jpdoc)               &
                     &        + 1018.  * trn(ji,jj,jk,jppoc)  ) * xstep    &
                     &        * xdiss(ji,jj,jk) * trn(ji,jj,jk,jpdoc)

# if defined key_degrad
                   zagg1   = zagg1   * facvol(ji,jj,jk)                 
                   zagg2   = zagg2   * facvol(ji,jj,jk)                 
                   zagg3   = zagg3   * facvol(ji,jj,jk)                 
                   zagg4   = zagg4   * facvol(ji,jj,jk)                 
                   zagg5   = zagg5   * facvol(ji,jj,jk)                 
                   zaggdoc = zaggdoc * facvol(ji,jj,jk)                 
# endif
                  zaggsh = ( zagg1 + zagg2 + zagg3 ) * rfact2 * xdiss(ji,jj,jk) / 1000.
                  zaggsi = ( zagg4 + zagg5 ) * xstep / 10.
                  zagg = 0.5 * xkr_stick * ( zaggsh + zaggsi )
                  !
                  znumdoc = trn(ji,jj,jk,jpnum) / ( trn(ji,jj,jk,jppoc) + rtrn )
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zaggdoc
                  tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + zaggdoc * znumdoc - zagg
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc

               ENDIF
            END DO
         END DO
      END DO

      IF( ln_diatrc ) THEN
         !
         ik1 = iksed + 1
         zrfact2 = 1.e3 * rfact2r
         IF( jnt == nrdttrc ) THEN
           CALL iom_put( "POCFlx"  , sinking (:,:,:)      * zrfact2 * tmask(:,:,:) )  ! POC export
           CALL iom_put( "NumFlx"  , sinking2 (:,:,:)     * zrfact2 * tmask(:,:,:) )  ! Num export
           CALL iom_put( "SiFlx"   , sinksil (:,:,:)      * zrfact2 * tmask(:,:,:) )  ! Silica export
           CALL iom_put( "CaCO3Flx", sinkcal (:,:,:)      * zrfact2 * tmask(:,:,:) )  ! Calcite export
           CALL iom_put( "xnum"    , znum3d  (:,:,:)                * tmask(:,:,:) )  ! Number of particles in aggregats
           CALL iom_put( "W1"      , wsbio3  (:,:,:)                * tmask(:,:,:) )  ! sinking speed of POC
           CALL iom_put( "W2"      , wsbio4  (:,:,:)                * tmask(:,:,:) )  ! sinking speed of aggregats
           CALL iom_put( "PMO"     , sinking (:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! POC export at 100m
           CALL iom_put( "PMO2"    , sinking2(:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! Num export at 100m
           CALL iom_put( "ExpFe1"  , sinkfer (:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! Export of iron at 100m
           CALL iom_put( "ExpSi"   , sinksil (:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! export of silica at 100m
           CALL iom_put( "ExpCaCO3", sinkcal (:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! export of calcite at 100m
         ENDIF
# if ! defined key_iomput
         trc2d(:,:  ,jp_pcs0_2d + 4)  = sinking (:,:,ik1)    * zrfact2 * tmask(:,:,1)
         trc2d(:,:  ,jp_pcs0_2d + 5)  = sinking2(:,:,ik1)    * zrfact2 * tmask(:,:,1)
         trc2d(:,:  ,jp_pcs0_2d + 6)  = sinkfer (:,:,ik1)    * zrfact2 * tmask(:,:,1)
         trc2d(:,:  ,jp_pcs0_2d + 7)  = sinksil (:,:,ik1)    * zrfact2 * tmask(:,:,1)
         trc2d(:,:  ,jp_pcs0_2d + 8)  = sinkcal (:,:,ik1)    * zrfact2 * tmask(:,:,1)
         trc3d(:,:,:,jp_pcs0_3d + 11) = sinking (:,:,:)      * zrfact2 * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 12) = sinking2(:,:,:)      * zrfact2 * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 13) = sinksil (:,:,:)      * zrfact2 * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 14) = sinkcal (:,:,:)      * zrfact2 * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 15) = znum3d  (:,:,:)                * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 16) = wsbio3  (:,:,:)                * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 17) = wsbio4  (:,:,:)                * tmask(:,:,:)
# endif
        !
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_alloc( jpi, jpj, jpk, znum3d )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!
      !! ** Purpose :   Initialization of sinking parameters
      !!                Kriest parameterization only
      !!
      !! ** Method  :   Read the nampiskrs namelist and check the parameters
      !!      called at the first timestep 
      !!
      !! ** input   :   Namelist nampiskrs
      !!----------------------------------------------------------------------
      INTEGER  ::   jk, jn, kiter
      REAL(wp) ::   znum, zdiv
      REAL(wp) ::   zws, zwr, zwl,wmax, znummax
      REAL(wp) ::   zmin, zmax, zl, zr, xacc
      !
      NAMELIST/nampiskrs/ xkr_sfact, xkr_stick ,  &
         &                xkr_nnano, xkr_ndiat, xkr_nmeso, xkr_naggr
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink_init')
      !
      REWIND( numnatp )                     ! read nampiskrs
      READ  ( numnatp, nampiskrs )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampiskrs'
         WRITE(numout,*) '    Sinking factor                           xkr_sfact    = ', xkr_sfact
         WRITE(numout,*) '    Stickiness                               xkr_stick    = ', xkr_stick
         WRITE(numout,*) '    Nbr of cell in nano size class           xkr_nnano    = ', xkr_nnano
         WRITE(numout,*) '    Nbr of cell in diatoms size class        xkr_ndiat    = ', xkr_ndiat
         WRITE(numout,*) '    Nbr of cell in mesozoo size class        xkr_nmeso    = ', xkr_nmeso
         WRITE(numout,*) '    Nbr of cell in aggregates size class     xkr_naggr    = ', xkr_naggr
      ENDIF


      ! max and min vertical particle speed
      xkr_wsbio_min = xkr_sfact * xkr_mass_min**xkr_eta
      xkr_wsbio_max = xkr_sfact * xkr_mass_max**xkr_eta
      WRITE(numout,*) ' max and min vertical particle speed ', xkr_wsbio_min, xkr_wsbio_max

      !
      !    effect of the sizes of the different living pools on particle numbers
      !    nano = 2um-20um -> mean size=6.32 um -> ws=2.596 -> xnum=xnnano=2.337
      !    diat and microzoo = 10um-200um -> 44.7 -> 8.732 -> xnum=xndiat=3.718
      !    mesozoo = 200um-2mm -> 632.45 -> 45.14 -> xnum=xnmeso=7.147
      !    aggregates = 200um-10mm -> 1414 -> 74.34 -> xnum=xnaggr=9.877
      !    doc aggregates = 1um
      ! ----------------------------------------------------------

      xkr_dnano = 1. / ( xkr_massp * xkr_nnano )
      xkr_ddiat = 1. / ( xkr_massp * xkr_ndiat )
      xkr_dmeso = 1. / ( xkr_massp * xkr_nmeso )
      xkr_daggr = 1. / ( xkr_massp * xkr_naggr )

      !!---------------------------------------------------------------------
      !!    'key_kriest'                                                  ???
      !!---------------------------------------------------------------------
      !  COMPUTATION OF THE VERTICAL PROFILE OF MAXIMUM SINKING SPEED
      !  Search of the maximum number of particles in aggregates for each k-level.
      !  Bissection Method
      !--------------------------------------------------------------------
      WRITE(numout,*)
      WRITE(numout,*)'    kriest : Compute maximum number of particles in aggregates'

      xacc     =  0.001_wp
      kiter    = 50
      zmin     =  1.10_wp
      zmax     = xkr_mass_max / xkr_mass_min
      xkr_frac = zmax

      DO jk = 1,jpk
         zl = zmin
         zr = zmax
         wmax = 0.5 * fse3t(1,1,jk) * rday / rfact2
         zdiv = xkr_zeta + xkr_eta - xkr_eta * zl
         znum = zl - 1.
         zwl =  xkr_wsbio_min * xkr_zeta / zdiv &
            & - ( xkr_wsbio_max * xkr_eta * znum * &
            &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
            & - wmax

         zdiv = xkr_zeta + xkr_eta - xkr_eta * zr
         znum = zr - 1.
         zwr =  xkr_wsbio_min * xkr_zeta / zdiv &
            & - ( xkr_wsbio_max * xkr_eta * znum * &
            &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
            & - wmax
iflag:   DO jn = 1, kiter
            IF    ( zwl == 0._wp ) THEN   ;   znummax = zl
            ELSEIF( zwr == 0._wp ) THEN   ;   znummax = zr
            ELSE
               znummax = ( zr + zl ) / 2.
               zdiv = xkr_zeta + xkr_eta - xkr_eta * znummax
               znum = znummax - 1.
               zws =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax
               IF( zws * zwl < 0. ) THEN   ;   zr = znummax
               ELSE                        ;   zl = znummax
               ENDIF
               zdiv = xkr_zeta + xkr_eta - xkr_eta * zl
               znum = zl - 1.
               zwl =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax

               zdiv = xkr_zeta + xkr_eta - xkr_eta * zr
               znum = zr - 1.
               zwr =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax
               !
               IF ( ABS ( zws )  <= xacc ) EXIT iflag
               !
            ENDIF
            !
         END DO iflag

         xnumm(jk) = znummax
         WRITE(numout,*) '       jk = ', jk, ' wmax = ', wmax,' xnum max = ', xnumm(jk)
         !
      END DO
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink_init')
      !
  END SUBROUTINE p4z_sink_init

#else

   SUBROUTINE p4z_sink ( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, jnt
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zagg1, zagg2, zagg3, zagg4
      REAL(wp) ::   zagg , zaggfe, zaggdoc, zaggdoc2, zaggdoc3
      REAL(wp) ::   zfact, zwsmax, zmax, zstep
      REAL(wp) ::   zrfact2
      INTEGER  ::   ik1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')
      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., fsdepw(ji,jj,jk+1) - zmax ) / 5000._wp
               wsbio4(ji,jj,jk) = wsbio2 + ( 200.- wsbio2 ) * zfact
            END DO
         END DO
      END DO

      ! limit the values of the sinking speeds to avoid numerical instabilities  
      wsbio3(:,:,:) = wsbio
      !
      ! OA Below, this is garbage. the ideal would be to find a time-splitting 
      ! OA algorithm that does not increase the computing cost by too much
      ! OA In ROMS, I have included a time-splitting procedure. But it is 
      ! OA too expensive as the loop is computed globally. Thus, a small e3t
      ! OA at one place determines the number of subtimesteps globally
      ! OA AWFULLY EXPENSIVE !! Not able to find a better approach. Damned !!

      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwsmax = 0.8 * fse3t(ji,jj,jk) / xstep
               wsbio4(ji,jj,jk) = MIN( wsbio4(ji,jj,jk), zwsmax )
               wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax )
            END DO
         END DO
      END DO

      wscal(:,:,:) = wsbio4(:,:,:)

      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------

      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0
      sinkfer2(:,:,:) = 0.e0

      !   Compute the sedimentation term using p4zsink2 for all the sinking particles
      !   -----------------------------------------------------

      CALL p4z_sink2( wsbio3, sinking , jppoc )
      CALL p4z_sink2( wsbio3, sinkfer , jpsfe )
      CALL p4z_sink2( wsbio4, sinking2, jpgoc )
      CALL p4z_sink2( wsbio4, sinkfer2, jpbfe )
      CALL p4z_sink2( wsbio4, sinksil , jpdsi )
      CALL p4z_sink2( wscal , sinkcal , jpcal )

      !  Exchange between organic matter compartments due to coagulation/disaggregation
      !  ---------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               zstep = xstep 
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               zfact = zstep * xdiss(ji,jj,jk)
               !  Part I : Coagulation dependent on turbulence
               zagg1 = 354.  * zfact * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jppoc)
               zagg2 = 4452. * zfact * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jpgoc)

               ! Part II : Differential settling

               !  Aggregation of small into large particles
               zagg3 =  4.7 * zstep * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jpgoc)
               zagg4 =  0.4 * zstep * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jppoc)

               zagg   = zagg1 + zagg2 + zagg3 + zagg4
               zaggfe = zagg * trn(ji,jj,jk,jpsfe) / ( trn(ji,jj,jk,jppoc) + rtrn )

               ! Aggregation of DOC to small particles
               zaggdoc  = ( 0.83 * trn(ji,jj,jk,jpdoc) + 271. * trn(ji,jj,jk,jppoc) ) * zfact * trn(ji,jj,jk,jpdoc)
               zaggdoc2 = 1.07e4 * zfact * trn(ji,jj,jk,jpgoc) * trn(ji,jj,jk,jpdoc)
               zaggdoc3 =   0.02 * ( 16706. * trn(ji,jj,jk,jppoc) + 231. * trn(ji,jj,jk,jpdoc) ) * zstep * trn(ji,jj,jk,jpdoc)

               !  Update the trends
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zagg + zaggdoc + zaggdoc3
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zagg + zaggdoc2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zaggfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zaggfe
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc - zaggdoc2 - zaggdoc3
               !
            END DO
         END DO
      END DO

      IF( ln_diatrc ) THEN
         zrfact2 = 1.e3 * rfact2r
         ik1  = iksed + 1
         IF( lk_iomput ) THEN
           IF( jnt == nrdttrc ) THEN
              CALL iom_put( "EPC100"  , ( sinking(:,:,ik1) + sinking2(:,:,ik1) ) * zrfact2 * tmask(:,:,1) ) ! Export of carbon at 100m
              CALL iom_put( "EPFE100" , ( sinkfer(:,:,ik1) + sinkfer2(:,:,ik1) ) * zrfact2 * tmask(:,:,1) ) ! Export of iron at 100m
              CALL iom_put( "EPCAL100",   sinkcal(:,:,ik1)                       * zrfact2 * tmask(:,:,1) ) ! Export of calcite  at 100m
              CALL iom_put( "EPSI100" ,   sinksil(:,:,ik1)                       * zrfact2 * tmask(:,:,1) ) ! Export of biogenic silica at 100m
           ENDIF
         ELSE
           trc2d(:,:,jp_pcs0_2d + 4) = sinking (:,:,ik1) * zrfact2 * tmask(:,:,1)
           trc2d(:,:,jp_pcs0_2d + 5) = sinking2(:,:,ik1) * zrfact2 * tmask(:,:,1)
           trc2d(:,:,jp_pcs0_2d + 6) = sinkfer (:,:,ik1) * zrfact2 * tmask(:,:,1)
           trc2d(:,:,jp_pcs0_2d + 7) = sinkfer2(:,:,ik1) * zrfact2 * tmask(:,:,1)
           trc2d(:,:,jp_pcs0_2d + 8) = sinksil (:,:,ik1) * zrfact2 * tmask(:,:,1)
           trc2d(:,:,jp_pcs0_2d + 9) = sinkcal (:,:,ik1) * zrfact2 * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink

   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------
   END SUBROUTINE p4z_sink_init

#endif



   SUBROUTINE p4z_sink2( pwsink, psinkflx, jp_tra )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index      
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztraz, zakz, zwsink2 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink2')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, ztraz, zakz, zwsink2 )

      zstep = rfact2 / 2.

      ztraz(:,:,:) = 0.e0
      zakz (:,:,:) = 0.e0

      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0
      IF( lk_degrad ) THEN
         zwsink2(:,:,:) = zwsink2(:,:,:) * facvol(:,:,:)
      ENDIF


      ! Vertical advective flux
      DO jn = 1, 2
         !  first guess of the slopes interior values
         DO jk = 2, jpkm1
            ztraz(:,:,jk) = ( trn(:,:,jk-1,jp_tra) - trn(:,:,jk,jp_tra) ) * tmask(:,:,jk)
         END DO
         ztraz(:,:,1  ) = 0.0
         ztraz(:,:,jpk) = 0.0

         ! slopes
         DO jk = 2, jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zign = 0.25 + SIGN( 0.25, ztraz(ji,jj,jk) * ztraz(ji,jj,jk+1) )
                  zakz(ji,jj,jk) = ( ztraz(ji,jj,jk) + ztraz(ji,jj,jk+1) ) * zign
               END DO
            END DO
         END DO
         
         ! Slopes limitation
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zakz(ji,jj,jk) = SIGN( 1., zakz(ji,jj,jk) ) *        &
                     &             MIN( ABS( zakz(ji,jj,jk) ), 2. * ABS(ztraz(ji,jj,jk+1)), 2. * ABS(ztraz(ji,jj,jk) ) )
               END DO
            END DO
         END DO
         
         ! vertical advective flux
         DO jk = 1, jpkm1
            DO jj = 1, jpj      
               DO ji = 1, jpi    
                  zigma = zwsink2(ji,jj,jk+1) * zstep / fse3w(ji,jj,jk+1)
                  zew   = zwsink2(ji,jj,jk+1)
                  psinkflx(ji,jj,jk+1) = -zew * ( trn(ji,jj,jk,jp_tra) - 0.5 * ( 1 + zigma ) * zakz(ji,jj,jk) ) * zstep
               END DO
            END DO
         END DO
         !
         ! Boundary conditions
         psinkflx(:,:,1  ) = 0.e0
         psinkflx(:,:,jpk) = 0.e0
         
         DO jk=1,jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
                  trn(ji,jj,jk,jp_tra) = trn(ji,jj,jk,jp_tra) + zflx
               END DO
            END DO
         END DO

      ENDDO

      DO jk=1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
               trb(ji,jj,jk,jp_tra) = trb(ji,jj,jk,jp_tra) + 2. * zflx
            END DO
         END DO
      END DO

      trn     (:,:,:,jp_tra) = trb(:,:,:,jp_tra)
      psinkflx(:,:,:)        = 2. * psinkflx(:,:,:)
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztraz, zakz, zwsink2 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink2')
      !
   END SUBROUTINE p4z_sink2


   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wsbio3 (jpi,jpj,jpk) , wsbio4  (jpi,jpj,jpk) , wscal(jpi,jpj,jpk) ,     &
         &      sinking(jpi,jpj,jpk) , sinking2(jpi,jpj,jpk)                      ,     &                
         &      sinkcal(jpi,jpj,jpk) , sinksil (jpi,jpj,jpk)                      ,     &                
#if defined key_kriest
         &      xnumm(jpk)                                                        ,     &                
#else
         &      sinkfer2(jpi,jpj,jpk)                                             ,     &                
#endif
         &      sinkfer(jpi,jpj,jpk)                                              , STAT=p4z_sink_alloc )                
         !
      IF( p4z_sink_alloc /= 0 ) CALL ctl_warn('p4z_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_sink_alloc
   
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sink                    ! Empty routine
   END SUBROUTINE p4z_sink
#endif 

   !!======================================================================
END MODULE  p4zsink
