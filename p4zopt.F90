MODULE p4zopt
   !!======================================================================
   !!                         ***  MODULE p4zopt  ***
   !! TOP - PISCES : Compute the light availability in the water column
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!----------------------------------------------------------------------
#if defined  key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_opt       : light availability in the water column
   !!----------------------------------------------------------------------
   USE trc            ! tracer variables
   USE oce_trc        ! tracer-ocean share variables
   USE sms_pisces     ! Source Minus Sink of PISCES
   USE iom            ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_opt        ! called in p4zbio.F90 module
   PUBLIC   p4z_opt_init   ! called in trcsms_pisces.F90 module
   PUBLIC   p4z_opt_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: etot, enano, ediat   !: PAR for phyto, nano and diat 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: emoy                 !: averaged PAR in the mixed layer

   INTEGER  ::   nksrp   ! levels below which the light cannot penetrate ( depth larger than 391 m)
   REAL(wp) ::   parlux = 0.43_wp / 3._wp

   REAL(wp), DIMENSION(3,61), PUBLIC ::   xkrgb   !: tabulated attenuation coefficients for RGB absorption
   
   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zopt.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_opt( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt  ***
      !!
      !! ** Purpose :   Compute the light availability in the water column
      !!              depending on the depth and the chlorophyll concentration
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   irgb
      REAL(wp) ::   zchl, zxsi0r
      REAL(wp) ::   zc0 , zc1 , zc2, zc3, z1_dep
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zdepmoy, zetmp, zetmp1, zetmp2
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zekg, zekr, zekb, ze0, ze1, ze2, ze3
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_opt')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj,      zdepmoy, zetmp, zetmp1, zetmp2       )
      CALL wrk_alloc( jpi, jpj, jpk, zekg, zekr, zekb, ze0, ze1, ze2, ze3 )

      !     Initialisation of variables used to compute PAR
      !     -----------------------------------------------
      ze1 (:,:,jpk) = 0._wp
      ze2 (:,:,jpk) = 0._wp
      ze3 (:,:,jpk) = 0._wp

      !                                        !* attenuation coef. function of Chlorophyll and wavelength (Red-Green-Blue)
      DO jk = 1, jpkm1                         !  --------------------------------------------------------
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zchl = ( trn(ji,jj,jk,jpnch) + trn(ji,jj,jk,jpdch) + rtrn ) * 1.e6
               zchl = MIN(  10. , MAX( 0.05, zchl )  )
               irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )
               !                                                         
               zekb(ji,jj,jk) = xkrgb(1,irgb) * fse3t(ji,jj,jk)
               zekg(ji,jj,jk) = xkrgb(2,irgb) * fse3t(ji,jj,jk)
               zekr(ji,jj,jk) = xkrgb(3,irgb) * fse3t(ji,jj,jk)
            END DO
         END DO
      END DO


      !                                        !* Photosynthetically Available Radiation (PAR)
      !                                        !  --------------------------------------
!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            zc1 = parlux * qsr(ji,jj) * EXP( -0.5 * zekb(ji,jj,1) )
            zc2 = parlux * qsr(ji,jj) * EXP( -0.5 * zekg(ji,jj,1) )
            zc3 = parlux * qsr(ji,jj) * EXP( -0.5 * zekr(ji,jj,1) )
            ze1  (ji,jj,1) = zc1
            ze2  (ji,jj,1) = zc2
            ze3  (ji,jj,1) = zc3
            etot (ji,jj,1) = (       zc1 +        zc2 +       zc3 )
            enano(ji,jj,1) = ( 2.1 * zc1 + 0.42 * zc2 + 0.4 * zc3 )
            ediat(ji,jj,1) = ( 1.6 * zc1 + 0.69 * zc2 + 0.7 * zc3 )
         END DO
      END DO

    
      DO jk = 2, nksrp      
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zc1 = ze1(ji,jj,jk-1) * EXP( -0.5 * ( zekb(ji,jj,jk-1) + zekb(ji,jj,jk) ) )
               zc2 = ze2(ji,jj,jk-1) * EXP( -0.5 * ( zekg(ji,jj,jk-1) + zekg(ji,jj,jk) ) )
               zc3 = ze3(ji,jj,jk-1) * EXP( -0.5 * ( zekr(ji,jj,jk-1) + zekr(ji,jj,jk) ) )
               ze1  (ji,jj,jk) = zc1
               ze2  (ji,jj,jk) = zc2
               ze3  (ji,jj,jk) = zc3
               etot (ji,jj,jk) = (       zc1 +        zc2 +       zc3 )
               enano(ji,jj,jk) = ( 2.1 * zc1 + 0.42 * zc2 + 0.4 * zc3 )
               ediat(ji,jj,jk) = ( 1.6 * zc1 + 0.69 * zc2 + 0.7 * zc3 )
            END DO
         END DO
      END DO

      IF( ln_qsr_bio ) THEN                    !* heat flux accros w-level (used in the dynamics)
         !                                     !  ------------------------
         zxsi0r = 1.e0 / rn_si0
         !
         ze0  (:,:,1) = rn_abs * qsr(:,:)
         ze1  (:,:,1) = parlux * qsr(:,:)             ! surface value : separation in R-G-B + near surface 
         ze2  (:,:,1) = parlux * qsr(:,:)
         ze3  (:,:,1) = parlux * qsr(:,:)
         etot3(:,:,1) =          qsr(:,:) * tmask(:,:,1)
         !
         DO jk = 2, nksrp + 1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  zc0 = ze0(ji,jj,jk-1) * EXP( -fse3t(ji,jj,jk-1) * zxsi0r )
                  zc1 = ze1(ji,jj,jk-1) * EXP( -zekb(ji,jj,jk-1 ) )
                  zc2 = ze2(ji,jj,jk-1) * EXP( -zekg(ji,jj,jk-1 ) )
                  zc3 = ze3(ji,jj,jk-1) * EXP( -zekr(ji,jj,jk-1 ) )
                  ze0(ji,jj,jk) = zc0
                  ze1(ji,jj,jk) = zc1
                  ze2(ji,jj,jk) = zc2
                  ze3(ji,jj,jk) = zc3
                  etot3(ji,jj,jk) = ( zc0 + zc1 + zc2 + zc3 ) * tmask(ji,jj,jk)
              END DO
              !
            END DO
            !
        END DO
        !
      ENDIF

      !                                        !* Euphotic depth and level
      neln(:,:) = 1                            !  ------------------------
      heup(:,:) = 300.

      DO jk = 2, nksrp
         DO jj = 1, jpj
           DO ji = 1, jpi
              IF( etot(ji,jj,jk) >= 0.0043 * qsr(ji,jj) )  THEN
                 neln(ji,jj) = jk+1                    ! Euphotic level : 1rst T-level strictly below Euphotic layer
                 !                                     ! nb: ensure the compatibility with nmld_trc definition in trd_mld_trc_zint
                 heup(ji,jj) = fsdepw(ji,jj,jk+1)      ! Euphotic layer depth
              ENDIF
           END DO
        END DO
      END DO
 
      heup(:,:) = MIN( 300., heup(:,:) )
!Joan
!      heup(:,:) = MIN( 100., heup(:,:) )

      !                                        !* mean light over the mixed layer
      zdepmoy(:,:)   = 0.e0                    !  -------------------------------
      zetmp  (:,:)   = 0.e0
      zetmp1 (:,:)   = 0.e0
      zetmp2 (:,:)   = 0.e0

      DO jk = 1, nksrp
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( fsdepw(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  zetmp  (ji,jj) = zetmp  (ji,jj) + etot (ji,jj,jk) * fse3t(ji,jj,jk)
                  zetmp1 (ji,jj) = zetmp1 (ji,jj) + enano(ji,jj,jk) * fse3t(ji,jj,jk)
                  zetmp2 (ji,jj) = zetmp2 (ji,jj) + ediat(ji,jj,jk) * fse3t(ji,jj,jk)
                  zdepmoy(ji,jj) = zdepmoy(ji,jj) + fse3t(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO
      !
      emoy(:,:,:) = etot(:,:,:)
      !
      DO jk = 1, nksrp
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( fsdepw(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
                  emoy (ji,jj,jk) = zetmp (ji,jj) * z1_dep
                  enano(ji,jj,jk) = zetmp1(ji,jj) * z1_dep
                  ediat(ji,jj,jk) = zetmp2(ji,jj) * z1_dep
               ENDIF
            END DO
         END DO
      END DO

      IF( ln_diatrc ) THEN        ! save output diagnostics
        !
        IF( lk_iomput ) THEN
           IF( jnt == nrdttrc ) THEN
              CALL iom_put( "Heup", heup(:,:  ) * tmask(:,:,1) )  ! euphotic layer deptht
              CALL iom_put( "PAR" , etot(:,:,:) * tmask(:,:,:) )  ! Photosynthetically Available Radiation
!!JOAN
              CALL iom_put( "PARmld", emoy(:,:,:) * tmask(:,:,:) )  ! PAR averaged in the ML
!!
	   ENDIF
        ELSE
           trc2d(:,:,  jp_pcs0_2d + 10) = heup(:,:  ) * tmask(:,:,1)  
           trc3d(:,:,:,jp_pcs0_3d + 3)  = etot(:,:,:) * tmask(:,:,:)
        ENDIF
        !
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj,      zdepmoy, zetmp, zetmp1, zetmp2 )
      CALL wrk_dealloc( jpi, jpj, jpk, zekg, zekr, zekb, ze0, ze1, ze2, ze3 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_opt')
      !
   END SUBROUTINE p4z_opt


   SUBROUTINE p4z_opt_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of tabulated attenuation coef
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_opt_init')
      !
      CALL trc_oce_rgb( xkrgb )                  ! tabulated attenuation coefficients
      nksrp = trc_oce_ext_lev( r_si2, 0.33e2 )     ! max level of light extinction (Blue Chl=0.01)
      !
      IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksrp, ' ref depth = ', gdepw_0(nksrp+1), ' m'
      !
                         etot (:,:,:) = 0._wp
                         enano(:,:,:) = 0._wp
                         ediat(:,:,:) = 0._wp
      IF( ln_qsr_bio )   etot3(:,:,:) = 0._wp
      ! 
      IF( nn_timing == 1 )  CALL timing_stop('p4z_opt_init')
      !
   END SUBROUTINE p4z_opt_init


   INTEGER FUNCTION p4z_opt_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( etot (jpi,jpj,jpk) , enano(jpi,jpj,jpk) ,     &
         &      ediat(jpi,jpj,jpk) , emoy (jpi,jpj,jpk) , STAT=p4z_opt_alloc )
         !
      IF( p4z_opt_alloc /= 0 ) CALL ctl_warn('p4z_opt_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_opt_alloc

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PISCES bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE p4z_opt                   ! Empty routine
   END SUBROUTINE p4z_opt
#endif 

   !!======================================================================
END MODULE  p4zopt
