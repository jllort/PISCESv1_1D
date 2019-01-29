MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/scavenging of organic compounds
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/scavenging of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zopt          !  optical model
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p4zint          !  interpolation and computation of various fields
   USE p4zlim
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   !! * Shared module variables
   REAL(wp), PUBLIC ::  xremik    = 0.3_wp     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xremip    = 0.025_wp   !: remineralisation rate of DOC
   REAL(wp), PUBLIC ::  nitrif    = 0.05_wp    !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::  xsirem    = 0.003_wp   !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xsiremlab = 0.025_wp   !: fast remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xsilab    = 0.31_wp    !: fraction of labile biogenic silica 
   REAL(wp), PUBLIC ::  xlam1     = 0.005_wp   !: scavenging rate of Iron 
   REAL(wp), PUBLIC ::  oxymin    = 1.e-6_wp   !: halk saturation constant for anoxia 
   REAL(wp), PUBLIC ::  ligand    = 0.6E-9_wp  !: ligand concentration in the ocean 


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr     !: denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitnh4   !: -    -    -    -   -


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zrem.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremip, zremik , zlam1b, zdepbac2
      REAL(wp) ::   zkeq  , zfeequi, zsiremin, zfesatur
      REAL(wp) ::   zsatur, zsatur2, znusil, zdep, zfactdep
      REAL(wp) ::   zbactfer, zorem, zorem2, zofer
      REAL(wp) ::   zosil, zdenom1, zscave, zaggdfe, zcoag
#if ! defined key_kriest
      REAL(wp) ::   zofer2, zdenom, zdenom2
#endif
      REAL(wp) ::   zlamfac, zonitr, zstep
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:  ) :: ztempbac 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zdepbac, zolimi, zolimi2
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_rem')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj,      ztempbac                 )
      CALL wrk_alloc( jpi, jpj, jpk, zdepbac, zolimi, zolimi2 )

       ! Initialisation of temprary arrys
       zdepbac (:,:,:) = 0._wp
       zolimi  (:,:,:) = 0._wp
       zolimi2 (:,:,:) = 0._wp
       ztempbac(:,:)   = 0._wp

      !  Computation of the mean phytoplankton concentration as
      !  a crude estimate of the bacterial biomass
      !   --------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep = MAX( hmld(ji,jj), heup(ji,jj) )
               IF( fsdept(ji,jj,jk) < zdep ) THEN
                  zdepbac(ji,jj,jk) = MIN( 0.7 * ( trn(ji,jj,jk,jpzoo) + 2.* trn(ji,jj,jk,jpmes) ), 4.e-6 )
                  ztempbac(ji,jj)   = zdepbac(ji,jj,jk)
               ELSE
                  zdepbac(ji,jj,jk) = MIN( 1., zdep / fsdept(ji,jj,jk) ) * ztempbac(ji,jj)
               ENDIF
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! denitrification factor computed from O2 levels
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.e-6  - trn(ji,jj,jk,jpoxy) )    &
                  &                                / ( oxymin + trn(ji,jj,jk,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               ! DOC ammonification. Depends on depth, phytoplankton biomass
               !     and a limitation term which is supposed to be a parameterization
               !     of the bacterial activity. 
               zremik = xremik * zstep / 1.e-6 * xlimbac(ji,jj,jk) * zdepbac(ji,jj,jk) 
               zremik = MAX( zremik, 2.e-4 * xstep )
               !     Ammonification in oxic waters with oxygen consumption
               !     -----------------------------------------------------
               zolimi (ji,jj,jk) = zremik * ( 1.- nitrfac(ji,jj,jk) ) * trn(ji,jj,jk,jpdoc) 
               zolimi2(ji,jj,jk) = MIN( ( trn(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimi(ji,jj,jk) ) 
               !     Ammonification in suboxic waters with denitrification
               !     -------------------------------------------------------
               denitr(ji,jj,jk)  = MIN(  ( trn(ji,jj,jk,jpno3) - rtrn ) / rdenit,   &
                  &                     zremik * nitrfac(ji,jj,jk) * trn(ji,jj,jk,jpdoc)  )
               !
               zolimi (ji,jj,jk) = MAX( 0.e0, zolimi (ji,jj,jk) )
               zolimi2(ji,jj,jk) = MAX( 0.e0, zolimi2(ji,jj,jk) )
               denitr (ji,jj,jk) = MAX( 0.e0, denitr (ji,jj,jk) )
               !
            END DO
         END DO
      END DO


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !    NH4 nitrification to NO3. Ceased for oxygen concentrations
               !    below 2 umol/L. Inhibited at strong light 
               !    ----------------------------------------------------------
               zonitr  =nitrif * zstep * trn(ji,jj,jk,jpnh4) / ( 1.+ emoy(ji,jj,jk) ) * ( 1.- nitrfac(ji,jj,jk) ) 
               denitnh4(ji,jj,jk) = nitrif * zstep * trn(ji,jj,jk,jpnh4) * nitrfac(ji,jj,jk) 
               !   Update of the tracers trends
               !   ----------------------------
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zonitr - denitnh4(ji,jj,jk)
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zonitr - rdenita * denitnh4(ji,jj,jk)
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2nit * zonitr
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2 * rno3 * zonitr + rno3 * ( rdenita - 1. ) * denitnh4(ji,jj,jk)
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               !    Bacterial uptake of iron. No iron is available in DOC. So
               !    Bacteries are obliged to take up iron from the water. Some
               !    studies (especially at Papa) have shown this uptake to be significant
               !    ----------------------------------------------------------
               zdepbac2 = zdepbac(ji,jj,jk) * zdepbac(ji,jj,jk)
               zbactfer = 20.e-6 *  rfact2 * prmax(ji,jj,jk)                                 &
                  &              * trn(ji,jj,jk,jpfer) / ( concfebac + trn(ji,jj,jk,jpfer) )    & 
                  &              * trn(ji,jj,jk,jpfer) / ( 5E-10 + trn(ji,jj,jk,jpfer) )    &
                  &              * zdepbac2 / ( xkgraz2 + zdepbac(ji,jj,jk) )               &
                  &              * ( 0.5 + SIGN( 0.5, trn(ji,jj,jk,jpfer) -2.e-11 )  )

               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zbactfer
#if defined key_kriest
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zbactfer
#else
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zbactfer
#endif
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !    POC disaggregation by turbulence and bacterial activity. 
               !    -------------------------------------------------------------
               zremip = xremip * zstep * tgfunc(ji,jj,jk) * ( 1.- 0.7 * nitrfac(ji,jj,jk) ) 

               !    POC disaggregation rate is reduced in anoxic zone as shown by
               !    sediment traps data. In oxic area, the exponent of the martin s
               !    law is around -0.87. In anoxic zone, it is around -0.35. This
               !    means a disaggregation constant about 0.5 the value in oxic zones
               !    -----------------------------------------------------------------
               zorem  = zremip * trn(ji,jj,jk,jppoc)
               zofer  = zremip * trn(ji,jj,jk,jpsfe)
#if ! defined key_kriest
               zorem2 = zremip * trn(ji,jj,jk,jpgoc)
               zofer2 = zremip * trn(ji,jj,jk,jpbfe)
#else
               zorem2 = zremip * trn(ji,jj,jk,jpnum)
#endif

               !  Update the appropriate tracers trends
               !  -------------------------------------

               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zorem
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zofer
#if defined key_kriest
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zorem
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) - zorem2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zofer
#else
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zorem2 - zorem
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) - zorem2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zofer2 - zofer
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) - zofer2
#endif

            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !     Remineralization rate of BSi depedant on T and saturation
               !     ---------------------------------------------------------
               zsatur   = ( sio3eq(ji,jj,jk) - trn(ji,jj,jk,jpsil) ) / ( sio3eq(ji,jj,jk) + rtrn )
               zsatur   = MAX( rtrn, zsatur )
               zsatur2  = zsatur * ( 1. + tsn(ji,jj,jk,jp_tem) / 400.)**4
               znusil   = 0.225  * ( 1. + tsn(ji,jj,jk,jp_tem) / 15.) * zsatur + 0.775 * zsatur2**9.25
               zdep     = MAX( hmld(ji,jj), heup(ji,jj) ) 
               zdep     = MAX( 0., fsdept(ji,jj,jk) - zdep )
               zfactdep = xsilab * EXP(-( xsiremlab - xsirem ) * zdep / wsbio2 )
               zsiremin = ( xsiremlab * zfactdep + xsirem * ( 1. - zfactdep ) ) * zstep * znusil
               zosil    = zsiremin * trn(ji,jj,jk,jpdsi)
               !
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zosil
               tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) + zosil
               !
            END DO
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem4')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      zfesatur = ligand
!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !  Compute de different ratios for scavenging of iron
               !  --------------------------------------------------

#if  defined key_kriest
               zdenom1 = trn(ji,jj,jk,jppoc) / &
           &           ( trn(ji,jj,jk,jppoc) + trn(ji,jj,jk,jpdsi) + trn(ji,jj,jk,jpcal) + rtrn )
#else
               zdenom = 1. / ( trn(ji,jj,jk,jppoc) + trn(ji,jj,jk,jpgoc) + trn(ji,jj,jk,jpdsi) + trn(ji,jj,jk,jpcal) + rtrn )
               zdenom1 = trn(ji,jj,jk,jppoc) * zdenom
               zdenom2 = trn(ji,jj,jk,jpgoc) * zdenom
#endif
               !  scavenging rate of iron. this scavenging rate depends on the load in particles
               !  on which they are adsorbed. The  parameterization has been taken from studies on Th
               !     ------------------------------------------------------------
               zkeq = fekeq(ji,jj,jk)
               zfeequi = ( -( 1. + zfesatur * zkeq - zkeq * trn(ji,jj,jk,jpfer) )               &
                  &        + SQRT( ( 1. + zfesatur * zkeq - zkeq * trn(ji,jj,jk,jpfer) )**2       &
                  &               + 4. * trn(ji,jj,jk,jpfer) * zkeq) ) / ( 2. * zkeq )

#if defined key_kriest
               zlam1b = 3.e-5 + xlam1 * (  trn(ji,jj,jk,jppoc)                   &
                  &                      + trn(ji,jj,jk,jpcal) + trn(ji,jj,jk,jpdsi)  ) * 1.e6
#else
               zlam1b = 3.e-5 + xlam1 * (  trn(ji,jj,jk,jppoc) + trn(ji,jj,jk,jpgoc)   &
                  &                      + trn(ji,jj,jk,jpcal) + trn(ji,jj,jk,jpdsi)  ) * 1.e6
#endif
               zscave = zfeequi * zlam1b * zstep

               !  Increased scavenging for very high iron concentrations
               !  found near the coasts due to increased lithogenic particles
               !  and let say it is unknown processes (precipitation, ...)
               !  -----------------------------------------------------------
               zlam1b  = xlam1 * MAX( 0.e0, ( trn(ji,jj,jk,jpfer) * 1.e9 - 1. ) )
               zcoag   = zfeequi * zlam1b * zstep
               zlamfac = MAX( 0.e0, ( gphit(ji,jj) + 55.) / 30. )
               zlamfac = MIN( 1.  , zlamfac )
               zdep    =  MIN(1., 1000. / fsdept(ji,jj,jk) )
#if ! defined key_kriest
               zlam1b = (  80.* ( trn(ji,jj,jk,jpdoc) + 35.e-6 )                           &
                  &     + 698.*   trn(ji,jj,jk,jppoc) + 1.05e4 * trn(ji,jj,jk,jpgoc)  )    &
                  &   * xdiss(ji,jj,jk) + 1E-4 * ( 1. - zlamfac ) * zdep
#else
               zlam1b = (  80.* (trn(ji,jj,jk,jpdoc) + 35E-6)              &
                  &     + 698.*  trn(ji,jj,jk,jppoc)  )                    &
                  &   * xdiss(ji,jj,jk) + 1E-4 * ( 1. - zlamfac ) * zdep
#endif
               zaggdfe = zlam1b * zstep * 0.5 * ( trn(ji,jj,jk,jpfer) - zfeequi )
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zscave - zaggdfe - zcoag
#if defined key_kriest
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * zdenom1
#else
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * zdenom1
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zscave * zdenom2
#endif
            END DO
         END DO
      END DO
      !

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem5')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------

      DO jk = 1, jpkm1
         tra(:,:,jk,jppo4) = tra(:,:,jk,jppo4) + zolimi (:,:,jk) + denitr(:,:,jk)
         tra(:,:,jk,jpnh4) = tra(:,:,jk,jpnh4) + zolimi (:,:,jk) + denitr(:,:,jk)
         tra(:,:,jk,jpno3) = tra(:,:,jk,jpno3) - denitr (:,:,jk) * rdenit
         tra(:,:,jk,jpdoc) = tra(:,:,jk,jpdoc) - zolimi (:,:,jk) - denitr(:,:,jk)
         tra(:,:,jk,jpoxy) = tra(:,:,jk,jpoxy) - zolimi2(:,:,jk) * o2ut
         tra(:,:,jk,jpdic) = tra(:,:,jk,jpdic) + zolimi (:,:,jk) + denitr(:,:,jk)
         tra(:,:,jk,jptal) = tra(:,:,jk,jptal) + rno3 * ( zolimi(:,:,jk) + ( rdenit + 1.) * denitr(:,:,jk) )
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem6')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj,      ztempbac                 )
      CALL wrk_dealloc( jpi, jpj, jpk, zdepbac, zolimi, zolimi2 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ xremik, xremip, nitrif, xsirem, xsiremlab, xsilab,   &
      &                   xlam1, oxymin, ligand 

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for remineralization, nampisrem'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    remineralisation rate of POC              xremip    =', xremip
         WRITE(numout,*) '    remineralization rate of DOC              xremik    =', xremik
         WRITE(numout,*) '    remineralization rate of Si               xsirem    =', xsirem
         WRITE(numout,*) '    fast remineralization rate of Si          xsiremlab =', xsiremlab
         WRITE(numout,*) '    fraction of labile biogenic silica        xsilab    =', xsilab
         WRITE(numout,*) '    scavenging rate of Iron                   xlam1     =', xlam1
         WRITE(numout,*) '    NH4 nitrification rate                    nitrif    =', nitrif
         WRITE(numout,*) '    halk saturation constant for anoxia       oxymin    =', oxymin
         WRITE(numout,*) '    ligand concentration in the ocean         ligand    =', ligand
      ENDIF
      !
      nitrfac (:,:,:) = 0._wp
      denitr  (:,:,:) = 0._wp
      denitnh4(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( denitr(jpi,jpj,jpk), denitnh4(jpi,jpj,jpk), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_warn('p4z_rem_alloc: failed to allocate arrays')
      !
   END FUNCTION p4z_rem_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_rem                    ! Empty routine
   END SUBROUTINE p4z_rem
#endif 

   !!======================================================================
END MODULE p4zrem
