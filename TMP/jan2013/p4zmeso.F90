MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_meso       :   Compute the sources/sinks for mesozooplankton
   !!   p4z_meso_init  :   Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zint          !  interpolation and computation of various fields
   USE p4zprod         !  production
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part2       = 0.5_wp     !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xprefc      = 1.0_wp     !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xprefp      = 0.3_wp     !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xprefz      = 1.0_wp     !: mesozoo preference for diatoms
   REAL(wp), PUBLIC ::  xprefpoc    = 0.3_wp     !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xthresh2zoo = 1E-8_wp    !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia = 1E-8_wp    !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy = 2E-7_wp    !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc = 1E-8_wp    !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2    = 0._wp      !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2     = 0.005_wp   !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2      = 0.04_wp    !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2    = 0.9_wp     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2     = 20E-6_wp   !: non assimilated fraction of P by mesozoo 
   REAL(wp), PUBLIC ::  unass2      = 0.3_wp     !: Efficicency of mesozoo growth 
   REAL(wp), PUBLIC ::  sigma2      = 0.6_wp     !: Fraction of mesozoo excretion as DOM 
   REAL(wp), PUBLIC ::  epsher2     = 0.3_wp     !: half sturation constant for grazing 2
   REAL(wp), PUBLIC ::  grazflux    = 3.E3_wp    !: mesozoo flux feeding rate

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmeso.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_meso( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam
      REAL(wp) :: zgraze2 , zdenom, zdenom2, zncratio
      REAL(wp) :: zfact   , zstep, zfood, zfoodlim
      REAL(wp) :: zepshert, zepsherv, zgrarsig, zgraztot, zgraztotf
      REAL(wp) :: zgrarem2, zgrafer2, zgrapoc2, zprcaca, zmortz2, zgrasrat
#if defined key_kriest
      REAL znumpoc
#endif
      REAL(wp) :: zrespz2, ztortz2, zgrazd, zgrazz, zgrazpof
      REAL(wp) :: zgrazn, zgrazpoc, zgraznf, zgrazf
      REAL(wp) :: zgrazfff, zgrazffe
      CHARACTER (len=25) :: charout
      REAL(wp) :: zrfact2
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_meso')
      !
         zrfact2 = 1.e3 * rfact2r
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompam   = MAX( ( trn(ji,jj,jk,jpmes) - 1.e-8 ), 0.e0 )
# if defined key_degrad
               zstep     = xstep * facvol(ji,jj,jk)
# else
               zstep     = xstep
# endif
               zfact     = zstep * tgfunc(ji,jj,jk) * zcompam

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz2   = resrat2 * zfact * trn(ji,jj,jk,jpmes) / ( xkmort + trn(ji,jj,jk,jpmes) )  &
                  &      + resrat2 * zfact * 3. * nitrfac(ji,jj,jk)

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation
               !  ---------------------------------------------------------------
               ztortz2   = mzrat2 * 1.e6 * zfact * trn(ji,jj,jk,jpmes)
               !

               zcompadi  = MAX( ( trn(ji,jj,jk,jpdia) - xthresh2dia ), 0.e0 )
               zcompaz   = MAX( ( trn(ji,jj,jk,jpzoo) - xthresh2zoo ), 0.e0 )
               zcompaph  = MAX( ( trn(ji,jj,jk,jpphy) - xthresh2phy ), 0.e0 )
               zcompapoc = MAX( ( trn(ji,jj,jk,jppoc) - xthresh2poc ), 0.e0 )

               zfood     = xprefc * zcompadi + xprefz * zcompaz + xprefp * zcompaph + xprefpoc * zcompapoc 
               zfoodlim  = MAX( 0., zfood - min(0.5*zfood,xthresh2) )
               zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
               zdenom2   = zdenom / ( zfood + rtrn )
               zgraze2   = grazrat2 * zstep * tgfunc2(ji,jj,jk) * trn(ji,jj,jk,jpmes) 

               zgrazd    = zgraze2  * xprefc   * zcompadi  * zdenom2 
               zgrazz    = zgraze2  * xprefz   * zcompaz   * zdenom2 
               zgrazn    = zgraze2  * xprefp   * zcompaph  * zdenom2 
               zgrazpoc  = zgraze2  * xprefpoc * zcompapoc * zdenom2 

               zgraznf   = zgrazn   * trn(ji,jj,jk,jpnfe) / ( trn(ji,jj,jk,jpphy) + rtrn)
               zgrazf    = zgrazd   * trn(ji,jj,jk,jpdfe) / ( trn(ji,jj,jk,jpdia) + rtrn)
               zgrazpof  = zgrazpoc * trn(ji,jj,jk,jpsfe) / ( trn(ji,jj,jk,jppoc) + rtrn)

               !  Mesozooplankton flux feeding on GOC
               !  ----------------------------------
# if ! defined key_kriest
               zgrazffe  = grazflux * zstep * wsbio4(ji,jj,jk)          &
                 &                 * tgfunc2(ji,jj,jk) * trn(ji,jj,jk,jpgoc) * trn(ji,jj,jk,jpmes)
               zgrazfff  = zgrazffe * trn(ji,jj,jk,jpbfe) / (trn(ji,jj,jk,jpgoc) + rtrn)
# else
               zgrazffe = grazflux * zstep * wsbio3(ji,jj,jk)     &
               zgrazfff   = zgrazffe * trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)
# endif
              !
              zgraztot   = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffe
              zgraztotf  = zgrazf + zgraznf + zgrazz * ferat3 + zgrazpof + zgrazfff 

              ! Total grazing ( grazing by microzoo is already computed in p4zmicro )
              grazing(ji,jj,jk) = grazing(ji,jj,jk) + zgraztot
              !    Mesozooplankton efficiency
              !    --------------------------
              zgrasrat   =  zgraztotf / ( zgraztot + rtrn )
              zncratio   = (  xprefc   * zcompadi * quotad(ji,jj,jk)  &
                  &         + xprefp   * zcompaph * quotan(ji,jj,jk)  &
                  &         + xprefz   * zcompaz                      &
                  &         + xprefpoc * zcompapoc   ) / ( zfood + rtrn )
               zepshert  = epsher2 * MIN( 1., zncratio )
               zepsherv  = zepshert * MIN( 1., zgrasrat / ferat3 )
               zgrarem2  = zgraztot * ( 1. - zepsherv - unass2 )
               zgrafer2  = zgraztot * MAX( 0. , ( 1. - unass2 ) * zgrasrat - ferat3 * zepshert ) 
               zgrapoc2  = zgraztot * unass2

               !   Update the arrays TRA which contain the biological sources and sinks
               zgrarsig  = zgrarem2 * sigma2
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarsig
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarsig
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem2 - zgrarsig
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarsig
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarsig
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgrarsig              
#if defined key_kriest
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc2
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + zgrapoc2 * xkr_dmeso
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zgraztotf * unass2
#else
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zgrapoc2
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zgraztotf * unass2
#endif
               zmortz2 = ztortz2 + zrespz2
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) - zmortz2 + zepsherv * zgraztot 
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazd
!Joan
		grazDM(ji,jj,jk) = zgrazd
!end
               diagraz(ji,jj,jk) = diagraz(ji,jj,jk) + ( zgrazd * zrfact2)
               nanograz(ji,jj,jk) = nanograz(ji,jj,jk) + ( zgrazn * zrfact2)
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zgrazz
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazn
!Joan
		grazNM(ji,jj,jk) = zgrazn
!end
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazn * trn(ji,jj,jk,jpnch) / ( trn(ji,jj,jk,jpphy) + rtrn )
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazd * trn(ji,jj,jk,jpdch) / ( trn(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpbsi) = tra(ji,jj,jk,jpbsi) - zgrazd * trn(ji,jj,jk,jpbsi) / ( trn(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zgrazd * trn(ji,jj,jk,jpbsi) / ( trn(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgraznf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazf

               zprcaca = xfracal(ji,jj,jk) * zgrazn
               ! calcite production
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zprcaca = part2 * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
#if defined key_kriest
               znumpoc = trn(ji,jj,jk,jpnum) / ( trn(ji,jj,jk,jppoc) + rtrn )
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortz2 - zgrazpoc - zgrazffe
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) - zgrazpoc * znumpoc &
                  &    + zmortz2  * xkr_dmeso - zgrazffe * znumpoc * wsbio4(ji,jj,jk) / ( wsbio3(ji,jj,jk) + rtrn )
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * zmortz2 - zgrazfff - zgrazpof
#else
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zgrazpoc
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zmortz2 - zgrazffe
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zgrazpof
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + ferat3 * zmortz2 - zgrazfff
#endif

            END DO
         END DO
      END DO
      !
      IF( ln_diatrc .AND. lk_iomput ) THEN
         zrfact2 = 1.e3 * rfact2r
         grazing(:,:,:) = grazing(:,:,:) * zrfact2 * tmask(:,:,:)   ! Total grazing of phyto by zoo
         prodcal(:,:,:) = prodcal(:,:,:) * zrfact2 * tmask(:,:,:)   ! Calcite production	
	 grazNM(:,:,:) = grazNM(:,:,:) * zrfact2 * tmask(:,:,:)
	 grazDM(:,:,:) = grazDM(:,:,:) * zrfact2 * tmask(:,:,:)
         IF( jnt == nrdttrc ) THEN
            CALL iom_put( "GRAZ" , grazing  )  ! Total grazing of phyto by zooplankton
            CALL iom_put( "PCAL" , prodcal  )  ! Calcite production
            CALL iom_put( "GRAZD", diagraz  )
	    CALL iom_put( "GRAZN", nanograz )
	    CALL iom_put( "GRAZNM", grazNM  )
	    CALL iom_put( "GRAZDM", grazDM  )
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_meso')
      !
   END SUBROUTINE p4z_meso

   SUBROUTINE p4z_meso_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismes/ part2, grazrat2, resrat2, mzrat2, xprefc, xprefp, xprefz,   &
         &                xprefpoc, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, &
         &                xthresh2, xkgraz2, epsher2, sigma2, unass2, grazflux

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampismes )


      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, nampismes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(numout,*) '    mesozoo preference for phyto                   xprefc       =', xprefc
         WRITE(numout,*) '    mesozoo preference for POC                     xprefp       =', xprefp
         WRITE(numout,*) '    mesozoo preference for zoo                     xprefz       =', xprefz
         WRITE(numout,*) '    mesozoo preference for poc                     xprefpoc     =', xprefpoc
         WRITE(numout,*) '    microzoo feeding threshold  for mesozoo        xthresh2zoo  =', xthresh2zoo
         WRITE(numout,*) '    diatoms feeding threshold  for mesozoo         xthresh2dia  =', xthresh2dia
         WRITE(numout,*) '    nanophyto feeding threshold for mesozoo        xthresh2phy  =', xthresh2phy
         WRITE(numout,*) '    poc feeding threshold for mesozoo              xthresh2poc  =', xthresh2poc
         WRITE(numout,*) '    feeding threshold for mesozooplankton          xthresh2     =', xthresh2
         WRITE(numout,*) '    exsudation rate of mesozooplankton             resrat2      =', resrat2
         WRITE(numout,*) '    mesozooplankton mortality rate                 mzrat2       =', mzrat2
         WRITE(numout,*) '    maximal mesozoo grazing rate                   grazrat2     =', grazrat2
         WRITE(numout,*) '    mesozoo flux feeding rate                      grazflux     =', grazflux
         WRITE(numout,*) '    non assimilated fraction of P by mesozoo       unass2       =', unass2
         WRITE(numout,*) '    Efficicency of Mesozoo growth                  epsher2      =', epsher2
         WRITE(numout,*) '    Fraction of mesozoo excretion as DOM           sigma2       =', sigma2
         WRITE(numout,*) '    half sturation constant for grazing 2          xkgraz2      =', xkgraz2
      ENDIF


   END SUBROUTINE p4z_meso_init


#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_meso                    ! Empty routine
   END SUBROUTINE p4z_meso
#endif 

   !!======================================================================
END MODULE  p4zmeso
