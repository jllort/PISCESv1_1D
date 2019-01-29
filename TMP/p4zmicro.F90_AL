MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_micro       :   Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init  :   Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim          !  Co-limitations
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zint          !  interpolation and computation of various fields
   USE p4zprod         !  production
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part       = 0.5_wp     !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::  xpref2c    = 0.2_wp     !: microzoo preference for POC 
   REAL(wp), PUBLIC ::  xpref2p    = 1.0_wp     !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::  xpref2d    = 0.6_wp     !: microzoo preference for diatoms
   REAL(wp), PUBLIC ::  xthreshdia = 1E-8_wp    !: diatoms feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshphy = 2E-7_wp    !: nanophyto threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshpoc = 1E-8_wp    !: poc threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthresh    = 0._wp      !: feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::  resrat     = 0.03_wp    !: exsudation rate of microzooplankton
   REAL(wp), PUBLIC ::  mzrat      = 0.0_wp     !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat    = 3.0_wp     !: maximal microzoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz     = 20E-6_wp   !: non assimilated fraction of P by microzoo 
   REAL(wp), PUBLIC ::  unass      = 0.3_wp     !: Efficicency of microzoo growth 
   REAL(wp), PUBLIC ::  sigma1     = 0.6_wp     !: Fraction of microzoo excretion as DOM 
   REAL(wp), PUBLIC ::  epsher     = 0.3_wp     !: half sturation constant for grazing 1 


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmicro.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_micro( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::  kt  ! ocean time step
      INTEGER, INTENT(in) ::  jnt 
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc
      REAL(wp) :: zgraze  , zdenom, zdenom2, zncratio
      REAL(wp) :: zfact   , zstep, zfood, zfoodlim
      REAL(wp) :: zepshert, zepsherv, zgrarsig, zgraztot, zgraztotf
      REAL(wp) :: zgrarem, zgrafer, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz, zgrasrat
      REAL(wp) :: zgrazp, zgrazm, zgrazsd
      REAL(wp) :: zgrazmf, zgrazsf, zgrazpf
      REAL(wp) :: zrfact2
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zgrazing
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_micro')
      !
      IF( ln_diatrc .AND. lk_iomput )  CALL wrk_alloc( jpi, jpj, jpk, zgrazing )
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaz = MAX( ( trn(ji,jj,jk,jpzoo) - 1.e-8 ), 0.e0 )
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               zfact   = zstep * tgfunc2(ji,jj,jk) * zcompaz

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz = resrat * zfact * trn(ji,jj,jk,jpzoo) / ( xkmort + trn(ji,jj,jk,jpzoo) )  &
                  &   + resrat * zfact * 3. * nitrfac(ji,jj,jk)

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation.
               !  ---------------------------------------------------------------
               ztortz = mzrat * 1.e6 * zfact * trn(ji,jj,jk,jpzoo)

               zcompadi  = MIN( MAX( ( trn(ji,jj,jk,jpdia) - xthreshdia ), 0.e0 ), xsizedia )
               zcompaph  = MAX( ( trn(ji,jj,jk,jpphy) - xthreshphy ), 0.e0 )
               zcompapoc = MAX( ( trn(ji,jj,jk,jppoc) - xthreshpoc ), 0.e0 )
               
               !     Microzooplankton grazing
               !     ------------------------
               zfood     = xpref2p * zcompaph + xpref2c * zcompapoc + xpref2d * zcompadi
               zfoodlim  = MAX( 0. , zfood - min(xthresh,0.5*zfood) )
               zdenom    = zfoodlim / ( xkgraz + zfoodlim )
               zdenom2   = zdenom / ( zfood + rtrn )
               zgraze    = grazrat * zstep * tgfunc2(ji,jj,jk) * trn(ji,jj,jk,jpzoo) 

               zgrazp    = zgraze  * xpref2p * zcompaph  * zdenom2 
               zgrazm    = zgraze  * xpref2c * zcompapoc * zdenom2 
               zgrazsd   = zgraze  * xpref2d * zcompadi  * zdenom2 

               zgrazpf   = zgrazp  * trn(ji,jj,jk,jpnfe) / (trn(ji,jj,jk,jpphy) + rtrn)
               zgrazmf   = zgrazm  * trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)
               zgrazsf   = zgrazsd * trn(ji,jj,jk,jpdfe) / (trn(ji,jj,jk,jpdia) + rtrn)
               !
               zgraztot  = zgrazp  + zgrazm  + zgrazsd 
               zgraztotf = zgrazpf + zgrazsf + zgrazmf 

               ! Grazing by microzooplankton
               IF( ln_diatrc .AND. lk_iomput )  zgrazing(ji,jj,jk) = zgraztot

               !    Various remineralization and excretion terms
               !    --------------------------------------------
               zgrasrat  = zgraztotf / ( zgraztot + rtrn )
               zncratio  = ( xpref2p * zcompaph * quotan(ji,jj,jk) &
                  &        + xpref2d * zcompadi * quotad(ji,jj,jk) + xpref2c * zcompapoc ) / ( zfood + rtrn )
               zepshert  = epsher * MIN( 1., zncratio )
               zepsherv  = zepshert * MIN( 1., zgrasrat / ferat3 )
               zgrafer   = zgraztot * MAX( 0. , ( 1. - unass ) * zgrasrat - ferat3 * zepsherv ) 
               zgrarem   = zgraztot * ( 1. - zepsherv - unass )
               zgrapoc   = zgraztot * unass

               !  Update of the TRA arrays
               !  ------------------------
               zgrarsig  = zgrarem * sigma1
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarsig
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarsig
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem - zgrarsig
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarsig
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zgraztotf * unass
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarsig
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgrarsig
#if defined key_kriest
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + zgrapoc * xkr_dmicro
#endif
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------
               zmortz = ztortz + zrespz
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zmortz + zepsherv * zgraztot 
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazp
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazsd
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazp  * trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)+rtrn)
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazsd * trn(ji,jj,jk,jpdch)/(trn(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazsd * trn(ji,jj,jk,jpdsi)/(trn(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazsd * trn(ji,jj,jk,jpdsi)/(trn(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgrazpf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazsf
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortz - zgrazm
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * zmortz - zgrazmf
               zprcaca = xfracal(ji,jj,jk) * zgrazp
               !
               ! calcite production
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zprcaca = part * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
#if defined key_kriest
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ( zmortz - zgrazm ) * xkr_dmicro
#endif
            END DO
         END DO
      END DO
      !
      IF( ln_diatrc .AND. lk_iomput .AND. jnt == nrdttrc ) THEN
         zrfact2 = 1.e3 * rfact2r
         CALL iom_put( "GRAZ1" , zgrazing(:,:,:) * zrfact2 * tmask(:,:,:) )  ! Total grazing of phyto by zooplankton
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_diatrc .AND. lk_iomput )  CALL wrk_dealloc( jpi, jpj, jpk, zgrazing )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_micro')
      !
   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampiszoo/ part, grazrat, resrat, mzrat, xpref2c, xpref2p, &
         &                xpref2d,  xthreshdia,  xthreshphy,  xthreshpoc, &
         &                xthresh, xkgraz, epsher, sigma1, unass

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampiszoo )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton, nampiszoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(numout,*) '    microzoo preference for POC                     xpref2c     =', xpref2c
         WRITE(numout,*) '    microzoo preference for nano                    xpref2p     =', xpref2p
         WRITE(numout,*) '    microzoo preference for diatoms                 xpref2d     =', xpref2d
         WRITE(numout,*) '    diatoms feeding threshold  for microzoo         xthreshdia  =', xthreshdia
         WRITE(numout,*) '    nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(numout,*) '    poc feeding threshold for microzoo              xthreshpoc  =', xthreshpoc
         WRITE(numout,*) '    feeding threshold for microzooplankton          xthresh     =', xthresh
         WRITE(numout,*) '    exsudation rate of microzooplankton             resrat      =', resrat
         WRITE(numout,*) '    microzooplankton mortality rate                 mzrat       =', mzrat
         WRITE(numout,*) '    maximal microzoo grazing rate                   grazrat     =', grazrat
         WRITE(numout,*) '    non assimilated fraction of P by microzoo       unass       =', unass
         WRITE(numout,*) '    Efficicency of microzoo growth                  epsher      =', epsher
         WRITE(numout,*) '    Fraction of microzoo excretion as DOM           sigma1      =', sigma1
         WRITE(numout,*) '    half sturation constant for grazing 1           xkgraz      =', xkgraz
      ENDIF

   END SUBROUTINE p4z_micro_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_micro                    ! Empty routine
   END SUBROUTINE p4z_micro
#endif 

   !!======================================================================
END MODULE  p4zmicro
