MODULE trcini_pisces
   !!======================================================================
   !!                         ***  MODULE trcini_pisces  ***
   !! TOP :   initialisation of the PISCES biochemical model
   !!======================================================================
   !! History :    -   !  1988-07  (E. Maier-Reiner) Original code
   !!              -   !  1999-10  (O. Aumont, C. Le Quere)
   !!              -   !  2002     (O. Aumont)  PISCES
   !!             1.0  !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcini.pisces.h90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_ini_pisces   : PISCES biochemical model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zche          !  Chemical model
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zrem          !  Remineralisation of organic matter
   USE p4zflx          !  Gas exchange
   USE p4zsed          !  Sedimentation
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmicro        !  Sources and sinks of microzooplankton
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p4zmort         !  Mortality terms for phytoplankton
   USE p4zlys          !  Calcite saturation
   USE p4zsed          !  Sedimentation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_pisces   ! called by trcini.F90 module

   REAL(wp) :: sco2   =  2.312e-3_wp
   REAL(wp) :: alka0  =  2.423e-3_wp
   REAL(wp) :: oxyg0  =  177.6e-6_wp 
   REAL(wp) :: po4    =  2.174e-6_wp 
   REAL(wp) :: bioma0 =  1.000e-8_wp  
   REAL(wp) :: silic1 =  91.65e-6_wp  
   REAL(wp) :: no3    =  31.04e-6_wp * 7.625_wp

#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_pisces.F90 3192 2011-12-05 13:55:12Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_pisces
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ini_pisces ***
      !!
      !! ** Purpose :   Initialisation of the PISCES biochemical model
      !!----------------------------------------------------------------------
      !
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_pisces :   PISCES biochemical model initialisation'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      CALL pisces_alloc()                          ! Allocate PISCES arrays

      !                                            ! Time-step
      rfact   = rdttrc(1)                          ! ---------
      rfactr  = 1. / rfact
      rfact2  = rfact / FLOAT( nrdttrc )
      rfact2r = 1. / rfact2

      IF(lwp) WRITE(numout,*) '    Passive Tracer  time step    rfact  = ', rfact, ' rdt = ', rdttra(1)
      IF(lwp) write(numout,*) '    PISCES  Biology time step    rfact2 = ', rfact2



      ! Set biological ratios
      ! ---------------------
      rno3    =  16._wp / 122._wp
      po4r    =   1._wp / 122._wp
      o2nit   =  32._wp / 122._wp
      rdenit  = 105._wp /  16._wp
      rdenita =   3._wp /  5._wp
      o2ut    = 131._wp / 122._wp

      CALL p4z_che        ! initialize the chemical constants

      ! Initialization of tracer concentration in case of  no restart 
      !--------------------------------------------------------------
      IF( .NOT. ln_rsttr ) THEN  
         
         trn(:,:,:,jpdic) = sco2
         trn(:,:,:,jpdoc) = bioma0
         trn(:,:,:,jptal) = alka0
         trn(:,:,:,jpoxy) = oxyg0
         trn(:,:,:,jpcal) = bioma0
         trn(:,:,:,jppo4) = po4 / po4r
         trn(:,:,:,jppoc) = bioma0
#  if ! defined key_kriest
         trn(:,:,:,jpgoc) = bioma0
         trn(:,:,:,jpbfe) = bioma0 * 5.e-6
#  else
         trn(:,:,:,jpnum) = bioma0 / ( 6. * xkr_massp )
#  endif
         trn(:,:,:,jpsil) = silic1
         trn(:,:,:,jpbsi) = bioma0 * 0.15
         trn(:,:,:,jpdsi) = bioma0 * 5.e-6
         trn(:,:,:,jpphy) = bioma0
         trn(:,:,:,jpdia) = bioma0
         trn(:,:,:,jpzoo) = bioma0
         trn(:,:,:,jpmes) = bioma0
         trn(:,:,:,jpfer) = 0.6E-9
         trn(:,:,:,jpsfe) = bioma0 * 5.e-6
         trn(:,:,:,jpdfe) = bioma0 * 5.e-6
         trn(:,:,:,jpnfe) = bioma0 * 5.e-6
         trn(:,:,:,jpnch) = bioma0 * 12. / 55.
         trn(:,:,:,jpdch) = bioma0 * 12. / 55.
         trn(:,:,:,jpno3) = no3
         trn(:,:,:,jpnh4) = bioma0

         ! initialize the half saturation constant for silicate
         ! ----------------------------------------------------
         xksi(:,:)    = 2.e-6
         xksimax(:,:) = xksi(:,:)

      ENDIF

      IF( .NOT. ln_rsttr ) THEN
         ! Initialization of chemical variables of the carbon cycle
         ! --------------------------------------------------------
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztmas   = tmask(ji,jj,jk)
                  ztmas1  = 1. - tmask(ji,jj,jk)
                  zcaralk = trn(ji,jj,jk,jptal) - borat(ji,jj,jk) / (  1. + 1.E-8 / ( rtrn + akb3(ji,jj,jk) )  )
                  zco3    = ( zcaralk - trn(ji,jj,jk,jpdic) ) * ztmas + 0.5e-3 * ztmas1
                  zbicarb = ( 2. * trn(ji,jj,jk,jpdic) - zcaralk )
                  hi(ji,jj,jk) = ( ak23(ji,jj,jk) * zbicarb / zco3 ) * ztmas + 1.e-9 * ztmas1
               END DO
            END DO
         END DO
         !
      END IF

      ! Time step duration for biology
      xstep = rfact2 / rday

      CALL p4z_sink_init      !  vertical flux of particulate organic matter
      CALL p4z_opt_init       !  Optic: PAR in the water column
      CALL p4z_lim_init       !  co-limitations by the various nutrients
      CALL p4z_prod_init      !  phytoplankton growth rate over the global ocean.
      CALL p4z_rem_init       !  remineralisation
      CALL p4z_mort_init      !  phytoplankton mortality 
      CALL p4z_micro_init     !  microzooplankton
      CALL p4z_meso_init      !  mesozooplankton
      CALL p4z_sed_init       !  sedimentation 
      CALL p4z_lys_init       !  calcite saturation
      CALL p4z_flx_init       !  gas exchange 

      ndayflxtr = 0

      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'Initialization of PISCES tracers done'
      IF(lwp) WRITE(numout,*) 
      !
   END SUBROUTINE trc_ini_pisces


   SUBROUTINE pisces_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE pisces_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of PISCES 
      !!----------------------------------------------------------------------
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =         sms_pisces_alloc()          ! Start of PISCES-related alloc routines...
      ierr = ierr +  p4z_che_alloc()
      ierr = ierr +  p4z_sink_alloc()
      ierr = ierr +  p4z_opt_alloc()
      ierr = ierr +  p4z_prod_alloc()
      ierr = ierr +  p4z_rem_alloc()
      ierr = ierr +  p4z_sed_alloc()
      ierr = ierr +  p4z_flx_alloc()
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'pisces_alloc: unable to allocate PISCES arrays' )
      !
   END SUBROUTINE pisces_alloc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No PISCES biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_pisces             ! Empty routine
   END SUBROUTINE trc_ini_pisces
#endif

   !!======================================================================
END MODULE trcini_pisces
