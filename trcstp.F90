MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================
   !! History :  1.0  !  2004-03  (C. Ethe)  Original
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE trc
   USE trctrp           ! passive tracers transport
   USE trcsms           ! passive tracers sources and sinks
   USE prtctl_trc       ! Print control for debbuging
   USE trcdia
   USE trcwri
   USE trcrst
   USE trdmod_trc_oce
   USE trdmld_trc
   USE iom
   USE in_out_manager
   USE trcsub
   !Sakina
   USE domain
   USE ioipsl, ONLY :   ymds2ju   ! for calendar
   !Joan (to output tpp)
   USE sms_pisces
   USE p4zopt


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_stp    ! called by step


   !Sakina
      REAL(wp)   :: valmax,timemax,propmax !: Value, timing, and proportion of diatoms at max Chl
      REAL(wp)   :: maxdev,timedev  !: Timing of the maximal derivative of Chl
      INTEGER,PARAMETER :: period=5 !: Period over which the derivative is calculated (in time step, should be an odd number)
      REAL(wp),DIMENSION(period) :: totchl,timetotchl  !: Estimation of the timing of the maximal derivative of Chl
      REAL(wp)   :: zjul ! Julian day of the begining of the run
      REAL(wp) :: parmlddev, ironCUM
   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcstp.F90 3192 2011-12-05 13:55:12Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_stp( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt      ! ocean time-step index
      INTEGER               ::  jk, jn  ! dummy loop indices
      REAL(wp)              ::  ztrai
      CHARACTER (len=25)    ::  charout
      !!-------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_stp')
      !
      IF( kt == nittrc000 ) THEN
                               CALL iom_close( numrtr )     ! close input  passive tracers restart file
         IF( lk_trdmld_trc  )  CALL trd_mld_trc_init        ! trends: Mixed-layer
      ENDIF
      !
      IF( lk_vvl ) THEN                              ! update ocean volume due to ssh temporal evolution
         DO jk = 1, jpk
            cvol(:,:,jk) = e1e2t(:,:) * fse3t(:,:,jk) * tmask(:,:,jk)
         END DO
         IF( lk_degrad )  cvol(:,:,:) = cvol(:,:,:) * facvol(:,:,:)      ! degrad option: reduction by facvol
         areatot         = glob_sum( cvol(:,:,:) )
      ENDIF
      !    
     IF( nn_dttrc /= 1 )   CALL trc_sub_stp( kt )  ! averaging physical variables for sub-stepping

     IF( MOD( kt , nn_dttrc ) == 0 ) THEN      ! only every nn_dttrc time step
         !
         IF(ln_ctl) THEN
            WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
            CALL prt_ctl_trc_info(charout)
         ENDIF
         !
         tra(:,:,:,:) = 0.e0
         !
                                   CALL trc_rst_opn  ( kt )       ! Open tracer restart file 
         IF( lk_iomput ) THEN  ;   CALL trc_wri      ( kt )       ! output of passive tracers with iom I/O manager
         ELSE                  ;   CALL trc_dia      ( kt )       ! output of passive tracers with old I/O manager
         ENDIF
                                   CALL trc_sms      ( kt )       ! tracers: sinks and sources
                                   CALL trc_trp      ( kt )       ! transport of passive tracers
         IF( lrst_trc )            CALL trc_rst_wri  ( kt )       ! write tracer restart file
         IF( lk_trdmld_trc  )      CALL trd_mld_trc  ( kt )       ! trends: Mixed-layer
         !
         IF( nn_dttrc /= 1   )     CALL trc_sub_reset( kt )       ! resetting physical variables when sub-stepping
         !
      ! Sakina
#if defined key_pisces
                                   CALL trc_bloom    ( kt )       ! Diagnostics of the phytoplanktonic bloom 
#endif
      ENDIF
      !
      ztrai = 0._wp                                                   !  content of all tracers
      DO jn = 1, jptra
         ztrai = ztrai + glob_sum( trn(:,:,:,jn) * cvol(:,:,:)   )
      END DO
      IF( lwp ) WRITE(numstr,9300) kt,  ztrai / areatot
9300  FORMAT(i10,e18.10)
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_stp')
      !
   END SUBROUTINE trc_stp

      !!-------------------------------------------------------------------
         SUBROUTINE trc_bloom( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_bloom  ***
      !!                      
      !! ** Purpose : Computation of various diagnostic of the phytoplanktonic bloom
      !! 
      !! ** Method  : 
      !!
      !! History :
      !!   1.0  !  03-May-2012  (Sakina Ayata)  Original code
      !!-------------------------------------------------------------------


      !! * Arguments
      INTEGER, INTENT(in) :: kt
      INTEGER             :: jk,jl,jj,ji,year,month,day
      REAL                :: val,dev

   ! Initialisation at the first time step
   ! -------------------------------------
   IF( kt  == nittrc000  ) THEN
      ironCUM   = 0._wp
      valmax    = 0._wp
      timemax   = 0._wp
      propmax   = 0._wp
      maxdev    = 0._wp
      timedev   = 0._wp
      totchl(:) = 0._wp
      timetotchl(:) = 0._wp
      year=NINT(nn_date0/10000.)
      month=NINT((nn_date0-year*10000)/100.)
      day=nn_date0-year*10000-month*100
      CALL ymds2ju(year,month,day,0.0,zjul)
   ENDIF

   ! Calculattion of diagnostics to characterize the phytoplanktonic bloom
   ! Value, timing, and proportion of diatom during the max of chl at the surface
   ! ------------------------------------------------------------------------------
   !DO jj = 1, jpj
   !   DO ji =1 ,jpi
   ! Code for 1D only

      ! Max chl
      val=tra(2,2,1,jpnch)+tra(2,2,1,jpdch)
      IF(val.GT.valmax) THEN
          valmax=val
          timemax=(kt*rn_rdt)/60/60/24+MOD(zjul+1,365.)   ! in days
          propmax=tra(2,2,1,jpdch)/val
      ENDIF

      ! Derivative of chl (estimated for each physical time step kt)
      if (kt .GT. nittrc000+period) then
         dev=(val-totchl(1))/(period-1)
         IF (dev.GT.maxdev) THEN
           maxdev=dev
           timedev=timetotchl((period+1)/2)  ! in days
           parmlddev=emoy(2,2,1)
	   ENDIF
      endif
      
      DO jj=1,period-1
          totchl(jj)=totchl(jj+1)
          timetotchl(jj)=timetotchl(jj+1)
      ENDDO
      totchl(period)=val
      timetotchl(period)=((kt-0.5)*rn_rdt)/60/60/24+MOD(zjul+1,365.)   ! in days

      ironCUM = MAX( 0., xtrad(2,2,nmln(2,2)-2,jpfer)) + ironCUM
   !  END DO
   !END DO

   ! Print results
   ! ---------------
   IF (kt ==nitend) THEN
     WRITE(*,*) 'Diagnostics_of_the_bloom_(at_the_surface) '
     WRITE(*,*) 'Max_Chl_(surface) =',valmax*1000000
     WRITE(*,*) 'Date_of_Max_Chl_(Julian_Day) =',timemax
     WRITE(*,*) 'Diatom_at_Max_Chl_(Proportion) =',propmax*100
     WRITE(*,*) 'Max_Dev_Chl_(surface) =',maxdev*1000000
     WRITE(*,*) 'Date_of_Max_Dev_Chl_(Julien_Day) =',timedev
     WRITE(*,*) 'Mean_PAR_in_MLD_tonet =', parmlddev
     WRITE(*,*) 'Total_primary_production_per_year =', tpp2 * 12./1.E12
     WRITE(*,*) 'Total_iron_injection_in_MLD_(umol/year) =', ironCUM*1.E3*86400.*365.
     WRITE(*,*) 'Warning:_Dates_from_(nn_date0) =',MOD(zjul+1,365.)
!Joan
     WRITE(*,'(8F12.6,20X)') valmax*1000000, timemax, propmax*100, maxdev*1000000, timedev, tpp2*12./1.E12 &
			&	, parmlddev, ironCUM*1.E3*86400.*365.
  ENDIF

   END SUBROUTINE trc_bloom

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp
#endif

   !!======================================================================
END MODULE trcstp
