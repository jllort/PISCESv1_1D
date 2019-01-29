MODULE p4zbio
   !!======================================================================
   !!                         ***  MODULE p4zbio  ***
   !! TOP :   PISCES bio-model
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_bio        :   computes the interactions between the different
   !!                      compartments of PISCES
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmort         !  Mortality terms for phytoplankton
   USE p4zmicro        !  Sources and sinks of microzooplankton
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p4zrem          !  Remineralisation of organic matter
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
  
   IMPLICIT NONE
   PRIVATE

   PUBLIC  p4z_bio    

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zbio.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_bio ( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_bio  ***
      !!
      !! ** Purpose :   Ecosystem model in the whole ocean: computes the
      !!              different interactions between the different compartments
      !!              of PISCES
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, jnt
      INTEGER  ::  ji, jj, jk, jn
      REAL(wp) ::  ztra
#if defined key_kriest
      REAL(wp) ::  zcoef1, zcoef2
#endif
      CHARACTER (len=25) :: charout

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_bio')
      !
      !     ASSIGN THE SHEAR RATE THAT IS USED FOR AGGREGATION
      !     OF PHYTOPLANKTON AND DETRITUS

      xdiss(:,:,:) = 1.
!!gm the use of nmld should be better here?
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( fsdepw(ji,jj,jk+1) > hmld(ji,jj) )   xdiss(ji,jj,jk) = 0.01
            END DO 
         END DO
      END DO

          
      CALL p4z_sink ( kt, jnt )     ! vertical flux of particulate organic matter
      CALL p4z_opt  ( kt, jnt )     ! Optic: PAR in the water column
      CALL p4z_lim  ( kt      )     ! co-limitations by the various nutrients
      CALL p4z_prod ( kt, jnt )     ! phytoplankton growth rate over the global ocean. 
      !                             ! (for each element : C, Si, Fe, Chl )
      CALL p4z_rem  ( kt      )     ! remineralization terms of organic matter+scavenging of Fe
      CALL p4z_mort ( kt      )     ! phytoplankton mortality
      !                             ! zooplankton sources/sinks routines 
      CALL p4z_micro( kt      )           ! microzooplankton
      CALL p4z_meso ( kt, jnt )           ! mesozooplankton

      !                             ! test if tracers concentrations fall below 0.
      xnegtr(:,:,:) = 1.e0
      DO jn = jp_pcs0, jp_pcs1
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( ( trn(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) < 0.e0 ) THEN 
                     ztra             = ABS(  ( trn(ji,jj,jk,jn) - rtrn ) &
                                            / ( tra(ji,jj,jk,jn) + rtrn ) )
                     xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
                  ENDIF
              END DO
            END DO
         END DO
      END DO
      !                                ! where at least 1 tracer concentration becomes negative
      !                                ! 
      DO jn = jp_pcs0, jp_pcs1
         trn(:,:,:,jn) = trn(:,:,:,jn) + xnegtr(:,:,:) * tra(:,:,:,jn)
      END DO


      tra(:,:,:,:) = 0.e0

#if defined key_kriest
      ! 
      zcoef1 = 1.e0 / xkr_massp 
      zcoef2 = 1.e0 / xkr_massp / 1.1
      DO jk = 1,jpkm1
         trn(:,:,jk,jpnum) = MAX(  trn(:,:,jk,jpnum), trn(:,:,jk,jppoc) * zcoef1 / xnumm(jk)  )
         trn(:,:,jk,jpnum) = MIN(  trn(:,:,jk,jpnum), trn(:,:,jk,jppoc) * zcoef2              )
      END DO
#endif

      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_bio')
      !
   END SUBROUTINE p4z_bio

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_bio                         ! Empty routine
   END SUBROUTINE p4z_bio
#endif 

   !!======================================================================
END MODULE  p4zbio
