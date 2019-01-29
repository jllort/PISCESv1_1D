MODULE trdtra
   !!======================================================================
   !!                       ***  MODULE  trdtra  ***
   !! Ocean diagnostics:  ocean tracers trends
   !!=====================================================================
   !! History :  1.0  !  2004-08  (C. Talandier) Original code
   !!            2.0  !  2005-04  (C. Deltel)    Add Asselin trend in the ML budget
   !!            3.3  !  2010-06  (C. Ethe) merge TRA-TRC 
   !!----------------------------------------------------------------------
#if  defined key_trdtra || defined key_trdmld || defined key_trdmld_trc || defined key_trdtrc
   !!----------------------------------------------------------------------
   !!   trd_tra      : Call the trend to be computed
   !!----------------------------------------------------------------------
   USE dom_oce          ! ocean domain 
   USE trdmod_oce       ! ocean active mixed layer tracers trends 
   USE trdmod           ! ocean active mixed layer tracers trends 
   USE trdmod_trc       ! ocean passive mixed layer tracers trends 
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! MPP library
   USE wrk_nemo        ! Memory allocation


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_tra          ! called by all  traXX modules
 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: trdtx, trdty, trdt  !:

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: trdtra.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trd_tra_alloc()
      !!----------------------------------------------------------------------------
      !!                  ***  FUNCTION trd_tra_alloc  ***
      !!----------------------------------------------------------------------------
      ALLOCATE( trdtx(jpi,jpj,jpk) , trdty(jpi,jpj,jpk) , trdt(jpi,jpj,jpk) , STAT= trd_tra_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( trd_tra_alloc )
      IF( trd_tra_alloc /= 0 )   CALL ctl_warn('trd_tra_alloc: failed to allocate arrays')
   END FUNCTION trd_tra_alloc


   SUBROUTINE trd_tra( kt, ctype, ktra, ktrd, ptrd, pun, ptra )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra  ***
      !! 
      !! ** Purpose : Dispatch all trends computation, e.g. vorticity, mld or 
      !!              integral constraints
      !!
      !! ** Method/usage : For the mixed-layer trend, the control surface can be either
      !!       a mixed layer depth (time varying) or a fixed surface (jk level or bowl). 
      !!      Choose control surface with nn_ctls in namelist NAMTRD :
      !!        nn_ctls = 0  : use mixed layer with density criterion 
      !!        nn_ctls = 1  : read index from file 'ctlsurf_idx'
      !!        nn_ctls > 1  : use fixed level surface jk = nn_ctls
      !!----------------------------------------------------------------------
      !
      INTEGER                         , INTENT(in)           ::  kt      ! time step
      CHARACTER(len=3)                , INTENT(in)           ::  ctype   ! tracers trends type 'TRA'/'TRC'
      INTEGER                         , INTENT(in)           ::  ktra    ! tracer index
      INTEGER                         , INTENT(in)           ::  ktrd    ! tracer trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)           ::  ptrd    ! tracer trend  or flux
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  pun     ! velocity 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  ptra    ! Tracer variablea
      !
      REAL(wp), POINTER, DIMENSION(:,:,:)  ::  ztrds
      !!----------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, jpk, ztrds )

      IF( .NOT. ALLOCATED( trdtx ) ) THEN       ! allocate trdtra arrays
         IF( trd_tra_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trd_tra : unable to allocate arrays' )
      ENDIF
      
      ! Control of optional arguments
      IF( ctype == 'TRA' .AND. ktra == jp_tem ) THEN 
         IF( PRESENT( ptra ) ) THEN    
            SELECT CASE( ktrd )            ! shift depending on the direction
            CASE( jptra_trd_xad )  ;  CALL trd_tra_adv( ptrd, pun, ptra, 'X', trdtx ) 
            CASE( jptra_trd_yad )  ;  CALL trd_tra_adv( ptrd, pun, ptra, 'Y', trdty ) 
            CASE( jptra_trd_zad )  ;  CALL trd_tra_adv( ptrd, pun, ptra, 'Z', trdt  ) 
            END SELECT
         ELSE
            trdt(:,:,:) = ptrd(:,:,:)
            IF( ktrd == jptra_trd_bbc .OR. ktrd == jptra_trd_qsr ) THEN
               ztrds(:,:,:) = 0.
               CALL trd_mod( trdt, ztrds, ktrd, ctype, kt )
            END IF
         END IF
      END IF

      IF( ctype == 'TRA' .AND. ktra == jp_sal ) THEN 
         IF( PRESENT( ptra ) ) THEN    
            SELECT CASE( ktrd )            ! shift depending on the direction
            CASE( jptra_trd_xad )  
                                CALL trd_tra_adv( ptrd, pun, ptra, 'X', ztrds ) 
                                CALL trd_mod( trdtx, ztrds, ktrd, ctype, kt   )
            CASE( jptra_trd_yad )  
                                CALL trd_tra_adv( ptrd, pun, ptra, 'Y', ztrds ) 
                                CALL trd_mod( trdty, ztrds, ktrd, ctype, kt   )
            CASE( jptra_trd_zad )    
                                CALL trd_tra_adv( ptrd, pun, ptra, 'Z', ztrds ) 
                                CALL trd_mod( trdt , ztrds, ktrd, ctype, kt   )
            END SELECT
         ELSE
            ztrds(:,:,:) = ptrd(:,:,:)
            CALL trd_mod( trdt, ztrds, ktrd, ctype, kt )  
         END IF
      END IF

      IF( ctype == 'TRC' ) THEN
         !
         IF( PRESENT( ptra ) ) THEN  
            SELECT CASE( ktrd )            ! shift depending on the direction
            CASE( jptra_trd_xad )  
                                CALL trd_tra_adv( ptrd, pun, ptra, 'X', ztrds ) 
                                CALL trd_mod_trc( ztrds, ktra, ktrd, kt       )
            CASE( jptra_trd_yad )  
                                CALL trd_tra_adv( ptrd, pun, ptra, 'Y', ztrds ) 
                                CALL trd_mod_trc( ztrds, ktra, ktrd, kt       )
            CASE( jptra_trd_zad )    
                                CALL trd_tra_adv( ptrd, pun, ptra, 'Z', ztrds ) 
                                CALL trd_mod_trc( ztrds, ktra, ktrd, kt       )
            END SELECT
         ELSE
            ztrds(:,:,:) = ptrd(:,:,:)
            CALL trd_mod_trc( ztrds, ktra, ktrd, kt       )  
         END IF
         !
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztrds )
      !
   END SUBROUTINE trd_tra


   SUBROUTINE trd_tra_adv( pf, pun, ptn, cdir, ptrd )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_adv  ***
      !! 
      !! ** Purpose :   transformed the i-, j- or k-advective flux into thes
      !!              i-, j- or k-advective trends, resp.
      !! ** Method  :   i-advective trends = -un. di-1[T] = -( di-1[fi] - tn di-1[un] )
      !!                k-advective trends = -un. di-1[T] = -( dj-1[fi] - tn dj-1[un] )
      !!                k-advective trends = -un. di+1[T] = -( dk+1[fi] - tn dk+1[un] )
      !!----------------------------------------------------------------------
      REAL(wp)        , INTENT(in ), DIMENSION(jpi,jpj,jpk) ::   pf      ! advective flux in one direction
      REAL(wp)        , INTENT(in ), DIMENSION(jpi,jpj,jpk) ::   pun     ! now velocity  in one direction
      REAL(wp)        , INTENT(in ), DIMENSION(jpi,jpj,jpk) ::   ptn     ! now or before tracer 
      CHARACTER(len=1), INTENT(in )                         ::   cdir    ! X/Y/Z direction
      REAL(wp)        , INTENT(out), DIMENSION(jpi,jpj,jpk) ::   ptrd    ! advective trend in one direction
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ii, ij, ik   ! index shift function of the direction
      REAL(wp) ::   zbtr         ! local scalar
      !!----------------------------------------------------------------------

      SELECT CASE( cdir )            ! shift depending on the direction
      CASE( 'X' )   ;   ii = 1   ; ij = 0   ;   ik = 0      ! i-advective trend
      CASE( 'Y' )   ;   ii = 0   ; ij = 1   ;   ik = 0      ! j-advective trend
      CASE( 'Z' )   ;   ii = 0   ; ij = 0   ;   ik =-1      ! k-advective trend
      END SELECT

      !                              ! set to zero uncomputed values
      ptrd(jpi,:,:) = 0.e0   ;   ptrd(1,:,:) = 0.e0
      ptrd(:,jpj,:) = 0.e0   ;   ptrd(:,1,:) = 0.e0
      ptrd(:,:,jpk) = 0.e0
      !
      !
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbtr    = 1.e0/ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
               ptrd(ji,jj,jk) = - zbtr * (      pf (ji,jj,jk) - pf (ji-ii,jj-ij,jk-ik)                    &
                 &                          - ( pun(ji,jj,jk) - pun(ji-ii,jj-ij,jk-ik) ) * ptn(ji,jj,jk)  )
            END DO
         END DO
      END DO
      !
   END SUBROUTINE trd_tra_adv

#   else
   !!----------------------------------------------------------------------
   !!   Default case :          Dummy module           No trend diagnostics
   !!----------------------------------------------------------------------
   USE par_oce      ! ocean variables trends
CONTAINS
   SUBROUTINE trd_tra( kt, ctype, ktra, ktrd, ptrd, pu, ptra )
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in)           ::  kt      ! time step
      CHARACTER(len=3)                , INTENT(in)           ::  ctype   ! tracers trends type 'TRA'/'TRC'
      INTEGER                         , INTENT(in)           ::  ktra    ! tracer index
      INTEGER                         , INTENT(in)           ::  ktrd    ! tracer trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)           ::  ptrd    ! tracer trend 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  pu      ! velocity 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  ptra    ! Tracer variable 
      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd(1,1,1), ptra(1,1,1), pu(1,1,1),   &
         &                                                               ktrd, ktra, ctype, kt
   END SUBROUTINE trd_tra
#   endif
   !!======================================================================
END MODULE trdtra
