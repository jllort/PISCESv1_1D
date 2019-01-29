MODULE trdmod_trc
   !!======================================================================
   !!                       ***  MODULE  trdmod_trc  ***
   !! Ocean diagnostics:  mixed layer passive tracer trends 
   !!======================================================================
   !! History :  3.0  !  2010-07  (C. Ethe)  Original code (from trdmod.F90)
   !!----------------------------------------------------------------------
#if   defined key_top && ( defined key_trdmld_trc   ||   defined key_trdtrc )
   !!----------------------------------------------------------------------
   !!   'key_trdmld_trc'                  mixed layer trend diagnostics
   !!   'key_trdtrc'                      3D trend diagnostics
   !!----------------------------------------------------------------------
   !!   trd_mod_trc      : passive tracer trends 
   !!----------------------------------------------------------------------
   USE trc               ! tracer definitions (trn, trb, tra, etc.)
   USE trcnam_trp
   USE trdmod_oce
   USE trdmod_trc_oce    ! definition of main arrays used for trends computations
   USE trdmld_trc        ! Mixed layer trends diag.
   USE iom               ! I/O library

   IMPLICIT NONE
   PRIVATE

   INTERFACE trd_mod_trc
      MODULE PROCEDURE trd_mod_trc_trp, trd_mod_trc_bio
   END INTERFACE

   PUBLIC trd_mod_trc

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Header:  $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trd_mod_trc_trp( ptrtrd, kjn, ktrd, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mod_trc  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   kjn                                 ! tracer index
      INTEGER, INTENT( in )  ::   ktrd                                ! tracer trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout )  ::   ptrtrd  ! Temperature or U trend
      CHARACTER (len=16) :: cltra
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
!         IF(lwp)WRITE(numout,*)
!         IF(lwp)WRITE(numout,*) 'trd_mod_trc:'
!         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Mixed layer trends for passive tracers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#if defined key_trdmld_trc 
      IF( lk_trdmld_trc .AND. ln_trdtrc( kjn ) ) THEN
         !
         SELECT CASE ( ktrd )
         CASE ( jptra_trd_xad     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_xad, '3D', kjn )
         CASE ( jptra_trd_yad     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_yad, '3D', kjn )
         CASE ( jptra_trd_zad     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_zad, '3D', kjn )
         CASE ( jptra_trd_ldf     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_ldf, '3D', kjn )
         CASE ( jptra_trd_bbl     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_bbl, '3D', kjn )
         CASE ( jptra_trd_zdf     )
            IF( ln_trcldf_iso ) THEN
               CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_ldf, '3D', kjn )
            ELSE
               CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_zdf, '3D', kjn )
            ENDIF
         CASE ( jptra_trd_dmp     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_dmp , '3D', kjn )
         CASE ( jptra_trd_nsr     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_sbc , '2D', kjn )
         CASE ( jptra_trd_sms     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_sms , '3D', kjn )
         CASE ( jptra_trd_radb    )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_radb, '3D', kjn )
         CASE ( jptra_trd_radn    )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_radn, '3D', kjn )
         CASE ( jptra_trd_atf     )   ;   CALL trd_mld_trc_zint( ptrtrd, jpmld_trc_atf , '3D', kjn )
         END SELECT
         !
      END IF
#endif 

      IF( lk_trdtrc .AND. ln_trdtrc( kjn ) ) THEN
         !
         SELECT CASE( ktrd )
         CASE( jptra_trd_xad  )       ;    WRITE (cltra,'("XAD_",4a)') 
         CASE( jptra_trd_yad  )       ;    WRITE (cltra,'("YAD_",4a)')
         CASE( jptra_trd_zad  )       ;    WRITE (cltra,'("ZAD_",4a)')
         CASE( jptra_trd_ldf  )       ;    WRITE (cltra,'("LDF_",4a)') 
         CASE( jptra_trd_bbl  )       ;    WRITE (cltra,'("BBL_",4a)')
         CASE( jptra_trd_nsr  )       ;    WRITE (cltra,'("FOR_",4a)')
         CASE( jptra_trd_zdf  )       ;    WRITE (cltra,'("ZDF_",4a)')
         CASE( jptra_trd_dmp  )       ;    WRITE (cltra,'("DMP_",4a)') 
         CASE( jptra_trd_sms  )       ;    WRITE (cltra,'("SMS_",4a)') 
         CASE( jptra_trd_atf  )       ;    WRITE (cltra,'("ATF_",4a)') 
         CASE( jptra_trd_radb )       ;    WRITE (cltra,'("RDB_",4a)') 
         CASE( jptra_trd_radn )       ;    WRITE (cltra,'("RDN_",4a)') 
         END SELECT
                                          cltra = TRIM(cltra)//TRIM(ctrcnm(kjn)) 
                                          CALL iom_put( cltra,  ptrtrd(:,:,:) )
         !
      END IF

   END SUBROUTINE trd_mod_trc_trp

   SUBROUTINE trd_mod_trc_bio( ptrbio, ktrd, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mod_bio  ***
      !!----------------------------------------------------------------------

      INTEGER, INTENT( in )  ::   kt                                  ! time step
      INTEGER, INTENT( in )  ::   ktrd                                ! bio trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout )  ::   ptrbio  ! Bio trend
      !!----------------------------------------------------------------------

#if defined key_trdmld_trc 
      CALL trd_mld_bio_zint( ptrbio, ktrd ) ! Verticaly integrated biological trends
#endif 

   END SUBROUTINE trd_mod_trc_bio
#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------

   INTERFACE trd_mod_trc
      MODULE PROCEDURE trd_mod_trc_trp, trd_mod_trc_bio
   END INTERFACE

CONTAINS

   SUBROUTINE trd_mod_trc_trp( ptrtrd, kjn, ktrd, kt )
      INTEGER               , INTENT( in )     ::   kt      ! time step
      INTEGER               , INTENT( in )     ::   kjn     ! tracer index
      INTEGER               , INTENT( in )     ::   ktrd    ! tracer trend index
      REAL, DIMENSION(:,:,:), INTENT( inout )  ::   ptrtrd  ! Temperature or U trend
      WRITE(*,*) 'trd_mod_trc_trp : You should not have seen this print! error?', ptrtrd(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kt
   END SUBROUTINE trd_mod_trc_trp

   SUBROUTINE trd_mod_trc_bio( ptrbio, ktrd, kt )
      INTEGER               , INTENT( in )     ::   kt      ! time step
      INTEGER               , INTENT( in )     ::   ktrd    ! tracer trend index
      REAL, DIMENSION(:,:,:), INTENT( inout )  ::   ptrbio  ! Temperature or U trend
      WRITE(*,*) 'trd_mod_trc_trp : You should not have seen this print! error?', ptrbio(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kt
   END SUBROUTINE trd_mod_trc_bio

#endif
   !!======================================================================
END MODULE trdmod_trc
