! Driver for RD_SN_DATA (new_main/subs/rd_sn_data.f).
!
! The alteration under test is the NISO=-1 header sentinel: an
! absent 'Number of isotopes:' line is now distinguished from an
! explicit 'Number of isotopes: 0'.  Upstream initialised NISO=0,
! so an explicit zero failed the presence check identically to a
! missing header.  With the sentinel, NISO<0 only fires when the
! header line is truly absent.
!
! The driver only exercises the header gate (first ~30 lines of the
! routine).  A file with NISO=0 passes the gate and prints the
! 'Read in SN_HYDRO_DATA data headers' banner before dying on the
! truncated body; a file missing the isotope header (or the age
! header) still STOPs at the gate.  Full SN_HYDRO_DATA ingestion
! needs the isotope/decay model state and stays an integration-level
! concern, not a unit test.
!
      PROGRAM DRV_RD_SN_DATA
      USE SET_KIND_MODULE
      USE MOD_CMFGEN
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=4, LU=10
      EXTERNAL RD_SN_DATA
      CALL RD_SN_DATA(ND,.TRUE.,LU)
      WRITE(*,'(A)')'RETURNED'
      END
