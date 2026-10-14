      subroutine get_edep_shock_power(r,v,density,nd)
	USE SET_KIND_MODULE
      use mod_cmfgen, only : clump_fac
      use shock_power_mod
      use control_variable_mod, only : prescribed_shock_power,
     1     vloc_shock_power, dvloc_shock_power, scl_pwr_by_fcl,
     2     rloc_shock_power, drloc_shock_power,
     3     do_nonmono_sob, do_nonmono_knop
      implicit none

      integer i, nd

      REAL(KIND=LDP) r(nd)
      REAL(KIND=LDP) v(nd)
      REAL(KIND=LDP) density(nd)
      REAL(KIND=LDP) wgt(nd)
      REAL(KIND=LDP) wrk(nd)

      REAL(KIND=LDP) t1

      ! RLOC_SHOCK_POWER > 0 selects deposition centred on a radius
      ! rather than a velocity.  The radius coordinate is single-valued
      ! on every grid, so it is the only well-defined shock location on
      ! a non-monotonic V field.
      logical :: dep_in_radius
      dep_in_radius = rloc_shock_power .gt. 0._LDP

      ! VLOC_SHOCK_POWER defines a velocity-space location.  A repeated
      ! velocity maps to more than one radial branch under the non-monotonic
      ! treatments, so require a single-valued velocity coordinate.
      ! DO_NONMONO_SOB and DO_NONMONO_KNOP default to F, so the legacy path
      ! is unchanged when both flags are F.  The guard is skipped for
      ! radius deposition, which does not use V at all.
      if ((do_nonmono_sob .or. do_nonmono_knop) .and.
     1    .not. dep_in_radius) then
         do i=1,nd
            if (v(i).le.0._LDP) then
               write(6,*) 'non-monotonic mode requires positive V in',
     1              ' get_edep_shock_power'
               write(6,*) i,v(i)
               stop
            endif
         enddo
         do i=1,nd-1
            if (v(i+1).ge.v(i)) then
               write(6,*) 'shock-power deposition requires decreasing V',
     1              ' under the non-monotonic flags; use',
     2              ' RLOC_SHOCK_POWER for a radius-centred deposit'
               write(6,*) i,v(i),v(i+1)
               stop
            endif
         enddo
      endif
!
      ! set up the deposition profile
      wgt(1:nd) = 0.0_LDP
      do i=1,nd
         write(6,*) 'In get_edep i, V, FCL: ',i,v(i),clump_fac(i)
         call flush(6)
         if (dep_in_radius) then
            t1 = (r(i) - rloc_shock_power) / drloc_shock_power
         else
            t1 = (v(i) - vloc_shock_power) / dvloc_shock_power
         endif
         if (abs(t1).lt.8._LDP) wgt(i) = exp(-t1*t1)
         ! dividing by clump_fac makes the shock edep max in clumped regions
         ! since CMFGEN scales the edep by clump_fac
         if (scl_pwr_by_fcl) then
            if (clump_fac(i).le.0._LDP) then
               write(6,*) 'clump_fac is zero in get_edep_shock_power'
               write(6,*) i,v(i),clump_fac(i)
               stop
            else
               wgt(i) = wgt(i) / clump_fac(i)
            endif
         endif
      enddo
      !do i=1,nd
      !   write(6,*) v(i),clump_fac(i),wgt(i)
      !enddo

      if (.not. allocated(shock_power)) allocate(shock_power(nd))
      shock_power(1:nd) = wgt(1:nd)

      ! Need to normalize with clumping factored in
      wrk = shock_power * r * r * clump_fac
      call lum_from_eta_v2(wrk,r,'LINMON',nd)
      if (sum(wrk).eq.0._LDP) then
         write(6,*) 'sum(wrk) is zero -- problem'
         stop
      endif

      shock_power = shock_power * prescribed_shock_power /
     1     (4._LDP*acos(-1._LDP)*sum(wrk)*1d30)

      ! Sanity check
      if (dep_in_radius) then
         write(6,*) 'shock deposition centred on R =',
     1        rloc_shock_power,' (10^10 cm), width',drloc_shock_power
      else
         write(6,*) 'shock deposition centred on V =',
     1        vloc_shock_power,' km/s, width',dvloc_shock_power
      endif
      write(6,*) 'prescribed_shock_power : ',prescribed_shock_power
      wrk = shock_power * r * r * clump_fac
      call lum_from_eta_v2(wrk,r,'LINMON',nd)
      write(6,*) 'prescribed_shock_power from CMFGEN :',
     1     4*acos(-1.0_LDP)*sum(wrk)*1.0E+30_LDP
      call flush(6)

      return
      end subroutine
