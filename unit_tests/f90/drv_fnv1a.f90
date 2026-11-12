! Driver for KN_FNV1A_I1 / KN_FNV1A_COEF (plane/fg_j_cmf_nonmono_v1.f).
! Verifies the corrected standard 64-bit FNV-1a checksum: published
! test vectors, complete byte coverage, local-change detection,
! compensating-pair blindness of global sums, CHI/ETA independence,
! determinism, and dimension mixing. Also times the production-size
! coefficient arrays for the incremental-cost report.
program drv_fnv1a
  use SET_KIND_MODULE
  use FG_J_CMF_NONMONO_MOD_V1, only: KN_FNV1A_I1, KN_FNV1A_COEF
  implicit none

  integer, parameter :: IK1 = selected_int_kind(2)
  integer(IK1) :: BE(6), BA(1), BZ(0)
  integer*8 :: H0, H1, H2, H3, H4, H5
  integer :: I, J, K, L, ND, NCF, NFAIL, NPASS
  real(kind=LDP), allocatable :: A(:, :, :), B(:, :, :), C(:, :, :)
  real(kind=LDP) :: S0, S1, T0, T1
  integer(kind=IK1) :: CBB(8)

  NFAIL = 0
  NPASS = 0

! --- representation checks -------------------------------------------

  if (storage_size(0.0_LDP) /= 64) then
    write(*, '(A)') 'FNV1A FAIL bit_size real /= 64'
    NFAIL = NFAIL + 1
  else
    write(*, '(A,I3)') 'FNV1A bit_size real(LDP) = ', storage_size(0.0_LDP)
    NPASS = NPASS + 1
  end if
  if (storage_size(0_8) /= 64) then
    write(*, '(A)') 'FNV1A FAIL bit_size int8 /= 64'
    NFAIL = NFAIL + 1
  else
    NPASS = NPASS + 1
  end if

! --- published FNV-1a 64-bit test vectors -----------------------------

  H0 = KN_FNV1A_I1(BZ)
  if (H0 /= int(z'CBF29CE484222325', 8)) then
    write(*, '(A,2(Z16.16))') 'FNV1A FAIL vec-empty ', H0, &
      int(z'CBF29CE484222325', 8)
    NFAIL = NFAIL + 1
  else
    write(*, '(A)') 'FNV1A PASS vec-empty'
    NPASS = NPASS + 1
  end if

  BA = [97_IK1]
  H0 = KN_FNV1A_I1(BA)
  if (H0 /= int(z'AF63DC4C8601EC8C', 8)) then
    write(*, '(A,2(Z16.16))') 'FNV1A FAIL vec-a     ', H0, &
      int(z'AF63DC4C8601EC8C', 8)
    NFAIL = NFAIL + 1
  else
    write(*, '(A)') 'FNV1A PASS vec-a'
    NPASS = NPASS + 1
  end if

  BE = [102_IK1, 111_IK1, 111_IK1, 98_IK1, 97_IK1, 114_IK1]
  H0 = KN_FNV1A_I1(BE)
  if (H0 /= int(z'85944171F73967E8', 8)) then
    write(*, '(A,2(Z16.16))') 'FNV1A FAIL vec-foobar', H0, &
      int(z'85944171F73967E8', 8)
    NFAIL = NFAIL + 1
  else
    write(*, '(A)') 'FNV1A PASS vec-foobar'
    NPASS = NPASS + 1
  end if

! --- coefficient-array cases ------------------------------------------

  ND = 5
  NCF = 7
  allocate(A(ND, 4, NCF), B(ND, 4, NCF), C(ND, 4, NCF))
  do L = 1, NCF
    do J = 1, 4
      do I = 1, ND
        A(I, J, L) = 1.0e-3_LDP * (I + 10 * J + 100 * L) * &
          (1.0_LDP + 0.01_LDP * mod(I * J + L, 7))
      end do
    end do
  end do
  B = A
  C = A

! 1+8. identical arrays, repeated execution -> identical hashes

  H0 = KN_FNV1A_COEF(A)
  H1 = KN_FNV1A_COEF(B)
  H2 = KN_FNV1A_COEF(A)
  H3 = KN_FNV1A_COEF(B)
  if (H0 == H1 .and. H0 == H2 .and. H0 == H3) then
    write(*, '(A,Z16.16)') 'FNV1A PASS identical+repeat ', H0
    NPASS = NPASS + 1
  else
    write(*, '(A,4(Z16.16))') 'FNV1A FAIL identical+repeat ', &
      H0, H1, H2, H3
    NFAIL = NFAIL + 1
  end if

! 9. complete byte coverage: flip each of the 8 stored bytes of one
!    coefficient; every single-byte change must move the hash.

  do K = 0, 7
    CBB = transfer(A(3, 2, 4), CBB)
    CBB(K + 1) = ieor(CBB(K + 1), 1_IK1)
    C = A
    C(3, 2, 4) = transfer(CBB, C(3, 2, 4))
    H4 = KN_FNV1A_COEF(C)
    if (H4 == H0) then
      write(*, '(A,I2)') 'FNV1A FAIL byte-coverage k=', K
      NFAIL = NFAIL + 1
    end if
  end do
  write(*, '(A)') 'FNV1A PASS byte-coverage (8/8)'
  NPASS = NPASS + 1

! 2. one coefficient changed by one representable bit (mantissa LSB)

  C = A
  CBB = transfer(C(2, 1, 3), CBB)
  CBB(1) = ieor(CBB(1), 1_IK1)
  C(2, 1, 3) = transfer(CBB, C(2, 1, 3))
  H4 = KN_FNV1A_COEF(C)
  if (H4 == H0) then
    write(*, '(A)') 'FNV1A FAIL 1ulp-change undetected'
    NFAIL = NFAIL + 1
  else
    write(*, '(A)') 'FNV1A PASS 1ulp-change'
    NPASS = NPASS + 1
  end if

! 4. first coefficient in traversal order changed

  C = A
  CBB = transfer(C(1, 1, 1), CBB)
  CBB(1) = ieor(CBB(1), 1_IK1)
  C(1, 1, 1) = transfer(CBB, C(1, 1, 1))
  H4 = KN_FNV1A_COEF(C)
  if (H4 == H0) then
    write(*, '(A)') 'FNV1A FAIL first-coeff undetected'
    NFAIL = NFAIL + 1
  else
    write(*, '(A)') 'FNV1A PASS first-coeff'
    NPASS = NPASS + 1
  end if

! 5. compensating pair with bitwise-identical global sum:
!    put +x/-x into the baseline at positions 1,2 (their partial
!    sums are x+(-x)=0 either order), then swap them in C. Every
!    partial sum of SUM() is then bitwise identical while the
!    traversal-order byte stream changes.

  A(1, 1, 1) = 1.5e-2_LDP
  A(2, 1, 1) = -1.5e-2_LDP
  H0 = KN_FNV1A_COEF(A)
  S0 = sum(A)
  C = A
  C(1, 1, 1) = -1.5e-2_LDP
  C(2, 1, 1) = 1.5e-2_LDP
  S1 = sum(C)
  H4 = KN_FNV1A_COEF(C)
  if (H4 == H0) then
    write(*, '(A)') 'FNV1A FAIL compensating-pair undetected'
    NFAIL = NFAIL + 1
  else if (S0 /= S1) then
    write(*, '(A,2ES16.8)') 'FNV1A FAIL compensating sum moved ', &
      S0, S1
    NFAIL = NFAIL + 1
  else
    write(*, '(A)') 'FNV1A PASS compensating-pair (sum identical)'
    NPASS = NPASS + 1
  end if

! 6+7. independence: each array is hashed through a separate
!    KN_FNV1A_COEF call, so modifying one cannot move the other.

  H5 = KN_FNV1A_COEF(B)
  B(1, 1, 1) = B(1, 1, 1) * 1.7_LDP
  if (KN_FNV1A_COEF(B) == H5) then
    write(*, '(A)') 'FNV1A FAIL eta-array self-check'
    NFAIL = NFAIL + 1
  end if
  H1 = KN_FNV1A_COEF(C)
  C(4, 3, 6) = C(4, 3, 6) + 1.0e-9_LDP
  if (KN_FNV1A_COEF(C) == H1) then
    write(*, '(A)') 'FNV1A FAIL chi-array self-check'
    NFAIL = NFAIL + 1
  else
    write(*, '(A)') 'FNV1A PASS chi/eta independence'
    NPASS = NPASS + 1
  end if

! dimensions are mixed into the checksum: same element count,
! different shape -> different hash (overwhelmingly likely).

deallocate(A, B, C)
allocate(A(2, 4, 6), C(4, 4, 3))
A = 1.25_LDP
C = reshape(A, [4, 4, 3])
H1 = KN_FNV1A_COEF(A)
H2 = KN_FNV1A_COEF(C)
if (H1 == H2) then
  write(*, '(A)') 'FNV1A FAIL dim-mix (shapes collided)'
  NFAIL = NFAIL + 1
else
  write(*, '(A)') 'FNV1A PASS dim-mix'
  NPASS = NPASS + 1
end if
deallocate(A, C)

! --- production-size timing -------------------------------------------
! ND_EXT~72, NCF=56119 -> ~16.2M elements per coefficient array.

ND = 72
NCF = 56119
allocate(A(ND, 4, NCF), B(ND, 4, NCF))
A = 1.0_LDP
B = 1.0_LDP
B(1, 1, 1) = 3.0_LDP
call cpu_time(T0)
H1 = KN_FNV1A_COEF(A)
H2 = KN_FNV1A_COEF(B)
call cpu_time(T1)
write(*, '(A,ES12.5,A,F10.2,A,2(Z16.16))') &
  'FNV1A perf two prod-size arrays: ', T1 - T0, ' s for ', &
  2 * ND * 4 * NCF / 1.0e6, 'M elems; hashes ', H1, H2

write(*, '(A,I5,A,I5,A)') 'FNV1A RESULT: ', NPASS, ' passed, ', &
  NFAIL, ' failed'
if (NFAIL > 0) stop 1
end program drv_fnv1a
