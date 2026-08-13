  PROGRAM test_quantize

!-----------------------------------------------------------------------
!  Standalone unit test for src/quantize.F (bitround + the digit tables).
!
!  Runs without MPI, netcdf, or a namelist -- the point is to prove the bit
!  arithmetic before it goes anywhere near the parallel writers, where a
!  wrong answer would show up as silently corrupted 1.4 TB of output.
!
!  Build + run:   cd shared/cm1_cd/tests && ./run_test_quantize.sh
!-----------------------------------------------------------------------

  use, intrinsic :: iso_fortran_env , only : int32 , real32
  use input , only : output_signif_digits , prcl_signif_digits
  use quantize_module , only : bitround , nsd_to_nsb , nsd_gridded , nsd_parcel

  implicit none

  integer :: nfail
  integer :: nsb
  real(real32) :: x , y , z
  real(real32), dimension(5) :: a , b
  integer(int32) :: i

  nfail = 0

  print *,'==================================================='
  print *,'  test_quantize'
  print *,'==================================================='

!-----------------------------------------------------------------------
!  nsd_to_nsb
!-----------------------------------------------------------------------

  call checki( 'nsd_to_nsb(0)  = 0  (lossless)' , nsd_to_nsb(0)  , 0  )
  call checki( 'nsd_to_nsb(-1) = 0  (lossless)' , nsd_to_nsb(-1) , 0  )
  call checki( 'nsd_to_nsb(1)  = 4'             , nsd_to_nsb(1)  , 4  )
  call checki( 'nsd_to_nsb(3)  = 10'            , nsd_to_nsb(3)  , 10 )
  call checki( 'nsd_to_nsb(4)  = 14'            , nsd_to_nsb(4)  , 14 )
  call checki( 'nsd_to_nsb(5)  = 17'            , nsd_to_nsb(5)  , 17 )
  call checki( 'nsd_to_nsb(6)  = 20'            , nsd_to_nsb(6)  , 20 )
  call checki( 'nsd_to_nsb(7)  = 23 (no-op)'    , nsd_to_nsb(7)  , 23 )

!-----------------------------------------------------------------------
!  no-ops: these must leave the bit pattern untouched
!-----------------------------------------------------------------------

  x = 3.14159265_real32
  y = x
  call bitround( y , 0 )
  call checkbits( 'nsb=0 is a no-op' , y , x )

  y = x
  call bitround( y , 23 )
  call checkbits( 'nsb=23 is a no-op' , y , x )

  y = x
  call bitround( y , 99 )
  call checkbits( 'nsb>23 is a no-op' , y , x )

!-----------------------------------------------------------------------
!  the discarded bits really are zero
!-----------------------------------------------------------------------

  x = 1.2345678_real32
  y = x
  call bitround( y , 10 )
  i = transfer( y , 0_int32 )
  call checki( 'low 13 mantissa bits cleared at nsb=10' ,   &
               iand( i , ishft(1_int32,23-10)-1_int32 ) , 0 )

!-----------------------------------------------------------------------
!  idempotence: rounding twice equals rounding once
!-----------------------------------------------------------------------

  x = 1234.56789_real32
  y = x
  call bitround( y , 10 )
  z = y
  call bitround( z , 10 )
  call checkbits( 'idempotent at nsb=10' , z , y )

!-----------------------------------------------------------------------
!  sign symmetry: -x must round to -(x rounded)
!-----------------------------------------------------------------------

  x = 1234.56789_real32
  y = x
  call bitround( y , 10 )
  z = -x
  call bitround( z , 10 )
  call checkbits( 'negative rounds symmetrically' , z , -y )

  x = 0.000173_real32
  y = x
  call bitround( y , 14 )
  z = -x
  call bitround( z , 14 )
  call checkbits( 'negative rounds symmetrically (small)' , z , -y )

!-----------------------------------------------------------------------
!  zero and signed zero survive
!-----------------------------------------------------------------------

  y = 0.0_real32
  call bitround( y , 10 )
  call checkbits( '+0 preserved' , y , 0.0_real32 )

  y = transfer( ishft(-1_int32,31) , 0.0_real32 )   ! -0.0
  z = y
  call bitround( z , 10 )
  call checkbits( '-0 preserved' , z , y )

!-----------------------------------------------------------------------
!  Inf / NaN pass through untouched
!-----------------------------------------------------------------------

  y = transfer( int(z'7F800000',int32) , 0.0_real32 )   ! +Inf
  z = y
  call bitround( z , 10 )
  call checkbits( '+Inf untouched' , z , y )

  y = transfer( int(z'FF800000',int32) , 0.0_real32 )   ! -Inf
  z = y
  call bitround( z , 10 )
  call checkbits( '-Inf untouched' , z , y )

  y = transfer( int(z'7FC00000',int32) , 0.0_real32 )   ! NaN
  z = y
  call bitround( z , 10 )
  call checkbits( 'NaN untouched' , z , y )

!-----------------------------------------------------------------------
!  huge finite value must not be rounded up into Inf
!-----------------------------------------------------------------------

  y = huge( 1.0_real32 )
  z = y
  call bitround( z , 10 )
  call checkl( 'huge() stays finite' , z .eq. y )

!-----------------------------------------------------------------------
!  relative error stays within the bound implied by the digit count.
!  With nsb kept bits the half-ulp bound is 2**-(nsb+1) relative.
!-----------------------------------------------------------------------

  call check_relerr(  3 )
  call check_relerr(  4 )
  call check_relerr(  5 )

!-----------------------------------------------------------------------
!  elemental: the same call works on an array, in place
!-----------------------------------------------------------------------

  a = (/ 1.5_real32 , -273.16_real32 , 1.0e-7_real32 , 98765.4_real32 , 0.0_real32 /)
  b = a
  call bitround( b , 10 )
  z = 0.0_real32
  do i = 1 , 5
    y = a(i)
    call bitround( y , 10 )
    if( transfer(y,0_int32) .ne. transfer(b(i),0_int32) ) z = 1.0_real32
  enddo
  call checkl( 'elemental array call matches scalar calls' , z .eq. 0.0_real32 )

!-----------------------------------------------------------------------
!  digit tables: off by default, per-variable extras when switched on
!-----------------------------------------------------------------------

  output_signif_digits = 0
  prcl_signif_digits   = 0
  call checki( 'gridded lossless when knob=0' , nsd_gridded('th') , 0 )
  call checki( 'parcel  lossless when knob=0' , nsd_parcel('x')   , 0 )

  prcl_signif_digits = 3
  call checki( 'parcel w  -> 3 (baseline)'      , nsd_parcel('w')     , 3 )
  call checki( 'parcel z  -> 4 (+1)'            , nsd_parcel('z')     , 4 )
  call checki( 'parcel x  -> 5 (+2)'            , nsd_parcel('x')     , 5 )
  call checki( 'parcel th -> 5 (+2)'            , nsd_parcel('th')    , 5 )
  call checki( 'parcel qv -> 4 (+1)'            , nsd_parcel('qv')    , 4 )
  call checki( 'parcel mtime stays lossless'    , nsd_parcel('mtime') , 0 )
  call checki( 'parcel name is trimmed/padded'  , nsd_parcel('x   ')  , 5 )

  output_signif_digits = 3
  call checki( 'gridded winterp -> 3 (baseline)' , nsd_gridded('winterp') , 3 )
  call checki( 'gridded th      -> 5 (+2)'       , nsd_gridded('th')      , 5 )
  call checki( 'gridded prs     -> 5 (+2)'       , nsd_gridded('prs')     , 5 )
  call checki( 'gridded qv      -> 4 (+1)'       , nsd_gridded('qv')      , 4 )
  call checki( 'gridded th0 stays lossless'      , nsd_gridded('th0')     , 0 )

  output_signif_digits = 6
  call checki( 'gridded th clamped to max_nsd=7' , nsd_gridded('th') , 7 )

!-----------------------------------------------------------------------

  print *,'==================================================='
  if( nfail.eq.0 )then
    print *,'  ALL TESTS PASSED'
  else
    print *,'  FAILURES: ',nfail
  endif
  print *,'==================================================='

  if( nfail.ne.0 ) stop 1

  CONTAINS

    subroutine checki( label , got , want )
    character(len=*), intent(in) :: label
    integer, intent(in) :: got , want
      if( got.eq.want )then
        print *,'  ok   : ',label
      else
        print *,'  FAIL : ',label,' got ',got,' want ',want
        nfail = nfail+1
      endif
    end subroutine checki

    subroutine checkl( label , cond )
    character(len=*), intent(in) :: label
    logical, intent(in) :: cond
      if( cond )then
        print *,'  ok   : ',label
      else
        print *,'  FAIL : ',label
        nfail = nfail+1
      endif
    end subroutine checkl

    ! compare bit patterns, not values, so that -0 vs +0 and NaN payloads
    ! are actually checked
    subroutine checkbits( label , got , want )
    character(len=*), intent(in) :: label
    real(real32), intent(in) :: got , want
      if( transfer(got,0_int32) .eq. transfer(want,0_int32) )then
        print *,'  ok   : ',label
      else
        print *,'  FAIL : ',label,' got ',got,' want ',want
        nfail = nfail+1
      endif
    end subroutine checkbits

    ! sweep a range of magnitudes and signs, and assert the relative error
    ! never exceeds the half-ulp bound for the kept bit count
    subroutine check_relerr( nsd )
    integer, intent(in) :: nsd
    integer :: k
    integer :: nb
    real(real32) :: v , q , err , worst , bound
    character(len=40) :: label
      nb = nsd_to_nsb( nsd )
      bound = 2.0_real32**(-(nb+1))
      worst = 0.0_real32
      do k = -20 , 20
        v = sign( 1.234567_real32 , real(k,real32) ) * (2.0_real32**k) * 1.7_real32
        q = v
        call bitround( q , nb )
        err = abs( (q-v)/v )
        worst = max( worst , err )
      enddo
      write(label,'(a,i1,a,i2,a)') 'relerr at nsd=',nsd,' (nsb=',nb,') within half-ulp'
      call checkl( trim(label) , worst.le.bound )
      print *,'         worst relative error = ',worst,'  bound = ',bound
    end subroutine check_relerr

  END PROGRAM test_quantize
