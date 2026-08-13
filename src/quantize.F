  MODULE quantize_module

  use input , only : output_signif_digits , prcl_signif_digits

  implicit none

  private
  public :: bitround , nsd_to_nsb , nsd_gridded , nsd_parcel

!-----------------------------------------------------------------------
!
!  Lossy "bit-rounding" (mantissa trimming) for compressed netcdf output.
!
!  WHAT IT DOES
!    Zeroes the low mantissa bits of a float32, keeping only the leading
!    nsb bits and rounding to nearest (ties to even).  The discarded bits
!    encode precision the model does not actually have, and once they are
!    zeroed the deflate filter -- which already runs on all compressed
!    output here -- packs the resulting long runs of zero bits far more
!    tightly.  Nothing else about the file changes: it is still an ordinary
!    deflated netcdf-4 variable, readable by any netcdf, with no plugin.
!
!  WHY IN THE MODEL, RATHER THAN nf90_def_var_quantize
!    This is exactly the transform netCDF's NC_QUANTIZE_BITROUND applies to
!    the data buffer before handing it to the deflate filter -- it is not an
!    HDF5 filter.  The netcdf-fortran on this machine is 4.5.3, which
!    predates nf90_def_var_quantize, and rebuilding netcdf against the
!    parallel HDF5 stack would put the working output_filetype=4 and
!    parallel-parcel writers at risk for no functional gain.  Doing the same
!    bit arithmetic here costs ~30 lines and no dependencies.
!
!  PRECISION IS RELATIVE, SO DIGITS ARE PER-VARIABLE
!    Bit-rounding preserves *relative* precision, which is the wrong knob
!    for fields carrying a large offset: at 3 significant digits a parcel x
!    position (up to 5e5 m) is wrong by 250 m -- a whole grid cell -- while
!    w or zvort at the same setting are over-resolved.  So the tables below
!    add extra digits to the offset-carrying fields, on top of the single
!    user-facing namelist knob.
!
!  See docs/claude/lossy_compression_plan.md for the measured error-vs-ratio
!  curves these defaults were chosen from.
!
!-----------------------------------------------------------------------

  ! IEEE-754 binary32 layout, used by bitround below.
  integer, parameter :: mant_bits = 23    ! explicit mantissa bits
  integer, parameter :: exp_all1  = 255   ! exponent field of Inf / NaN

  ! log2(10): decimal significant digits -> mantissa bits.
  real, parameter :: bits_per_digit = 3.3219281

  ! Beyond this many decimal digits, bit-rounding a float32 is a no-op
  ! (23 explicit mantissa bits ~ 7.2 decimal digits).
  integer, parameter :: max_nsd = 7

  CONTAINS

!-----------------------------------------------------------------------
!  Keep the leading nsb mantissa bits of x; zero the rest, rounding to
!  nearest with ties to even.  nsb <= 0 or nsb >= 23 is a no-op, so
!  "lossless" is the natural default rather than a special case.
!
!  Elemental, so the same call works on the 1d parcel buffers and the 2d/3d
!  gridded buffers, in place, with no temporary.
!
!  Note the dummy is real(real32), not real: if this fork is ever built in
!  double precision (-r8), that mismatch is a compile error rather than a
!  silently wrong bit pattern.  Today OPTS carries no -r8, and the netcdf
!  variables are nf90_float, so single precision is correct.
!-----------------------------------------------------------------------

    elemental subroutine bitround( x , nsb )
    use, intrinsic :: iso_fortran_env , only : int32 , real32
    implicit none

    real(real32), intent(inout) :: x
    integer, intent(in) :: nsb

    integer(int32) :: i , j , shift , half , mask

      ! nothing to do: lossless, or asking for more bits than exist
      if( nsb.le.0 .or. nsb.ge.mant_bits ) return

      i = transfer( x , 0_int32 )

      ! leave Inf and NaN exactly as they are (exponent field all ones).
      ! ishft is a logical shift, so this is correct for negative x too.
      if( iand( ishft(i,-mant_bits) , exp_all1 ) .eq. exp_all1 ) return

      shift = mant_bits - nsb                 ! number of bits to discard
      half  = ishft( 1_int32 , shift-1 )      ! weight of the top discarded bit
      mask  = ishft( -1_int32 , shift )       ! ones above the discarded bits

      ! Round to nearest, ties to even: add half an ulp-of-the-kept-precision,
      ! minus one, plus the lowest kept bit (which breaks exact ties toward
      ! the even value), then mask.  The sign bit is the top bit and the add
      ! only touches the mantissa/exponent below it, so this is correct for
      ! both signs -- the bit pattern of a negative float is sign-magnitude,
      ! and increasing the magnitude field increases |x| either way.
      j = i + ( half - 1_int32 ) + iand( ishft(i,-shift) , 1_int32 )
      j = iand( j , mask )

      ! guard the one pathological case: a value within rounding distance of
      ! the largest finite float32 would carry out of the exponent field and
      ! become Inf (or flip the sign bit).  Leave such values untouched.
      if( iand( ishft(j,-mant_bits) , exp_all1 ) .eq. exp_all1 ) return

      x = transfer( j , 0.0_real32 )

    end subroutine bitround

!-----------------------------------------------------------------------
!  Decimal significant digits -> mantissa bits to keep.
!  Returns 0 ("lossless") for nsd <= 0, which is the namelist default.
!-----------------------------------------------------------------------

    pure integer function nsd_to_nsb( nsd )
    implicit none

    integer, intent(in) :: nsd

      if( nsd.le.0 )then
        nsd_to_nsb = 0
      else
        nsd_to_nsb = min( mant_bits , max( 1 , ceiling( nsd*bits_per_digit ) ) )
      endif

    end function nsd_to_nsb

!-----------------------------------------------------------------------
!  Significant digits for one gridded (cm1out_*.nc) variable.
!
!  The namelist knob output_signif_digits sets the baseline; the table adds
!  digits back for fields whose values carry a large offset relative to the
!  signal (see the module header).  0 means this variable stays lossless.
!-----------------------------------------------------------------------

    integer function nsd_gridded( name )
    implicit none

    character(len=*), intent(in) :: name

    integer :: extra

      if( output_signif_digits.le.0 )then
        nsd_gridded = 0
        return
      endif

      select case( trim(adjustl(name)) )

        ! base state: written once, negligible size, and every perturbation
        ! diagnostic downstream is computed against it.  Keep it exact.
        case( 'th0' , 'prs0' , 'pi0' , 'rho0' , 'qv0' , 'u0' , 'v0' )
          nsd_gridded = 0
          return

        ! large offset, small signal: absolute values ~300 K, ~1e5 Pa, etc.,
        ! while the physically interesting part is the perturbation.
        case( 'th' , 't' , 'prs' , 'pi' , 'rho' , 'psfc' )
          extra = 2

        ! positive-definite with a narrow working range
        case( 'qv' )
          extra = 1

        ! everything else -- perturbations, velocities, hydrometeors, number
        ! concentrations, vorticity, surface fluxes -- is zero-centred or
        ! spans decades, so relative precision is the right thing.
        case default
          extra = 0

      end select

      nsd_gridded = min( max_nsd , output_signif_digits+extra )

    end function nsd_gridded

!-----------------------------------------------------------------------
!  Significant digits for one parcel (cm1out_pdata_*.nc) variable.
!  Same scheme as nsd_gridded, with the position fields called out: on the
!  sc_lt21 domain x reaches 5e5 m against dx = 250 m, so positions need the
!  most digits of anything in the file.
!-----------------------------------------------------------------------

    integer function nsd_parcel( name )
    implicit none

    character(len=*), intent(in) :: name

    integer :: extra

      if( prcl_signif_digits.le.0 )then
        nsd_parcel = 0
        return
      endif

      select case( trim(adjustl(name)) )

        ! model time, replicated across every parcel: constant, so it costs
        ! nothing after deflate either way.  Keep it exact.
        case( 'mtime' )
          nsd_parcel = 0
          return

        ! horizontal position: largest magnitudes in the file, and trajectory
        ! analysis is the whole point of the parcel output.
        case( 'x' , 'y' )
          extra = 2

        ! height: smaller magnitude than x/y, but dz_bot is only 25 m
        case( 'z' )
          extra = 1

        ! offset-carrying thermodynamics (see nsd_gridded)
        case( 'th' , 't' , 'prs' )
          extra = 2

        case( 'qv' , 'rho' , 'qsl' , 'qsi' )
          extra = 1

        ! velocities, buoyancy, vertical pressure gradient, vorticity,
        ! hydrometeors, number concentrations, surface diagnostics
        case default
          extra = 0

      end select

      nsd_parcel = min( max_nsd , prcl_signif_digits+extra )

    end function nsd_parcel

  END MODULE quantize_module
