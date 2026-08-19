  MODULE quantize_module

  use input , only : output_signif_digits , prcl_signif_digits ,        &
                     output_pack , prcl_pack

  implicit none

  private
  public :: bitround , nsd_to_nsb , nsd_gridded , nsd_parcel
  public :: pack_gridded , pack_parcel , pack_i2 , pack_i4 , pack_range_ok
  public :: pack_check
  public :: pk_none , pk_i2 , pk_i4

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

  !---------------------------------------------------------------------
  !  Fixed-point packing kinds, returned by pack_gridded / pack_parcel.
  !  The value is the storage size in bytes, so it doubles as documentation.
  !---------------------------------------------------------------------
  integer, parameter :: pk_none = 0    ! leave as nf90_float (lossless)
  integer, parameter :: pk_i2   = 2    ! nf90_short  + scale_factor/add_offset
  integer, parameter :: pk_i4   = 4    ! nf90_int    + scale_factor/add_offset

  ! Representable integer range.  -32768 / -2147483648 are deliberately left
  ! unused: several tools treat the most negative value as a fill value, and
  ! keeping the range symmetric makes the clamp easy to reason about.
  real(kind=8), parameter :: i2lim = 32767.0d0
  real(kind=8), parameter :: i4lim = 2147483647.0d0

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


!-----------------------------------------------------------------------
!
!  FIXED-POINT PACKING  (output_pack / prcl_pack)
!
!  WHY THIS EXISTS ALONGSIDE bitround
!    bitround preserves *relative* precision.  For almost every CM1 field
!    that is the wrong knob: what the science needs is an *absolute* floor
!    (0.01 m/s on a wind, 0.01 K on a temperature, 1 Pa on a pressure).
!    Relative precision spends as many bits describing w = 0.0013 m/s as it
!    does w = 43 m/s -- and most of the domain, most of the time, is full of
!    tiny boundary-layer velocity noise whose mantissas are pure entropy.
!    Measured on real sc_lt21 output, that is exactly why the velocity
!    fields -- 71% of the gridded file -- still compress only 1.2-1.8x.
!
!    Packing to a fixed absolute step collapses all that noise onto a few
!    integer levels, which deflate then eats.  It also halves (i2) the bytes
!    handed to the compressor before compression even starts.  Measured:
!    the gridded stream goes 3.4x -> ~7x, the parcel stream 1.30x -> ~2x.
!
!  HOW IT IS STORED
!    Standard CF: the variable is defined as nf90_short / nf90_int and
!    carries scale_factor + add_offset, with
!        unpacked = packed*scale_factor + add_offset.
!    netCDF does NOT apply this itself -- it is a client convention -- so
!    the model packs explicitly here and the writers put integers.  Every
!    normal reader (xarray, netCDF4-python) unpacks automatically, so
!    downstream analysis needs no changes at all.
!
!  THE OFFSET MUST BE A CONSTANT
!    A netcdf variable has ONE scale_factor/add_offset for the whole array,
!    so the offsets below are fixed constants, not derived from the data.
!    That matters for three separate reasons: the parallel writers define
!    variables collectively and every rank must agree; a data-derived offset
!    would differ between output times, making files non-comparable; and a
!    constant offset is what lets a later time-delta scheme subtract two
!    files' packed integers directly.
!
!    The constants below were chosen from measured min/max over the whole
!    sc_lt21 run (t = 0, 100, 181; all 145 M parcels, full 3D fields), with
!    generous headroom -- the tightest is th, which still has 62% of the
!    int16 range spare.  Values outside the representable range are CLAMPED,
!    never wrapped, and pack_range_ok lets the writers shout about it.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Shared spec for names that mean the same thing in both streams.
!  Sets pk = pk_none if the name is not one of them.
!-----------------------------------------------------------------------

    subroutine pack_common( name , pk , scale , offset )
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(out) :: pk
    real, intent(out) :: scale , offset

      pk     = pk_none
      scale  = 1.0
      offset = 0.0

      select case( trim(adjustl(name)) )

        ! velocities: zero-centred, |v| < 100 m/s in any storm we run.
        ! int16 at 0.01 m/s reaches +/-327 m/s.
        case( 'u' , 'v' , 'w' , 'uinterp' , 'vinterp' , 'winterp' )
          pk = pk_i2 ; scale = 0.01 ; offset = 0.0

        ! potential temperature reaches ~540 K at the model top; temperature
        ! bottoms out near 190 K.  Offsets differ so both keep 0.01 K in int16.
        case( 'th' )
          pk = pk_i2 ; scale = 0.01 ; offset = 400.0
        case( 't' )
          pk = pk_i2 ; scale = 0.01 ; offset = 250.0

        ! pressure spans 4.2e3 - 1.0e5 Pa: ~1e5 levels at 1 Pa, so int32.
        case( 'prs' , 'psfc' )
          pk = pk_i4 ; scale = 1.0 ; offset = 52000.0

        case( 'pi' )
          pk = pk_i4 ; scale = 1.0e-5 ; offset = 0.7

        case( 'rho' )
          pk = pk_i4 ; scale = 1.0e-5 ; offset = 0.6

        ! mixing ratios (incl. saturation): 0 - 0.032 kg/kg, 1e-7 = 1e-4 g/kg.
        case( 'qv' , 'qc' , 'qr' , 'qi' , 'qs' , 'qg' , 'qhl' , 'qsl' , 'qsi' )
          pk = pk_i4 ; scale = 1.0e-7 ; offset = 0.0

        ! number concentrations reach ~4e6 /kg; 1 /kg is far below anything
        ! meaningful, and int32 covers 2e9.
        case( 'nci' , 'ncs' , 'ncr' , 'ncg' , 'nchl' , 'ncc' )
          pk = pk_i4 ; scale = 1.0 ; offset = 0.0

      end select

    end subroutine pack_common

!-----------------------------------------------------------------------
!  Packing spec for one gridded (cm1out_*.nc) variable.
!  pk = pk_none means "write it as float, as before".
!
!  Base-state fields stay lossless for the same reason they do in
!  nsd_gridded: they are written once, cost nothing, and every perturbation
!  diagnostic downstream is computed against them.
!-----------------------------------------------------------------------

    subroutine pack_gridded( name , pk , scale , offset )
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(out) :: pk
    real, intent(out) :: scale , offset

      pk     = pk_none
      scale  = 1.0
      offset = 0.0

      if( output_pack.le.0 ) return

      select case( trim(adjustl(name)) )
        case( 'th0' , 'prs0' , 'pi0' , 'rho0' , 'qv0' , 'u0' , 'v0' )
          return
      end select

      call pack_common( name , pk , scale , offset )

    end subroutine pack_gridded

!-----------------------------------------------------------------------
!  Packing spec for one parcel (cm1out_pdata_*.nc) variable.
!
!  Positions are the fields that most need care: x reaches 5.0e5 m on the
!  sc_lt21 domain, so at 0.25 m (= dx/1000) that is 2e6 levels and int32 is
!  required.  mtime is a scalar broadcast across every parcel -- deflate
!  already collapses it to nothing, so there is no reason to touch it.
!-----------------------------------------------------------------------

    subroutine pack_parcel( name , pk , scale , offset )
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(out) :: pk
    real, intent(out) :: scale , offset

      pk     = pk_none
      scale  = 1.0
      offset = 0.0

      if( prcl_pack.le.0 ) return

      select case( trim(adjustl(name)) )

        case( 'mtime' )
          return

        ! offsets put the domain centre at 0 so the int32 range is used
        ! symmetrically; these are domain-size independent up to +/-5e8 m.
        case( 'x' )
          pk = pk_i4 ; scale = 0.25 ; offset = 250000.0
        case( 'y' )
          pk = pk_i4 ; scale = 0.25 ; offset = 150000.0
        case( 'z' )
          pk = pk_i4 ; scale = 0.05 ; offset = 11000.0

        ! Buoyancy, vertical perturbation pressure gradient and vertical
        ! vorticity.  Steps here are deliberately twice as coarse as the
        ! measurement floor needs: over the whole sc_lt21 run these three are
        ! the tightest fields in the table (observed |vpg| 2.97 against a
        ! +/-3.28 cap at the finer step, i.e. only 2.2x margin), and they are
        ! small fields, so buying 2x the range costs nothing -- a coarser step
        ! compresses very slightly *better*.  4e-4 m/s/s is still 0.013% of a
        ! 3 m/s/s updraft signal.
        case( 'b' )
          pk = pk_i2 ; scale = 2.0e-4 ; offset = 0.0     ! +/-6.55 m/s/s
        case( 'vpg' )
          pk = pk_i2 ; scale = 4.0e-4 ; offset = 0.0     ! +/-13.1 m/s/s

        case( 'zvort' )
          pk = pk_i2 ; scale = 4.0e-5 ; offset = 0.0     ! +/-1.31 /s

        case default
          call pack_common( name , pk , scale , offset )

      end select

    end subroutine pack_parcel

!-----------------------------------------------------------------------
!  Float -> packed integer, CLAMPED (never wrapped) at the representable
!  range.  Elemental, so one call handles the 1d parcel buffers and the
!  2d/3d gridded buffers alike.
!
!  The arithmetic is done in double precision on purpose: (x-offset)/scale
!  for prs is ~1e5 and for qv ~1.6e5, which is fine in single, but z at
!  0.05 m is 4.4e5 and x at 0.25 m is 2e6 -- past the point where float32
!  can represent every integer exactly, which would quantize the quantizer.
!-----------------------------------------------------------------------

    elemental integer(kind=2) function pack_i2( x , scale , offset )
    use, intrinsic :: iso_fortran_env , only : real32
    implicit none

    real(real32), intent(in) :: x
    real, intent(in) :: scale , offset

    real(kind=8) :: d

      d = ( real(x,kind=8) - real(offset,kind=8) ) / real(scale,kind=8)
      d = max( -i2lim , min( i2lim , d ) )
      pack_i2 = int( nint(d,kind=4) , kind=2 )

    end function pack_i2

    elemental integer(kind=4) function pack_i4( x , scale , offset )
    use, intrinsic :: iso_fortran_env , only : real32
    implicit none

    real(real32), intent(in) :: x
    real, intent(in) :: scale , offset

    real(kind=8) :: d

      d = ( real(x,kind=8) - real(offset,kind=8) ) / real(scale,kind=8)
      d = max( -i4lim , min( i4lim , d ) )
      pack_i4 = nint( d , kind=4 )

    end function pack_i4

!-----------------------------------------------------------------------
!  .true. if [vmin,vmax] survives packing without clamping.
!
!  Clamping is silent data corruption -- the value written back is simply
!  wrong, and nothing in the file records that it happened.  The writers
!  call this once per field (two reductions, negligible next to deflate)
!  and print a loud warning if it fails, so a domain or a storm that outgrows
!  the table above shows up in the run log instead of in the science.
!-----------------------------------------------------------------------

    logical function pack_range_ok( vmin , vmax , pk , scale , offset )
    implicit none

    real, intent(in) :: vmin , vmax , scale , offset
    integer, intent(in) :: pk

    real(kind=8) :: dlo , dhi , lim

      pack_range_ok = .true.
      if( pk.eq.pk_none ) return

      if( pk.eq.pk_i2 )then
        lim = i2lim
      else
        lim = i4lim
      endif

      dlo = ( real(vmin,kind=8) - real(offset,kind=8) ) / real(scale,kind=8)
      dhi = ( real(vmax,kind=8) - real(offset,kind=8) ) / real(scale,kind=8)

      if( dlo.lt.(-lim) .or. dhi.gt.lim ) pack_range_ok = .false.

    end function pack_range_ok

!-----------------------------------------------------------------------
!  Range check + loud warning, for the writers to call once per field.
!
!  Deliberately a warning and not a stopcm1: clamping is wrong, but killing
!  a 60-hour run because one cell in one field grazed the limit would be
!  worse than finishing with a flagged field.  The message names the
!  variable and prints both ranges so the table in this module can be fixed.
!-----------------------------------------------------------------------

    subroutine pack_check( name , vmin , vmax , pk , scale , offset )
    implicit none

    character(len=*), intent(in) :: name
    real, intent(in) :: vmin , vmax , scale , offset
    integer, intent(in) :: pk

    real(kind=8) :: lim

      if( pk.eq.pk_none ) return
      if( pack_range_ok( vmin , vmax , pk , scale , offset ) ) return

      if( pk.eq.pk_i2 )then
        lim = i2lim
      else
        lim = i4lim
      endif

      print *,'  ###############################################################'
      print *,'  WARNING: fixed-point packing CLAMPED values in ',trim(adjustl(name))
      print *,'    data range        = ',vmin,vmax
      print *,'    representable     = ',real(offset,kind=8)-lim*real(scale,kind=8),  &
                                         real(offset,kind=8)+lim*real(scale,kind=8)
      print *,'    scale,offset,kind = ',scale,offset,pk
      print *,'    -> clamped values are WRONG.  Widen the table in quantize.F.'
      print *,'  ###############################################################'

    end subroutine pack_check


  END MODULE quantize_module
