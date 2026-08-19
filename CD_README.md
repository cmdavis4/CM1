This readme describes the changes I have made to the version of CM1 attached to this readme.

FLEXVARS:
* Initial forcing
    * `var1`: Location of initial forcing in x direction. Implemented for warm bubble (`iinit=1`) and cold pool (`iinit=2`) forcings. For the warm bubble, sets the center of the bubble in the x direction. For the cold pool, sets the eastern edge of the cold pool.
    * `var2`: Center of bubble in y direction. Only implemented for the warm bubble forcing.
    * `var3`: Vertical location of initial forcing
        * Warm bubble: Height of the center of the bubble above the ground (m)
        * Cold pool: Depth of cold pool (m)
    * `var4`: Horizontal radius of bubble (m). Warm bubble only
    * `var5`: Vertical radius of bubble (m). Warm bubble only
    * `var6`: Max potential temperature perturbation
        * Warm bubble: Max potential temperature perturbation at bubble center (K)
        * Cold pool: Max potential temperature perturbation at surface (K)
    * `var7`: Flag controlling whether $q_v$ or the relative humidity should be kept constant in the warm bubble forcing. `var7=0` keeps $q_v$ constant, `var7=1` keeps RH constant. Warm bubble only
* Initial profile
    * `var8`: Mixed-layer $q_{v0}$ value for the Weisman-Klemp sounding. I think this was 11e-3 (ie 11 g/kg) in the default CM1 setup. The values in the original Weisman-Klemp 1982 paper were 11-14 g/kg.

OUTPUT COMPRESSION (`&param9`):

Three real namelist variables (not flex vars) control compression of the netcdf
output. All three default to today's behavior, so an existing `namelist.input`
that does not mention them produces byte-identical output.

* `output_deflate_level` (default `1`): zlib level, 1-9, for every compressed
  netcdf variable — gridded and parcel, serial and parallel writers. Was
  hardcoded to 1 in three places before. **Measured on real `sc_lt21` output,
  raising it is not worth it**: level 6 buys 1.7% (parcels) to 3.2% (gridded)
  for 2.3-2.8x the compression CPU, and level 9 buys ~2-4% for 10-16x. Leave
  it at 1 unless you have a reason.
* `output_signif_digits` (default `0` = lossless): decimal significant digits
  to keep in gridded `cm1out_*.nc` variables.
* `prcl_signif_digits` (default `0` = lossless): same, for parcel
  `cm1out_pdata_*.nc` variables.

The two `*_signif_digits` knobs switch on **lossy bit-rounding**: the low
mantissa bits of each float32 are zeroed (round-to-nearest, ties to even) so
that the deflate filter — which already runs on all this output — packs them
much more tightly. This is exactly the transform netCDF >= 4.9's
`nf90_def_var_quantize` applies internally; it is done in Fortran here because
the netcdf-fortran available on the run machine is 4.5.3. Implementation:
`src/quantize.F`, applied at the write sites in `writeout.F` (`writeo`,
`writeosp`, `writeocomm3`) and `writeout_nc.F` (`writepdata_nc`,
`writepdata_nc_par`). Quantized variables carry
`_QuantizeBitRoundNumberOfSignificantBits` and `quantization` attributes
recording what was done — lossy output is not reversible, so this is the only
record of a file's real precision.

Bit-rounding preserves *relative* precision, which is the wrong knob for
fields carrying a large offset (parcel `x` reaches 5e5 m against dx = 250 m;
`th` is ~300 K with O(1 K) perturbations). So `src/quantize.F` holds a small
per-variable table that adds digits back for those fields on top of the global
knob: positions and `th`/`t`/`prs` get +2, `z`/`qv`/`rho`/`qs*` get +1,
everything zero-centred gets the baseline, and `mtime` plus the base-state
fields stay lossless. `prcl_signif_digits = 3` therefore means 3 digits for
`w`, 4 for `z`, 5 for `x` and `th`.

FIXED-POINT PACKING (`output_pack`, `prcl_pack`):

Two more `&param9` integers, both default `0` (off):

* `output_pack` (default `0`): gridded `cm1out_*.nc` fixed-point packing
* `prcl_pack`   (default `0`): parcel `cm1out_pdata_*.nc` fixed-point packing

Set to `1` and the listed variables are stored as `nf90_short` / `nf90_int`
with CF `scale_factor` + `add_offset` instead of `nf90_float`. **Downstream
analysis needs no changes** -- xarray and netCDF4-python both unpack
`scale_factor`/`add_offset` automatically and hand back floats.

This is the knob that actually matters for size, and it supersedes
`*_signif_digits` for most fields. Bit-rounding keeps *relative* precision,
which is the wrong thing for CM1: it spends as many bits on `w = 0.0013 m/s`
as on `w = 43 m/s`, and most of the domain is full of tiny boundary-layer
velocity noise. That is why the velocity fields -- 71% of the gridded file --
still compressed only 1.2-1.8x. Packing to an absolute step collapses that
noise onto a few integer levels. **Measured on real `sc_lt21` output: the
gridded stream goes 3.4x -> ~7x.**

Precision is per-variable, from a table in `src/quantize.F`:

| fields | step | storage |
|---|---|---|
| `u v w uinterp vinterp winterp` | 0.01 m/s | int16 |
| `th` (offset 400), `t` (offset 250) | 0.01 K | int16 |
| `prs`, `psfc` | 1 Pa | int32 |
| `pi`, `rho` | 1e-5 | int32 |
| all mixing ratios incl. `qsl`/`qsi` | 1e-7 kg/kg | int32 |
| `nc*` number concentrations | 1 /kg | int32 |
| parcel `x`, `y` | 0.25 m | int32 |
| parcel `z` | 0.05 m | int32 |
| parcel `b` / `vpg` | 1e-4 / 2e-4 m/s/s | int16 |
| parcel `zvort` | 2e-5 /s | int16 |
| base state, `mtime`, 2d surface fields, anything untabled | -- | float (lossless) |

The offsets are **constants**, not derived from the data: a netcdf variable
has one `add_offset` for the whole array, all ranks must agree (the writers
define collectively), and a data-derived offset would differ between output
times. They were chosen from measured min/max over the whole `sc_lt21` run
with generous headroom -- the tightest is `th`, which still has 62% of the
int16 range spare.

Out-of-range values are **clamped, never wrapped**, and every write path calls
`pack_check`, which prints a loud `WARNING: fixed-point packing CLAMPED` block
naming the variable and both ranges. If you change domain size or run a much
more violent storm, watch for that in the run log.

One real (tiny) caveat: reading a packed variable as float32 gives an error of
`step/2` **plus** one float32 rounding of `packed*scale + offset` at the
variable's magnitude. For `th` that is 0.005 K + ~3e-5 K. Decoding in float64
is exactly within `step/2`. Both are asserted in `tests/lossy_verify`.

PARCEL ID COORDINATE:

`xh` in `cm1out_pdata_*.nc` is just the ramp `1..nparcels`, but it was written
into **every** per-time file, contiguous and uncompressed: 581 MB of every
8.6 GB file at 1e8 parcels, 6.8% of the whole parcel stream, ~105 GB per run.
It is now chunked + deflated (under 5 MB), and the serial writer builds it in
one buffered `put_var` instead of `nparcels` single-element calls. In the
parallel writer the `xh` write is now COLLECTIVE, since it is a compressed
variable. Parcel output is therefore no longer byte-identical to pre-2026-08
files -- the *contents* are unchanged and verified, only the encoding differs.

Unit test (no MPI, netcdf, or namelist needed, runs in seconds):

```
cd tests && ./run_test_quantize.sh
```

End-to-end (11 tiny 4-rank runs; checks that the defaults still reproduce
pre-feature output bit for bit, that the parallel and serial writers agree,
and every packing error bound):

```
cd tests/lossy_verify && ./setup.sh && sbatch <root>/verify.slurm
# then: python compare.py
```

**`src/Makefile` is gitignored** (it is machine-specific), so adding
`quantize.F` to it does not propagate through git. On any other machine, after
pulling, add to `src/Makefile`:

```
SRC = ... \
	quantize.F \
	...

quantize.o: input.o
writeout.o: ... quantize.o
writeout_nc.o: constants.o input.o quantize.o
```