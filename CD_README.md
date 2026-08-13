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

Unit test (no MPI, netcdf, or namelist needed, runs in seconds):

```
cd tests && ./run_test_quantize.sh
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