"""Check the lossy-compression verification runs (see setup.sh).

Four things are asserted, in decreasing order of how much they matter:

  1. REGRESSION.  With the compression knobs at their defaults, the new binary
     must reproduce the pre-feature reference output in _ftype4_verify/ bit for
     bit.  If this fails, the feature is not actually off by default and
     nothing else here is worth reading.

  2. WRITE-PATH AGREEMENT.  The parallel collective writers (filetype=4,
     var16=1) and the serial ones (filetype=2, var16=0) must produce the same
     bits with quantization on, exactly as they do with it off.  This is the
     "does quantization survive the collective path" check -- bit-rounding is
     a pre-filter transform on each rank's own buffer, so it should be
     orthogonal to nf90_collective, but that is an assumption worth testing.

  3. ERROR BOUND.  Every quantized variable must be within the relative
     tolerance implied by its significant-digit count, and every variable the
     digit table declares lossless must be EXACTLY unchanged.

  4. RATIO.  Report what the size actually bought, per stream.

Run after the SLURM job finishes:
    ~/programs/miniforge3/envs/common/bin/python compare.py
"""
import glob
import os
import sys

import numpy as np
import netCDF4 as nc

ROOT = "/gpfs/wolf2/arm/atm124/proj-shared/cmdavis4/projects/dc_entrainment/sc/_lossy_verify"
REF = "/gpfs/wolf2/arm/atm124/proj-shared/cmdavis4/projects/dc_entrainment/sc/_ftype4_verify"

BITS_PER_DIGIT = np.log2(10.0)

# must mirror the tables in src/quantize.F
GRID_EXTRA = {"th": 2, "t": 2, "prs": 2, "pi": 2, "rho": 2, "psfc": 2, "qv": 1}
GRID_LOSSLESS = {"th0", "prs0", "pi0", "rho0", "qv0", "u0", "v0"}
PRCL_EXTRA = {"x": 2, "y": 2, "th": 2, "t": 2, "prs": 2,
              "z": 1, "qv": 1, "rho": 1, "qsl": 1, "qsi": 1}
PRCL_LOSSLESS = {"mtime"}
MAX_NSD = 7

# must mirror the pack tables in src/quantize.F: name -> (step, offset, bytes)
_VEL = {v: (0.01, 0.0, 2) for v in ("u", "v", "w", "uinterp", "vinterp", "winterp")}
_Q = {v: (1e-7, 0.0, 4) for v in ("qv", "qc", "qr", "qi", "qs", "qg", "qhl", "qsl", "qsi")}
_NC = {v: (1.0, 0.0, 4) for v in ("nci", "ncs", "ncr", "ncg", "nchl", "ncc")}
PACK_COMMON = dict(th=(0.01, 400.0, 2), t=(0.01, 250.0, 2),
                   prs=(1.0, 52000.0, 4), psfc=(1.0, 52000.0, 4),
                   pi=(1e-5, 0.7, 4), rho=(1e-5, 0.6, 4), **_VEL, **_Q, **_NC)
PACK_GRID = {k: v for k, v in PACK_COMMON.items()
             if k not in ("th0", "prs0", "pi0", "rho0", "qv0", "u0", "v0")}
PACK_PRCL = dict(PACK_COMMON,
                 x=(0.25, 250000.0, 4), y=(0.25, 150000.0, 4), z=(0.05, 11000.0, 4),
                 b=(2e-4, 0.0, 2), vpg=(4e-4, 0.0, 2), zvort=(4e-5, 0.0, 2))
PACK_PRCL.pop("mtime", None)

nfail = 0


def report(ok, label, detail=""):
    global nfail
    if not ok:
        nfail += 1
    print(f"  {'ok  ' if ok else 'FAIL'} : {label}{('   ' + detail) if detail else ''}")


def nsd_for(name, baseline, extra_tab, lossless_set):
    if baseline <= 0 or name in lossless_set:
        return 0
    return min(MAX_NSD, baseline + extra_tab.get(name, 0))


def tolerance(nsd):
    """Half-ulp relative bound for keeping nsd decimal digits."""
    nsb = min(23, max(1, int(np.ceil(nsd * BITS_PER_DIGIT))))
    return 2.0 ** (-(nsb + 1))


def files(case, pattern):
    return sorted(os.path.basename(f) for f in glob.glob(f"{ROOT}/{case}/{pattern}"))


def bitwise_identical(fa, fb):
    """Compare every data variable's raw bits (packed ints compared raw)."""
    a, b = nc.Dataset(fa), nc.Dataset(fb)
    a.set_auto_maskandscale(False)
    b.set_auto_maskandscale(False)
    try:
        if set(a.variables) != set(b.variables):
            return False, f"variable sets differ: {set(a.variables) ^ set(b.variables)}"
        for v in a.variables:
            da = np.asarray(a[v][:])
            db = np.asarray(b[v][:])
            if da.shape != db.shape:
                return False, f"{v} shape {da.shape} vs {db.shape}"
            if da.dtype != db.dtype:
                return False, f"{v} dtype {da.dtype} vs {db.dtype}"
            if da.dtype.kind == "f":
                ia = da.astype(da.dtype).view(f"u{da.dtype.itemsize}")
                ib = db.astype(db.dtype).view(f"u{db.dtype.itemsize}")
                if not np.array_equal(ia, ib):
                    n = int((ia != ib).sum())
                    return False, f"{v}: {n} of {ia.size} values differ in bits"
            elif not np.array_equal(da, db):
                return False, f"{v} differs"
        return True, ""
    finally:
        a.close()
        b.close()


def check_bitwise(case, refdir, pattern, label):
    fs = files(case, pattern)
    if not fs:
        report(False, label, "NO OUTPUT FILES (run failed?)")
        return
    for fn in fs:
        rp = f"{refdir}/{fn}"
        if not os.path.exists(rp):
            report(False, f"{label} [{fn}]", "reference missing")
            continue
        ok, why = bitwise_identical(f"{ROOT}/{case}/{fn}", rp)
        report(ok, f"{label} [{fn}]", why)


def check_lossy(ll_case, lossy_case, pattern, baseline, extra_tab, lossless_set, label):
    fs = files(ll_case, pattern)
    if not fs:
        report(False, label, "NO OUTPUT FILES (run failed?)")
        return
    size_ll = sum(os.path.getsize(f"{ROOT}/{ll_case}/{f}") for f in fs)
    size_lo = sum(os.path.getsize(f"{ROOT}/{lossy_case}/{f}") for f in files(lossy_case, pattern))
    worst = {}
    for fn in fs:
        a = nc.Dataset(f"{ROOT}/{ll_case}/{fn}")
        b = nc.Dataset(f"{ROOT}/{lossy_case}/{fn}")
        try:
            for v in a.variables:
                if v not in b.variables:
                    report(False, f"{label} [{fn}]", f"{v} missing from lossy run")
                    continue
                da = np.asarray(a[v][:], dtype="f8")
                db = np.asarray(b[v][:], dtype="f8")
                if da.shape != db.shape:
                    report(False, f"{label} [{fn}]", f"{v} shape mismatch")
                    continue
                nsd = nsd_for(v, baseline, extra_tab, lossless_set)
                nz = (da != 0) & np.isfinite(da) & np.isfinite(db)
                rel = np.max(np.abs((db[nz] - da[nz]) / da[nz])) if nz.any() else 0.0
                prev = worst.get(v, (0.0, nsd))[0]
                worst[v] = (max(prev, rel), nsd)
        finally:
            a.close()
            b.close()

    for v, (rel, nsd) in sorted(worst.items()):
        if nsd == 0:
            report(rel == 0.0, f"{label}: {v} declared lossless is exact",
                   "" if rel == 0.0 else f"rel err {rel:.3e}")
        else:
            tol = tolerance(nsd)
            report(rel <= tol, f"{label}: {v} within {nsd}-digit tolerance",
                   f"rel err {rel:.2e} <= {tol:.2e}")

    if size_lo:
        print(f"         size: {size_ll/1e6:.2f} MB lossless -> {size_lo/1e6:.2f} MB lossy "
              f"= {size_ll/size_lo:.2f}x")


print("=" * 72)
print("  1. REGRESSION: defaults must reproduce the pre-feature output exactly")
print("=" * 72)
check_bitwise("regress_ft2", f"{REF}/even_ft2", "cm1out_0*.nc", "ft2 gridded vs reference")
check_bitwise("regress_ft4", f"{REF}/even_ft4", "cm1out_0*.nc", "ft4 gridded vs reference")

print()
print("=" * 72)
print("  2. WRITE-PATH AGREEMENT: parallel vs serial writers, quantization on")
print("=" * 72)
for fn in files("prcl_ft4_lossy", "cm1out_0*.nc"):
    ok, why = bitwise_identical(f"{ROOT}/prcl_ft4_lossy/{fn}", f"{ROOT}/prcl_ft2_lossy/{fn}")
    report(ok, f"gridded ft4 == ft2 under quantization [{fn}]", why)
for fn in files("prcl_ft4_lossy", "cm1out_pdata_*.nc"):
    ok, why = bitwise_identical(f"{ROOT}/prcl_ft4_lossy/{fn}", f"{ROOT}/prcl_ft2_lossy/{fn}")
    report(ok, f"parcel parallel == serial under quantization [{fn}]", why)

print()
print("=" * 72)
print("  3+4. ERROR BOUNDS AND RATIOS")
print("=" * 72)
for ll, lo, tag in [("prcl_ft4_lossless", "prcl_ft4_lossy", "ft4 even"),
                    ("uneven_ft4_lossless", "uneven_ft4_lossy", "ft4 uneven")]:
    print(f"--- gridded, {tag} ---")
    check_lossy(ll, lo, "cm1out_0*.nc", 3, GRID_EXTRA, GRID_LOSSLESS, f"gridded {tag}")
    print(f"--- parcels, {tag} ---")
    check_lossy(ll, lo, "cm1out_pdata_*.nc", 3, PRCL_EXTRA, PRCL_LOSSLESS, f"parcel {tag}")

# ---------------------------------------------------------------------------
#  5. FIXED-POINT PACKING (output_pack / prcl_pack)
# ---------------------------------------------------------------------------

def check_pack(ll_case, pk_case, pattern, table, label):
    """Packed values must round-trip to within half a step, the on-disk type
    must be what the table asks for, and anything the table does not name must
    stay float and stay EXACTLY unchanged.

    Two separate error measures, because they mean different things:

      f64 : decode the raw integers in float64.  This is what the MODEL did,
            and it must be <= step/2, full stop.
      f32 : what a reader that decodes into float32 actually sees.  That adds
            one float32 rounding of (packed*scale + offset), so the bound is
            step/2 + eps32(|value|).  For th (~300 K) that extra term is
            ~3e-5 K -- negligible, but real, and worth measuring rather than
            hiding behind a fudged tolerance.
    """
    fs = files(ll_case, pattern)
    if not fs:
        report(False, label, "NO OUTPUT FILES (run failed?)")
        return
    size_ll = sum(os.path.getsize(f"{ROOT}/{ll_case}/{f}") for f in fs)
    size_pk = sum(os.path.getsize(f"{ROOT}/{pk_case}/{f}") for f in files(pk_case, pattern))
    e64, e32, types, mag = {}, {}, {}, {}
    for fn in fs:
        a = nc.Dataset(f"{ROOT}/{ll_case}/{fn}")
        b = nc.Dataset(f"{ROOT}/{pk_case}/{fn}")
        b.set_auto_maskandscale(False)
        try:
            for v in a.variables:
                if v not in b.variables:
                    report(False, f"{label} [{fn}]", f"{v} missing from packed run")
                    continue
                ref = np.asarray(a[v][:], dtype="f8")
                raw = np.asarray(b[v][:])
                types[v] = raw.dtype
                if ref.shape != raw.shape:
                    report(False, f"{label} [{fn}]", f"{v} shape mismatch")
                    continue
                if v in table:
                    sc = np.float64(b[v].scale_factor)
                    off = np.float64(b[v].add_offset)
                    d64 = raw.astype("f8") * sc + off
                    d32 = (raw.astype("f4") * np.float32(sc)
                           + np.float32(off)).astype("f8")
                else:
                    d64 = d32 = ref if raw.dtype.kind == "f" else raw.astype("f8")
                    if raw.dtype.kind == "f":
                        d64 = d32 = raw.astype("f8")
                e64[v] = max(e64.get(v, 0.0), float(np.max(np.abs(d64 - ref))) if ref.size else 0.0)
                e32[v] = max(e32.get(v, 0.0), float(np.max(np.abs(d32 - ref))) if ref.size else 0.0)
                mag[v] = max(mag.get(v, 0.0), float(np.max(np.abs(ref))) if ref.size else 0.0)
        finally:
            a.close(); b.close()

    for v in sorted(e64):
        if v in table:
            step, _off, nbytes = table[v]
            want = np.dtype("i2") if nbytes == 2 else np.dtype("i4")
            report(types[v] == want, f"{label}: {v} stored as {want}", f"got {types[v]}")
            report(e64[v] <= step / 2,
                   f"{label}: {v} f64 decode within half-step ({step:g})",
                   f"err {e64[v]:.4e} <= {step/2:.4e}")
            # float32 readers additionally carry one rounding at this magnitude
            tol32 = step / 2 + np.spacing(np.float32(mag[v])) * 2
            report(e32[v] <= tol32,
                   f"{label}: {v} f32 decode within half-step + eps32",
                   f"err {e32[v]:.4e} <= {tol32:.4e}")
        else:
            report(types[v].kind == "f", f"{label}: {v} left as float", f"got {types[v]}")
            report(e64[v] == 0.0, f"{label}: {v} unpacked var is exact",
                   "" if e64[v] == 0.0 else f"max abs err {e64[v]:.3e}")
    if size_pk:
        print(f"         size: {size_ll/1e6:.2f} MB float -> {size_pk/1e6:.2f} MB packed "
              f"= {size_ll/size_pk:.2f}x")


def check_parcel_ids(case):
    """Item 2: xh is compressed now -- it must still be exactly 1..nparcels."""
    for fn in files(case, "cm1out_pdata_*.nc"):
        d = nc.Dataset(f"{ROOT}/{case}/{fn}")
        try:
            xh = np.asarray(d["xh"][:], dtype="f8")
            n = xh.size
            ok = np.array_equal(xh, np.arange(1, n + 1, dtype="f8"))
            filt = d["xh"].filters() or {}
            report(ok, f"xh is exactly 1..{n} [{fn}]")
            report(bool(filt.get("zlib")), f"xh is deflated [{fn}]", str(filt))
        finally:
            d.close()


print()
print("=" * 72)
print("  5. FIXED-POINT PACKING")
print("=" * 72)
for fn in files("pack_ft4", "cm1out_0*.nc"):
    ok, why = bitwise_identical(f"{ROOT}/pack_ft4/{fn}", f"{ROOT}/pack_ft2/{fn}")
    report(ok, f"gridded ft4 == ft2 when packed [{fn}]", why)
for fn in files("pack_ft4", "cm1out_pdata_*.nc"):
    ok, why = bitwise_identical(f"{ROOT}/pack_ft4/{fn}", f"{ROOT}/pack_ft2/{fn}")
    report(ok, f"parcel parallel == serial when packed [{fn}]", why)

print("--- gridded, even ---")
check_pack("prcl_ft4_lossless", "pack_ft4", "cm1out_0*.nc", PACK_GRID, "gridded pack")
print("--- parcels, even ---")
check_pack("prcl_ft4_lossless", "pack_ft4", "cm1out_pdata_*.nc", PACK_PRCL, "parcel pack")
print("--- gridded, uneven decomposition ---")
check_pack("uneven_ft4_lossless", "pack_uneven_ft4", "cm1out_0*.nc", PACK_GRID, "gridded pack uneven")
print("--- parcel id coordinate (item 2) ---")
check_parcel_ids("pack_ft4")
check_parcel_ids("prcl_ft2_lossless")

print()
print("=" * 72)
print("  ALL CHECKS PASSED" if nfail == 0 else f"  FAILURES: {nfail}")
print("=" * 72)
sys.exit(1 if nfail else 0)
