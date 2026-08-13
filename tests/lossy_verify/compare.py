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
    """Compare every data variable's raw bits."""
    a, b = nc.Dataset(fa), nc.Dataset(fb)
    try:
        if set(a.variables) != set(b.variables):
            return False, f"variable sets differ: {set(a.variables) ^ set(b.variables)}"
        for v in a.variables:
            da = np.asarray(a[v][:])
            db = np.asarray(b[v][:])
            if da.shape != db.shape:
                return False, f"{v} shape {da.shape} vs {db.shape}"
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

print()
print("=" * 72)
print("  ALL CHECKS PASSED" if nfail == 0 else f"  FAILURES: {nfail}")
print("=" * 72)
sys.exit(1 if nfail else 0)
