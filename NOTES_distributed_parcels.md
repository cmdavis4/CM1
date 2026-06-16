# Distributed parcels for CM1 (scaling to ~1e8 parcels)

## Problem
The stock parcel scheme **replicates `pdata(nparcels,npvals)` on every MPI rank** and does a
global `MPI_ALLREDUCE` over `3*nparcels` every timestep ([parcel.F:627](src/parcel.F#L627)).
Memory per rank = `nparcels*npvals*4` bytes *regardless of rank count* → at 1e8 parcels,
npvals=28 that's ~11 GB/rank × ppnode → OOM hang (observed at 1e8). Per-step global comm is a
second wall.

## Goal
Each rank stores **only the parcels physically inside its subdomain**, advects them locally, and
**migrates** parcels to a neighbor rank when they cross a subdomain boundary. No global parcel
array, no per-step global all-reduce. Target: run 1e7–1e8 parcels.

This is **parcel-internal**: it only *reads* the Eulerian wind/density fields (unchanged). It
*frees* memory currently consumed by replicated parcel arrays — net positive for the rest of the
model. No correctness/perf repercussion for the Eulerian solver.

## User decisions (locked)
- **Output:** keep current architecture; add a flag to switch **serial gather→rank0 write**
  (existing `writepdata_nc`, preserves format+deflate) vs **parallel collective write** (their
  HDF5 has parallel compression). Flag = flex var.
- **Restart:** keep working — gather to rank0 on write (same restart format), scatter on read.

## Flags (flex vars; default 0 = unchanged legacy behavior)
- `var17` — parcel architecture: 0 = legacy replicate-all (current), 1 = distributed.
- `var16` — output mode when distributed: 0 = serial gather→rank0, 1 = parallel collective.
- Legacy path is kept as the **default and as a verification oracle** (distributed must
  reproduce legacy trajectories on small cases).

## Key existing facts
- Ownership predicate (reuse verbatim): parcel (x,y) is local iff it maps to cell
  `iflag∈[1,ni], jflag∈[1,nj]` within `[xf(1),xf(ni+1)]×[yf(1),yf(nj+1)]`, with the
  boundary tie-break at [parcel.F:186-196](src/parcel.F#L186-L196) (`myi/myj/nodex/nodey`).
- 2D decomposition: `nodex,nodey,myi,myj`. Neighbor ranks incl. corners:
  `mywest,myeast,mysouth,mynorth,mysw,mynw,myne,myse` (in `input` module).
- Per-step update: [parcel_driver](src/parcel.F#L11) (called from [solve3.F:768](src/solve3.F#L768)).
- Diagnostics at write: [parcel_interp](src/parcel.F#L657) → `MPI_REDUCE` to rank0
  ([parcel.F:1611](src/parcel.F#L1611)).
- Serial writer: [parcel_write](src/parcel.F#L1627)/[writepdata_nc](src/writeout_nc.F#L1810).
- Restart parcels: [restart_write.F:1777](src/restart_write.F#L1777), restart_read mirrors.

## Data-structure design
New module state in `parcel_module` (active only when `var17=1`):
- `np_loc` — number of parcels currently owned by this rank.
- `np_max` — local capacity (allocated). Start `np_max ≈ headroom * nparcels/numprocs`
  (headroom ~4), **grow by realloc** if exceeded (handle clustering gracefully).
- `pdata(np_max,npvals)` — local parcels only (reuse name; legacy dims to `nparcels`).
- `pid(np_max)` — int global parcel ID (persistent identity for output ordering/trajectories).
- `nparcels` stays the **global total** (output dim, ID range, restart count).

## Phases
1. **Foundation** (buildable, inert): add flags, module state, and helpers
   `parcel_owner_rank(x,y)->rank` (global xfref/yfref → global cell → owner) and
   `parcel_is_local(x,y)` (reuse predicate). [this commit]
2. **Init** (`init3d`): when `var17=1`, each rank generates the same deterministic parcel set
   (mass-weighted, same seed) and **keeps only locals** (filter by `parcel_is_local`), assigning
   `pid = generation index`. No communication. Legacy path unchanged.
3. **Driver** (`parcel_driver`): loop over `1..np_loc` (all owned), advect (existing RK2);
   after the step, **migrate** parcels that left the subdomain via two-pass face exchange
   (x then y) using `mywest/myeast/mysouth/mynorth` (two passes handle diagonals; <1 cell/step
   by CFL). Pack (npvals+pid), Isend/Irecv counts then data, unpack/compact. Replace the global
   all-reduce. Preserve periodic wrap + 2D (nx=1/ny=1/axisymm) special cases.
4. **Diagnostics** (`parcel_interp`): compute for `1..np_loc` only (each local parcel is in
   this rank's subdomain so fields are valid); drop sentinel/reduce.
5. **Output**:
   - serial (`var16=0`): indexed `MPI_Gatherv` one variable-slab at a time into rank0 in `pid`
     order (rank0 holds ≤1 slab ≈ `nparcels*4` ≈ 400 MB @1e8), then existing `writepdata_nc`.
   - parallel (`var16=1`): redistribute parcels to contiguous global-ID blocks (Alltoallv by
     `pid/blocksize`), open file with `nf90_mpiio`+collective, chunk `(blocksize,1)` + deflate,
     each rank writes its hyperslab.
6. **Restart**: write = serial gather→rank0 (same format). read = rank0 reads all, scatters by
   `parcel_owner_rank` (Alltoallv).
7. **Verify**: small interactive `mpirun -np 4` run, tiny domain, compare distributed vs legacy
   parcel trajectories (should match to round-off); check global parcel-count conservation each
   step (sum np_loc == nparcels); then scale test.

## Verification invariants to assert in code (debug)
- `SUM(np_loc) == nparcels` every step (no lost/duplicated parcels).
- every local parcel passes `parcel_is_local` after migration.
