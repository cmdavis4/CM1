#!/bin/bash
# Build the lossy-compression verification run directories on scratch.
#
# Eight tiny 4-rank runs, in four groups:
#
#   regress_ft2 / regress_ft4     iprcl=0, lossless.  Must reproduce the
#                                 pre-feature reference output in
#                                 _ftype4_verify/ BIT FOR BIT.  This is the
#                                 test that matters most: it proves the
#                                 feature is genuinely off by default.
#
#   prcl_ft4_{lossless,lossy}     parcels via the PARALLEL collective writer
#                                 (var16=1) + gridded filetype=4
#
#   prcl_ft2_{lossless,lossy}     parcels via the SERIAL writer (var16=0) +
#                                 gridded filetype=2 (rank-0 gather).  Pairs
#                                 with the above to show quantization gives
#                                 the same answer on both write paths.
#
#   uneven_ft4_{lossless,lossy}   nx=ny=41 on 4 ranks, so tiles are ragged
#                                 and the staggered-grid edge handling in the
#                                 filetype=4 writer is exercised.
#
# Reference for the harness layout: sc/_ftype4_verify (the filetype=4
# bitwise-identical verification this mirrors).

set -e

ROOT=/gpfs/wolf2/arm/atm124/proj-shared/cmdavis4/projects/dc_entrainment/sc/_lossy_verify
REF=/gpfs/wolf2/arm/atm124/proj-shared/cmdavis4/projects/dc_entrainment/sc/_ftype4_verify
SRCDIR="$(cd "$(dirname "$0")/../../src" && pwd)"
EXE="$SRCDIR/../run/cm1.exe"

if [ ! -x "$EXE" ]; then
  echo "no cm1.exe at $EXE -- build it first (cd src && make)" >&2
  exit 1
fi

# parcel settings for the parcel-bearing cases.  var18=1 (mass-weighted by
# COUNT) rather than the production var18=2 (by density): the count is then
# exactly nparcels, which keeps the two runs in a lossless/lossy pair
# trivially comparable.
NPARCELS=200000

mkdir -p "$ROOT"

make_case () {
  local name=$1 filetype=$2 nxy=$3 iprcl=$4 var16=$5 gdig=$6 pdig=$7
  local d="$ROOT/$name"

  rm -rf "$d"
  mkdir -p "$d"
  cp "$REF/even_ft2/input_sounding" "$d/"
  cp "$EXE" "$d/"
  cp "$REF/even_ft2/namelist.input" "$d/namelist.input"

  local n="$d/namelist.input"

  sed -i "s/^ nx           =.*/ nx           =      $nxy,/" "$n"
  sed -i "s/^ ny           =.*/ ny           =      $nxy,/" "$n"
  sed -i "s/^ output_filetype  =.*/ output_filetype  = $filetype,/" "$n"

  # parcels
  sed -i "s/^ iprcl     =.*/ iprcl     =  $iprcl,/" "$n"
  sed -i "s/^ nparcels  =.*/ nparcels  =  $NPARCELS,/" "$n"
  sed -i "s/^ prclfrq =.*/ prclfrq = 1,/" "$n"
  sed -i "s/^ var16     =.*/ var16     =   $var16,/" "$n"
  sed -i "s/^ var17     =.*/ var17     =   1,/" "$n"
  sed -i "s/^ var18     =.*/ var18     =   1,/" "$n"

  # the compression knobs are new, so they are not in the reference namelist:
  # insert them into &param9 rather than substituting.
  sed -i "/^ output_filetype  =/a\\ output_deflate_level = 1,\\n output_signif_digits = $gdig,\\n prcl_signif_digits   = $pdig," "$n"

  echo "  $name: filetype=$filetype nx=ny=$nxy iprcl=$iprcl var16=$var16 gdig=$gdig pdig=$pdig"
}

echo "building cases under $ROOT"

#          name                 ftype nxy iprcl var16 gdig pdig
make_case  regress_ft2              2  40     0     1    0    0
make_case  regress_ft4              4  40     0     1    0    0
make_case  prcl_ft4_lossless        4  40     1     1    0    0
make_case  prcl_ft4_lossy           4  40     1     1    3    3
make_case  prcl_ft2_lossless        2  40     1     0    0    0
make_case  prcl_ft2_lossy           2  40     1     0    3    3
make_case  uneven_ft4_lossless      4  41     1     1    0    0
make_case  uneven_ft4_lossy         4  41     1     1    3    3

cp "$(dirname "$0")/verify.slurm" "$ROOT/"
cp "$(dirname "$0")/compare.py"   "$ROOT/"

echo
echo "submit with:  sbatch $ROOT/verify.slurm"
