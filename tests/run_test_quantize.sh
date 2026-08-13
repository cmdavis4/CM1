#!/bin/bash
# Build and run the standalone unit test for src/quantize.F.
#
# Deliberately does NOT use src/Makefile: this compiles only input.F and
# quantize.F, so it runs in seconds and needs no MPI, netcdf, or namelist.
# Matches the production build's compiler and flags, minus -DMPI/-DNETCDF.

set -e

cd "$(dirname "$0")"
SRC=../src
BUILD=./build

module load intel/20.0.4 2>/dev/null || true

FC=${FC:-ifort}
OPTS="-O3 -ip -assume byterecl -fp-model precise -ftz -no-fma"

rm -rf "$BUILD"
mkdir -p "$BUILD"

# .F files are preprocessed to .f90 the same way src/Makefile does it
for f in input quantize ; do
  cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding "$SRC/$f.F" > "$BUILD/$f.f90"
done

cd "$BUILD"
$FC $OPTS -c input.f90
$FC $OPTS -c quantize.f90
$FC $OPTS -c ../test_quantize.F90 -o test_quantize.o
$FC $OPTS input.o quantize.o test_quantize.o -o test_quantize

./test_quantize
