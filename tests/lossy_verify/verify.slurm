#!/bin/bash
#SBATCH -A atm124
#SBATCH -p batch_short
#SBATCH -N 1
#SBATCH -J lossy_verify
#SBATCH -t 00:30:00
#SBATCH --exclusive
#SBATCH --output=/gpfs/wolf2/arm/atm124/proj-shared/cmdavis4/projects/dc_entrainment/sc/_lossy_verify/verify.slurm.out
#SBATCH --error=/gpfs/wolf2/arm/atm124/proj-shared/cmdavis4/projects/dc_entrainment/sc/_lossy_verify/verify.slurm.err

module load intel/20.0.4
module load openmpi/4.0.4
module load netcdf-c
module load netcdf-fortran

# required for the var16=1 parallel collective HDF5 writes (see docs/claude/parcels.md)
export OMPI_MCA_io=romio321
export OMPI_MCA_coll_hcoll_enable=0

ROOT=/gpfs/wolf2/arm/atm124/proj-shared/cmdavis4/projects/dc_entrainment/sc/_lossy_verify

for d in regress_ft2 regress_ft4 \
         prcl_ft4_lossless prcl_ft4_lossy \
         prcl_ft2_lossless prcl_ft2_lossy \
         uneven_ft4_lossless uneven_ft4_lossy \
         pack_ft4 pack_ft2 pack_uneven_ft4 ; do
  echo "############### RUN $d ###############"
  cd "$ROOT/$d"
  srun -n 4 ./cm1.exe > cm1.out 2>&1
  echo "EXIT($d)=$?"
  ls -la cm1out_*.nc 2>/dev/null
done
echo "ALL_DONE"
