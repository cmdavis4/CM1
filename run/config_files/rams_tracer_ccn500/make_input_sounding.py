"""Write a CM1 input_sounding from the horizontal-mean RAMS t=0 state.

Usage: python make_input_sounding.py <a-A t=0 .h5> <-head.txt> <output path>
"""
import sys

import h5py
import numpy as np

analysis_path, header_path, output_path = sys.argv[1:4]

# RAMS thermodynamic constants (rconstants.f90)
rams_cp = 1004.0
rams_rgas = 287.0
reference_pressure_hpa = 1000.0

# Read the scalar-level heights from the header
header_lines = open(header_path).read().split("\n")
ztn_start = header_lines.index("__ztn01")
n_levels = int(header_lines[ztn_start + 1])
level_heights = np.array(
    [float(value) for value in header_lines[ztn_start + 2 : ztn_start + 2 + n_levels]]
)

# Horizontal means over the interior (drop the boundary rows)
with h5py.File(analysis_path) as analysis_file:
    interior = (slice(None), slice(1, -1), slice(1, -1))
    theta_mean = analysis_file["THETA"][interior].mean(axis=(1, 2))
    vapor_mean = analysis_file["RV"][interior].mean(axis=(1, 2))
    exner_mean = analysis_file["PI"][interior].mean(axis=(1, 2))
    u_mean = analysis_file["UP"][interior].mean(axis=(1, 2))
    v_mean = analysis_file["VP"][interior].mean(axis=(1, 2))

# Pressure from the RAMS Exner function (PI = cp * (p/p00)^(R/cp))
pressure_hpa = reference_pressure_hpa * (exner_mean / rams_cp) ** (rams_cp / rams_rgas)

# Surface values: midpoint of the below-ground level (k=1) and the first real level (k=2)
surface_pressure_hpa = 0.5 * (pressure_hpa[0] + pressure_hpa[1])
surface_theta = 0.5 * (theta_mean[0] + theta_mean[1])
surface_vapor = 0.5 * (vapor_mean[0] + vapor_mean[1])

# Write the CM1 format: surface line, then z(m) theta(K) qv(g/kg) u v
with open(output_path, "w") as sounding_file:
    sounding_file.write(
        f"{surface_pressure_hpa:10.3f} {surface_theta:10.4f} {surface_vapor * 1e3:10.5f}\n"
    )
    # Skip the below-ground level k=1
    for k in range(1, n_levels):
        sounding_file.write(
            f"{level_heights[k]:10.2f} {theta_mean[k]:10.4f} {vapor_mean[k] * 1e3:10.5f}"
            f" {u_mean[k]:8.3f} {v_mean[k]:8.3f}\n"
        )

print("surface p, theta, qv:", surface_pressure_hpa, surface_theta, surface_vapor * 1e3)
print("levels:", n_levels - 1, "top z:", level_heights[-1])
print("max |u|,|v|:", np.abs(u_mean).max(), np.abs(v_mean).max())
