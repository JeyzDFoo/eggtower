"""
Configuration constants for the Egg Tower structural analysis.
"""
from materials import BFRP, ARAMID_CABLE, CABLE_MATERIALS

# Tower geometry constraints
TOTAL_HEIGHT = 1600  # m
D_BASE = 30  # m - diameter at bottom
D_TOP = 5    # m - minimum diameter at top
ASPECT_RATIO = 5  # h/d ratio for each egg
THICKNESS_FACTOR = 1.40  # Shell thickness multiplier (1.40 needed for SF ≥ 2.0 with guy forces)

# Wind parameters
RHO_AIR = 1.225  # kg/m³ (air density at sea level)
V_REF = 75  # m/s (reference wind speed at 10m height) = 270 km/h
CD = 1.2  # drag coefficient for ellipsoid
CD_CABLE = 1.0  # drag coefficient for cables (cylinder)
CL_VORTEX = 0.4  # lift coefficient for vortex shedding oscillations
STROUHAL = 0.2  # Strouhal number for ellipsoid

# Cable distribution
N_CABLES_PER_LEVEL = 20  # number of cable attachment points around circumference

# Default materials
SHELL_MATERIAL = BFRP
CABLE_MATERIAL = ARAMID_CABLE  # Kevlar 49 Parafil - no creep, dielectric, proven for guy wires

# Legacy compatibility - extract values from BFRP material
RHO = BFRP.density
SIGMA_MAX = BFRP.compressive_strength
TENSILE_STRENGTH = BFRP.tensile_strength
E_MODULUS = BFRP.youngs_modulus
POISSON = BFRP.poisson_ratio
