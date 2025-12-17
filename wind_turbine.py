"""
Wind Turbine Integration for Egg Tower.

Two turbine configurations are analyzed:

1. JUNCTION TURBINES (Horizontal Disk)
   - Placed at egg-to-egg junctions
   - Horizontal spinning disk captures vertical updraft
   - Limited by flow deflection losses

2. INTERNAL DUCT TURBINES (Egg-Integrated) *** PREFERRED ***
   - Uses the hollow egg interior as a wind concentrator
   - Like a reversed jet engine / Diffuser-Augmented Wind Turbine (DAWT)
   - Inlets on windward shell → accelerated flow → turbine → exhaust leeward

Internal Duct Physics:
---------------------
The egg shell acts as a natural convergent-divergent duct:

    [Inlet]  →  [Convergent Section]  →  [Throat/Turbine]  →  [Divergent/Exhaust]
    
    Wind enters     Shell curves        Minimum area,       Pressure recovery,
    at shell       inward toward       max velocity,       flow exits
    surface        centerline          turbine here        leeward side

By Bernoulli: v_throat = v_inlet × √(A_inlet / A_throat)
With proper design, can achieve 2-3x velocity amplification = 8-27x power!

Gyroscopic Stabilization:
------------------------
Turbines spinning around the vertical axis (core of egg) provide
angular momentum that resists lateral tower sway.
"""
import json
import os
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from config import RHO_AIR, V_REF


# Path to optimized tower configuration
OPTIMIZED_TOWER_PATH = os.path.join(os.path.dirname(__file__), 'optimized_tower.json')


def load_optimized_tower() -> dict:
    """
    Load the optimized tower configuration from JSON file.
    
    Returns:
        Dictionary with tower configuration including:
        - d_base, d_top, target_height
        - n_eggs
        - eggs: list of egg geometry dictionaries
    
    Raises:
        FileNotFoundError: If optimized_tower.json doesn't exist
    """
    if not os.path.exists(OPTIMIZED_TOWER_PATH):
        raise FileNotFoundError(
            f"Optimized tower file not found: {OPTIMIZED_TOWER_PATH}\n"
            f"Run main.py first to generate the optimized tower configuration."
        )
    
    with open(OPTIMIZED_TOWER_PATH, 'r') as f:
        config = json.load(f)
    
    return config


def get_optimized_eggs() -> List[dict]:
    """
    Get the list of eggs from the optimized tower configuration.
    
    Returns:
        List of egg dictionaries with geometry and stress data.
    """
    config = load_optimized_tower()
    return config['eggs']


@dataclass
class WindTurbine:
    """
    Specifications for a wind turbine.
    
    Can be configured for:
    - Internal duct mounting (vertical axis, inside egg core)
    - Junction mounting (horizontal disk at egg junction)
    """
    name: str
    rotor_diameter: float       # m - diameter of the rotor
    rated_power_kW: float       # kW - power at rated wind speed
    rated_wind_speed: float     # m/s - wind speed at rated power
    cut_in_speed: float         # m/s - minimum operational wind speed
    cut_out_speed: float        # m/s - maximum operational wind speed
    rotor_mass_kg: float        # kg - mass of rotor + hub
    nacelle_mass_kg: float      # kg - mass of nacelle (generator, gearbox)
    rpm_rated: float            # rpm - rotor speed at rated power
    efficiency: float = 0.40    # Cp (power coefficient, theoretical max = 0.593 Betz limit)
    

# Horizontal-disk turbine models (vertical axis, captures updrafts)
# These are sized to fit within the junction throat between eggs
TURBINE_MODELS = {
    'micro': WindTurbine(
        name="Micro Disk (5m)",
        rotor_diameter=5.0,
        rated_power_kW=15,
        rated_wind_speed=12,
        cut_in_speed=2.0,       # Lower cut-in for vertical flow
        cut_out_speed=40,       # Lower cut-out (updrafts less extreme)
        rotor_mass_kg=200,      # Heavier disk for gyroscopic effect
        nacelle_mass_kg=400,
        rpm_rated=150,
        efficiency=0.35
    ),
    'small': WindTurbine(
        name="Small Disk (10m)", 
        rotor_diameter=10.0,
        rated_power_kW=80,
        rated_wind_speed=12,
        cut_in_speed=2.0,
        cut_out_speed=40,
        rotor_mass_kg=1200,     # Heavier disk
        nacelle_mass_kg=3000,
        rpm_rated=100,
        efficiency=0.38
    ),
    'medium': WindTurbine(
        name="Medium Disk (20m)",
        rotor_diameter=20.0,
        rated_power_kW=400,
        rated_wind_speed=12,
        cut_in_speed=2.0,
        cut_out_speed=45,
        rotor_mass_kg=6000,
        nacelle_mass_kg=20000,
        rpm_rated=50,
        efficiency=0.40
    ),
    'large': WindTurbine(
        name="Large Disk (30m)",
        rotor_diameter=30.0,
        rated_power_kW=1200,
        rated_wind_speed=12,
        cut_in_speed=2.0,
        cut_out_speed=45,
        rotor_mass_kg=12000,
        nacelle_mass_kg=45000,
        rpm_rated=40,
        efficiency=0.42
    ),
}


# ============================================================================
# INTERNAL DUCT TURBINE SYSTEM (Egg as Wind Concentrator)
# ============================================================================

@dataclass
class InternalDuctTurbine:
    """
    Turbine designed to sit at the core of an egg, capturing accelerated duct flow.
    
    The egg shell acts as a natural convergent-divergent duct:
    - Inlets on windward side of shell
    - Flow converges toward central turbine
    - Turbine captures high-velocity throat flow
    - Flow diverges and exits leeward side
    """
    name: str
    rotor_diameter: float       # m - must fit inside egg core
    rated_power_kW: float       # kW - power at rated (accelerated) wind speed
    rated_wind_speed: float     # m/s - throat wind speed at rated power
    cut_in_speed: float         # m/s - minimum operational wind speed
    cut_out_speed: float        # m/s - maximum operational wind speed
    rotor_mass_kg: float        # kg - mass of rotor
    generator_mass_kg: float    # kg - mass of generator (mounted below rotor)
    rpm_rated: float            # rpm - rotor speed at rated power
    efficiency: float = 0.45    # Higher efficiency possible with ducted flow


# Internal duct turbines - sized for egg cores
INTERNAL_TURBINE_MODELS = {
    'micro_duct': InternalDuctTurbine(
        name="Micro Duct Turbine (1.5m)",
        rotor_diameter=1.5,
        rated_power_kW=15,
        rated_wind_speed=25,      # Expects accelerated flow
        cut_in_speed=5.0,
        cut_out_speed=60,
        rotor_mass_kg=80,
        generator_mass_kg=200,
        rpm_rated=500,
        efficiency=0.42
    ),
    'mini_duct': InternalDuctTurbine(
        name="Mini Duct Turbine (2m)",
        rotor_diameter=2.0,
        rated_power_kW=25,
        rated_wind_speed=25,
        cut_in_speed=5.0,
        cut_out_speed=60,
        rotor_mass_kg=120,
        generator_mass_kg=350,
        rpm_rated=400,
        efficiency=0.43
    ),
    'small_duct': InternalDuctTurbine(
        name="Small Duct Turbine (3m)",
        rotor_diameter=3.0,
        rated_power_kW=50,
        rated_wind_speed=25,      # Expects accelerated flow
        cut_in_speed=5.0,
        cut_out_speed=60,
        rotor_mass_kg=200,
        generator_mass_kg=500,
        rpm_rated=300,
        efficiency=0.45
    ),
    'medium_duct': InternalDuctTurbine(
        name="Medium Duct Turbine (5m)",
        rotor_diameter=5.0,
        rated_power_kW=150,
        rated_wind_speed=25,
        cut_in_speed=5.0,
        cut_out_speed=60,
        rotor_mass_kg=500,
        generator_mass_kg=1200,
        rpm_rated=200,
        efficiency=0.45
    ),
    'large_duct': InternalDuctTurbine(
        name="Large Duct Turbine (8m)",
        rotor_diameter=8.0,
        rated_power_kW=400,
        rated_wind_speed=25,
        cut_in_speed=5.0,
        cut_out_speed=60,
        rotor_mass_kg=1200,
        generator_mass_kg=3000,
        rpm_rated=120,
        efficiency=0.45
    ),
    'xlarge_duct': InternalDuctTurbine(
        name="XL Duct Turbine (12m)",
        rotor_diameter=12.0,
        rated_power_kW=800,
        rated_wind_speed=25,
        cut_in_speed=5.0,
        cut_out_speed=60,
        rotor_mass_kg=2500,
        generator_mass_kg=6000,
        rpm_rated=80,
        efficiency=0.45
    ),
}


def calculate_egg_duct_geometry(egg: dict, inlet_fraction: float = 0.3) -> dict:
    """
    Calculate the internal duct geometry for an egg used as a wind concentrator.
    
    The egg shell has:
    - Inlets on windward side (fraction of surface area)
    - Convergent section toward center
    - Throat at the core (where turbine sits)
    - Divergent section to leeward exhaust
    
    Parameters:
        egg: egg properties dictionary
        inlet_fraction: fraction of windward surface used as inlet (0.3 = 30%)
    
    Returns:
        Dictionary with duct geometry
    """
    d_outer = egg['d']  # max outer diameter
    h = egg['h']  # height
    t_avg = egg.get('t_avg', d_outer / 250)  # shell thickness
    
    # Internal diameter at equator (inside the shell)
    d_internal = d_outer - 2 * t_avg
    r_internal = d_internal / 2
    
    # Egg internal cross-sectional area at equator (widest point)
    A_equator = np.pi * r_internal**2
    
    # Inlet area (windward half of shell, fraction open)
    # Approximate as fraction of the projected internal frontal area
    A_inlet = inlet_fraction * A_equator
    
    # Throat diameter (where turbine sits, inside the egg core)
    # The throat ratio scales inversely with egg size:
    # - Large eggs (>50m): 10% throat - plenty of room for big turbines
    # - Medium eggs (20-50m): 12% throat - need proportionally larger core
    # - Small eggs (<20m): 15% throat - maximize turbine size in limited space
    if d_internal > 50:
        throat_ratio = 0.10
    elif d_internal > 20:
        throat_ratio = 0.12
    else:
        throat_ratio = 0.15
    
    d_throat = d_internal * throat_ratio
    A_throat = np.pi * (d_throat / 2)**2
    
    # Area ratio determines velocity amplification
    # v_throat / v_inlet = A_inlet / A_throat (continuity)
    area_ratio = A_inlet / A_throat if A_throat > 0 else 1.0
    
    # Velocity amplification (with duct efficiency losses)
    # Real ducts achieve 60-80% of theoretical due to friction, separation
    duct_efficiency = 0.70
    v_amplification = np.sqrt(area_ratio) * duct_efficiency  # sqrt because A ∝ v for ducts
    
    # For a well-designed DAWT (Diffuser Augmented Wind Turbine)
    # amplification of 1.5-2.5x is achievable
    v_amplification = min(v_amplification, 2.5)  # Cap at realistic limit
    v_amplification = max(v_amplification, 1.0)  # Floor at 1.0
    
    return {
        'egg_d_outer': d_outer,
        'egg_d_internal': d_internal,
        'egg_h': h,
        't_shell': t_avg,
        'A_equator_m2': A_equator,
        'inlet_fraction': inlet_fraction,
        'A_inlet_m2': A_inlet,
        'd_throat_m': d_throat,
        'A_throat_m2': A_throat,
        'area_ratio': area_ratio,
        'duct_efficiency': duct_efficiency,
        'v_amplification': v_amplification,
        'power_amplification': v_amplification**3,
    }


def calculate_internal_turbine_power(
    turbine: InternalDuctTurbine,
    wind_speed_ambient: float,
    v_amplification: float
) -> dict:
    """
    Calculate power output for an internal duct turbine.
    
    The turbine receives accelerated flow from the egg duct.
    """
    # Throat wind speed
    v_throat = wind_speed_ambient * v_amplification
    
    # Check operational limits
    if v_throat < turbine.cut_in_speed:
        return {
            'status': 'below_cut_in',
            'power_kW': 0,
            'v_ambient': wind_speed_ambient,
            'v_throat': v_throat,
        }
    elif v_throat > turbine.cut_out_speed:
        return {
            'status': 'above_cut_out',
            'power_kW': 0,
            'v_ambient': wind_speed_ambient,
            'v_throat': v_throat,
        }
    
    # Swept area
    A = np.pi * (turbine.rotor_diameter / 2)**2
    
    # Available wind power at throat
    P_wind = 0.5 * RHO_AIR * A * v_throat**3
    
    # Extracted power (with efficiency)
    P_extract = P_wind * turbine.efficiency
    
    # Cap at rated power
    if P_extract > turbine.rated_power_kW * 1000:
        P_extract = turbine.rated_power_kW * 1000
        status = 'rated'
    else:
        status = 'variable'
    
    return {
        'status': status,
        'power_W': P_extract,
        'power_kW': P_extract / 1000,
        'v_ambient': wind_speed_ambient,
        'v_throat': v_throat,
        'P_wind_kW': P_wind / 1000,
        'capacity_factor': P_extract / (turbine.rated_power_kW * 1000),
    }


def design_internal_duct_system(
    egg: dict,
    inlet_fraction: float = 0.25
) -> dict:
    """
    Design an internal duct turbine system for a single egg.
    
    The egg acts as a wind concentrator with:
    - Shell openings on windward side (inlets)
    - Internal convergent duct toward core
    - Turbine at the throat (core)
    - Divergent section and leeward exhaust
    
    Parameters:
        egg: egg properties dictionary
        inlet_fraction: fraction of windward shell used as inlet
    
    Returns:
        Dictionary with internal duct system design
    """
    # Calculate duct geometry
    duct = calculate_egg_duct_geometry(egg, inlet_fraction)
    
    # Select appropriate turbine size for this egg
    # Turbine diameter should be ~80% of throat diameter
    target_turbine_d = duct['d_throat_m'] * 0.8
    
    # Find best matching turbine
    best_turbine = None
    best_match = float('inf')
    for key, turbine in INTERNAL_TURBINE_MODELS.items():
        match = abs(turbine.rotor_diameter - target_turbine_d)
        if match < best_match and turbine.rotor_diameter <= duct['d_throat_m']:
            best_match = match
            best_turbine = (key, turbine)
    
    if best_turbine is None:
        # Throat too small for any turbine
        return {
            'feasible': False,
            'reason': f"Throat diameter {duct['d_throat_m']:.1f}m too small for turbines",
            'egg_index': egg.get('index', -1),
            'egg_d': egg['d'],
            'duct': duct,
        }
    
    turbine_key, turbine = best_turbine
    
    # Wind speed at egg center height
    z_center = egg['z_base'] + egg['h'] / 2
    v_ambient = wind_speed_at_height(z_center)
    
    # Calculate power with duct amplification
    power = calculate_internal_turbine_power(turbine, v_ambient, duct['v_amplification'])
    
    # Gyroscopic moment (turbine spins around vertical axis in core)
    r = turbine.rotor_diameter / 2
    omega = turbine.rpm_rated * 2 * np.pi / 60
    I_rotor = 0.5 * turbine.rotor_mass_kg * r**2
    L = I_rotor * omega
    sway_rate = 0.01  # rad/s reference
    tau_gyro = sway_rate * L
    
    # Total mass
    total_mass = turbine.rotor_mass_kg + turbine.generator_mass_kg
    
    return {
        'feasible': True,
        'egg_index': egg.get('index', -1),
        'egg_d': egg['d'],
        'egg_h': egg['h'],
        'z_center': z_center,
        'turbine_model': turbine.name,
        'turbine_diameter': turbine.rotor_diameter,
        'd_throat': duct['d_throat_m'],
        'inlet_fraction': inlet_fraction,
        'v_ambient': v_ambient,
        'v_amplification': duct['v_amplification'],
        'v_throat': v_ambient * duct['v_amplification'],
        'power_kW': power['power_kW'],
        'power_status': power['status'],
        'mass_kg': total_mass,
        'gyro_torque_kNm': tau_gyro / 1000,
        'duct': duct,
        'power_detail': power,
    }


def analyze_internal_duct_integration(
    eggs: List[dict],
    inlet_fraction: float = 0.25,
    min_egg_diameter: float = 10.0
) -> dict:
    """
    Analyze internal duct turbine integration for all eggs in the tower.
    
    Each egg that's large enough gets an internal duct turbine system.
    
    Parameters:
        eggs: list of analyzed eggs
        inlet_fraction: fraction of windward shell used as inlet
        min_egg_diameter: minimum egg diameter for turbine installation (default 10m)
    
    Returns:
        Comprehensive analysis dictionary
    """
    installations = []
    total_power_kW = 0
    total_mass_kg = 0
    total_gyro_torque_kNm = 0
    
    for egg in eggs:
        # Check if egg is large enough
        if egg['d'] < min_egg_diameter:
            continue
        
        # Design internal duct system for this egg
        design = design_internal_duct_system(egg, inlet_fraction)
        
        if not design['feasible']:
            continue
        
        installations.append(design)
        total_power_kW += design['power_kW']
        total_mass_kg += design['mass_kg']
        total_gyro_torque_kNm += design['gyro_torque_kNm']
    
    # Annual energy estimate (35% capacity factor - higher due to duct amplification)
    hours_per_year = 8760
    capacity_factor = 0.35
    annual_energy_MWh = total_power_kW * hours_per_year * capacity_factor / 1000
    
    # Revenue estimate ($50/MWh wholesale)
    revenue_per_MWh = 50
    annual_revenue_usd = annual_energy_MWh * revenue_per_MWh
    
    return {
        'system_type': 'Internal Duct (Egg as Wind Concentrator)',
        'n_installations': len(installations),
        'total_power_kW': total_power_kW,
        'total_power_MW': total_power_kW / 1000,
        'total_mass_kg': total_mass_kg,
        'total_mass_tonnes': total_mass_kg / 1000,
        'total_gyro_torque_kNm': total_gyro_torque_kNm,
        'annual_energy_MWh': annual_energy_MWh,
        'annual_revenue_usd': annual_revenue_usd,
        'capacity_factor': capacity_factor,
        'installations': installations,
    }


def print_internal_duct_analysis(analysis: dict, detailed: bool = True):
    """Print internal duct turbine analysis."""
    print(f"\n{'='*90}")
    print("INTERNAL DUCT TURBINE ANALYSIS (Egg as Wind Concentrator)")
    print("Each egg acts like a jet engine intake - accelerating wind to central turbine")
    print(f"{'='*90}")
    
    print(f"\n📊 SUMMARY")
    print(f"   System Type:          {analysis['system_type']}")
    print(f"   Eggs with Turbines:   {analysis['n_installations']}")
    print(f"   Total Rated Power:    {analysis['total_power_MW']:.2f} MW")
    print(f"   Added Mass:           {analysis['total_mass_tonnes']:.1f} tonnes")
    
    print(f"\n⚡ ENERGY PRODUCTION")
    print(f"   Capacity Factor:      {analysis['capacity_factor']*100:.0f}%")
    print(f"   Annual Energy:        {analysis['annual_energy_MWh']:.0f} MWh")
    print(f"   Annual Revenue:       ${analysis['annual_revenue_usd']:,.0f}")
    
    print(f"\n🌀 GYROSCOPIC STABILIZATION")
    print(f"   Total Gyro Torque:    {analysis['total_gyro_torque_kNm']:.1f} kN·m")
    print(f"   (turbines spin around vertical axis in egg cores)")
    
    if detailed and analysis['installations']:
        print(f"\n{'─'*90}")
        print("EGG-BY-EGG INSTALLATION DETAILS")
        print(f"{'─'*90}")
        print(f"{'Egg':<5} {'Height':<10} {'Egg Ø':<10} {'Throat':<10} {'Turb Ø':<10} {'V_amb':<8} {'V_thr':<8} {'Power':<10} {'Status'}")
        print(f"{'#':<5} {'(m)':<10} {'(m)':<10} {'(m)':<10} {'(m)':<10} {'(m/s)':<8} {'(m/s)':<8} {'(kW)':<10}")
        print(f"{'─'*90}")
        
        for inst in analysis['installations']:
            print(f"{inst['egg_index']:<5} "
                  f"{inst['z_center']:<10.0f} "
                  f"{inst['egg_d']:<10.1f} "
                  f"{inst['d_throat']:<10.1f} "
                  f"{inst['turbine_diameter']:<10.1f} "
                  f"{inst['v_ambient']:<8.1f} "
                  f"{inst['v_throat']:<8.1f} "
                  f"{inst['power_kW']:<10.1f} "
                  f"{inst['power_status']}")


def wind_speed_at_height(z: float, v_ref: float = None, z_ref: float = 10) -> float:
    """
    Wind speed varies with height according to power law.
    v(z) = v_ref × (z/z_ref)^α where α ≈ 0.16 for open terrain
    
    Note: For power generation, we use AVERAGE operating wind speed (~10 m/s at 10m),
    not the extreme design wind speed (75 m/s) used for structural analysis.
    """
    if v_ref is None:
        # Use average operating wind for power calculations, not extreme design wind
        v_ref = 10.0  # m/s - typical average wind for power generation
    alpha = 0.16
    return v_ref * (z / z_ref) ** alpha


def calculate_junction_geometry(egg_lower: dict, egg_upper: dict) -> dict:
    """
    Calculate the geometry at the junction between two eggs.
    
    At the junction point, both eggs taper to smaller diameters:
    - Lower egg tapers toward its top (pole)
    - Upper egg tapers toward its base (pole)
    
    The gap between them is where wind can flow through.
    
    Parameters:
        egg_lower: properties of the lower egg
        egg_upper: properties of the upper egg
    
    Returns:
        Dictionary with junction geometry
    """
    # For an ellipsoidal egg, the diameter at the poles is effectively 0
    # but we need the clearance between outer shell and any central structure
    
    # Junction height (where the two eggs meet)
    z_junction = egg_lower['z_top']
    
    # The eggs are stacked, so at the junction:
    # - Lower egg has diameter approaching 0 at its top pole
    # - Upper egg has diameter approaching 0 at its bottom pole
    
    # The "funnel" effect comes from the curved surfaces on either side
    # Think of it as wind flowing around the equators and converging
    
    # Approximate the "throat" area as an annular ring
    # between the egg shells at their narrowest point
    
    # The effective flow area at the junction is the gap
    # For stacked eggs, this is essentially the open space minus shell thickness
    
    # Max diameter of lower egg (at its equator, ~60% up)
    d_lower = egg_lower['d']
    
    # Max diameter of upper egg  
    d_upper = egg_upper['d']
    
    # The junction is between the tip of the lower egg and base of upper
    # Approximate flow area at equator (reference)
    h_lower = egg_lower['h']
    h_upper = egg_upper['h']
    
    # For the lower egg, diameter at height z above its base:
    # Using ellipse equation: (z-a)²/a² + r²/b² = 1
    # At the top (z = h), r = 0
    # At equator (z ~ 0.6h for egg shape), r = d/2
    
    # Height of equator relative to junction (looking down from junction)
    z_equator_lower = 0.6 * h_lower  # equator position
    dist_to_equator_lower = h_lower - z_equator_lower  # distance from top to equator
    
    # Upper egg equator (looking up from junction)
    z_equator_upper = 0.6 * h_upper
    dist_to_equator_upper = z_equator_upper  # distance from base to equator
    
    return {
        'z_junction': z_junction,
        'd_lower': d_lower,
        'd_upper': d_upper,
        'h_lower': h_lower,
        'h_upper': h_upper,
        'd_junction_lower_equator': d_lower,  # reference diameter
        'd_junction_upper_equator': d_upper,
        'funnel_height_lower': dist_to_equator_lower,
        'funnel_height_upper': dist_to_equator_upper,
    }


def calculate_venturi_acceleration(
    d_outer: float, 
    d_throat: float,
    n_turbines: int = 1,
    turbine_diameter: float = 10.0,
    blockage_fraction: float = 0.7
) -> dict:
    """
    Calculate vertical updraft velocity through the egg junction throat.
    
    Physics: Horizontal wind hits the curved egg surface and is deflected.
    Some flow goes around the egg, some is forced upward through the junction.
    The constriction at the junction accelerates this vertical flow.
    
    Flow model:
    1. Horizontal wind approaches at velocity v_h
    2. Stagnation pressure converts to upward flow: v_up ≈ 0.3-0.5 × v_h
    3. Venturi acceleration through throat: v_throat = v_up × (A_junction/A_throat)
    
    For horizontal-disk turbines, we care about VERTICAL velocity through throat.
    
    Parameters:
        d_outer: egg diameter at the junction (m)
        d_throat: effective throat diameter where turbine sits (m)
        n_turbines: number of turbines (typically 1 large central disk)
        turbine_diameter: diameter of turbine rotor disk (m)
        blockage_fraction: fraction of throat area occupied by turbine (0.7 = 70%)
    
    Returns:
        Dictionary with venturi analysis for vertical flow
    """
    # Junction throat area (annular ring at the tip-to-base connection)
    # The throat is where the lower egg's tip meets the upper egg's base
    # Approximate as a circular opening with diameter proportional to egg size
    A_throat = np.pi * (d_throat / 2) ** 2
    
    # Turbine swept area (horizontal disk)
    A_turbine = np.pi * (turbine_diameter / 2) ** 2
    
    # Capture area - portion of approaching horizontal flow that gets deflected upward
    # This depends on the egg's frontal area and how much flow is redirected
    # Approximately half the egg's frontal area contributes to updraft
    A_capture = 0.5 * np.pi * (d_outer / 2) ** 2
    
    # Flow fraction that goes through throat (vs around the egg)
    # Depends on geometry - larger throat = more flow through
    throat_flow_fraction = min(0.4, A_throat / A_capture)
    
    # Vertical velocity as fraction of horizontal wind
    # Stagnation point flow + Venturi acceleration
    # Typical: 20-40% of horizontal wind converts to updraft through throat
    v_vertical_fraction = throat_flow_fraction * 0.8  # 80% conversion efficiency
    
    # Additional Venturi acceleration through the throat
    # If turbine blocks part of throat, remaining flow accelerates
    effective_throat = A_throat * (1 - blockage_fraction * 0.5)  # partial blockage
    area_ratio = A_throat / effective_throat if effective_throat > 0 else 1.0
    
    # Total velocity amplification (vertical wind / horizontal wind)
    v_amplification = v_vertical_fraction * min(area_ratio, 1.5)
    v_amplification = max(v_amplification, 0.15)  # Floor at 15% of horizontal wind
    
    return {
        'd_outer': d_outer,
        'd_throat': d_throat,
        'n_turbines': n_turbines,
        'turbine_diameter': turbine_diameter,
        'A_capture_m2': A_capture,
        'A_throat_m2': A_throat,
        'A_turbine_m2': A_turbine,
        'throat_flow_fraction': throat_flow_fraction,
        'area_ratio': area_ratio,
        'v_amplification': v_amplification,  # v_vertical / v_horizontal
        'power_factor': v_amplification ** 3,  # power scales with v³
    }


def calculate_turbine_power(
    turbine: WindTurbine,
    wind_speed: float,
    v_amplification: float = 1.0
) -> dict:
    """
    Calculate power output for a single turbine.
    
    P = 0.5 × ρ × A × v³ × Cp
    
    Where:
        ρ = air density
        A = rotor swept area
        v = wind speed
        Cp = power coefficient (efficiency)
    """
    # Effective wind speed with venturi amplification
    v_eff = wind_speed * v_amplification
    
    # Check operational limits
    if v_eff < turbine.cut_in_speed:
        return {
            'status': 'below_cut_in',
            'power_kW': 0,
            'wind_speed': wind_speed,
            'v_effective': v_eff,
        }
    elif v_eff > turbine.cut_out_speed:
        return {
            'status': 'above_cut_out',
            'power_kW': 0,
            'wind_speed': wind_speed,
            'v_effective': v_eff,
        }
    
    # Swept area
    A = np.pi * (turbine.rotor_diameter / 2) ** 2
    
    # Available wind power
    P_wind = 0.5 * RHO_AIR * A * v_eff ** 3
    
    # Extracted power (with efficiency)
    P_extract = P_wind * turbine.efficiency
    
    # Cap at rated power
    if P_extract > turbine.rated_power_kW * 1000:
        P_extract = turbine.rated_power_kW * 1000
        status = 'rated'
    else:
        status = 'variable'
    
    return {
        'status': status,
        'power_W': P_extract,
        'power_kW': P_extract / 1000,
        'wind_speed': wind_speed,
        'v_effective': v_eff,
        'P_wind_kW': P_wind / 1000,
        'capacity_factor': P_extract / (turbine.rated_power_kW * 1000),
    }


def calculate_gyroscopic_moment(
    turbine: WindTurbine,
    omega_rad_s: float = None,
    sway_rate_rad_s: float = 0.01
) -> dict:
    """
    Calculate gyroscopic stabilization from spinning turbine rotor.
    
    A spinning rotor acts as a gyroscope, resisting changes to its orientation.
    When the tower sways, the gyroscopic effect creates a restoring torque.
    
    Gyroscopic torque: τ = ω × L = ω × I × Ω
    Where:
        τ = gyroscopic torque (N·m)
        I = rotor moment of inertia (kg·m²)
        Ω = rotor angular velocity (rad/s)
        ω = precession rate (tower sway rate, rad/s)
    
    Parameters:
        turbine: turbine specifications
        omega_rad_s: rotor angular velocity (rad/s), defaults from rpm_rated
        sway_rate_rad_s: tower sway rate (rad/s), typical ~0.01-0.1
    
    Returns:
        Dictionary with gyroscopic analysis
    """
    # Rotor angular velocity
    if omega_rad_s is None:
        omega_rad_s = turbine.rpm_rated * 2 * np.pi / 60
    
    # Approximate moment of inertia for the rotor
    # Model as thin ring at blade tip (simplified)
    # I = m × r² for point masses at radius r
    # For blades, use ~0.5 × m × r² as approximation
    r = turbine.rotor_diameter / 2
    I_rotor = 0.5 * turbine.rotor_mass_kg * r ** 2
    
    # Angular momentum
    L = I_rotor * omega_rad_s
    
    # Gyroscopic torque from tower sway
    # τ = ω_sway × L (perpendicular to both)
    tau_gyro = sway_rate_rad_s * L
    
    # For n turbines on opposite sides, some cancel, some add
    # Depends on rotation direction and placement
    # If all rotate same direction: net effect on sway perpendicular to tower axis
    
    return {
        'turbine': turbine.name,
        'I_rotor_kg_m2': I_rotor,
        'omega_rad_s': omega_rad_s,
        'omega_rpm': omega_rad_s * 60 / (2 * np.pi),
        'L_angular_momentum_kg_m2_s': L,
        'sway_rate_rad_s': sway_rate_rad_s,
        'tau_gyro_Nm': tau_gyro,
        'tau_gyro_kNm': tau_gyro / 1000,
    }


def design_turbine_ring(
    egg_lower: dict,
    egg_upper: dict,
    turbine_type: str = 'small',
    target_coverage: float = 0.7
) -> dict:
    """
    Design a horizontal-disk turbine for the junction between two eggs.
    
    Configuration: Single large horizontal disk turbine at the throat,
    capturing vertical updraft flow accelerated by the Venturi effect.
    
    The disk spins around a vertical axis (like a helicopter rotor),
    providing both power generation and gyroscopic stabilization.
    
    Parameters:
        egg_lower: lower egg properties
        egg_upper: upper egg properties  
        turbine_type: key for TURBINE_MODELS
        target_coverage: fraction of throat area covered by turbine (0.7 = 70%)
    
    Returns:
        Dictionary with turbine design
    """
    turbine = TURBINE_MODELS[turbine_type]
    
    # Junction geometry
    junction = calculate_junction_geometry(egg_lower, egg_upper)
    z = junction['z_junction']
    
    # Use the smaller of the two egg diameters at the junction
    # This is the throat diameter where the turbine sits
    d_throat = min(junction['d_lower'], junction['d_upper'])
    d_outer = max(junction['d_lower'], junction['d_upper'])
    
    # Single central turbine sized to fit the throat
    # Turbine diameter should be less than throat diameter
    effective_turbine_d = min(turbine.rotor_diameter, d_throat * 0.8)
    n_turbines = 1  # Single central disk
    
    # Horizontal wind speed at this height
    v_horizontal = wind_speed_at_height(z)
    
    # Venturi analysis for vertical updraft
    venturi = calculate_venturi_acceleration(
        d_outer=d_outer,
        d_throat=d_throat,
        n_turbines=n_turbines,
        turbine_diameter=effective_turbine_d,
        blockage_fraction=target_coverage
    )
    
    # Vertical wind speed through turbine
    v_vertical = v_horizontal * venturi['v_amplification']
    
    # Power calculation using vertical wind speed
    power_result = calculate_turbine_power(turbine, v_vertical, v_amplification=1.0)
    
    # Total power (single turbine)
    total_power_kW = power_result['power_kW']
    
    # Gyroscopic effect - horizontal disk with vertical axis
    # This provides stabilization against lateral tower sway
    gyro = calculate_gyroscopic_moment(turbine)
    total_gyro_torque_kNm = gyro['tau_gyro_kNm']
    
    # Total mass
    total_mass_kg = turbine.rotor_mass_kg + turbine.nacelle_mass_kg
    
    return {
        'junction_z': z,
        'd_throat': d_throat,
        'd_outer': d_outer,
        'turbine': turbine.name,
        'effective_turbine_d': effective_turbine_d,
        'n_turbines': n_turbines,
        'wind_speed_horizontal_m_s': v_horizontal,
        'v_amplification': venturi['v_amplification'],
        'wind_speed_vertical_m_s': v_vertical,
        'power_per_turbine_kW': power_result['power_kW'],
        'total_power_kW': total_power_kW,
        'total_mass_kg': total_mass_kg,
        'total_mass_tonnes': total_mass_kg / 1000,
        'gyro_torque_kNm': total_gyro_torque_kNm,
        'total_gyro_torque_kNm': total_gyro_torque_kNm,
        'venturi': venturi,
        'power_result': power_result,
        'gyro_result': gyro,
    }


def analyze_tower_turbine_integration(
    eggs: List[dict],
    turbine_type: str = 'small',
    min_diameter: float = 15.0,
    max_junctions: int = None
) -> dict:
    """
    Analyze wind turbine integration for the entire tower.
    
    Places a single horizontal-disk turbine at each egg junction where:
    - Throat diameter is large enough to fit the turbine disk
    - Wind speeds are sufficient for power generation
    
    Parameters:
        eggs: list of analyzed eggs from tower analysis
        turbine_type: key for TURBINE_MODELS
        min_diameter: minimum throat diameter for turbine placement
        max_junctions: limit on number of junctions (None = all valid)
    
    Returns:
        Comprehensive analysis dictionary
    """
    turbine = TURBINE_MODELS[turbine_type]
    
    junction_rings = []
    total_power_kW = 0
    total_mass_kg = 0
    total_gyro_torque_kNm = 0
    
    for i in range(len(eggs) - 1):
        egg_lower = eggs[i]
        egg_upper = eggs[i + 1]
        
        # Throat diameter is the smaller of the two eggs at junction
        d_throat = min(egg_lower['d'], egg_upper['d'])
        
        # Check if throat is large enough for the turbine
        # Turbine needs to fit inside with some clearance
        if d_throat < turbine.rotor_diameter * 1.2:
            continue
            
        # Check minimum diameter threshold
        if d_throat < min_diameter:
            continue
        
        # Design turbine for this junction
        ring = design_turbine_ring(egg_lower, egg_upper, turbine_type)
        ring['junction_index'] = i
        ring['egg_indices'] = (i, i + 1)
        
        junction_rings.append(ring)
        total_power_kW += ring['total_power_kW']
        total_mass_kg += ring['total_mass_kg']
        total_gyro_torque_kNm += ring['total_gyro_torque_kNm']
        
        if max_junctions and len(junction_rings) >= max_junctions:
            break
    
    # Annual energy estimate (assume 30% capacity factor average)
    hours_per_year = 8760
    capacity_factor = 0.30
    annual_energy_MWh = total_power_kW * hours_per_year * capacity_factor / 1000
    
    # Revenue estimate ($50/MWh wholesale)
    revenue_per_MWh = 50
    annual_revenue_usd = annual_energy_MWh * revenue_per_MWh
    
    return {
        'turbine_model': turbine.name,
        'n_junction_rings': len(junction_rings),
        'total_turbines': sum(r['n_turbines'] for r in junction_rings),
        'total_power_kW': total_power_kW,
        'total_power_MW': total_power_kW / 1000,
        'total_mass_kg': total_mass_kg,
        'total_mass_tonnes': total_mass_kg / 1000,
        'total_gyro_torque_kNm': total_gyro_torque_kNm,
        'annual_energy_MWh': annual_energy_MWh,
        'annual_revenue_usd': annual_revenue_usd,
        'junction_rings': junction_rings,
    }


def estimate_sway_damping(
    total_gyro_torque_kNm: float,
    tower_tip_mass_kg: float,
    tower_height_m: float,
    natural_period_s: float = 10.0
) -> dict:
    """
    Estimate the gyroscopic damping effect on tower sway.
    
    The gyroscopic torque opposes the rate of angular change (sway).
    This is similar to adding rotational damping.
    
    Parameters:
        total_gyro_torque_kNm: total gyroscopic torque from all turbines
        tower_tip_mass_kg: effective mass at tower tip
        tower_height_m: tower height
        natural_period_s: natural sway period of tower
    
    Returns:
        Dictionary with damping analysis
    """
    # Natural frequency
    omega_n = 2 * np.pi / natural_period_s
    
    # Maximum sway rate (assuming ~1m tip displacement amplitude)
    sway_amplitude_m = 1.0
    sway_amplitude_rad = sway_amplitude_m / tower_height_m
    max_sway_rate = sway_amplitude_rad * omega_n
    
    # Gyroscopic restoring moment at max sway rate
    # (scale from design sway rate to actual max)
    design_sway_rate = 0.01  # rad/s (used in gyro calculation)
    gyro_moment_kNm = total_gyro_torque_kNm * (max_sway_rate / design_sway_rate)
    
    # Compare to aerodynamic overturning moment
    # Rough estimate: wind load creates moment = F_wind × h/2
    # At 70 m/s wind on a 30m diameter egg: F ~ 0.5 × 1.225 × 70² × 0.5 × π × 15 × 50 ≈ 7 MN
    typical_wind_moment_kNm = 50000  # rough estimate
    
    gyro_effectiveness = gyro_moment_kNm / typical_wind_moment_kNm
    
    # Equivalent damping ratio contribution
    # ζ_gyro ≈ τ_gyro / (2 × m × ω_n × h × v_max)
    # Simplified estimate
    zeta_gyro = gyro_moment_kNm * 1000 / (2 * tower_tip_mass_kg * omega_n * tower_height_m * max_sway_rate * tower_height_m)
    
    return {
        'natural_period_s': natural_period_s,
        'omega_n_rad_s': omega_n,
        'sway_amplitude_m': sway_amplitude_m,
        'max_sway_rate_rad_s': max_sway_rate,
        'gyro_moment_at_max_sway_kNm': gyro_moment_kNm,
        'typical_wind_moment_kNm': typical_wind_moment_kNm,
        'gyro_effectiveness_percent': gyro_effectiveness * 100,
        'estimated_damping_ratio_contribution': zeta_gyro,
    }


def print_turbine_analysis(analysis: dict, detailed: bool = False):
    """Print wind turbine integration analysis."""
    print(f"\n{'='*80}")
    print("HORIZONTAL-DISK WIND TURBINE INTEGRATION")
    print("(Vertical axis, captures updraft through junction throat)")
    print(f"{'='*80}")
    
    print(f"\n📊 SUMMARY")
    print(f"   Turbine Model:        {analysis['turbine_model']}")
    print(f"   Junctions with Turbines: {analysis['n_junction_rings']}")
    print(f"   Total Turbines:       {analysis['total_turbines']} (1 per junction)")
    print(f"   Total Rated Power:    {analysis['total_power_MW']:.2f} MW")
    print(f"   Added Mass:           {analysis['total_mass_tonnes']:.1f} tonnes")
    
    print(f"\n⚡ ENERGY PRODUCTION")
    print(f"   Annual Energy:        {analysis['annual_energy_MWh']:.0f} MWh")
    print(f"   Annual Revenue:       ${analysis['annual_revenue_usd']:,.0f}")
    
    print(f"\n🌀 GYROSCOPIC STABILIZATION (vertical axis)")
    print(f"   Total Gyro Torque:    {analysis['total_gyro_torque_kNm']:.1f} kN·m")
    print(f"   (resists lateral tower sway)")
    
    if detailed and analysis['junction_rings']:
        print(f"\n{'─'*90}")
        print("JUNCTION DETAILS (Horizontal disk capturing vertical updraft)")
        print(f"{'─'*90}")
        print(f"{'Jct':<5} {'Height':<10} {'Throat':<10} {'Disk Ø':<10} {'V_horiz':<10} {'V_vert':<10} {'Power':<10} {'Status':<12}")
        print(f"{'#':<5} {'(m)':<10} {'(m)':<10} {'(m)':<10} {'(m/s)':<10} {'(m/s)':<10} {'(kW)':<10} {'':<12}")
        print(f"{'─'*90}")
        
        for ring in analysis['junction_rings'][:10]:  # Limit output
            status = ring['power_result'].get('status', 'unknown')
            print(f"{ring['junction_index']:<5} "
                  f"{ring['junction_z']:<10.1f} "
                  f"{ring.get('d_throat', ring.get('d_ring', 0)):<10.1f} "
                  f"{ring.get('effective_turbine_d', 0):<10.1f} "
                  f"{ring.get('wind_speed_horizontal_m_s', ring.get('wind_speed_m_s', 0)):<10.1f} "
                  f"{ring.get('wind_speed_vertical_m_s', ring.get('v_effective_m_s', 0)):<10.1f} "
                  f"{ring['total_power_kW']:<10.1f} "
                  f"{status:<12}")
        
        if len(analysis['junction_rings']) > 10:
            print(f"   ... and {len(analysis['junction_rings']) - 10} more junctions")


def evaluate_integration_scenarios(eggs: List[dict]) -> dict:
    """
    Evaluate multiple turbine integration scenarios.
    
    Compare different turbine sizes and coverage options.
    """
    scenarios = []
    
    for turbine_type in ['micro', 'small', 'medium']:
        for coverage in [0.2, 0.3, 0.5]:
            # Modify design_turbine_ring call by using different targets
            analysis = analyze_tower_turbine_integration(
                eggs, 
                turbine_type=turbine_type,
                min_diameter=TURBINE_MODELS[turbine_type].rotor_diameter * 2
            )
            
            scenarios.append({
                'turbine_type': turbine_type,
                'coverage': coverage,
                'power_MW': analysis['total_power_MW'],
                'mass_tonnes': analysis['total_mass_tonnes'],
                'annual_MWh': analysis['annual_energy_MWh'],
                'revenue_usd': analysis['annual_revenue_usd'],
                'power_per_mass': analysis['total_power_kW'] / max(1, analysis['total_mass_kg']) * 1000,  # W/kg
            })
    
    return {
        'scenarios': scenarios,
        'best_power': max(scenarios, key=lambda x: x['power_MW']),
        'best_efficiency': max(scenarios, key=lambda x: x['power_per_mass']),
    }


# ============================================================================
# FEASIBILITY ASSESSMENT
# ============================================================================

def assess_feasibility(eggs: List[dict]) -> dict:
    """
    Comprehensive feasibility assessment for wind turbine integration.
    
    Evaluates INTERNAL DUCT system (preferred) vs Junction turbines:
    1. Geometric compatibility (do turbines fit?)
    2. Structural impact (added mass and loads)
    3. Energy production potential
    4. Gyroscopic stabilization benefit
    5. Economic viability
    """
    # Run INTERNAL DUCT analysis (preferred system)
    # Use min_egg_diameter=10 to include all eggs that can fit micro turbines
    duct_analysis = analyze_internal_duct_integration(eggs, inlet_fraction=0.25, min_egg_diameter=10)
    
    # Also run junction analysis for comparison
    junction_analysis = analyze_tower_turbine_integration(eggs, turbine_type='small')
    
    # Use internal duct as primary (better power/mass ratio)
    analysis = duct_analysis
    
    # Calculate structural impact (handle both 'mass' and 'shell_mass_kg' keys)
    def get_egg_mass(e):
        return e.get('mass', e.get('shell_mass_kg', 0))
    
    tower_mass = sum(get_egg_mass(e) for e in eggs)
    mass_increase_percent = (analysis['total_mass_kg'] / tower_mass) * 100
    
    # Gyroscopic effectiveness
    damping = estimate_sway_damping(
        analysis['total_gyro_torque_kNm'],
        tower_tip_mass_kg=get_egg_mass(eggs[-1]),
        tower_height_m=eggs[-1]['z_top']
    )
    
    # Economic metrics for internal duct
    # CAPEX scales with turbine size: ~$2000/kW for small turbines
    # Plus shell modification costs (~$50k per egg for inlet/outlet openings)
    capex_per_kW = 2000  # $/kW installed
    shell_mod_per_egg = 50_000  # Shell cutting and reinforcement
    total_capex = (
        analysis['total_power_kW'] * capex_per_kW + 
        analysis['n_installations'] * shell_mod_per_egg
    )
    simple_payback_years = total_capex / max(1, analysis['annual_revenue_usd'])
    
    return {
        'system_type': 'Internal Duct (Egg as Wind Concentrator)',
        'geometric_feasibility': {
            'suitable_eggs': analysis['n_installations'],
            'total_turbines': analysis['n_installations'],
            'assessment': 'FEASIBLE' if analysis['n_installations'] >= 3 else 'LIMITED',
        },
        'structural_impact': {
            'added_mass_tonnes': analysis['total_mass_tonnes'],
            'mass_increase_percent': mass_increase_percent,
            'assessment': 'MINIMAL' if mass_increase_percent < 0.1 else 'ACCEPTABLE',
        },
        'energy_production': {
            'rated_power_MW': analysis['total_power_MW'],
            'annual_energy_MWh': analysis['annual_energy_MWh'],
            'assessment': 'GOOD' if analysis['annual_energy_MWh'] > 1000 else 'MODEST',
        },
        'gyroscopic_stabilization': {
            'gyro_effectiveness_percent': damping['gyro_effectiveness_percent'],
            'damping_contribution': damping['estimated_damping_ratio_contribution'],
            'assessment': 'SUPPLEMENTARY',  # Small turbines = modest gyro effect
        },
        'economics': {
            'estimated_capex_usd': total_capex,
            'annual_revenue_usd': analysis['annual_revenue_usd'],
            'simple_payback_years': simple_payback_years,
            'assessment': 'ATTRACTIVE' if simple_payback_years < 15 else 'REASONABLE' if simple_payback_years < 25 else 'MARGINAL',
        },
        'overall_recommendation': _generate_duct_recommendation(analysis, damping, mass_increase_percent, simple_payback_years),
        'detailed_analysis': analysis,
        'junction_comparison': junction_analysis,
        'damping_analysis': damping,
    }


def _generate_duct_recommendation(analysis: dict, damping: dict, mass_increase: float, payback: float) -> str:
    """Generate overall recommendation for internal duct system."""
    
    pros = []
    cons = []
    
    if analysis['total_power_MW'] > 0.5:
        pros.append(f"Meaningful power generation ({analysis['total_power_MW']:.2f} MW)")
    
    if analysis['annual_energy_MWh'] > 1000:
        pros.append(f"Significant annual energy ({analysis['annual_energy_MWh']:.0f} MWh)")
    
    if mass_increase < 0.1:
        pros.append(f"Negligible mass impact ({mass_increase:.3f}%)")
    
    if payback < 15:
        pros.append(f"Attractive payback period ({payback:.1f} years)")
    elif payback < 25:
        pros.append(f"Reasonable payback period ({payback:.1f} years)")
    else:
        cons.append(f"Long payback period ({payback:.1f} years)")
    
    pros.append("Uses existing hollow structure - minimal modification")
    pros.append("High power-to-mass ratio (duct acceleration)")
    
    if len(pros) >= 4:
        return "RECOMMENDED: Internal duct turbines are highly feasible with excellent power-to-mass ratio."
    elif len(pros) >= 2:
        return "RECOMMENDED WITH STUDY: Good potential - detailed CFD analysis of duct flow recommended."
    else:
        return "CONDITIONAL: May require design optimization for economic viability."


def _generate_recommendation(analysis: dict, damping: dict, mass_increase: float) -> str:
    """Generate overall recommendation text."""
    
    pros = []
    cons = []
    
    if analysis['total_power_MW'] > 1:
        pros.append(f"Substantial power generation ({analysis['total_power_MW']:.1f} MW)")
    
    if analysis['annual_energy_MWh'] > 5000:
        pros.append(f"Meaningful annual energy production ({analysis['annual_energy_MWh']:.0f} MWh)")
    
    if damping['gyro_effectiveness_percent'] > 1:
        pros.append("Gyroscopic effect provides supplementary stabilization")
    
    if mass_increase > 3:
        cons.append(f"Adds {mass_increase:.1f}% to tower mass")
    
    if len(pros) >= 2 and len(cons) <= 1:
        return "RECOMMENDED: Wind turbine integration is feasible and provides meaningful benefits."
    elif len(pros) >= 1:
        return "CONDITIONAL: Integration provides some benefits but requires detailed engineering study."
    else:
        return "NOT RECOMMENDED: Benefits do not justify complexity and added mass."


def print_feasibility_report(feasibility: dict):
    """Print comprehensive feasibility report."""
    print(f"\n{'='*90}")
    print("WIND TURBINE INTEGRATION FEASIBILITY ASSESSMENT")
    print(f"System: {feasibility.get('system_type', 'Unknown')}")
    print(f"{'='*90}")
    
    print(f"\n📐 GEOMETRIC FEASIBILITY")
    geo = feasibility['geometric_feasibility']
    n_items = geo.get('suitable_eggs', geo.get('suitable_junctions', 0))
    print(f"   Suitable Locations:   {n_items}")
    print(f"   Total Turbines:       {geo['total_turbines']}")
    print(f"   Assessment:           {geo['assessment']}")
    
    print(f"\n🏗️ STRUCTURAL IMPACT")
    struct = feasibility['structural_impact']
    print(f"   Added Mass:           {struct['added_mass_tonnes']:.1f} tonnes")
    print(f"   Mass Increase:        {struct['mass_increase_percent']:.3f}%")
    print(f"   Assessment:           {struct['assessment']}")
    
    print(f"\n⚡ ENERGY PRODUCTION")
    energy = feasibility['energy_production']
    print(f"   Rated Power:          {energy['rated_power_MW']:.2f} MW")
    print(f"   Annual Energy:        {energy['annual_energy_MWh']:.0f} MWh")
    print(f"   Assessment:           {energy['assessment']}")
    
    print(f"\n🌀 GYROSCOPIC STABILIZATION")
    gyro = feasibility['gyroscopic_stabilization']
    print(f"   Effectiveness:        {gyro['gyro_effectiveness_percent']:.2f}%")
    print(f"   Damping Contribution: {gyro['damping_contribution']:.4f}")
    print(f"   Assessment:           {gyro['assessment']}")
    
    print(f"\n💰 ECONOMICS")
    econ = feasibility['economics']
    print(f"   Estimated CAPEX:      ${econ['estimated_capex_usd']:,.0f}")
    print(f"   Annual Revenue:       ${econ['annual_revenue_usd']:,.0f}")
    print(f"   Simple Payback:       {econ['simple_payback_years']:.1f} years")
    print(f"   Assessment:           {econ['assessment']}")
    
    print(f"\n{'─'*90}")
    print(f"📋 OVERALL RECOMMENDATION")
    print(f"   {feasibility['overall_recommendation']}")
    print(f"{'─'*90}")


if __name__ == "__main__":
    # Test with a sample tower
    from geometry import build_tower
    from structural import analyze_tower
    
    print("Building test tower...")
    eggs, ratio = build_tower(d_base=100, d_top=5, target_height=1600)
    analyzed_eggs = analyze_tower(eggs)
    
    print(f"Tower has {len(analyzed_eggs)} eggs")
    
    # Run feasibility assessment
    feasibility = assess_feasibility(analyzed_eggs)
    print_feasibility_report(feasibility)
    
    # Detailed turbine analysis
    analysis = analyze_tower_turbine_integration(analyzed_eggs, turbine_type='small')
    print_turbine_analysis(analysis, detailed=True)
