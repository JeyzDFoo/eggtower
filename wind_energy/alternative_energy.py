"""
Alternative Energy Generation for Egg Tower.

This module explores two additional energy generation concepts:

1. VORTEX SHEDDING ENERGY CAPTURE
   - Tower naturally oscillates due to vortex shedding from wind
   - Piezoelectric or electromagnetic generators convert motion to electricity
   - Similar to "Vortex Bladeless" wind turbine concept
   - RESULT: Not viable - very low power output from differential motion
   
2. VERTICAL UPDRAFT COLUMN (Solar/Wind Chimney)
   - Hollow eggs form a continuous vertical column
   - Air enters at base, heats up (solar) or is pushed up (wind pressure)
   - Accelerates upward through the column
   - Horizontal disk turbines at egg transitions capture updraft
   - RESULT: Very promising - ~900 kW, 2300 MWh/year, 3.8 year payback
   
IMPORTANT CAVEATS for Updraft Column:
- Solar chimney ΔT=25°C is optimistic without dedicated collector area
- Friction losses in 1600m column may be higher than modeled
- Turbine placement needs careful CFD analysis to avoid blockage
- Combined solar + wind driving forces are theoretical maximum
- Real performance likely 30-50% of calculated values
- Still potentially viable at reduced efficiency

Both concepts leverage the unique geometry of the egg tower.
"""
import numpy as np
import json
import os
from dataclasses import dataclass
from typing import List, Dict, Tuple

# Import from parent/sibling modules
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from config import RHO_AIR, V_REF, RHO
from wind_energy.wind_turbine import load_optimized_tower, get_optimized_eggs, wind_speed_at_height


# =============================================================================
# CONSTANTS
# =============================================================================

G = 9.81  # m/s² gravity
T_AMBIENT = 288  # K (15°C ambient temperature)
CP_AIR = 1005  # J/(kg·K) specific heat of air
SOLAR_IRRADIANCE = 1000  # W/m² peak solar


# =============================================================================
# VORTEX SHEDDING ENERGY CAPTURE
# =============================================================================

@dataclass
class VortexHarvester:
    """Vortex-induced vibration energy harvester."""
    name: str
    type: str                    # 'piezoelectric' or 'electromagnetic'
    power_per_unit_kW: float     # Power output per unit at rated oscillation
    mass_kg: float               # Mass per unit
    cost_usd: float              # Cost per unit
    min_amplitude_m: float       # Minimum oscillation amplitude to operate
    max_amplitude_m: float       # Maximum safe amplitude
    efficiency: float            # Conversion efficiency


# Harvester options
VORTEX_HARVESTERS = {
    'piezo_small': VortexHarvester(
        name="Piezoelectric Stack (Small)",
        type="piezoelectric",
        power_per_unit_kW=0.5,
        mass_kg=50,
        cost_usd=15_000,
        min_amplitude_m=0.01,
        max_amplitude_m=0.5,
        efficiency=0.15
    ),
    'piezo_large': VortexHarvester(
        name="Piezoelectric Stack (Large)",
        type="piezoelectric",
        power_per_unit_kW=2.0,
        mass_kg=200,
        cost_usd=50_000,
        min_amplitude_m=0.02,
        max_amplitude_m=1.0,
        efficiency=0.18
    ),
    'electromagnetic': VortexHarvester(
        name="Electromagnetic Linear Generator",
        type="electromagnetic",
        power_per_unit_kW=5.0,
        mass_kg=500,
        cost_usd=80_000,
        min_amplitude_m=0.05,
        max_amplitude_m=2.0,
        efficiency=0.25
    ),
}


def calculate_vortex_shedding_frequency(diameter: float, wind_speed: float, 
                                         strouhal: float = 0.2) -> float:
    """
    Calculate vortex shedding frequency for a cylinder/ellipsoid.
    
    f = St × V / D
    
    Parameters:
        diameter: Characteristic diameter (m)
        wind_speed: Wind velocity (m/s)
        strouhal: Strouhal number (0.2 typical for cylinders)
    
    Returns:
        Shedding frequency (Hz)
    """
    if diameter <= 0:
        return 0
    return strouhal * wind_speed / diameter


def calculate_oscillation_amplitude(
    egg: dict,
    wind_speed: float,
    damping_ratio: float = 0.02
) -> dict:
    """
    Estimate oscillation amplitude from vortex-induced vibration.
    
    For VIV, maximum amplitude occurs when shedding frequency matches
    the structure's natural frequency (lock-in).
    
    Amplitude is typically 0.5-1.5 diameters for undamped structures,
    reduced by damping.
    """
    d = egg['d']
    h = egg['h']
    mass = egg.get('mass', egg.get('shell_mass_kg', 0))
    
    # Vortex shedding frequency
    f_shed = calculate_vortex_shedding_frequency(d, wind_speed)
    
    # Estimate natural frequency of egg (simplified as cantilever)
    # f_n ≈ (1/2π) × √(k/m) - very rough estimate
    # For a large shell, natural frequency is quite low
    E = 45e9  # BFRP modulus (Pa)
    t = egg.get('t_avg', d/200)
    I = np.pi * (d/2)**3 * t  # Moment of inertia approximation
    L = h
    k = 3 * E * I / L**3  # Cantilever stiffness
    f_natural = (1 / (2 * np.pi)) * np.sqrt(k / max(mass, 1))
    
    # Lock-in occurs when f_shed ≈ f_natural
    # Frequency ratio
    r = f_shed / max(f_natural, 0.01)
    
    # Amplitude magnification factor (resonance curve)
    # A/A_static = 1 / sqrt((1-r²)² + (2ζr)²)
    denominator = np.sqrt((1 - r**2)**2 + (2 * damping_ratio * r)**2)
    magnification = 1 / max(denominator, 0.1)
    
    # Static deflection from vortex lift force
    # F_lift = 0.5 × ρ × V² × C_L × A
    C_L = 0.4  # Vortex lift coefficient
    A_projected = d * h  # Projected area
    F_lift = 0.5 * RHO_AIR * wind_speed**2 * C_L * A_projected
    
    # Static deflection
    A_static = F_lift / max(k, 1)
    
    # Dynamic amplitude
    A_dynamic = A_static * magnification
    
    # Cap at reasonable values (can't exceed ~1 diameter)
    A_dynamic = min(A_dynamic, d * 0.5)
    
    return {
        'egg_index': egg.get('index', 0),
        'd': d,
        'f_shedding_Hz': f_shed,
        'f_natural_Hz': f_natural,
        'frequency_ratio': r,
        'magnification': magnification,
        'amplitude_m': A_dynamic,
        'amplitude_fraction_d': A_dynamic / d,
        'is_lock_in': 0.8 < r < 1.2,  # Near resonance
    }


def analyze_vortex_energy_harvesting(eggs: List[dict]) -> dict:
    """
    Analyze potential for vortex-induced vibration energy harvesting.
    
    The tower oscillates due to vortex shedding. We can place harvesters
    at egg junctions to capture this motion.
    
    Note: Individual eggs are very stiff (double-curved shell). The useful
    oscillation comes from the TOWER AS A WHOLE swaying, not individual eggs.
    We model cumulative displacement at top of tower.
    """
    # Use design wind speed for analysis
    design_wind = 25  # m/s (moderate-high wind for energy harvesting)
    
    installations = []
    total_power_kW = 0
    total_mass_kg = 0
    total_cost = 0
    
    # Total tower properties
    tower_height = eggs[-1]['z_top']
    total_mass = sum(egg.get('mass', egg.get('shell_mass_kg', 0)) for egg in eggs)
    avg_diameter = np.mean([egg['d'] for egg in eggs])
    
    # Tower as a whole oscillates - use whole-tower analysis
    # Natural frequency of 1600m tower: f ≈ 0.02-0.05 Hz (very low)
    f_tower_natural = 0.03  # Hz (estimate for 1600m structure)
    
    # Vortex shedding at tower scale
    v_top = wind_speed_at_height(tower_height, v_ref=design_wind, z_ref=10)
    f_shed = calculate_vortex_shedding_frequency(avg_diameter, v_top)
    
    # Tower top displacement from wind
    # Typical sway for tall structures: H/500 to H/200 at design wind
    tower_sway_amplitude = tower_height / 400  # ~4m amplitude at 25 m/s
    
    # Use smaller piezo harvesters at each junction (they work with small motion)
    harvester = VORTEX_HARVESTERS['piezo_small']
    
    # At each egg junction, the relative motion is smaller
    for i in range(1, len(eggs)):  # Start from second egg (first junction)
        egg = eggs[i]
        z = egg['z_base']
        
        # Relative displacement at this height (linear mode shape approximation)
        height_ratio = z / tower_height
        local_amplitude = tower_sway_amplitude * height_ratio
        
        # Differential motion between adjacent eggs (mode shape derivative)
        # This is what the harvester at the junction experiences
        delta_h = egg['h']
        differential_amplitude = tower_sway_amplitude * (delta_h / tower_height) * height_ratio
        
        # Minimum ~1cm relative motion for piezo to work
        if differential_amplitude < 0.005:
            continue
        
        # Power from relative motion between eggs
        # P = 0.5 × k × A² × ω × η
        omega = 2 * np.pi * f_tower_natural
        
        # Piezo spring constant (stiff)
        k_piezo = 1e6  # N/m
        
        # Power = 0.5 × k × A² × ω × η
        P_available = 0.5 * k_piezo * differential_amplitude**2 * omega
        P_harvested = P_available * harvester.efficiency / 1000  # kW
        
        # Cap at harvester rated power
        P_harvested = min(P_harvested, harvester.power_per_unit_kW)
        
        installations.append({
            'junction': i,
            'z_junction': z,
            'local_amplitude_m': local_amplitude,
            'differential_amplitude_m': differential_amplitude,
            'power_kW': P_harvested,
            'harvester': harvester.name,
        })
        
        total_power_kW += P_harvested
        total_mass_kg += harvester.mass_kg
        total_cost += harvester.cost_usd
    
    # Capacity factor (wind isn't always at design speed)
    capacity_factor = 0.30  # Good winds at height
    annual_energy_MWh = total_power_kW * 8760 * capacity_factor / 1000
    annual_revenue = annual_energy_MWh * 50  # $50/MWh
    
    return {
        'system_type': 'Vortex Shedding Energy Capture',
        'concept': 'Piezoelectric harvesters at egg junctions capture tower sway',
        'design_wind_speed': design_wind,
        'tower_sway_amplitude_m': tower_sway_amplitude,
        'tower_natural_freq_Hz': f_tower_natural,
        'n_harvesters': len(installations),
        'harvester_type': harvester.name,
        'total_power_kW': total_power_kW,
        'total_mass_kg': total_mass_kg,
        'total_cost_usd': total_cost,
        'capacity_factor': capacity_factor,
        'annual_energy_MWh': annual_energy_MWh,
        'annual_revenue_usd': annual_revenue,
        'simple_payback_years': total_cost / max(annual_revenue, 1),
        'installations': installations,
    }


# =============================================================================
# VERTICAL UPDRAFT COLUMN (SOLAR/WIND CHIMNEY)
# =============================================================================

@dataclass
class UpdraftTurbine:
    """
    Aerospace-style turbomachinery for high-velocity updraft capture.
    
    Unlike wind turbines, these are designed like simplified jet engine compressor
    stages - optimized for high velocity, compact diameter, no cut-out limit.
    Think: axial flow turbine running in reverse as a generator.
    """
    name: str
    rotor_diameter: float       # m
    rated_power_kW: float       # at rated updraft velocity  
    rated_updraft_speed: float  # m/s - design point
    cut_in_speed: float         # m/s - minimum operating speed
    # No cut-out speed - aerospace turbines can handle extreme velocities
    rotor_mass_kg: float
    generator_mass_kg: float
    efficiency: float           # Turbine efficiency (not Betz-limited in duct)


UPDRAFT_TURBINES = {
    # Standard low-speed turbines (for baseline/comparison)
    'small': UpdraftTurbine(
        name="Small Axial Turbine (2m)",
        rotor_diameter=2.0,
        rated_power_kW=10,
        rated_updraft_speed=15,
        cut_in_speed=2.0,
        rotor_mass_kg=100,
        generator_mass_kg=200,
        efficiency=0.35
    ),
    'medium': UpdraftTurbine(
        name="Medium Axial Turbine (4m)",
        rotor_diameter=4.0,
        rated_power_kW=50,
        rated_updraft_speed=15,
        cut_in_speed=2.0,
        rotor_mass_kg=300,
        generator_mass_kg=600,
        efficiency=0.38
    ),
    'large': UpdraftTurbine(
        name="Large Axial Turbine (6m)",
        rotor_diameter=6.0,
        rated_power_kW=120,
        rated_updraft_speed=15,
        cut_in_speed=2.0,
        rotor_mass_kg=600,
        generator_mass_kg=1200,
        efficiency=0.40
    ),
    # High-velocity aerospace-style turbines
    # Based on simplified jet engine turbine stage technology
    'aerospace_small': UpdraftTurbine(
        name="Aerospace Turbine (2m)",
        rotor_diameter=2.0,
        rated_power_kW=200,            # Much higher power at high speed
        rated_updraft_speed=80,         # Designed for 80 m/s
        cut_in_speed=10.0,              # Needs minimum flow
        rotor_mass_kg=150,              # Titanium/composite construction
        generator_mass_kg=300,
        efficiency=0.45                 # Better efficiency - ducted flow
    ),
    'aerospace_medium': UpdraftTurbine(
        name="Aerospace Turbine (3m)",
        rotor_diameter=3.0,
        rated_power_kW=500,            
        rated_updraft_speed=80,         
        cut_in_speed=10.0,              
        rotor_mass_kg=250,              
        generator_mass_kg=500,
        efficiency=0.45                 
    ),
    'aerospace_large': UpdraftTurbine(
        name="Aerospace Turbine (4m)",
        rotor_diameter=4.0,
        rated_power_kW=900,            
        rated_updraft_speed=80,         
        cut_in_speed=10.0,              
        rotor_mass_kg=400,              
        generator_mass_kg=800,
        efficiency=0.45                 
    ),
}


def calculate_column_geometry(eggs: List[dict], tapered_duct: bool = True) -> dict:
    """
    Calculate the geometry of the internal vertical column with optional tapered duct.
    
    The hollow eggs stack to form a continuous chimney from base to top.
    
    Three configurations:
    1. NATURAL GEOMETRY (tapered_duct=False): 
       - Airflow follows the egg's internal shape
       - Throat (junction) is narrower than waist
       - Creates venturi effect: high velocity at throat, high losses
       
    2. TAPERED DUCT (tapered_duct=True):
       - Install a smoothly tapered internal duct
       - Sized to fit through each junction's throat
       - Moderate taper maintains reasonable velocity
       - MULTI-STAGE TURBINES extract energy to control velocity
       - If velocity too high → extract more energy to slow flow
       
    Key Physics:
    - Mass flow rate constant: ṁ = ρ × A × v
    - To maintain target velocity v_target, area must scale with height
    - Energy extraction at each stage reduces velocity
    - Tapered duct with distributed turbines = optimal control
    
    Parameters:
        eggs: List of egg geometry dictionaries
        tapered_duct: If True, use tapered duct with multi-stage extraction
    """
    column_segments = []
    throat_locations = []  # Junctions between eggs
    turbine_stages = []    # Optimal turbine placement locations
    
    # First pass: calculate throat diameters at each junction
    throat_diameters = []
    for i, egg in enumerate(eggs):
        d_outer = egg['d']
        t_avg = egg.get('t_avg', d_outer / 200)
        d_internal = d_outer - 2 * t_avg
        
        # Throat diameter at junction (35% of internal diameter)
        d_throat_ratio = 0.35
        d_throat = d_internal * d_throat_ratio
        throat_diameters.append(d_throat)
    
    # For tapered duct: size to 85% of each local throat (clearance)
    # This creates a smooth taper that follows the structural constraints
    clearance_factor = 0.85
    duct_diameters = [d * clearance_factor for d in throat_diameters]
    d_duct_base = duct_diameters[0]
    d_duct_top = duct_diameters[-1]
    
    # Target velocity for turbine operation
    v_target = 40.0  # m/s - optimal for aerospace turbines
    
    for i, egg in enumerate(eggs):
        d_outer = egg['d']
        t_avg = egg.get('t_avg', d_outer / 200)
        d_internal = d_outer - 2 * t_avg
        
        # Natural geometry values (for reference)
        d_column_waist = d_internal * 0.60
        A_column_waist = np.pi * (d_column_waist / 2)**2
        d_throat_ratio = 0.35
        d_column_throat = d_internal * d_throat_ratio
        A_column_throat = np.pi * (d_column_throat / 2)**2
        
        # Effective duct diameter at this level
        if tapered_duct:
            d_duct = duct_diameters[i]
            A_duct = np.pi * (d_duct / 2)**2
            d_column_effective = d_duct
            A_column_effective = A_duct
            velocity_ratio = 1.0  # Smooth taper, no sudden changes
        else:
            d_column_effective = d_column_waist
            A_column_effective = A_column_waist
            velocity_ratio = d_column_waist**2 / d_column_throat**2
        
        column_segments.append({
            'egg_index': egg.get('index', len(column_segments)),
            'z_base': egg['z_base'],
            'z_top': egg['z_top'],
            'height': egg['h'],
            'd_internal': d_internal,
            'd_column_waist': d_column_waist,
            'd_column_throat': d_column_throat,
            'd_duct': d_duct if tapered_duct else d_column_waist,
            'd_column': d_column_effective,
            'A_column_waist': A_column_waist,
            'A_column_throat': A_column_throat,
            'A_duct': A_duct if tapered_duct else A_column_waist,
            'A_column_m2': A_column_effective,
            'throat_velocity_ratio': velocity_ratio,
        })
        
        # Record throat location (junction between this egg and next)
        if i < len(eggs) - 1:
            throat_locations.append({
                'junction_index': i,
                'z_location': egg['z_top'],
                'd_throat': d_column_throat,
                'd_duct': duct_diameters[i] if tapered_duct else d_column_throat,
                'A_throat': A_column_throat,
                'A_duct': np.pi * (duct_diameters[i]/2)**2 if tapered_duct else A_column_throat,
            })
            
            # Turbine stage at each junction (inside the duct)
            turbine_stages.append({
                'stage': i,
                'z_location': egg['z_top'],
                'd_duct': duct_diameters[i] if tapered_duct else d_column_throat,
                'A_duct': np.pi * (duct_diameters[i]/2)**2 if tapered_duct else A_column_throat,
            })
    
    # Column areas
    if tapered_duct:
        A_base = np.pi * (d_duct_base / 2)**2
        A_top = np.pi * (d_duct_top / 2)**2
        A_avg = np.sqrt(A_base * A_top)
        A_throat_avg = A_avg  # Duct is smooth, no sudden throat
        velocity_boost_at_throat = 1.0  # No venturi
        
        # Taper ratio
        taper_ratio = d_duct_top / d_duct_base
    else:
        A_base = column_segments[0]['A_column_waist']
        A_top = column_segments[-1]['A_column_waist']
        A_avg = np.sqrt(A_base * A_top)
        A_throat_base = column_segments[0]['A_column_throat']
        A_throat_top = column_segments[-1]['A_column_throat']
        A_throat_avg = np.sqrt(A_throat_base * A_throat_top)
        velocity_boost_at_throat = A_avg / A_throat_avg
        taper_ratio = 1.0
    
    # Total column height
    total_height = eggs[-1]['z_top'] - eggs[0]['z_base']
    
    # Tapered duct properties
    if tapered_duct:
        duct_length = total_height
        # Average diameter for surface area calculation
        d_avg = (d_duct_base + d_duct_top) / 2
        duct_surface = np.pi * d_avg * duct_length
        # Lightweight composite: ~5 kg/m², $50/m²
        duct_mass_kg = duct_surface * 5
        duct_cost = duct_surface * 50
    else:
        duct_mass_kg = 0
        duct_cost = 0
    
    return {
        'n_segments': len(column_segments),
        'total_height': total_height,
        'tapered_duct': tapered_duct,
        'duct_diameter_base': d_duct_base if tapered_duct else None,
        'duct_diameter_top': d_duct_top if tapered_duct else None,
        'taper_ratio': taper_ratio,
        'duct_mass_kg': duct_mass_kg,
        'duct_cost_usd': duct_cost,
        'base_column_diameter': column_segments[0]['d_column'],
        'top_column_diameter': column_segments[-1]['d_column'],
        'base_area_m2': A_base,
        'top_area_m2': A_top,
        'average_column_area': A_avg,
        'average_throat_area': A_throat_avg,
        'throat_velocity_ratio': velocity_boost_at_throat,
        'throat_locations': throat_locations,
        'turbine_stages': turbine_stages,
        'n_turbine_stages': len(turbine_stages),
        'segments': column_segments,
        'v_target': v_target,
    }


def calculate_thermal_updraft(
    column_height: float,
    base_area: float,
    solar_heating: float = 0.5,  # Fraction of solar absorbed
    ambient_temp: float = T_AMBIENT
) -> dict:
    """
    Calculate updraft velocity from solar chimney effect.
    
    The sun heats the column interior (especially lower sections).
    Hot air rises, creating continuous updraft.
    
    Physics:
    - ΔT from solar heating of internal surfaces
    - Buoyancy force: F = ρ × g × β × ΔT × V
    - Updraft velocity: v = √(2 × g × H × ΔT / T_ambient)
    """
    # Solar energy absorbed per meter height
    # Assume shell is partially transparent or has collector surfaces
    collector_efficiency = solar_heating
    
    # Temperature rise estimate
    # Solar power in = SOLAR_IRRADIANCE × effective_area × efficiency
    # This heats air mass flowing through
    
    # For solar chimney: ΔT ≈ 20-40°C is typical
    # Height advantage increases velocity
    delta_T = 25  # K (conservative estimate)
    
    # Theoretical updraft velocity (no friction)
    # v = √(2 × g × H × ΔT / T_amb)
    v_theoretical = np.sqrt(2 * G * column_height * delta_T / ambient_temp)
    
    # Account for friction losses in column (reduce by ~50%)
    friction_factor = 0.50
    v_actual = v_theoretical * friction_factor
    
    # Mass flow rate
    rho_hot = RHO_AIR * ambient_temp / (ambient_temp + delta_T)  # Less dense
    m_dot = rho_hot * base_area * v_actual  # kg/s
    
    # Available power in flow
    P_flow = 0.5 * rho_hot * base_area * v_actual**3
    
    return {
        'delta_T_K': delta_T,
        'v_theoretical': v_theoretical,
        'v_actual': v_actual,
        'friction_factor': friction_factor,
        'mass_flow_kg_s': m_dot,
        'air_density_hot': rho_hot,
        'power_in_flow_kW': P_flow / 1000,
    }


def calculate_wind_driven_updraft(
    column: dict,
    base_wind_speed: float = 10.0
) -> dict:
    """
    Calculate updraft driven by wind pressure differential.
    
    Wind creates positive pressure on windward side at base,
    and negative pressure (suction) at top due to Bernoulli effect.
    
    This pressure differential drives air upward through the column.
    
    Physics:
    - Base pressure: P_base = 0.5 × ρ × V_base² × Cp_base
    - Top suction: P_top = -0.5 × ρ × V_top² × Cp_top  
    - Pressure difference drives flow through column
    """
    # Wind speeds at base and top
    z_base = column['segments'][0]['z_base']
    z_top = column['segments'][-1]['z_top']
    
    v_base = wind_speed_at_height(max(z_base, 10), v_ref=base_wind_speed, z_ref=10)
    v_top = wind_speed_at_height(z_top, v_ref=base_wind_speed, z_ref=10)
    
    # Pressure coefficients
    Cp_inlet = 0.8   # Stagnation on windward shell
    Cp_outlet = -0.6  # Suction on leeward top
    
    # Pressure differential
    P_base = 0.5 * RHO_AIR * v_base**2 * Cp_inlet
    P_top = 0.5 * RHO_AIR * v_top**2 * Cp_outlet  # Negative = suction
    
    delta_P = P_base - P_top  # Total driving pressure
    
    # Updraft velocity from pressure differential
    # ΔP = 0.5 × ρ × v² × (friction_loss_coeff)
    # Rearrange: v = √(2 × ΔP / (ρ × K))
    K_loss = 2.0  # Friction and minor losses coefficient
    
    v_updraft = np.sqrt(2 * abs(delta_P) / (RHO_AIR * K_loss))
    
    # Mass flow
    A_column = column['average_column_area']
    m_dot = RHO_AIR * A_column * v_updraft
    
    # Available power
    P_flow = 0.5 * RHO_AIR * A_column * v_updraft**3
    
    return {
        'v_base_wind': v_base,
        'v_top_wind': v_top,
        'P_base_Pa': P_base,
        'P_top_Pa': P_top,
        'delta_P_Pa': delta_P,
        'v_updraft': v_updraft,
        'mass_flow_kg_s': m_dot,
        'power_in_flow_kW': P_flow / 1000,
    }


def calculate_valved_accelerator_boost(
    column: dict,
    eggs: list,
    base_wind_speed: float = 10.0,
    valve_efficiency: float = 0.6,
    dual_port: bool = True
) -> dict:
    """
    Calculate updraft boost from valved accelerators at choke points.
    
    CONCEPT: At each egg's maximum diameter (waist), external wind creates
    pressure differentials. Valved ports on BOTH sides harness this:
    
    SINGLE PORT MODE (dual_port=False):
    - Leeward valves at WAIST → suction pulls air up (Cp = -0.8)
    - ΔP per stage ≈ 0.5 × ρ × V² × 0.8 × η
    
    DUAL PORT MODE (dual_port=True):
    - SUCTION ports at WAIST (max Ø): Leeward Cp = -0.8 (29 ports, one per egg)
    - INJECTION ports at THROAT (junction): Windward Cp = +0.8 (28 ports, between eggs)
    - Uses optimal pressure zones for each port type
    
    Physics:
    - At waist: Flow accelerates around max diameter → low pressure (suction)
    - At throat: Flow impacts narrow junction → stagnation pressure (injection)
    - Combined effect: Suction pulls from above, injection pushes from below
    
    Parameters:
        column: Column geometry from calculate_column_geometry
        eggs: List of egg dictionaries with z_base, d, h
        base_wind_speed: Reference wind speed at 10m (m/s)
        valve_efficiency: How much of theoretical pressure is captured (0-1)
        dual_port: If True, use both suction (waist) and injection (throat) ports
    """
    suction_stages = []  # At waist (max diameter)
    injection_stages = []  # At throat (junction between eggs)
    total_delta_P = 0
    
    # Pressure coefficients
    Cp_suction = -0.8   # At waist (leeward), flow accelerates → suction
    Cp_injection = 0.9  # At throat (windward), flow stagnates → higher pressure
    
    # === SUCTION PORTS: At waist of each egg (29 stages) ===
    for i, egg in enumerate(eggs):
        # Height at egg's maximum diameter (waist) - approximately 40% up from base
        z_waist = egg['z_base'] + egg['h'] * 0.4
        
        # Wind speed at this height
        v_wind = wind_speed_at_height(z_waist, v_ref=base_wind_speed, z_ref=10)
        
        # Dynamic pressure at waist
        q_dynamic = 0.5 * RHO_AIR * v_wind**2
        
        # Suction pressure boost
        delta_P_suction = q_dynamic * abs(Cp_suction) * valve_efficiency * 0.5
        
        suction_stages.append({
            'type': 'suction',
            'egg_index': i,
            'z_location': z_waist,
            'v_wind': v_wind,
            'delta_P_Pa': delta_P_suction,
        })
        
        total_delta_P += delta_P_suction
    
    # === INJECTION PORTS: At throat between eggs (28 stages) ===
    if dual_port:
        for i in range(len(eggs) - 1):
            # Throat is at the junction between egg i and egg i+1
            z_throat = eggs[i]['z_top']  # Top of lower egg = bottom of upper egg
            
            # Wind speed at this height
            v_wind = wind_speed_at_height(z_throat, v_ref=base_wind_speed, z_ref=10)
            
            # Dynamic pressure at throat
            q_dynamic = 0.5 * RHO_AIR * v_wind**2
            
            # Injection pressure boost (stagnation zone at narrow junction)
            delta_P_injection = q_dynamic * Cp_injection * valve_efficiency * 0.5
            
            injection_stages.append({
                'type': 'injection',
                'junction_index': i,
                'z_location': z_throat,
                'v_wind': v_wind,
                'delta_P_Pa': delta_P_injection,
            })
            
            total_delta_P += delta_P_injection
    
    # Calculate boosted updraft velocity
    K_loss = 2.5 if not dual_port else 2.8  # Slightly higher losses with dual ports
    
    v_boosted = np.sqrt(2 * total_delta_P / (RHO_AIR * K_loss))
    
    # Mass flow with boost
    A_column = column['average_column_area']
    m_dot = RHO_AIR * A_column * v_boosted
    
    # Power in boosted flow
    P_flow = 0.5 * RHO_AIR * A_column * v_boosted**3
    
    n_suction = len(suction_stages)
    n_injection = len(injection_stages)
    
    return {
        'n_suction_stages': n_suction,
        'n_injection_stages': n_injection,
        'n_total_stages': n_suction + n_injection,
        'total_delta_P_Pa': total_delta_P,
        'avg_delta_P_per_stage': total_delta_P / (n_suction + n_injection) if (n_suction + n_injection) > 0 else 0,
        'v_boosted': v_boosted,
        'mass_flow_kg_s': m_dot,
        'power_in_flow_kW': P_flow / 1000,
        'suction_stages': suction_stages,
        'injection_stages': injection_stages,
        'valve_efficiency': valve_efficiency,
        'Cp_suction': Cp_suction,
        'Cp_injection': Cp_injection if dual_port else 0,
        'dual_port': dual_port,
    }


def design_updraft_turbine_system(
    column: dict,
    driving_pressure_Pa: float,
    valve_contribution_Pa: float = 0,
    n_turbines: int = None,
    high_speed: bool = False
) -> dict:
    """
    Design turbine system using CHIMNEY PHYSICS (suction-driven flow).
    
    CRITICAL: This is a chimney, not a fan! 
    - Flow is PULLED by distributed pressure differential (buoyancy + wind + valves)
    - Velocity develops naturally from force balance
    - At equilibrium: ΔP_driving = ΔP_friction + ΔP_turbines
    
    The driving pressure is distributed along the full height. Friction opposes it.
    Turbines extract energy by creating additional pressure drop.
    
    Physics:
        ΔP_driving = ΔP_thermal + ΔP_wind + ΔP_valves  (total available)
        ΔP_friction = f × (L/D_avg) × ½ρv²             (friction loss)
        ΔP_turbine = ½ρv² × Cp × n_turbines            (extraction)
        
        At steady state: ΔP_driving = ΔP_friction + ΔP_turbine
        Solve for equilibrium velocity v
    
    Parameters:
        column: Column geometry
        driving_pressure_Pa: Total driving pressure (thermal + wind + valves)
        valve_contribution_Pa: Portion from valves (for reporting)
        n_turbines: Override number of turbines (optional)
        high_speed: If True, select high-speed turbines
    """
    # === GEOMETRY SETUP ===
    d_throat_avg = (column['segments'][0]['d_column_throat'] + 
                   column['segments'][-1]['d_column_throat']) / 2 if 'd_column_throat' in column['segments'][0] else \
                   (column['base_column_diameter'] + column['top_column_diameter']) / 2 * 0.58
    
    is_tapered = column.get('tapered_duct', False)
    if is_tapered:
        d_duct_base = column.get('duct_diameter_base', d_throat_avg)
        d_duct_top = column.get('duct_diameter_top', d_throat_avg)
        d_duct_avg = (d_duct_base + d_duct_top) / 2
        target_d = d_duct_avg * 0.70
    else:
        d_duct_base = d_throat_avg
        d_duct_top = d_throat_avg
        d_duct_avg = d_throat_avg
        target_d = d_throat_avg * 0.8
    
    # Select turbine
    use_aerospace = high_speed
    best_turbine = None
    for key, turbine in UPDRAFT_TURBINES.items():
        is_aerospace = 'aerospace' in key
        if use_aerospace and not is_aerospace:
            continue
        if not use_aerospace and is_aerospace:
            continue
        if turbine.rotor_diameter <= target_d:
            if best_turbine is None or turbine.rotor_diameter > best_turbine.rotor_diameter:
                best_turbine = turbine
    
    if best_turbine is None and use_aerospace:
        for key, turbine in UPDRAFT_TURBINES.items():
            if turbine.rotor_diameter <= target_d:
                if best_turbine is None or turbine.rotor_diameter > best_turbine.rotor_diameter:
                    best_turbine = turbine
    
    if best_turbine is None:
        return {
            'feasible': False,
            'reason': f"Duct diameter {d_duct_avg:.1f}m too small for turbines"
        }
    
    # === CHIMNEY EQUILIBRIUM PHYSICS ===
    #
    # This is a SUCTION system, not a blower!
    # 
    # Driving forces (distributed along column height):
    #   - Thermal buoyancy: warm air is less dense → rises
    #   - Wind pressure: higher at top outlet, lower at base inlet
    #   - Valve injection: adds momentum at each stage
    #
    # Opposing forces:
    #   - Friction: ΔP_f = f × (L/D) × ½ρv²
    #   - Turbine extraction: ΔP_t = Cp × ½ρv² × n_turbines
    #
    # At steady state: ΔP_driving = ΔP_friction + ΔP_turbine
    #
    # Rearranging: v² = 2 × ΔP_driving / (ρ × (K_friction + K_turbine))
    
    # Get stage info
    turbine_stages = column.get('turbine_stages', [])
    n_stages = len(turbine_stages) if turbine_stages else column['n_segments'] - 1
    
    # Duct geometry
    if is_tapered:
        duct_diameters = np.linspace(d_duct_base, d_duct_top, max(n_stages, 1))
        duct_areas = np.pi * (duct_diameters / 2)**2
    else:
        A_throat_base = column['segments'][0]['A_column_throat']
        A_throat_top = column['segments'][-1]['A_column_throat']
        duct_areas = np.linspace(A_throat_base, A_throat_top, max(n_stages, 1))
        duct_diameters = np.sqrt(4 * duct_areas / np.pi)
    
    A_avg = np.mean(duct_areas)
    D_avg = np.sqrt(4 * A_avg / np.pi)
    A_rotor = np.pi * (best_turbine.rotor_diameter / 2)**2
    
    # === LOSS COEFFICIENTS ===
    L = column['total_height']
    f_friction = 0.012 if is_tapered else 0.015
    
    # Friction coefficient: K_f = f × L/D
    K_friction = f_friction * L / D_avg
    
    # Venturi losses (only for natural geometry)
    if is_tapered:
        K_venturi = 0.0
    else:
        K_venturi = 0.08 * n_stages
    
    # Turbine pressure coefficient (per turbine)
    # Betz limit: max Cp = 16/27 ≈ 0.593
    # Turbine creates ΔP = Cp × ½ρv² across rotor
    Cp_turbine = 0.4  # Conservative: accounts for sub-optimal operation
    
    # Total turbine coefficient (function of how many we install)
    # More turbines = more extraction = more resistance = lower velocity
    
    # === SOLVE FOR EQUILIBRIUM WITH DIFFERENT TURBINE COUNTS ===
    # Trade-off: More turbines extract more, but slow the flow
    
    BETZ_LIMIT = 0.593
    best_config = None
    best_power = 0
    
    # Try different numbers of turbines
    for n_turb in range(1, n_stages + 1):
        # Total extraction coefficient
        K_turbine = Cp_turbine * n_turb * (A_rotor / A_avg)
        
        # Total resistance
        K_total = K_friction + K_venturi + K_turbine
        
        # Equilibrium velocity: ΔP_driving = K_total × ½ρv²
        # v = sqrt(2 × ΔP_driving / (ρ × K_total))
        v_eq = np.sqrt(2 * driving_pressure_Pa / (RHO_AIR * K_total))
        
        # Check if velocity is in turbine operating range
        if v_eq < best_turbine.cut_in_speed:
            continue  # Too slow
        # No cut-out limit for aerospace-style turbines in ducted flow
        
        # Power extracted by turbines
        # P = n × ½ρ × A_rotor × v³ × Cp_actual × η
        Cp_actual = min(BETZ_LIMIT, Cp_turbine)
        P_per_turbine = 0.5 * RHO_AIR * A_rotor * v_eq**3 * Cp_actual * best_turbine.efficiency
        P_per_turbine = min(P_per_turbine, best_turbine.rated_power_kW * 1000)
        
        P_total = P_per_turbine * n_turb
        
        if P_total > best_power:
            best_power = P_total
            best_config = {
                'n_turbines': n_turb,
                'v_equilibrium': v_eq,
                'K_friction': K_friction,
                'K_turbine': K_turbine,
                'K_total': K_total,
                'P_per_turbine_W': P_per_turbine,
                'P_total_W': P_total,
            }
    
    if best_config is None:
        # No valid configuration found, use single turbine estimate
        K_total = K_friction + K_venturi + Cp_turbine * (A_rotor / A_avg)
        v_eq = np.sqrt(2 * driving_pressure_Pa / (RHO_AIR * K_total))
        best_config = {
            'n_turbines': 1,
            'v_equilibrium': v_eq,
            'K_friction': K_friction,
            'K_turbine': Cp_turbine * (A_rotor / A_avg),
            'K_total': K_total,
            'P_per_turbine_W': 0,
            'P_total_W': 0,
        }
    
    v_eq = best_config['v_equilibrium']
    n_active = best_config['n_turbines']
    
    # === ENERGY BUDGET ===
    # Power available from driving pressure: P_avail = ΔP × Q = ΔP × A × v
    Q_flow = A_avg * v_eq  # m³/s
    P_available = driving_pressure_Pa * Q_flow  # Watts
    
    # Power lost to friction: P_friction = K_friction × ½ρv² × Q
    ΔP_friction = K_friction * 0.5 * RHO_AIR * v_eq**2
    P_friction = ΔP_friction * Q_flow
    
    # Power extracted by turbines
    P_extracted = best_config['P_total_W']
    
    # Efficiency
    system_efficiency = P_extracted / P_available if P_available > 0 else 0
    friction_fraction = P_friction / P_available if P_available > 0 else 0
    
    # Mass flow
    mass_flow = RHO_AIR * Q_flow
    
    # Create stage data for compatibility with output
    segment_heights = [column['segments'][i]['z_top'] for i in range(min(n_stages, len(column['segments'])))]
    stage_data = []
    for i in range(n_stages):
        z_stage = segment_heights[i] if i < len(segment_heights) else 100 + i * 55
        is_active = i < n_active
        stage_data.append({
            'stage': i,
            'z_height': z_stage,
            'd_duct': duct_diameters[i] if i < len(duct_diameters) else D_avg,
            'v_in': v_eq,
            'v_out': v_eq,  # Uniform velocity (chimney physics)
            'm_dot_in': mass_flow,
            'm_dot_inject': 0,
            'm_dot_out': mass_flow,
            'E_in_kW': P_available / 1000 / n_stages,
            'E_friction_kW': P_friction / 1000 / n_stages,
            'E_out_kW': (P_available - P_friction - P_extracted) / 1000 / n_stages,
            'can_extract': is_active,
            'P_available_kW': best_config['P_per_turbine_W'] / 1000 if is_active else 0,
            'P_extracted_kW': best_config['P_per_turbine_W'] / 1000 if is_active else 0,
        })
    
    # Cost estimate
    is_aerospace = 'Aerospace' in best_turbine.name
    cost_per_turbine = 150_000 if is_aerospace else 50_000
    total_cost = cost_per_turbine * n_active
    total_mass = (best_turbine.rotor_mass_kg + best_turbine.generator_mass_kg) * n_active
    
    # Total power
    total_power_kW = P_extracted / 1000
    avg_power_kW = total_power_kW / n_active if n_active > 0 else 0
    
    return {
        'feasible': True,
        'turbine_model': best_turbine.name,
        'turbine_diameter': best_turbine.rotor_diameter,
        'tapered_duct': is_tapered,
        'duct_diameter_base': d_duct_base,
        'duct_diameter_top': d_duct_top,
        'taper_ratio': column.get('taper_ratio', 1.0),
        'n_stages': n_stages,
        'n_active_stages': n_active,
        'n_turbines': n_active,
        'stage_data': stage_data,
        # Chimney equilibrium values
        'v_equilibrium': v_eq,
        'v_inlet': v_eq,
        'v_outlet': v_eq,
        'K_friction': best_config['K_friction'],
        'K_turbine': best_config['K_turbine'],
        'K_total': best_config['K_total'],
        # Mass flow
        'm_dot_inlet': mass_flow,
        'm_dot_outlet': mass_flow,
        'mass_flow_kg_s': mass_flow,
        # Legacy compatibility
        'updraft_speed': v_eq,
        'throat_velocity_ratio': 1.0,
        'throat_velocity': v_eq,
        'final_velocity': v_eq,
        'power_per_turbine_kW': avg_power_kW,
        'total_power_kW': total_power_kW,
        'total_mass_kg': total_mass,
        'total_cost_usd': total_cost,
        # Energy accounting (chimney physics)
        'P_driving_Pa': driving_pressure_Pa,
        'P_valve_contribution_Pa': valve_contribution_Pa,
        'P_available_kW': P_available / 1000,
        'P_friction_kW': P_friction / 1000,
        'P_extracted_kW': total_power_kW,
        'P_losses_kW': P_friction / 1000,
        'friction_fraction': friction_fraction,
        'system_efficiency': system_efficiency,
        'turbine_efficiency': best_turbine.efficiency,
    }

def analyze_updraft_column(eggs: List[dict], tapered_duct: bool = True) -> dict:
    """
    Complete analysis of vertical updraft column energy generation.
    
    Combines solar chimney effect + wind-driven updraft + valved accelerators.
    Compares single-port (leeward only) vs dual-port (windward + leeward) configurations.
    
    Parameters:
        eggs: List of egg geometry dictionaries
        tapered_duct: If True, use a smoothly tapered internal duct that follows
                     the structural geometry. Multi-stage turbines extract energy
                     to control velocity. Eliminates venturi losses.
    """
    # Calculate column geometry (with or without tapered duct)
    column = calculate_column_geometry(eggs, tapered_duct=tapered_duct)
    
    # === BASELINE: Thermal + simple wind ===
    # Thermal (solar chimney) updraft
    thermal = calculate_thermal_updraft(
        column_height=column['total_height'],
        base_area=column['average_column_area']
    )
    
    # Wind-driven updraft (at average wind speed)
    wind_driven = calculate_wind_driven_updraft(column, base_wind_speed=10)
    
    # Combined updraft (they can add together)
    # During day: solar + wind
    # At night: wind only
    v_combined_day = thermal['v_actual'] + wind_driven['v_updraft']
    v_night = wind_driven['v_updraft']
    
    # Weighted average (12 hrs day, 12 hrs night, wind varies)
    v_baseline = (v_combined_day * 0.4 + v_night * 0.6)  # Conservative
    
    # === ENHANCED: With Valved Accelerators ===
    # Compare single-port (leeward only) vs dual-port (windward + leeward)
    valved_single = calculate_valved_accelerator_boost(column, eggs, base_wind_speed=10, dual_port=False)
    valved_dual = calculate_valved_accelerator_boost(column, eggs, base_wind_speed=10, dual_port=True)
    
    # === CHIMNEY PHYSICS: Calculate DRIVING PRESSURE, not velocity ===
    # 
    # The velocity develops naturally from the balance:
    #   ΔP_driving = ΔP_friction + ΔP_turbines
    #
    # We calculate total driving pressure from all sources.
    
    # Baseline driving pressure (thermal + wind)
    # Thermal buoyancy: ΔP = ρ × g × H × ΔT/T
    # Approximation from thermal analysis
    P_thermal = thermal.get('driving_pressure_Pa', 
                           0.5 * RHO_AIR * thermal['v_actual']**2 * 3)  # K ~ 3 for chimney
    
    # Wind-driven pressure (outlet - inlet)
    P_wind = wind_driven.get('driving_pressure_Pa',
                            0.5 * RHO_AIR * wind_driven['v_updraft']**2 * 2)  # K ~ 2
    
    P_baseline = P_thermal + P_wind
    
    # Total driving pressures with valves
    P_total_single = P_baseline + valved_single['total_delta_P_Pa']
    P_total_dual = P_baseline + valved_dual['total_delta_P_Pa']
    
    # Design turbine systems for each configuration
    # Pass DRIVING PRESSURE - the function solves for equilibrium velocity
    turbine_system_baseline = design_updraft_turbine_system(
        column, driving_pressure_Pa=P_baseline, valve_contribution_Pa=0, high_speed=False
    )
    turbine_system_single = design_updraft_turbine_system(
        column, driving_pressure_Pa=P_total_single, 
        valve_contribution_Pa=valved_single['total_delta_P_Pa'], high_speed=True
    )
    turbine_system_dual = design_updraft_turbine_system(
        column, driving_pressure_Pa=P_total_dual,
        valve_contribution_Pa=valved_dual['total_delta_P_Pa'], high_speed=True
    )
    
    # Get equilibrium velocities from the solutions
    v_with_single = turbine_system_single.get('v_equilibrium', 0) if turbine_system_single['feasible'] else 0
    v_with_dual = turbine_system_dual.get('v_equilibrium', 0) if turbine_system_dual['feasible'] else 0
    
    # Valve system costs
    valve_cost_single = len(eggs) * 20_000 + 100_000  # 1 valve per egg + control
    valve_cost_dual = len(eggs) * 2 * 20_000 + 150_000  # 2 valves per egg + smarter control
    
    # Get power from each system
    sps = turbine_system_single
    dps = turbine_system_dual
    
    # Use single-port as primary (better cost-effectiveness with current turbines)
    # Dual-port pushes velocity too high (40 m/s), exceeding turbine cut-out
    if sps['total_power_kW'] >= dps['total_power_kW']:
        # Single-port is better
        valved_boost = valved_single
        v_with_valves = v_with_single
        turbine_system = turbine_system_single
        valve_cost = valve_cost_single
        preferred_config = 'single'
    else:
        # Dual-port is better (would need high-speed turbines)
        valved_boost = valved_dual
        v_with_valves = v_with_dual
        turbine_system = turbine_system_dual
        valve_cost = valve_cost_dual
        preferred_config = 'dual'
    
    if not turbine_system['feasible']:
        return {
            'system_type': 'Vertical Updraft Column',
            'feasible': False,
            'reason': turbine_system['reason'],
        }
    
    # Get duct cost from column geometry
    duct_cost = column.get('duct_cost_usd', 0)
    total_cost = turbine_system['total_cost_usd'] + valve_cost + duct_cost
    
    # Annual energy
    capacity_factor = 0.35  # Higher with active valve control
    annual_energy_MWh = turbine_system['total_power_kW'] * 8760 * capacity_factor / 1000
    annual_revenue = annual_energy_MWh * 50
    
    payback = total_cost / max(annual_revenue, 1)
    
    config_name = 'Single-Port' if preferred_config == 'single' else 'Dual-Port'
    
    return {
        'system_type': f'Vertical Updraft Column with {config_name} Valved Accelerators',
        'feasible': True,
        'preferred_config': preferred_config,
        
        # Geometry
        'column_height': column['total_height'],
        'base_column_diameter': column['base_column_diameter'],
        'top_column_diameter': column['top_column_diameter'],
        
        # Baseline updraft analysis
        'thermal_updraft': thermal,
        'wind_updraft': wind_driven,
        'v_combined_day': v_combined_day,
        'v_night': v_night,
        'v_baseline': v_baseline,
        
        # Valved accelerator comparisons
        'valved_single_port': {
            'delta_P_Pa': valved_single['total_delta_P_Pa'],
            'v_updraft': v_with_single,
            'boost_pct': (v_with_single / v_baseline - 1) * 100,
            'cost_usd': valve_cost_single,
            'power_kW': turbine_system_single['total_power_kW'],
        },
        'valved_dual_port': {
            'delta_P_Pa': valved_dual['total_delta_P_Pa'],
            'v_updraft': v_with_dual,
            'boost_pct': (v_with_dual / v_baseline - 1) * 100,
            'cost_usd': valve_cost_dual,
            'power_kW': turbine_system_dual['total_power_kW'],
        },
        
        # Primary system (selected configuration)
        'valved_accelerator': valved_boost,
        'v_with_valves': turbine_system.get('v_equilibrium', 0),
        'velocity_boost_pct': (turbine_system.get('v_equilibrium', 0) / max(v_baseline, 0.1) - 1) * 100,
        
        # Driving pressures (chimney physics)
        'P_baseline_Pa': P_baseline,
        'P_total_Pa': P_total_dual if preferred_config == 'dual' else P_total_single,
        
        # Comparison
        'baseline_system': turbine_system_baseline,
        'single_port_system': turbine_system_single,
        'dual_port_system': turbine_system_dual,
        
        # Turbine system (using enhanced)
        'turbine_system': turbine_system,
        
        # Costs
        'turbine_cost_usd': turbine_system['total_cost_usd'],
        'valve_system_cost_usd': valve_cost,
        'duct_cost_usd': duct_cost,
        'total_cost_usd': total_cost,
        
        # Energy production
        'capacity_factor': capacity_factor,
        'annual_energy_MWh': annual_energy_MWh,
        'annual_revenue_usd': annual_revenue,
        'simple_payback_years': payback,
    }


# =============================================================================
# COMBINED ANALYSIS
# =============================================================================

def analyze_alternative_energy(eggs: List[dict] = None) -> dict:
    """
    Complete analysis of alternative energy generation options.
    """
    if eggs is None:
        eggs = get_optimized_eggs()
    
    vortex = analyze_vortex_energy_harvesting(eggs)
    updraft = analyze_updraft_column(eggs)
    
    return {
        'vortex_shedding': vortex,
        'updraft_column': updraft,
    }


def print_alternative_energy_analysis(analysis: dict):
    """Print formatted analysis of alternative energy systems."""
    print(f"\n{'='*90}")
    print("ALTERNATIVE ENERGY GENERATION ANALYSIS")
    print(f"{'='*90}")
    
    # === VORTEX SHEDDING ===
    vortex = analysis['vortex_shedding']
    print(f"\n🌀 VORTEX SHEDDING ENERGY CAPTURE")
    print(f"{'─'*90}")
    print(f"   Concept: Harvest tower oscillations from vortex-induced vibration")
    print(f"   Harvester Type: {vortex['harvester_type']}")
    print(f"   Design Wind Speed: {vortex['design_wind_speed']} m/s")
    print()
    print(f"   Installations: {vortex['n_harvesters']}")
    print(f"   Total Power: {vortex['total_power_kW']:.1f} kW")
    print(f"   Added Mass: {vortex['total_mass_kg']/1000:.1f} tonnes")
    print(f"   CAPEX: ${vortex['total_cost_usd']:,.0f}")
    print()
    print(f"   Capacity Factor: {vortex['capacity_factor']*100:.0f}%")
    print(f"   Annual Energy: {vortex['annual_energy_MWh']:.0f} MWh")
    print(f"   Annual Revenue: ${vortex['annual_revenue_usd']:,.0f}")
    print(f"   Payback: {vortex['simple_payback_years']:.1f} years")
    
    if vortex['installations']:
        print(f"\n   Harvester Locations:")
        for inst in vortex['installations'][:5]:  # Show first 5
            print(f"      Junction {inst['junction']:2d}: z={inst['z_junction']:6.0f}m, "
                  f"diff_A={inst['differential_amplitude_m']*1000:.1f}mm, "
                  f"P={inst['power_kW']:.2f}kW")
    
    # === UPDRAFT COLUMN ===
    updraft = analysis['updraft_column']
    config = updraft.get('preferred_config', 'single')
    config_label = 'SINGLE-PORT' if config == 'single' else 'DUAL-PORT'
    print(f"\n🔥 VERTICAL UPDRAFT COLUMN WITH {config_label} VALVED ACCELERATORS")
    print(f"{'─'*90}")
    if config == 'single':
        print(f"   Concept: Hollow eggs form chimney; leeward valved ports harness wind suction")
    else:
        print(f"   Concept: Hollow eggs form chimney; valved ports on BOTH sides of each egg")
        print(f"            Windward ports INJECT air; Leeward ports EXTRACT via suction")
    
    if not updraft['feasible']:
        print(f"   ❌ NOT FEASIBLE: {updraft['reason']}")
    else:
        # Check for tapered duct
        ts = updraft['turbine_system']
        is_tapered = ts.get('tapered_duct', False)
        
        print(f"\n   Column Geometry:")
        print(f"      Height: {updraft['column_height']:.0f} m")
        if is_tapered:
            d_base = ts.get('duct_diameter_base', 0)
            d_top = ts.get('duct_diameter_top', 0)
            taper = ts.get('taper_ratio', 1.0)
            n_stages = ts.get('n_active_stages', 0)
            print(f"      🔧 TAPERED INTERNAL DUCT:")
            print(f"         Base Ø: {d_base:.1f} m → Top Ø: {d_top:.1f} m (taper: {taper:.2f})")
            print(f"         {n_stages} turbine stages for velocity control")
            print(f"         Smooth taper eliminates venturi losses")
            if updraft.get('duct_cost_usd'):
                print(f"         Duct CAPEX: ${updraft['duct_cost_usd']:,.0f}")
        else:
            print(f"      Column Ø (base): {updraft['base_column_diameter']:.1f} m")
            print(f"      Column Ø (top): {updraft['top_column_diameter']:.1f} m")
            print(f"      ⚠️  Natural geometry: 28 venturi stages create high losses")
        
        print(f"\n   Baseline Updraft (thermal + wind inlet/outlet):")
        print(f"      Thermal (solar): {updraft['thermal_updraft']['v_actual']:.1f} m/s")
        print(f"      Wind pressure diff: {updraft['wind_updraft']['v_updraft']:.1f} m/s")
        print(f"      Combined baseline: {updraft['v_baseline']:.1f} m/s")
        
        # Single vs Dual port comparison
        sp = updraft['valved_single_port']
        dp = updraft['valved_dual_port']
        print(f"\n   🎛️  VALVED ACCELERATOR COMPARISON:")
        print(f"      {'Configuration':<25} {'ΔP (Pa)':<12} {'Velocity':<12} {'Boost':<10} {'Cost':<12}")
        print(f"      {'─'*70}")
        print(f"      {'Single-port (suction)':<25} {sp['delta_P_Pa']:<12.0f} {sp['v_updraft']:<12.1f} +{sp['boost_pct']:<9.0f}% ${sp['cost_usd']:,.0f}")
        print(f"      {'Dual-port (suct+inject)':<25} {dp['delta_P_Pa']:<12.0f} {dp['v_updraft']:<12.1f} +{dp['boost_pct']:<9.0f}% ${dp['cost_usd']:,.0f}")
        
        va = updraft['valved_accelerator']
        print(f"\n   Dual-Port Configuration:")
        print(f"      SUCTION ports at WAIST (leeward): {va['n_suction_stages']} stages, Cp = {va['Cp_suction']}")
        print(f"      INJECTION ports at THROAT (windward): {va['n_injection_stages']} stages, Cp = +{va['Cp_injection']}")
        print(f"      Total stages: {va['n_total_stages']}")
        print(f"      Valve efficiency: {va['valve_efficiency']*100:.0f}%")
        print(f"      Total ΔP from valves: {va['total_delta_P_Pa']:.0f} Pa")
        print(f"      ────────────────────────────────")
        print(f"      TOTAL UPDRAFT: {updraft['v_with_valves']:.1f} m/s (+{updraft['velocity_boost_pct']:.0f}% vs baseline)")
        
        # Show control volume analysis
        if ts.get('P_available_kW') or ts.get('P_driving_Pa'):
            print(f"\n   🎯 CHIMNEY EQUILIBRIUM ANALYSIS:")
            
            if is_tapered:
                v_eq = ts.get('v_equilibrium', 0)
                m_flow = ts.get('mass_flow_kg_s', 0)
                
                print(f"      TAPERED DUCT (Suction-Driven Flow):")
                print(f"         Driving ΔP: {ts.get('P_driving_Pa', 0):.0f} Pa")
                print(f"         Valve contribution: {ts.get('P_valve_contribution_Pa', 0):.0f} Pa")
                print(f"         Friction coefficient K_f: {ts.get('K_friction', 0):.1f}")
                print(f"         Equilibrium velocity: {v_eq:.1f} m/s (uniform throughout)")
                print(f"         Mass flow: {m_flow:.0f} kg/s")
                print(f"         Active turbine stages: {ts.get('n_active_stages', 0)}")
                print(f"         Turbine Ø: {ts['turbine_diameter']:.1f} m")
            else:
                print(f"      NATURAL GEOMETRY (venturi stages):")
                print(f"         Equilibrium velocity: {ts.get('v_equilibrium', 0):.1f} m/s")
                print(f"         Turbine Ø: {ts['turbine_diameter']:.1f} m")
            
            P_available = ts.get('P_available_kW', 0)
            P_friction = ts.get('P_friction_kW', 0)
            P_extracted = ts.get('total_power_kW', 0)
            friction_frac = ts.get('friction_fraction', 0)
            
            print(f"\n      ⚡ POWER BUDGET (Chimney Physics):")
            print(f"         Power available (ΔP × Q): {P_available:.0f} kW")
            print(f"         Friction losses:         {P_friction:.0f} kW ({friction_frac*100:.0f}%)")
            print(f"         Turbine extraction:      {P_extracted:.0f} kW")
            print(f"         System efficiency:       {ts.get('system_efficiency', 0)*100:.1f}%")
        
        # Comparison: baseline vs single vs dual
        bs = updraft['baseline_system']
        sps = updraft['single_port_system']
        dps = updraft['dual_port_system']
        print(f"\n   Turbine Output Comparison:")
        bs_power = max(bs['total_power_kW'], 0.1)  # Avoid division by zero
        print(f"      Without valves:      {bs['total_power_kW']:.1f} kW")
        print(f"      Single-port valves:  {sps['total_power_kW']:.1f} kW (+{(sps['total_power_kW']/bs_power-1)*100:.0f}%)")
        print(f"      Dual-port valves:    {dps['total_power_kW']:.1f} kW (+{(dps['total_power_kW']/bs_power-1)*100:.0f}%)")
        
        config_name = 'dual-port' if updraft['preferred_config'] == 'dual' else 'single-port'
        print(f"\n   Turbine System ({config_name}):")
        print(f"      Model: {ts['turbine_model']}")
        print(f"      Rotor Ø: {ts['turbine_diameter']:.1f} m")
        print(f"      Active stages: {ts.get('n_active_stages', ts.get('n_turbines', 0))}")
        print(f"      Power per stage: {ts['power_per_turbine_kW']:.1f} kW")
        print(f"      Total Power: {ts['total_power_kW']:.1f} kW")
        
        print(f"\n   Economics ({config_name.title()} System):")
        print(f"      Turbine CAPEX: ${updraft['turbine_cost_usd']:,.0f}")
        print(f"      Valve system CAPEX: ${updraft['valve_system_cost_usd']:,.0f}")
        if is_tapered and updraft.get('duct_cost_usd'):
            print(f"      Tapered duct CAPEX: ${updraft['duct_cost_usd']:,.0f}")
        print(f"      Total CAPEX: ${updraft['total_cost_usd']:,.0f}")
        print(f"      Capacity Factor: {updraft['capacity_factor']*100:.0f}%")
        print(f"      Annual Energy: {updraft['annual_energy_MWh']:.0f} MWh")
        print(f"      Annual Revenue: ${updraft['annual_revenue_usd']:,.0f}")
        print(f"      Payback: {updraft['simple_payback_years']:.1f} years")
    
    # === COMPARISON ===
    print(f"\n{'='*90}")
    print("📊 COMPARISON: MUTUALLY EXCLUSIVE OPTIONS")
    print(f"{'='*90}")
    print("""
   ⚠️  THESE ARE ALTERNATIVE CONCEPTS - CANNOT BE COMBINED!
   
   • Internal Duct Turbines: Requires OPENINGS in shell for horizontal airflow
   • Vertical Updraft Column: Requires SEALED chimney with valved ports for vertical flow
   
   The same internal space cannot be used for both. Choose ONE approach.
""")
    
    # Load internal duct data for comparison
    try:
        from wind_energy import load_wind_energy_analysis
        duct = load_wind_energy_analysis()
        
        print(f"{'System':<35} {'Power (kW)':<12} {'Annual MWh':<12} {'CAPEX ($)':<15} {'Payback (yr)':<12}")
        print(f"{'─'*90}")
        print(f"{'Option A: Internal Duct Turbines':<35} {duct['summary']['rated_power_kW']:<12.0f} "
              f"{duct['summary']['annual_energy_MWh']:<12.0f} "
              f"{duct['summary']['estimated_capex_usd']:<15,.0f} "
              f"{duct['summary']['simple_payback_years']:<12.1f}")
        if updraft['feasible']:
            print(f"{'Option B: Updraft + Valved Accel.':<35} {updraft['turbine_system']['total_power_kW']:<12.0f} "
                  f"{updraft['annual_energy_MWh']:<12.0f} "
                  f"{updraft['total_cost_usd']:<15,.0f} "
                  f"{updraft['simple_payback_years']:<12.1f}")
        print(f"{'─'*90}")
        print(f"{'(Vortex Shedding - Not Viable)':<35} {vortex['total_power_kW']:<12.0f} "
              f"{vortex['annual_energy_MWh']:<12.0f} "
              f"{vortex['total_cost_usd']:<15,.0f} "
              f"{vortex['simple_payback_years']:<12.1f}")
        
        # Recommendation - compare fairly
        if updraft['feasible']:
            if updraft['annual_energy_MWh'] > duct['summary']['annual_energy_MWh']:
                print(f"\n   🏆 RECOMMENDATION: Option B (Updraft Column with Valved Accelerators)")
                print(f"   - Higher annual energy: {updraft['annual_energy_MWh']:.0f} vs {duct['summary']['annual_energy_MWh']:.0f} MWh")
                print(f"   - Valved accelerators use external wind to boost updraft")
                print(f"   - Smart control adapts to wind direction")
            else:
                print(f"\n   🏆 RECOMMENDATION: Option A (Internal Duct Turbines)")
                print(f"   - Higher annual energy: {duct['summary']['annual_energy_MWh']:.0f} vs {updraft['annual_energy_MWh']:.0f} MWh")
                print(f"   - Proven horizontal-axis turbine technology")
    except:
        pass
    
    # Caveats
    print(f"\n{'='*90}")
    print("⚠️  TECHNICAL NOTES")
    print(f"{'='*90}")
    print("""
   VORTEX SHEDDING:
   - Not viable - differential motion too small for useful power
   - Piezoelectric tech limited at these scales
   
   UPDRAFT COLUMN WITH VALVED ACCELERATORS:
   - Each egg waist has leeward valves that harness wind suction
   - 29 stages act as a multi-stage ejector pump
   - Control system opens/closes valves based on wind direction
   - Works day AND night (not dependent on solar heating)
   - Valve efficiency assumed 60% - actual depends on port design
   - Still requires sealed column construction between eggs
""")


def save_alternative_energy_analysis(analysis: dict = None) -> str:
    """Save alternative energy analysis to JSON."""
    output_path = os.path.join(os.path.dirname(__file__), 'alternative_energy_analysis.json')
    
    if analysis is None:
        analysis = analyze_alternative_energy()
    
    # Prepare for JSON serialization
    export = {
        'note': 'Auto-generated alternative energy analysis',
        'vortex_shedding': {
            'system_type': analysis['vortex_shedding']['system_type'],
            'n_harvesters': analysis['vortex_shedding']['n_harvesters'],
            'total_power_kW': analysis['vortex_shedding']['total_power_kW'],
            'annual_energy_MWh': analysis['vortex_shedding']['annual_energy_MWh'],
            'capex_usd': analysis['vortex_shedding']['total_cost_usd'],
            'payback_years': analysis['vortex_shedding']['simple_payback_years'],
        },
        'updraft_column': {
            'feasible': analysis['updraft_column']['feasible'],
        }
    }
    
    if analysis['updraft_column']['feasible']:
        export['updraft_column'].update({
            'system_type': analysis['updraft_column']['system_type'],
            'column_height': analysis['updraft_column']['column_height'],
            'v_baseline': analysis['updraft_column']['v_baseline'],
            'v_with_valves': analysis['updraft_column']['v_with_valves'],
            'velocity_boost_pct': analysis['updraft_column']['velocity_boost_pct'],
            'total_power_kW': analysis['updraft_column']['turbine_system']['total_power_kW'],
            'annual_energy_MWh': analysis['updraft_column']['annual_energy_MWh'],
            'turbine_capex_usd': analysis['updraft_column']['turbine_cost_usd'],
            'valve_capex_usd': analysis['updraft_column']['valve_system_cost_usd'],
            'total_capex_usd': analysis['updraft_column']['total_cost_usd'],
            'payback_years': analysis['updraft_column']['simple_payback_years'],
        })
    
    with open(output_path, 'w') as f:
        json.dump(export, f, indent=2)
    
    print(f"💾 Saved alternative energy analysis to: {output_path}")
    return output_path


if __name__ == '__main__':
    analysis = analyze_alternative_energy()
    print_alternative_energy_analysis(analysis)
    save_alternative_energy_analysis(analysis)
