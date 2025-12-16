"""
Wind loading and cable force calculations including vortex shedding.
"""
import numpy as np
from config import RHO_AIR, V_REF, CD, CD_CABLE, CL_VORTEX, STROUHAL, CABLE_MATERIALS


def get_guy_indices(n_eggs, n_guy_levels=3):
    """Get the egg indices where guy cables are attached."""
    if n_guy_levels >= n_eggs:
        return list(range(n_eggs))
    return [int(i * (n_eggs - 1) / (n_guy_levels - 1)) for i in range(n_guy_levels)]


def wind_speed_at_height(z, v_ref=V_REF, z_ref=10):
    """
    Wind speed varies with height according to power law.
    v(z) = v_ref × (z/z_ref)^α where α ≈ 0.16 for open terrain
    """
    alpha = 0.16
    return v_ref * (z / z_ref) ** alpha


def vortex_shedding_frequency(v, D, St=STROUHAL):
    """
    Calculate vortex shedding frequency.
    
    f_s = St × v / D
    
    Parameters:
        v: wind speed (m/s)
        D: characteristic diameter (m)
        St: Strouhal number (≈0.2 for cylinders/ellipsoids)
    
    Returns:
        Shedding frequency (Hz)
    """
    return St * v / D


def vortex_lift_force(egg, rho_air=RHO_AIR, Cl=CL_VORTEX):
    """
    Calculate peak cross-wind lift force due to vortex shedding.
    
    F_lift = 0.5 × ρ × v² × Cl × A_projected
    
    This is the oscillating force perpendicular to wind direction.
    Actual force varies sinusoidally at shedding frequency.
    """
    z_center = egg['z_base'] + egg['h'] / 2
    v = wind_speed_at_height(z_center)
    
    # Dynamic pressure
    q = 0.5 * rho_air * v**2
    
    # Projected area (side view)
    a = egg['h'] / 2
    b = egg['d'] / 2
    A_proj = np.pi * a * b
    
    # Peak lift force (oscillating)
    F_lift = q * Cl * A_proj
    
    # Shedding frequency
    f_s = vortex_shedding_frequency(v, egg['d'])
    
    return {
        'F_lift': F_lift,
        'shedding_freq': f_s,
        'wind_speed': v
    }


def cable_wind_drag(n_strands, d_strand, cable_length, z_center, rho_air=RHO_AIR, Cd_cable=CD_CABLE):
    """
    Calculate wind drag on cables between two eggs.
    
    Cables are distributed around the circumference. The projected area
    depends on each cable's angle relative to the wind direction:
    
    - Cables at θ=0° (windward/leeward): parallel to wind, minimal projected area
    - Cables at θ=90° (sides): perpendicular to wind, full projected area
    
    For a cable at angle θ from the wind direction:
        Projected width = d_strand × |sin(θ)|
    
    Integrating over the full circumference (0 to π for one side):
        Effective area factor = (1/π) × ∫|sin(θ)|dθ from 0 to π = 2/π ≈ 0.637
    
    This means the effective projected area is about 64% of the total
    frontal area if all cables were perpendicular to the wind.
    
    Parameters:
        n_strands: total number of strands around circumference
        d_strand: diameter of each strand (m)
        cable_length: length of cable segment (m)
        z_center: height at center of cable segment (m)
    
    Returns:
        Dictionary with cable drag force and breakdown
    """
    # Wind speed at cable height
    v = wind_speed_at_height(z_center)
    q = 0.5 * rho_air * v**2
    
    # Total frontal area if all cables were perpendicular
    A_total = n_strands * d_strand * cable_length
    
    # Effective area factor: cables are distributed around circumference
    # Average of |sin(θ)| over 0 to π = 2/π ≈ 0.637
    ANGULAR_FACTOR = 2 / np.pi  # ≈ 0.637
    
    A_effective = A_total * ANGULAR_FACTOR
    
    # Wind drag force on cables
    F_cable_drag = q * Cd_cable * A_effective
    
    return {
        'n_strands': n_strands,
        'd_strand': d_strand,
        'cable_length': cable_length,
        'A_total': A_total,
        'A_effective': A_effective,
        'angular_factor': ANGULAR_FACTOR,
        'F_cable_drag': F_cable_drag,
        'F_cable_drag_kN': F_cable_drag / 1000
    }


def estimate_tower_natural_frequency(eggs, E=45e9, rho_shell=2000):
    """
    Estimate the fundamental natural frequency of the tower.
    
    Simplified cantilever beam approximation:
    f_n = (1.875²/2π) × √(EI/ρAL⁴)
    
    For a tapered structure, we use average properties.
    """
    n = len(eggs)
    total_height = eggs[-1]['z_base'] + eggs[-1]['h']
    
    # Average diameter and moment of inertia
    d_avg = np.mean([e['d'] for e in eggs])
    t_avg = d_avg / 1000  # approximate shell thickness
    
    # Moment of inertia for thin-walled cylinder
    r = d_avg / 2
    I = np.pi * r**3 * t_avg
    
    # Mass per unit length
    A_shell = 2 * np.pi * r * t_avg
    m_per_length = A_shell * rho_shell
    
    # Cantilever fundamental frequency
    lambda_1 = 1.875  # first mode eigenvalue
    f_n = (lambda_1**2 / (2 * np.pi)) * np.sqrt(E * I / (m_per_length * total_height**4))
    
    return f_n


def check_vortex_lock_in(f_shedding, f_natural, lock_in_range=(0.8, 1.2)):
    """
    Check if vortex shedding frequency is in lock-in range with natural frequency.
    
    Lock-in occurs when: 0.8 < f_s/f_n < 1.2 (approximately)
    This causes resonance and greatly amplified oscillations.
    
    Returns:
        dict with lock-in status and frequency ratio
    """
    ratio = f_shedding / f_natural if f_natural > 0 else float('inf')
    in_lock_in = lock_in_range[0] < ratio < lock_in_range[1]
    
    return {
        'f_shedding': f_shedding,
        'f_natural': f_natural,
        'ratio': ratio,
        'in_lock_in': in_lock_in,
        'risk_level': 'CRITICAL' if in_lock_in else ('HIGH' if 0.5 < ratio < 1.5 else 'LOW')
    }


def analyze_vortex_shedding(eggs):
    """
    Comprehensive vortex shedding analysis for the tower.
    
    Returns analysis for each egg level and overall assessment.
    """
    # Estimate tower natural frequency
    f_natural = estimate_tower_natural_frequency(eggs)
    
    results = []
    max_lift_force = 0
    critical_eggs = []
    
    for i, egg in enumerate(eggs):
        vortex = vortex_lift_force(egg)
        lock_in = check_vortex_lock_in(vortex['shedding_freq'], f_natural)
        
        result = {
            'egg_index': i,
            'z': egg['z_base'] + egg['h'] / 2,
            'd': egg['d'],
            'wind_speed': vortex['wind_speed'],
            'F_lift_kN': vortex['F_lift'] / 1000,
            'shedding_freq': vortex['shedding_freq'],
            'lock_in': lock_in
        }
        results.append(result)
        
        max_lift_force = max(max_lift_force, vortex['F_lift'])
        if lock_in['in_lock_in']:
            critical_eggs.append(i)
    
    return {
        'f_natural': f_natural,
        'eggs': results,
        'max_lift_force_kN': max_lift_force / 1000,
        'critical_eggs': critical_eggs,
        'has_lock_in_risk': len(critical_eggs) > 0
    }


def wind_force_on_egg(egg, rho_air=RHO_AIR, Cd=CD, Cl=CL_VORTEX, include_vortex=True):
    """
    Calculate wind force on a single egg including vortex shedding effects.
    
    Drag force (along-wind): F_drag = 0.5 × ρ × v² × Cd × A_projected
    Lift force (cross-wind): F_lift = 0.5 × ρ × v² × Cl × A_projected (oscillating)
    
    For structural design, we consider the vector sum of drag and peak lift.
    """
    z_center = egg['z_base'] + egg['h'] / 2
    v = wind_speed_at_height(z_center)
    
    # Dynamic pressure
    q = 0.5 * rho_air * v**2
    
    # Projected area (side view of ellipsoid)
    a = egg['h'] / 2
    b = egg['d'] / 2
    A_proj = np.pi * a * b
    
    # Drag force (along-wind)
    F_drag = q * Cd * A_proj
    
    # Vortex-induced lift force (cross-wind, oscillating)
    F_lift = q * Cl * A_proj if include_vortex else 0
    
    # Combined force for cable sizing (SRSS - Square Root Sum of Squares)
    # This is conservative: assumes peak lift occurs with full drag
    F_combined = np.sqrt(F_drag**2 + F_lift**2)
    
    # Vortex shedding frequency
    f_s = vortex_shedding_frequency(v, egg['d'])
    
    return {
        'z': z_center,
        'wind_speed': v,
        'dynamic_pressure': q,
        'projected_area': A_proj,
        'wind_force': F_combined,  # Use combined force for cable sizing
        'F_drag': F_drag,
        'F_lift': F_lift,
        'shedding_freq': f_s,
        'drag_lift_ratio': F_lift / F_drag if F_drag > 0 else 0
    }


def cable_tension_equilibrium(eggs, level_idx, F_eggs_cumulative, cable_length, sin_from_vert, 
                                z_center, n_cables_per_side, theta_angles, cos_thetas, sum_cos,
                                F_cable_drag_above=0):
    """
    Find equilibrium cable tension at a level, accounting for cable self-drag.
    
    The cable tension T creates a cable size, which has drag, which increases T.
    We solve: T = f(T) where f includes cable drag contribution.
    
    Uses scipy.optimize.brentq for robust root finding.
    """
    from scipy.optimize import brentq
    from config import CABLE_MATERIAL
    
    def tension_residual(T_guess):
        """
        Residual function: T_computed - T_guess = 0 at equilibrium.
        """
        # Given T_guess, calculate total strands needed
        total_strands = 0
        for j, theta in enumerate(theta_angles):
            # Tension at this angular position
            T_at_angle = T_guess * cos_thetas[j] / cos_thetas[0]  # Scale from max
            total_strands += CABLE_MATERIAL.fiber_count(T_at_angle) * 2  # Both sides
        
        # Cable drag for these strands
        drag_result = cable_wind_drag(
            total_strands,
            CABLE_MATERIAL.fiber_diameter,
            cable_length,
            z_center
        )
        F_cable_drag = drag_result['F_cable_drag']
        
        # Total horizontal force including cable drag
        F_total = F_eggs_cumulative + F_cable_drag_above + F_cable_drag
        
        # Resulting max tension (at θ=0)
        F_horiz_max = F_total * cos_thetas[0] / sum_cos
        T_computed = F_horiz_max / sin_from_vert if sin_from_vert > 0.01 else F_horiz_max * 100
        
        return T_computed - T_guess
    
    # Find bounds: minimum is eggs-only, maximum is some multiple
    F_horiz_max_eggs = (F_eggs_cumulative + F_cable_drag_above) * cos_thetas[0] / sum_cos
    T_min = F_horiz_max_eggs / sin_from_vert if sin_from_vert > 0.01 else F_horiz_max_eggs * 100
    T_max = T_min * 20  # Cable drag shouldn't increase tension by more than 20x
    
    # Check if solution exists in range
    r_min = tension_residual(T_min)
    r_max = tension_residual(T_max)
    
    if r_min * r_max > 0:
        # No sign change - either no solution or need wider bounds
        # If r_min > 0, even minimum tension gives higher computed T (runaway)
        # In this case, the system is unstable - return the max as a warning
        if r_min > 0:
            return T_max, T_max, True  # Unstable flag
        else:
            return T_min, 0, False
    
    # Find equilibrium
    T_equilibrium = brentq(tension_residual, T_min, T_max, xtol=1000)  # 1 kN tolerance
    
    # Calculate the cable drag at equilibrium
    total_strands = 0
    for j, theta in enumerate(theta_angles):
        T_at_angle = T_equilibrium * cos_thetas[j] / cos_thetas[0]
        total_strands += CABLE_MATERIAL.fiber_count(T_at_angle) * 2
    
    drag_result = cable_wind_drag(total_strands, CABLE_MATERIAL.fiber_diameter, cable_length, z_center)
    
    return T_equilibrium, drag_result['F_cable_drag'], False


def calculate_cable_forces(eggs, n_cables_per_side=8, include_cable_drag=False):
    """
    Calculate tensile force in each cable segment with distributed load.
    
    Cables are distributed around the egg circumference. Wind load is shared
    based on each cable's angular position relative to the wind direction.
    
    For n cables distributed from 0° to 90° on one side (mirrored on other):
    - Cable at θ=0° (windward) takes full horizontal component
    - Cable at θ=90° (side) takes zero horizontal component
    - Force on cable at angle θ: F_cable(θ) = F_total × cos(θ) / Σcos(θ)
    
    NOTE on cable drag: With near-vertical cables (θ < 2°), cable self-drag
    creates an unstable feedback loop. The current egg-to-egg cable geometry
    has angles of ~0.5°, making cable drag feedback unstable (k = 5.0 > 1).
    
    For stability, cable angles need sin(θ) > drag_ratio (~0.04), i.e., θ > 2.3°.
    Real tall structures use ground-anchored guys at steeper angles (30-45°)
    for horizontal wind resistance, with egg-to-egg cables only for continuity.
    
    Parameters:
        eggs: list of egg dictionaries
        n_cables_per_side: number of cable attachment points per quarter (0° to 90°)
        include_cable_drag: whether to include wind drag on cables (unstable, off by default)
    """
    from config import CABLE_MATERIAL
    
    n = len(eggs)
    cable_forces = []
    
    # Calculate wind force on each egg (drag + vortex lift)
    for egg in eggs:
        wind_data = wind_force_on_egg(egg)
        egg['wind_force'] = wind_data['wind_force']
        egg['wind_speed'] = wind_data['wind_speed']
        egg['F_egg_only'] = wind_data['wind_force']
    
    # Angular positions for cables (0° = windward, 90° = perpendicular)
    theta_angles = np.linspace(0, np.pi/2, n_cables_per_side)
    cos_thetas = np.cos(theta_angles)
    sum_cos = np.sum(cos_thetas) * 2  # both sides of windward
    
    # Track cumulative cable drag from levels above
    cable_drag_above = 0
    any_unstable = False
    
    # Process from top down
    for i in range(n - 2, -1, -1):  # Start from second-to-top, go to bottom
        # Egg forces from this level and above
        F_eggs_cumulative = sum(e['wind_force'] for e in eggs[i+1:])
        F_eggs_cumulative += eggs[i]['wind_force'] / 2
        
        current_egg = eggs[i]
        next_egg = eggs[i + 1]
        
        z1 = current_egg['z_base'] + current_egg['h'] / 2
        z2 = next_egg['z_base'] + next_egg['h'] / 2
        x1 = current_egg['d'] / 2
        x2 = next_egg['d'] / 2
        
        dz = z2 - z1
        dx = x1 - x2
        cable_length = np.sqrt(dz**2 + dx**2)
        sin_from_vert = dx / cable_length if cable_length > 0 else 0.01
        z_center = (z1 + z2) / 2
        
        if include_cable_drag:
            T_max, F_cable_drag, unstable = cable_tension_equilibrium(
                eggs, i, F_eggs_cumulative, cable_length, sin_from_vert,
                z_center, n_cables_per_side, theta_angles, cos_thetas, sum_cos,
                cable_drag_above
            )
            if unstable:
                any_unstable = True
            cable_drag_above += F_cable_drag
        else:
            F_horiz_max = F_eggs_cumulative * cos_thetas[0] / sum_cos
            T_max = F_horiz_max / sin_from_vert if sin_from_vert > 0.01 else F_horiz_max * 100
            F_cable_drag = 0
        
        # Calculate tensions at all angular positions
        F_horizontal_total = F_eggs_cumulative + (cable_drag_above if include_cable_drag else 0)
        cable_tensions_by_angle = []
        for j, theta in enumerate(theta_angles):
            F_horiz_cable = F_horizontal_total * cos_thetas[j] / sum_cos
            T_cable = T_max * cos_thetas[j] / cos_thetas[0]  # Scale from max
            cable_tensions_by_angle.append({
                'theta_deg': np.degrees(theta),
                'F_horizontal': F_horiz_cable,
                'T_cable': T_cable,
                'T_cable_kN': T_cable / 1000
            })
        
        cable_forces.append({
            'from_egg': i,
            'to_egg': i + 1,
            'z1': z1,
            'z2': z2,
            'd1': current_egg['d'],
            'cable_length': cable_length,
            'cable_angle_deg': np.degrees(np.arctan2(dz, dx)),
            'F_horizontal_total': F_horizontal_total,
            'F_cable_drag': F_cable_drag,
            'n_cables_per_side': n_cables_per_side,
            'n_cables_total': n_cables_per_side * 2,
            'cables_by_angle': cable_tensions_by_angle,
            'T_cable': T_max,
            'T_cable_kN': T_max / 1000,
            'unstable': unstable if include_cable_drag else False
        })
    
    # Reverse to get bottom-to-top order
    cable_forces = cable_forces[::-1]
    
    # Ground cables (no cable drag iteration needed - fixed geometry)
    F_horizontal_ground = sum(e['wind_force'] for e in eggs) + cable_drag_above
    
    bottom_egg = eggs[0]
    z1 = bottom_egg['h'] / 2
    z2 = 0
    x1 = bottom_egg['d'] / 2
    x2 = bottom_egg['d'] * 1.2
    
    dz = z1 - z2
    dx = x2 - x1
    cable_length = np.sqrt(dz**2 + dx**2)
    sin_from_vert = dx / cable_length if cable_length > 0 else 0.01
    
    ground_tensions_by_angle = []
    for j, theta in enumerate(theta_angles):
        F_horiz_cable = F_horizontal_ground * cos_thetas[j] / sum_cos
        T_cable = F_horiz_cable / sin_from_vert if sin_from_vert > 0.01 else F_horiz_cable * 100
        ground_tensions_by_angle.append({
            'theta_deg': np.degrees(theta),
            'F_horizontal': F_horiz_cable,
            'T_cable': T_cable,
            'T_cable_kN': T_cable / 1000
        })
    
    T_ground_max = ground_tensions_by_angle[0]['T_cable']
    
    cable_forces.append({
        'from_egg': 0,
        'to_egg': 'ground',
        'z1': z1,
        'z2': z2,
        'd1': bottom_egg['d'],
        'cable_length': cable_length,
        'cable_angle_deg': np.degrees(np.arctan2(dz, dx)),
        'F_horizontal_total': F_horizontal_ground,
        'F_cable_drag': 0,
        'n_cables_per_side': n_cables_per_side,
        'n_cables_total': n_cables_per_side * 2,
        'cables_by_angle': ground_tensions_by_angle,
        'T_cable': T_ground_max,
        'T_cable_kN': T_ground_max / 1000,
        'unstable': False
    })
    
    if any_unstable:
        print("⚠️  WARNING: Cable drag feedback is unstable at some levels!")
    
    return cable_forces


def calculate_alignment_cable_forces(eggs, guy_indices, n_cables_per_side=20):
    """
    Calculate alignment cable forces based on section-local loads.
    
    Alignment cables transfer wind loads from intermediate eggs to guy-anchored
    levels. They do NOT carry the full cumulative tower wind load - that's the
    job of the guy cables.
    
    Physics model:
    - Eggs between guy levels experience wind load
    - This load accumulates toward the nearer guy level (load path)
    - Maximum shear occurs at section midpoint
    - Alignment cables prevent buckling by transferring these section-local loads
    
    Parameters:
        eggs: list of egg dictionaries with wind_force calculated
        guy_indices: list of egg indices where guy cables attach
        n_cables_per_side: number of cable attachment points per quarter
    
    Returns:
        List of cable force dictionaries for each egg-to-egg connection
    """
    from config import CABLE_MATERIAL
    
    n = len(eggs)
    cable_forces = []
    
    # Ensure wind forces are calculated
    for egg in eggs:
        if 'wind_force' not in egg:
            wind_data = wind_force_on_egg(egg)
            egg['wind_force'] = wind_data['wind_force']
    
    # Create sections between guy levels
    sections = []
    for i in range(len(guy_indices)):
        start_idx = guy_indices[i]
        end_idx = guy_indices[i + 1] if i < len(guy_indices) - 1 else n
        sections.append((start_idx, end_idx))
    
    # Angular distribution
    theta_angles = np.linspace(0, np.pi/2, n_cables_per_side)
    cos_thetas = np.cos(theta_angles)
    sum_cos = np.sum(cos_thetas) * 2
    
    # Calculate alignment cable forces for each egg-to-egg connection
    for i in range(n - 1):
        current_egg = eggs[i]
        next_egg = eggs[i + 1]
        
        # Determine which section this connection is in
        section_idx = None
        for s_idx, (start, end) in enumerate(sections):
            if start <= i < end:
                section_idx = s_idx
                break
        
        if section_idx is None:
            section_idx = len(sections) - 1
        
        start_idx, end_idx = sections[section_idx]
        
        # Alignment cables provide lateral stability, not horizontal force equilibrium
        # They must resist:
        # 1. Small differential lateral forces between adjacent eggs
        # 2. Provide pretension for structural integrity
        # 3. Prevent buckling by maintaining egg alignment
        #
        # The GUY CABLES resist the full horizontal wind load.
        # Alignment cables only need to handle local stability requirements.
        #
        # Design load: a fraction of the local egg's wind force
        # This represents the force to maintain alignment during wind events
        ALIGNMENT_FACTOR = 0.1  # 10% of local wind force for stability
        
        F_local = current_egg['wind_force']
        F_alignment = F_local * ALIGNMENT_FACTOR
        
        # Geometry
        z1 = current_egg['z_base'] + current_egg['h'] / 2
        z2 = next_egg['z_base'] + next_egg['h'] / 2
        x1 = current_egg['d'] / 2
        x2 = next_egg['d'] / 2
        
        dz = z2 - z1
        dx = x1 - x2
        cable_length = np.sqrt(dz**2 + dx**2)
        sin_from_vert = dx / cable_length if cable_length > 0 else 0.01
        
        # For alignment cables, tension is based on maintaining structural integrity
        # NOT on resisting the full horizontal force component
        # Use direct force (not divided by sin_from_vert since they're not the load path)
        T_max = F_alignment / max(sin_from_vert, 0.05)  # Reasonable pretension level
        
        # Tensions at all angular positions
        cable_tensions_by_angle = []
        for j, theta in enumerate(theta_angles):
            F_horiz_cable = F_alignment * cos_thetas[j] / sum_cos
            T_cable = T_max * cos_thetas[j] / cos_thetas[0]
            cable_tensions_by_angle.append({
                'theta_deg': np.degrees(theta),
                'F_horizontal': F_horiz_cable,
                'T_cable': T_cable,
                'T_cable_kN': T_cable / 1000
            })
        
        cable_forces.append({
            'from_egg': i,
            'to_egg': i + 1,
            'z1': z1,
            'z2': z2,
            'd1': current_egg['d'],
            'cable_length': cable_length,
            'cable_angle_deg': np.degrees(np.arctan2(dz, dx)),
            'sin_from_vert': sin_from_vert,
            'section_idx': section_idx,
            'F_local': F_local,
            'F_alignment': F_alignment,
            'F_shear': F_alignment,  # For backwards compatibility
            'F_shear_kN': F_alignment / 1000,
            'F_horizontal_total': F_alignment,
            'F_cable_drag': 0,
            'n_cables_per_side': n_cables_per_side,
            'n_cables_total': n_cables_per_side * 2,
            'cables_by_angle': cable_tensions_by_angle,
            'T_cable': T_max,
            'T_cable_kN': T_max / 1000,
            'is_guy_level': i in guy_indices or (i + 1) in guy_indices
        })
    
    return cable_forces


def get_configured_material_key():
    """
    Get the material key for the configured CABLE_MATERIAL.
    """
    from config import CABLE_MATERIAL
    for key, mat in CABLE_MATERIALS.items():
        if mat.name == CABLE_MATERIAL.name:
            return key
    return 'aramid'  # fallback


def size_cable(T_cable, cable_length, material_key=None):
    """
    Size a cable for a given tension and material.
    Uses configured CABLE_MATERIAL if material_key not specified.
    """
    if material_key is None:
        material_key = get_configured_material_key()
    
    mat = CABLE_MATERIALS[material_key]
    
    # Use CableMaterial class methods
    d_cable = mat.cable_diameter(T_cable)
    A_gross = np.pi * (d_cable / 2)**2
    A_eff = A_gross * mat.fill_factor
    sigma_allow = mat.tensile_strength / mat.safety_factor
    mass = A_gross * cable_length * mat.density
    
    result = {
        'material': mat.name,
        'T_cable_kN': T_cable / 1000,
        'sigma_allow_MPa': sigma_allow / 1e6,
        'A_eff_m2': A_eff,
        'A_gross_m2': A_gross,
        'd_cable_m': d_cable,
        'd_cable_mm': d_cable * 1000,
        'cable_length_m': cable_length,
        'mass_kg': mass,
        'mass_tonnes': mass / 1000
    }
    
    # Laced configuration for materials with fiber_diameter
    if mat.fiber_diameter is not None:
        n_fibers = mat.fiber_count(T_cable)
        A_fiber = np.pi * (mat.fiber_diameter / 2)**2 * mat.fill_factor
        mass_per_fiber = A_fiber * cable_length * mat.density / mat.fill_factor
        result['laced'] = {
            'd_fiber_mm': mat.fiber_diameter * 1000,
            'n_fibers': n_fibers,
            'mass_per_fiber_kg': mass_per_fiber,
            'total_laced_mass_kg': n_fibers * mass_per_fiber
        }
    
    return result


def analyze_cable_system(cable_forces, material_key=None):
    """
    Analyze the full cable system for a given material.
    Uses configured CABLE_MATERIAL if material_key not specified.
    Accounts for all cables around the circumference at each angular position.
    Returns strand density (strands per meter of circumference) for each level.
    """
    if material_key is None:
        material_key = get_configured_material_key()
    
    mat = CABLE_MATERIALS[material_key]
    results = []
    total_mass = 0
    total_strands = 0
    
    for cable in cable_forces:
        cable_length = cable['cable_length']
        diameter = cable.get('d1', 20.0)  # diameter at this level
        circumference = np.pi * diameter
        
        # Analyze each angular position
        level_strands = 0
        level_mass = 0
        cables_detail = []
        
        if 'cables_by_angle' in cable:
            for cab in cable['cables_by_angle']:
                sized = size_cable(cab['T_cable'], cable_length, material_key)
                n_strands = sized['laced']['n_fibers'] if 'laced' in sized else 1
                strand_mass = sized['laced']['total_laced_mass_kg'] if 'laced' in sized else sized['mass_kg']
                
                # Multiply by 2 for both sides (mirror symmetry)
                level_strands += n_strands * 2
                level_mass += strand_mass * 2
                
                cables_detail.append({
                    'theta_deg': cab['theta_deg'],
                    'T_cable_kN': cab['T_cable_kN'],
                    'n_strands': n_strands,
                    'mass_kg': strand_mass
                })
        else:
            # Fallback for old format
            sized = size_cable(cable['T_cable'], cable_length, material_key)
            n_strands = sized['laced']['n_fibers'] if 'laced' in sized else 1
            level_strands = n_strands * 2
            level_mass = sized['mass_kg'] * 2
        
        # Calculate strand density (strands per meter of circumference)
        strand_density = level_strands / circumference
        
        results.append({
            'from_egg': cable['from_egg'],
            'to_egg': cable['to_egg'],
            'z1': cable['z1'],
            'diameter': diameter,
            'circumference': circumference,
            'cable_length': cable_length,
            'T_max_kN': cable['T_cable_kN'],
            'n_cables_total': cable.get('n_cables_total', 2),
            'total_strands_at_level': level_strands,
            'strand_density': strand_density,  # strands per meter circumference
            'total_mass_at_level_kg': level_mass,
            'cables_detail': cables_detail
        })
        
        total_strands += level_strands
        total_mass += level_mass
    
    return results, total_mass, total_strands
