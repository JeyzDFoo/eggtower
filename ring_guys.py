"""
Ring structures and ground-anchored guy wire system.

This module implements a hybrid cable system:
1. Small egg-to-egg cables for structural continuity
2. Ring structures at key levels extending beyond egg diameter
3. Ground-anchored guy wires from rings for wind resistance

This design solves the cable drag stability problem by using
steeper cable angles (to ground) rather than nearly-vertical
egg-to-egg cables.

Guy Wire Placement Algorithm:
----------------------------
Starting from the top of the tower, we calculate how many eggs down
the alignment wires can handle before they reach their tension capacity.
When the accumulated horizontal wind force exceeds what the alignment
cables can resist, we need a guy wire level to transfer that force to ground.
"""
import numpy as np
from wind_cables import wind_speed_at_height, cable_wind_drag, wind_force_on_egg
from config import RHO_AIR, CD, CD_CABLE, CABLE_MATERIAL, N_CABLES_PER_LEVEL
from materials import CABLE_MATERIALS


def get_configured_material_key():
    """Get the material key for the configured CABLE_MATERIAL."""
    for key, mat in CABLE_MATERIALS.items():
        if mat.name == CABLE_MATERIAL.name:
            return key
    return 'aramid'  # fallback


def calculate_guy_placement_from_alignment_capacity(eggs, n_alignment_cables=9, 
                                                      arm_radius_factor=10.0,
                                                      cable_material=None, safety_factor=1.5,
                                                      n_strands_per_cable=50):
    """
    Calculate where guy wires are needed based on alignment wire capacity.
    
    GEOMETRY - Alignment cables with horizontal arms (proportional extension):
    ==========================================================================
    
    Cables are attached at the egg/egg interface (tip contact point) via
    horizontal arms that extend proportionally to the egg radius:
    
            ___
           /   \  egg i+1  (radius r2)
           \___/
             │    ← tip contact point
        ─────┼─────  ← arm extends to r2 × arm_radius_factor
             │         cables attach at arm ends
            ___
           /   \  egg i  (radius r1 > r2)
           \___/
             │
        ─────┼─────  ← arm extends to r1 × arm_radius_factor
             │
    
    The MECHANICAL ADVANTAGE of proportional scaling:
    - Tower tapers: r1 > r2 (lower eggs are bigger)
    - Arm radii: r_arm_lower = r1 × factor, r_arm_upper = r2 × factor
    - Horizontal distance: dr = factor × (r1 - r2) = factor × taper
    - The factor AMPLIFIES the cable angle!
    
    Example with factor = 5 (for 5:1 aspect ratio egg):
    - Natural taper: dr = r1 - r2 ≈ 0.5m per egg → θ ≈ 0.5°
    - With factor 5: dr = 5 × 0.5m = 2.5m per egg → θ ≈ 2.6°
    - 5× steeper angle → 5× more horizontal force capacity!
    
    Parameters:
        eggs: list of egg dictionaries from build_tower()
        n_alignment_cables: number of alignment cables around circumference (default 9)
        arm_radius_factor: arm extends to egg_radius × this factor (default 5.0)
        cable_material: material key for cable sizing (uses config default if None)
        safety_factor: additional safety factor for cable capacity
        n_strands_per_cable: number of strands per cable position
    
    Returns:
        Dictionary with guy placement analysis
    """
    if cable_material is None:
        cable_material = get_configured_material_key()
    
    mat = CABLE_MATERIALS[cable_material]
    n = len(eggs)
    
    # Calculate wind force on each egg
    for egg in eggs:
        if 'wind_force' not in egg:
            wind_data = wind_force_on_egg(egg)
            egg['wind_force'] = wind_data['wind_force']
    
    # =========================================================================
    # ALIGNMENT CABLE GEOMETRY WITH HORIZONTAL ARMS
    # =========================================================================
    #
    #         ___
    #        /   \  egg i+1  (smaller, at top)
    #        \___/
    #          │    
    #     ─────●─────  arm at level i+1, radius = r_arm_upper
    #          │         
    #          │     cable runs from arm to arm
    #          │         
    #     ─────●─────  arm at level i, radius = r_arm_lower  
    #          │
    #         ___
    #        /   \  egg i  (larger, at bottom)
    #        \___/
    #
    # Cable geometry:
    #   - Vertical distance: dz = distance between arm levels
    #   - Horizontal distance: dr = r_arm_lower - r_arm_upper (tower tapers)
    #   - Cable length: L = sqrt(dz² + dr²)
    #   - Angle from vertical: sin(θ) = dr / L
    #
    # For horizontal force equilibrium:
    #   F_horizontal = Σ T_cable × sin(θ) × cos(φ)
    #   where φ is the angular position around the circumference
    #
    # =========================================================================
    
    # Angular distribution of cables around 360°
    # e.g., 9 cables at 0°, 40°, 80°, 120°, 160°, 200°, 240°, 280°, 320°
    phi_angles = np.linspace(0, 2*np.pi, n_alignment_cables, endpoint=False)
    
    # For wind in the 0° direction, contribution is cos(φ)
    # Only cables with cos(φ) > 0 contribute to resisting wind (front half)
    # Sum of cos(φ) for cables in front half
    cos_phis = np.cos(phi_angles)
    sum_cos_positive = np.sum(np.maximum(cos_phis, 0))  # Only count positive contributions
    
    # Work from top down, accumulating horizontal shear
    guy_levels = []  # Indices where guy wires are needed
    section_analysis = []  # Detailed analysis of each section
    
    accumulated_shear = 0   # Horizontal shear force accumulated from above (N)
    section_start_idx = n - 1  # Start at top egg
    
    for i in range(n - 2, -1, -1):  # From second-to-top down to bottom
        current_egg = eggs[i]
        upper_egg = eggs[i + 1]
        
        # Add wind force from the egg above this connection
        F_wind_upper = upper_egg['wind_force']
        accumulated_shear += F_wind_upper
        
        # =====================================================================
        # ARM GEOMETRY AT THIS LEVEL - PROPORTIONAL EXTENSION
        # =====================================================================
        # Each arm extends from the tower axis to: egg_radius × arm_radius_factor
        # 
        # Lower egg (current_egg): r_arm = r_egg × factor
        # Upper egg (upper_egg): r_arm = r_egg × factor
        #
        # This AMPLIFIES the tower taper!
        #   dr = r_lower × factor - r_upper × factor = factor × (r_lower - r_upper)
        #
        # Example with factor = 5:
        #   - Tower taper per egg: r_lower - r_upper ≈ 0.27m
        #   - With factor: dr = 5 × 0.27m = 1.35m
        #   - For egg height 55m: angle = arctan(1.35/55) = 1.4°
        #   - This is 5× steeper than without factor!
        #
        # Physical meaning:
        #   - For 5:1 aspect ratio egg with d=11m, arm extends to r_arm = 5 × 5.5m = 27.5m
        #   - This is 5× the egg diameter from center axis
        #   - Just like the user suggested!
        # =====================================================================
        
        # Arm radii at each level (proportional to egg radius)
        r_egg_lower = current_egg['d'] / 2
        r_egg_upper = upper_egg['d'] / 2
        
        r_arm_lower = r_egg_lower * arm_radius_factor
        r_arm_upper = r_egg_upper * arm_radius_factor
        
        # Vertical distance between cable attachment points
        # Cable runs from arm at top of egg i to arm at top of egg i+1
        z_lower = current_egg['z_base'] + current_egg['h']  # top of egg i
        z_upper = upper_egg['z_base'] + upper_egg['h']      # top of egg i+1
        dz = z_upper - z_lower  # vertical distance
        
        # Horizontal distance (due to taper in egg radii)
        # Note: arm_extension cancels out
        dr = r_arm_lower - r_arm_upper  # = r_egg_lower - r_egg_upper
        
        # Cable geometry
        cable_length = np.sqrt(dz**2 + dr**2)
        sin_from_vert = dr / cable_length if cable_length > 0 else 0.001
        angle_from_vert = np.degrees(np.arcsin(min(sin_from_vert, 1.0)))
        
        # =====================================================================
        # FORCE EQUILIBRIUM
        # =====================================================================
        # Horizontal force must be resisted by cable tension:
        #   F_shear = Σ T × sin(θ) × cos(φ)
        #
        # For uniform cable tension:
        #   F_shear = T × sin(θ) × sum_cos_positive
        #   T = F_shear / (sin(θ) × sum_cos_positive)
        # =====================================================================
        
        T_required = accumulated_shear / (sin_from_vert * sum_cos_positive) if sin_from_vert > 0.01 else accumulated_shear * 100
        
        # Cable capacity
        T_capacity = mat.fiber_working_load() * n_strands_per_cable / safety_factor
        
        # Utilization
        utilization = T_required / T_capacity
        can_handle = utilization <= 1.0
        
        section_analysis.append({
            'egg_index': i,
            'z': z_lower,
            'cable_to_egg': i + 1,
            'r_arm_lower': r_arm_lower,
            'r_arm_upper': r_arm_upper,
            'dz': dz,
            'dr': dr,
            'cable_length': cable_length,
            'sin_from_vert': sin_from_vert,
            'angle_from_vert_deg': angle_from_vert,
            'F_wind_egg_kN': F_wind_upper / 1000,
            'accumulated_shear_kN': accumulated_shear / 1000,
            'T_required_kN': T_required / 1000,
            'T_capacity_kN': T_capacity / 1000,
            'utilization': utilization,
            'can_handle': can_handle,
            'section_start': section_start_idx
        })
        
        if not can_handle:
            # =====================================================================
            # GUY WIRE NEEDED
            # =====================================================================
            # When alignment cables reach capacity, we need a guy wire.
            # The guy wire takes the accumulated horizontal shear to ground.
            # This resets the accumulator for the next section.
            # =====================================================================
            guy_level_idx = i + 1
            guy_levels.append({
                'egg_index': guy_level_idx,
                'z': eggs[guy_level_idx]['z_base'] + eggs[guy_level_idx]['h'],
                'accumulated_shear_kN': accumulated_shear / 1000,
                'eggs_in_section': section_start_idx - guy_level_idx + 1,
                'section_top_idx': section_start_idx,
                'section_bottom_idx': guy_level_idx
            })
            
            # Reset accumulator for next section
            accumulated_shear = 0
            section_start_idx = i
    
    # =========================================================================
    # GROUND LEVEL - Final guy to foundation
    # =========================================================================
    # Ground always takes whatever remains
    if accumulated_shear > 0 or len(guy_levels) == 0:
        # Add bottom egg's contribution
        accumulated_shear += eggs[0]['wind_force']
        
        guy_levels.append({
            'egg_index': 0,
            'z': eggs[0]['h'] / 2,
            'accumulated_shear_kN': accumulated_shear / 1000,
            'eggs_in_section': section_start_idx + 1,
            'section_top_idx': section_start_idx,
            'section_bottom_idx': 0
        })
    
    # Reverse to get bottom-to-top order
    guy_levels = guy_levels[::-1]
    section_analysis = section_analysis[::-1]
    
    # Summary statistics
    guy_indices = [g['egg_index'] for g in guy_levels]
    n_guy_levels = len(guy_levels)
    
    # Calculate eggs per section
    eggs_per_section = []
    for g in guy_levels:
        eggs_per_section.append(g['eggs_in_section'])
    
    # Total wind shear on tower
    total_wind_shear = sum(e['wind_force'] for e in eggs) / 1000  # kN
    
    # Get typical cable angle for reporting
    avg_angle = np.mean([s['angle_from_vert_deg'] for s in section_analysis]) if section_analysis else 0
    
    return {
        'guy_levels': guy_levels,
        'guy_indices': guy_indices,
        'n_guy_levels': n_guy_levels,
        'section_analysis': section_analysis,
        'eggs_per_section': eggs_per_section,
        'avg_eggs_per_section': np.mean(eggs_per_section) if eggs_per_section else 0,
        'cable_material': mat.name,
        'T_capacity_kN': mat.fiber_working_load() * n_strands_per_cable / 1000 / safety_factor,
        'n_alignment_cables': n_alignment_cables,
        'n_strands_per_cable': n_strands_per_cable,
        'arm_radius_factor': arm_radius_factor,
        'sum_cos_positive': sum_cos_positive,
        'avg_cable_angle_deg': avg_angle,
        'total_wind_shear_kN': total_wind_shear,
        'physics_model': 'proportional_arm_shear_transfer'
    }


def print_guy_placement_analysis(analysis):
    """Print the guy wire placement analysis based on proportional arm geometry."""
    arm_factor = analysis.get('arm_radius_factor', 1.0)
    print()
    print("═" * 80)
    print("GUY WIRE PLACEMENT ANALYSIS (Proportional Arm Model)")
    print("═" * 80)
    print()
    print("  GEOMETRY (Mechanical Advantage via Proportional Arms):")
    print("  ─────────────────────────────────────────────────────────────────────")
    print(f"  • Arms extend to: egg_radius × {arm_factor:.1f} (proportional scaling)")
    print(f"  • This amplifies cable angle by factor of {arm_factor:.1f}!")
    print("  • Cable runs from arm-to-arm between adjacent levels")
    print("  • Steeper angle → more horizontal force capacity")
    print()
    print("  PHYSICS:")
    print("  ─────────────────────────────────────────────────────────────────────")
    print("  • Cables resist horizontal shear via: T × sin(θ) × cos(φ)")
    print(f"  • θ = {arm_factor:.1f} × tower_taper_angle (amplified!)")
    print("  • φ = angular position around circumference")
    print("  • Guy wires take shear to ground when cables reach capacity")
    print()
    print(f"  CABLE PARAMETERS:")
    print(f"    Material: {analysis['cable_material']}")
    print(f"    Cables around circumference: {analysis['n_alignment_cables']}")
    print(f"    Strands per cable: {analysis['n_strands_per_cable']}")
    print(f"    Cable capacity: {analysis['T_capacity_kN']:.0f} kN")
    print(f"    Average cable angle: {analysis['avg_cable_angle_deg']:.1f}°")
    print(f"    Effective cos sum: {analysis['sum_cos_positive']:.2f}")
    print()
    print(f"  RESULTS:")
    print(f"    Guy wire levels required: {analysis['n_guy_levels']}")
    print(f"    Average eggs per section: {analysis['avg_eggs_per_section']:.1f}")
    print(f"    Total wind shear on tower: {analysis['total_wind_shear_kN']:.0f} kN")
    print()
    
    print(f"  {'Level':<6} {'Egg#':<6} {'Height':<10} {'Shear':<12} {'Eggs in':<10}")
    print(f"  {'':<6} {'':<6} {'(m)':<10} {'(kN)':<12} {'Section':<10}")
    print("  " + "-" * 50)
    
    for i, level in enumerate(analysis['guy_levels']):
        shear = level.get('accumulated_shear_kN', 0)
        print(f"  {i+1:<6} {level['egg_index']:<6} {level['z']:<10.0f} "
              f"{shear:<12.0f} {level['eggs_in_section']:<10}")
    
    print()


def calculate_ring_guy_system(eggs, n_guy_levels=3, anchor_radius=400, 
                               n_guys_per_ring=3, cable_material=None,
                               include_cable_drag=True, max_iterations=10,
                               guy_indices=None):
    """
    Design a ring + ground-anchored guy wire system.
    
    Parameters:
        eggs: list of egg dictionaries from build_tower()
        n_guy_levels: number of levels with guy wires to ground (ignored if guy_indices provided)
        anchor_radius: distance from tower base to ground anchors (m)
        n_guys_per_ring: number of guy wires per ring (distributed around)
        cable_material: material key for sizing cables (uses config default if None)
        include_cable_drag: whether to iteratively include cable self-drag
        max_iterations: max iterations for cable drag convergence
        guy_indices: specific egg indices for guy wire placement (overrides n_guy_levels)
    
    Returns:
        Dictionary with guy wire analysis and ring specifications
    """
    n = len(eggs)
    
    # Use configured material if not specified
    if cable_material is None:
        cable_material = get_configured_material_key()
    
    mat = CABLE_MATERIALS[cable_material]
    
    # If guy_indices provided, use them directly
    if guy_indices is not None:
        guy_indices = sorted(guy_indices)
        n_guy_levels = len(guy_indices)
    else:
        # Optimize guy level placement based on height for good cable angles
        # Avoid placing guys too low (steep angle = inefficient)
        # Minimum height for reasonable cable angle: h_min where θ >= 25°
        # tan(25°) ≈ 0.47, so h_min ≈ anchor_radius * tan(25°) ≈ 0.47 * anchor_radius
        min_height = anchor_radius * 0.47  # ~25° angle minimum
        
        # Find eggs above minimum height
        valid_eggs = [i for i, egg in enumerate(eggs) 
                      if egg['z_base'] + egg['h']/2 >= min_height]
        
        if len(valid_eggs) < n_guy_levels:
            # Fall back to using all eggs if not enough above min height
            valid_eggs = list(range(n))
        
        # Distribute guy levels among valid eggs
        # Weight placement toward upper portion where wind loads are higher
        if n_guy_levels >= len(valid_eggs):
            guy_indices = valid_eggs
        else:
            # Place guys at optimal heights: emphasize upper levels
            # Use square root distribution to weight toward top
            guy_indices = []
            for i in range(n_guy_levels):
                # Fraction from 0 to 1, weighted toward higher values
                frac = (i / (n_guy_levels - 1)) ** 0.7 if n_guy_levels > 1 else 0.5
                idx = int(frac * (len(valid_eggs) - 1))
                egg_idx = valid_eggs[idx]
                if egg_idx not in guy_indices:
                    guy_indices.append(egg_idx)
            
            # Ensure we have exactly n_guy_levels
            while len(guy_indices) < n_guy_levels:
                # Add missing levels
                for idx in valid_eggs:
                    if idx not in guy_indices:
                        guy_indices.append(idx)
                        break
        
        guy_indices = sorted(guy_indices)[:n_guy_levels]
    
    # Calculate wind force on each egg
    for egg in eggs:
        z_center = egg['z_base'] + egg['h'] / 2
        v = wind_speed_at_height(z_center)
        q = 0.5 * RHO_AIR * v**2
        
        # Projected area
        a = egg['h'] / 2
        b = egg['d'] / 2
        A_proj = np.pi * a * b  # Ellipse area
        
        # Combined drag + vortex lift (SRSS)
        Cl = 0.4
        F_drag = q * CD * A_proj
        F_lift = q * Cl * A_proj
        egg['wind_force'] = np.sqrt(F_drag**2 + F_lift**2)
    
    # Iterative calculation for cable drag
    # Each level's guy wires resist: (1) egg forces above, (2) their own cable drag
    # Cable drag at level i is resisted by level i's guys (not accumulated)
    
    for iteration in range(max_iterations if include_cable_drag else 1):
        guy_levels = []
        total_guy_mass = 0
        total_ring_mass = 0
        converged = True
        
        for i, guy_idx in enumerate(guy_indices):
            egg = eggs[guy_idx]
            z_attach = egg['z_base'] + egg['h'] / 2
            r_egg = egg['d'] / 2
            
            # Wind force from eggs this guy level must resist
            if i < len(guy_indices) - 1:
                next_guy_idx = guy_indices[i + 1]
                F_wind_eggs = sum(eggs[j]['wind_force'] for j in range(guy_idx, next_guy_idx))
            else:
                F_wind_eggs = sum(eggs[j]['wind_force'] for j in range(guy_idx, n))
            
            # Guy wire geometry
            cable_length = np.sqrt(z_attach**2 + anchor_radius**2)
            sin_theta = anchor_radius / cable_length
            theta_deg = np.degrees(np.arcsin(sin_theta))
            
            # Stability factor
            k_factor = 0.04 / sin_theta
            stable = k_factor < 1.0
            
            if include_cable_drag and stable:
                # Solve for equilibrium: T = (F_eggs + F_drag(T)) / (sin(θ) * n_guys/2)
                # F_drag(T) = k_factor * T * sin(θ) * n_guys/2  (approximately)
                # T = F_eggs / (sin(θ) * n_guys/2) / (1 - k_factor)
                # This only works if k_factor < 1
                amplification = 1 / (1 - k_factor)
                T_per_guy = F_wind_eggs / (sin_theta * n_guys_per_ring / 2) * amplification
                F_cable_drag = F_wind_eggs * (amplification - 1)
            else:
                T_per_guy = F_wind_eggs / (sin_theta * n_guys_per_ring / 2)
                F_cable_drag = 0
            
            F_wind_total = F_wind_eggs + F_cable_drag
            
            # Size the cables
            n_strands_per_guy = mat.fiber_count(T_per_guy)
            cable_mass_per_guy = n_strands_per_guy * mat.fiber_mass_per_meter() * cable_length
            total_guy_mass_at_level = cable_mass_per_guy * n_guys_per_ring
            total_guy_mass += total_guy_mass_at_level
            
            # Ring structure
            ring_radius = r_egg + 5
            ring_circumference = 2 * np.pi * ring_radius
            ring_mass = ring_circumference * 50
            total_ring_mass += ring_mass
            
            guy_levels.append({
                'egg_index': guy_idx,
                'z_attach': z_attach,
                'r_egg': r_egg,
                'ring_radius': ring_radius,
                'anchor_radius': anchor_radius,
                'cable_length': cable_length,
                'theta_deg': theta_deg,
                'sin_theta': sin_theta,
                'k_factor': k_factor,
                'stable': stable,
                'amplification': amplification if (include_cable_drag and stable) else 1.0,
                'F_wind_eggs_kN': F_wind_eggs / 1000,
                'F_cable_drag_kN': F_cable_drag / 1000,
                'F_wind_total_kN': F_wind_total / 1000,
                'T_per_guy_kN': T_per_guy / 1000,
                'n_guys': n_guys_per_ring,
                'n_strands_per_guy': n_strands_per_guy,
                'total_strands': n_strands_per_guy * n_guys_per_ring,
                'cable_mass_tonnes': total_guy_mass_at_level / 1000,
                'ring_mass_tonnes': ring_mass / 1000
            })
        
        # No iteration needed with analytical solution
        break
    
    # Summary
    return {
        'n_guy_levels': n_guy_levels,
        'anchor_radius': anchor_radius,
        'n_guys_per_ring': n_guys_per_ring,
        'levels': guy_levels,
        'total_guy_cable_mass_tonnes': total_guy_mass / 1000,
        'total_ring_mass_tonnes': total_ring_mass / 1000,
        'total_system_mass_tonnes': (total_guy_mass + total_ring_mass) / 1000,
        'all_stable': all(level['stable'] for level in guy_levels),
        'max_k_factor': max(level['k_factor'] for level in guy_levels),
        'max_tension_MN': max(level['T_per_guy_kN'] for level in guy_levels) / 1000,
        'iterations': iteration + 1 if include_cable_drag else 1,
        'include_cable_drag': include_cable_drag
    }


def print_ring_guy_analysis(analysis):
    """Print formatted analysis of ring + guy wire system."""
    print()
    print("═" * 80)
    print("RING + GROUND-ANCHORED GUY WIRE SYSTEM")
    print("═" * 80)
    print(f"  Configuration:")
    print(f"    Guy wire levels:     {analysis['n_guy_levels']}")
    print(f"    Anchor radius:       {analysis['anchor_radius']} m from tower base")
    print(f"    Guys per ring:       {analysis['n_guys_per_ring']}")
    if analysis.get('include_cable_drag'):
        print(f"    Cable drag:          Included (converged in {analysis.get('iterations', '?')} iterations)")
    print()
    
    stable_str = "✓ ALL STABLE" if analysis['all_stable'] else "⚠️ SOME UNSTABLE"
    print(f"  Stability: {stable_str} (max k = {analysis['max_k_factor']:.2f})")
    print()
    
    print(f"  {'Level':<6} {'Height':<8} {'θ':<8} {'k':<6} {'Amp':<6} {'F_eggs':<10} {'F_drag':<10} {'T/guy':<10} {'Mass':<8}")
    print(f"  {'':6} {'(m)':<8} {'(deg)':<8} {'':<6} {'':<6} {'(kN)':<10} {'(kN)':<10} {'(MN)':<10} {'(t)':<8}")
    print("  " + "-" * 80)
    
    for level in analysis['levels']:
        stable_mark = "✓" if level['stable'] else "✗"
        amp = level.get('amplification', 1.0)
        print(f"  {level['egg_index']:<6} {level['z_attach']:<8.0f} {level['theta_deg']:<8.1f} "
              f"{level['k_factor']:<5.2f}{stable_mark} {amp:<6.2f} {level['F_wind_eggs_kN']:<10.0f} "
              f"{level['F_cable_drag_kN']:<10.0f} "
              f"{level['T_per_guy_kN']/1000:<10.2f} "
              f"{level['cable_mass_tonnes']:<8.1f}")
    
    print("  " + "-" * 76)
    print()
    print(f"  MASS SUMMARY:")
    print(f"    Guy cable mass:      {analysis['total_guy_cable_mass_tonnes']:.0f} tonnes")
    print(f"    Ring structure mass: {analysis['total_ring_mass_tonnes']:.0f} tonnes")
    print(f"    Total system mass:   {analysis['total_system_mass_tonnes']:.0f} tonnes")
    print()
    
    # Calculate total forces
    total_egg_force = sum(l['F_wind_eggs_kN'] for l in analysis['levels'])
    total_cable_drag = sum(l['F_cable_drag_kN'] for l in analysis['levels'])
    print(f"  FORCE SUMMARY:")
    print(f"    Total egg wind force:  {total_egg_force:.0f} kN")
    print(f"    Total cable drag:      {total_cable_drag:.0f} kN ({100*total_cable_drag/total_egg_force:.1f}% of egg force)")
    print(f"    Combined force:        {total_egg_force + total_cable_drag:.0f} kN")


def optimize_anchor_radius(eggs, n_guy_levels=3, target_k=0.2):
    """
    Find optimal anchor radius for a target k factor.
    
    Parameters:
        eggs: list of egg dictionaries
        n_guy_levels: number of guy levels
        target_k: target stability factor (lower = more stable)
    
    Returns:
        Optimal anchor radius and analysis
    """
    # k = 0.04 / sin(θ) = target_k
    # sin(θ) = 0.04 / target_k
    # sin(θ) = r / sqrt(z² + r²)
    
    # For the lowest guy level (most critical due to longest cable)
    # We want sin(θ_min) = 0.04 / target_k
    
    sin_theta_min = 0.04 / target_k
    
    # At height z, for anchor at r:
    # r / sqrt(z² + r²) = sin_theta_min
    # r² = sin²(z² + r²)
    # r² = sin² × z² + sin² × r²
    # r²(1 - sin²) = sin² × z²
    # r² × cos² = sin² × z²
    # r = z × tan(θ)
    
    theta = np.arcsin(sin_theta_min)
    
    # Use the top guy level (highest z) to determine anchor radius
    # This ensures all lower levels have steeper angles (smaller k)
    top_egg_idx = len(eggs) - 1
    z_top = eggs[top_egg_idx]['z_base'] + eggs[top_egg_idx]['h'] / 2
    
    r_optimal = z_top * np.tan(theta)
    
    # Run analysis with this radius
    analysis = calculate_ring_guy_system(eggs, n_guy_levels=n_guy_levels, 
                                         anchor_radius=r_optimal)
    
    return r_optimal, analysis
