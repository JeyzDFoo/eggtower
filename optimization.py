"""
Optimization functions for egg tower design.
"""
from config import (TOTAL_HEIGHT, D_BASE, D_TOP, ASPECT_RATIO, 
                    N_CABLES_PER_LEVEL, V_REF)
from geometry import build_tower
from structural import analyze_tower
from wind_cables import calculate_alignment_cable_forces, analyze_cable_system, get_guy_indices, wind_force_on_egg
from ring_guys import calculate_ring_guy_system


def evaluate_diameter(d_base):
    """
    Evaluate a single base diameter and return key metrics.
    """
    eggs, ratio = build_tower(d_base=d_base, d_top=D_TOP, target_height=TOTAL_HEIGHT, aspect_ratio=ASPECT_RATIO)
    analyzed = analyze_tower(eggs)
    
    # Calculate wind forces on each egg
    for egg in analyzed:
        wind_data = wind_force_on_egg(egg)
        egg['wind_force'] = wind_data['wind_force']
    
    # Guy system analysis (determines guy_indices)
    guy_system = calculate_ring_guy_system(analyzed)
    guy_indices = [level['egg_index'] for level in guy_system['levels']]
    guy_cable_mass = guy_system['total_system_mass_tonnes']
    
    # Alignment cables (section-local loads, not full cumulative)
    cable_forces = calculate_alignment_cable_forces(analyzed, guy_indices, n_cables_per_side=N_CABLES_PER_LEVEL)
    cable_results, alignment_cable_mass, total_strands = analyze_cable_system(cable_forces)
    
    n_eggs = len(eggs)
    total_height = analyzed[-1]['z_base'] + analyzed[-1]['h']
    shell_mass = sum(e['mass'] for e in analyzed) / 1000  # tonnes
    max_tension = max(c['T_cable'] for c in cable_forces) / 1e6  # MN
    max_strand_density = max(r['strand_density'] for r in cable_results) if cable_results else 0
    alignment_mass_tonnes = alignment_cable_mass / 1000
    total_mass = shell_mass + alignment_mass_tonnes + guy_cable_mass
    
    return {
        'd_base': d_base,
        'n_eggs': n_eggs,
        'height': total_height,
        'shell_mass': shell_mass,
        'max_tension': max_tension,
        'max_strand_density': max_strand_density,
        'alignment_cable_mass': alignment_mass_tonnes,
        'guy_cable_mass': guy_cable_mass,
        'total_mass': total_mass,
        'eggs': analyzed,
        'cable_forces': cable_forces,
        'guy_indices': guy_indices
    }


def binary_search_optimal_diameter(objective='total_mass', d_min=20, d_max=300, tol=0.5):
    """
    Binary search to find optimal base diameter.
    
    Parameters:
        objective: 'total_mass', 'max_ropes', or 'max_tension'
        d_min, d_max: search bounds (meters)
        tol: tolerance for convergence (meters)
    
    Returns:
        Optimal result dictionary
    """
    print(f"\n{'═'*100}")
    print(f"BINARY SEARCH FOR OPTIMAL BASE DIAMETER")
    print(f"{'═'*100}")
    print(f"Objective: Minimize {objective}")
    print(f"Search range: {d_min}m to {d_max}m | Tolerance: {tol}m")
    print(f"{'═'*100}")
    
    iteration = 0
    history = []
    
    print(f"\n{'Iter':<6} {'D_low':<10} {'D_mid':<10} {'D_high':<10} {'Objective':<16} {'# Eggs':<8} {'Total Mass':<14}")
    print(f"{'─'*80}")
    
    while (d_max - d_min) > tol:
        iteration += 1
        
        d_mid = (d_min + d_max) / 2
        d_left = (d_min + d_mid) / 2
        d_right = (d_mid + d_max) / 2
        
        r_left = evaluate_diameter(d_left)
        r_mid = evaluate_diameter(d_mid)
        r_right = evaluate_diameter(d_right)
        
        obj_left = r_left[objective]
        obj_mid = r_mid[objective]
        obj_right = r_right[objective]
        
        if obj_left <= obj_mid and obj_left <= obj_right:
            d_max = d_mid
            best = r_left
        elif obj_right <= obj_mid and obj_right <= obj_left:
            d_min = d_mid
            best = r_right
        else:
            d_min = d_left
            d_max = d_right
            best = r_mid
        
        history.append(best)
        print(f"{iteration:<6} {d_min:<10.1f} {(d_min+d_max)/2:<10.1f} {d_max:<10.1f} "
              f"{best[objective]:<16,.2f} {best['n_eggs']:<8} {best['total_mass']:<14,.0f}")
    
    d_optimal = (d_min + d_max) / 2
    optimal = evaluate_diameter(d_optimal)
    
    print(f"{'─'*80}")
    print(f"\n  CONVERGED after {iteration} iterations")
    print(f"  {'─'*60}")
    print(f"  Optimal D_base:     {optimal['d_base']:.1f} m")
    print(f"  Number of eggs:     {optimal['n_eggs']}")
    print(f"  Total height:       {optimal['height']:.1f} m")
    print(f"  Shell mass:         {optimal['shell_mass']:,.0f} tonnes")
    print(f"  Alignment cables:   {optimal['alignment_cable_mass']:,.0f} tonnes")
    print(f"  Guy cables:         {optimal['guy_cable_mass']:,.0f} tonnes")
    print(f"  TOTAL MASS:         {optimal['total_mass']:,.0f} tonnes")
    print(f"  Max cable tension:  {optimal['max_tension']:.2f} MN")
    print(f"  Max strand density: {optimal['max_strand_density']:.1f} strands/m")
    
    return optimal


def run_diameter_sweep():
    """
    Sweep through different base diameters to find optimal configuration.
    """
    print(f"\n{'═'*100}")
    print("BASE DIAMETER PARAMETER SWEEP")
    print(f"{'═'*100}")
    print(f"Target height: {TOTAL_HEIGHT}m | Top diameter: {D_TOP}m | Aspect ratio: {ASPECT_RATIO}")
    print(f"Wind: {V_REF} m/s at 10m height | Material: Dyneema SK78 (50mm laced)")
    print(f"{'═'*100}")
    
    d_base_values = [30, 50, 75, 100, 125, 150, 200]
    
    print(f"\n{'D_base':<10} {'# Eggs':<10} {'Height':<12} {'Shell Mass':<14} {'Max Tension':<14} "
          f"{'Strand Den.':<14} {'Cable Mass':<14} {'Total Mass':<14}")
    print(f"{'(m)':<10} {'':<10} {'(m)':<12} {'(tonnes)':<14} {'(MN)':<14} "
          f"{'(per m)':<14} {'(tonnes)':<14} {'(tonnes)':<14}")
    print(f"{'─'*110}")
    
    results = []
    
    for d_base in d_base_values:
        result = evaluate_diameter(d_base)
        print(f"{d_base:<10} {result['n_eggs']:<10} {result['height']:<12.1f} {result['shell_mass']:<14,.0f} "
              f"{result['max_tension']:<14.2f} {result['max_strand_density']:<14.1f} {result['cable_mass']:<14,.0f} {result['total_mass']:<14,.0f}")
        results.append(result)
    
    print(f"{'─'*110}")
    
    min_tension = min(results, key=lambda x: x['max_tension'])
    min_mass = min(results, key=lambda x: x['total_mass'])
    min_density = min(results, key=lambda x: x['max_strand_density'])
    
    print(f"\n  OPTIMAL CONFIGURATIONS:")
    print(f"  {'─'*60}")
    print(f"  Minimum cable tension:  D_base = {min_tension['d_base']}m → {min_tension['max_tension']:.2f} MN")
    print(f"  Minimum total mass:     D_base = {min_mass['d_base']}m → {min_mass['total_mass']:,.0f} tonnes")
    print(f"  Minimum strand density: D_base = {min_density['d_base']}m → {min_density['max_strand_density']:.1f} strands/m")
    
    return results
