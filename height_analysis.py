"""
Analyze maximum tower height with fixed number of guy wire levels.

This script determines how tall the tower can be while keeping the same
number of guy wire levels (9) and maintaining SF ≥ 2.0 on all eggs.
"""
import numpy as np
from config import D_TOP, ASPECT_RATIO, THICKNESS_FACTOR, V_REF
from geometry import build_tower
from structural import analyze_tower
from wind_cables import wind_force_on_egg
from ring_guys import (calculate_guy_placement_from_alignment_capacity, 
                       calculate_ring_guy_system)


def find_min_base_diameter(target_height, min_sf=2.0, d_min=15, d_max=100, tol=0.5):
    """
    Find minimum base diameter that achieves SF >= min_sf for a given height.
    Larger base = lower stress, so we search for the smallest diameter that works.
    """
    # Binary search for minimum viable diameter
    while (d_max - d_min) > tol:
        d_mid = (d_min + d_max) / 2
        
        eggs, ratio = build_tower(
            d_base=d_mid,
            d_top=D_TOP,
            target_height=target_height,
            aspect_ratio=ASPECT_RATIO,
            thickness_factor=THICKNESS_FACTOR
        )
        
        # Calculate wind forces
        for egg in eggs:
            wind_data = wind_force_on_egg(egg)
            egg['wind_force'] = wind_data['wind_force']
        
        # Guy placement
        placement = calculate_guy_placement_from_alignment_capacity(
            eggs, n_alignment_cables=9, arm_radius_factor=10.0, n_strands_per_cable=50
        )
        guy_indices = placement['guy_indices']
        guy_system = calculate_ring_guy_system(eggs, guy_indices=guy_indices)
        
        # Re-analyze with guy forces
        analyzed = analyze_tower(eggs, guy_system=guy_system)
        
        # Find min SF
        min_sf_found = float('inf')
        for i, egg in enumerate(analyzed):
            if i < len(analyzed) - 1:
                if egg['safety_factor'] < min_sf_found:
                    min_sf_found = egg['safety_factor']
        
        if min_sf_found >= min_sf:
            # This diameter works, try smaller
            d_max = d_mid
        else:
            # Need larger diameter
            d_min = d_mid
    
    return d_max  # Return the safe value


def analyze_height(target_height, max_guy_levels=9, verbose=False):
    """
    Analyze if a given tower height can work with limited guy levels.
    
    Returns dict with:
    - feasible: True if SF ≥ 2.0 for all eggs and guy_levels ≤ max
    - n_eggs: number of eggs
    - guy_levels_needed: number of guy levels required by alignment capacity
    - min_sf: minimum safety factor
    - critical_egg: egg index with minimum SF
    - total_mass: total tower mass in tonnes
    """
    # Find minimum base diameter that achieves SF >= 2.0
    d_base = find_min_base_diameter(target_height, min_sf=2.0)
    
    # Build tower with optimized diameter
    eggs, ratio = build_tower(
        d_base=d_base,
        d_top=D_TOP,
        target_height=target_height,
        aspect_ratio=ASPECT_RATIO,
        thickness_factor=THICKNESS_FACTOR
    )
    
    # Calculate wind forces
    for egg in eggs:
        wind_data = wind_force_on_egg(egg)
        egg['wind_force'] = wind_data['wind_force']
    
    # Calculate guy placement based on alignment capacity
    placement = calculate_guy_placement_from_alignment_capacity(
        eggs,
        n_alignment_cables=9,
        arm_radius_factor=10.0,
        n_strands_per_cable=50
    )
    
    guy_levels_needed = placement['n_guy_levels']
    guy_indices = placement['guy_indices']
    
    # Calculate guy system with exact placement
    guy_system = calculate_ring_guy_system(eggs, guy_indices=guy_indices)
    
    # Re-analyze tower with guy forces
    analyzed = analyze_tower(eggs, guy_system=guy_system)
    
    # Find minimum SF and check feasibility
    min_sf = float('inf')
    critical_egg = 0
    all_pass = True
    
    for i, egg in enumerate(analyzed):
        if i < len(analyzed) - 1:  # Skip top egg (no load above)
            sf = egg['safety_factor']
            if sf < min_sf:
                min_sf = sf
                critical_egg = i
            if sf < 2.0:
                all_pass = False
    
    # Calculate masses
    shell_mass = sum(e['mass'] for e in analyzed) / 1000
    guy_mass = guy_system['total_system_mass_tonnes']
    
    feasible = all_pass and (guy_levels_needed <= max_guy_levels)
    
    result = {
        'height': target_height,
        'n_eggs': len(eggs),
        'd_base': d_base,
        'guy_levels_needed': guy_levels_needed,
        'min_sf': min_sf,
        'critical_egg': critical_egg,
        'shell_mass': shell_mass,
        'guy_mass': guy_mass,
        'total_mass': shell_mass + guy_mass,
        'feasible': feasible,
        'sf_ok': all_pass,
        'guy_ok': guy_levels_needed <= max_guy_levels
    }
    
    if verbose:
        status = "✓ FEASIBLE" if feasible else "✗ INFEASIBLE"
        print(f"Height {target_height:5.0f}m: {len(eggs):2d} eggs, "
              f"{guy_levels_needed:2d} guy levels, "
              f"min SF={min_sf:.2f}, "
              f"mass={shell_mass + guy_mass:.0f}t  {status}")
    
    return result


def find_max_height(max_guy_levels=9, h_min=1600, h_max=3000, tol=25):
    """
    Binary search to find maximum feasible tower height.
    """
    print(f"\n{'═'*80}")
    print(f"MAXIMUM HEIGHT ANALYSIS (with {max_guy_levels} guy levels)")
    print(f"{'═'*80}")
    print(f"Constraints: SF ≥ 2.0, guy levels ≤ {max_guy_levels}")
    print(f"Search range: {h_min}m to {h_max}m")
    print(f"{'═'*80}\n")
    
    # First check baseline
    print("Baseline check at 1600m:")
    baseline = analyze_height(1600, max_guy_levels, verbose=True)
    if not baseline['feasible']:
        print("ERROR: Even baseline 1600m is not feasible!")
        return baseline
    
    print("\nSearching for maximum height...")
    print("-" * 70)
    
    # Binary search for max feasible height
    h_low = h_min
    h_high = h_max
    best_result = baseline
    
    while (h_high - h_low) > tol:
        h_mid = (h_low + h_high) / 2
        
        result = analyze_height(h_mid, max_guy_levels, verbose=True)
        
        if result['feasible']:
            h_low = h_mid
            best_result = result
        else:
            h_high = h_mid
    
    print("-" * 70)
    print(f"\n{'═'*80}")
    print(f"MAXIMUM FEASIBLE HEIGHT: {best_result['height']:.0f}m")
    print(f"{'═'*80}")
    print(f"  Number of eggs: {best_result['n_eggs']}")
    print(f"  Base diameter: {best_result['d_base']:.1f}m")
    print(f"  Guy levels: {best_result['guy_levels_needed']}")
    print(f"  Minimum SF: {best_result['min_sf']:.2f} (egg #{best_result['critical_egg']})")
    print(f"  Total mass: {best_result['total_mass']:.0f} tonnes")
    print(f"{'═'*80}")
    
    return best_result


def sweep_heights(h_values, max_guy_levels=9):
    """
    Analyze multiple tower heights and show comparison.
    """
    print(f"\n{'═'*90}")
    print(f"TOWER HEIGHT COMPARISON (with {max_guy_levels} guy levels max)")
    print(f"{'═'*90}")
    
    print(f"\n{'Height':>8} {'Eggs':>6} {'D_base':>8} {'Guy Lvls':>10} {'Min SF':>10} "
          f"{'Mass(t)':>10} {'Status':>15}")
    print("-" * 90)
    
    results = []
    for h in h_values:
        r = analyze_height(h, max_guy_levels)
        results.append(r)
        
        sf_status = "SF OK" if r['sf_ok'] else f"SF FAIL"
        guy_status = "Guy OK" if r['guy_ok'] else f"Need {r['guy_levels_needed']}"
        status = "✓ OK" if r['feasible'] else f"✗ {sf_status if not r['sf_ok'] else guy_status}"
        
        print(f"{h:>8.0f} {r['n_eggs']:>6d} {r['d_base']:>8.1f} {r['guy_levels_needed']:>10d} "
              f"{r['min_sf']:>10.2f} {r['total_mass']:>10.0f} {status:>15}")
    
    return results


if __name__ == "__main__":
    # First sweep a range of heights
    heights = [1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000]
    sweep_heights(heights, max_guy_levels=9)
    
    # Then find exact maximum
    find_max_height(max_guy_levels=9, h_min=1600, h_max=3000, tol=25)
