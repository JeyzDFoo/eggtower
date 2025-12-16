"""
Cost analysis for egg tower construction.

Material cost estimates (2024 USD):
- BFRP (Basalt Fiber Reinforced Polymer): $5-15/kg depending on quality
- Aramid (Kevlar 49) fiber: $20-40/kg for raw fiber
- Steel (for comparison): $1-2/kg

Using conservative mid-range estimates.
"""
from height_analysis import analyze_height, find_min_base_diameter
from config import D_TOP, ASPECT_RATIO, THICKNESS_FACTOR
from geometry import build_tower
from structural import analyze_tower
from wind_cables import wind_force_on_egg
from ring_guys import calculate_guy_placement_from_alignment_capacity, calculate_ring_guy_system


# Material costs in USD per kg
COST_BFRP_PER_KG = 10.0       # Mid-range BFRP
COST_ARAMID_PER_KG = 30.0     # Kevlar 49 fiber
COST_STEEL_PER_KG = 1.50      # Structural steel (for reference)

# Installation multipliers (labor, equipment, logistics)
INSTALL_MULTIPLIER_SHELL = 3.0    # Complex shell fabrication and lifting
INSTALL_MULTIPLIER_CABLE = 2.0    # Cable installation is simpler


def estimate_costs(height, d_base=None, verbose=True):
    """
    Estimate construction costs for a tower configuration.
    
    Returns cost breakdown in USD.
    """
    # Find optimal diameter if not provided
    if d_base is None:
        d_base = find_min_base_diameter(height, min_sf=2.0)
    
    # Build tower
    eggs, ratio = build_tower(
        d_base=d_base,
        d_top=D_TOP,
        target_height=height,
        aspect_ratio=ASPECT_RATIO,
        thickness_factor=THICKNESS_FACTOR
    )
    
    # Calculate wind forces and guy system
    for egg in eggs:
        wind_data = wind_force_on_egg(egg)
        egg['wind_force'] = wind_data['wind_force']
    
    placement = calculate_guy_placement_from_alignment_capacity(
        eggs, n_alignment_cables=9, arm_radius_factor=10.0, n_strands_per_cable=50
    )
    guy_indices = placement['guy_indices']
    guy_system = calculate_ring_guy_system(eggs, guy_indices=guy_indices)
    
    # Re-analyze with guy forces
    analyzed = analyze_tower(eggs, guy_system=guy_system)
    
    # Mass breakdown (kg)
    shell_mass_kg = sum(e['mass'] for e in analyzed)
    guy_cable_mass_kg = guy_system['total_system_mass_tonnes'] * 1000
    
    # Estimate alignment cable mass (roughly 20% of guy cable mass based on prior runs)
    # This is approximate - could calculate more precisely
    alignment_cable_mass_kg = guy_cable_mass_kg * 0.3
    
    total_cable_mass_kg = guy_cable_mass_kg + alignment_cable_mass_kg
    
    # Material costs
    shell_material_cost = shell_mass_kg * COST_BFRP_PER_KG
    cable_material_cost = total_cable_mass_kg * COST_ARAMID_PER_KG
    
    # Installed costs (includes fabrication, transport, installation)
    shell_installed_cost = shell_material_cost * INSTALL_MULTIPLIER_SHELL
    cable_installed_cost = cable_material_cost * INSTALL_MULTIPLIER_CABLE
    
    total_material_cost = shell_material_cost + cable_material_cost
    total_installed_cost = shell_installed_cost + cable_installed_cost
    
    # Cost per meter of height
    cost_per_meter = total_installed_cost / height
    
    result = {
        'height': height,
        'd_base': d_base,
        'n_eggs': len(eggs),
        'guy_levels': len(guy_indices),
        
        # Masses (tonnes)
        'shell_mass_t': shell_mass_kg / 1000,
        'cable_mass_t': total_cable_mass_kg / 1000,
        'total_mass_t': (shell_mass_kg + total_cable_mass_kg) / 1000,
        
        # Costs (USD millions)
        'shell_material_M': shell_material_cost / 1e6,
        'cable_material_M': cable_material_cost / 1e6,
        'total_material_M': total_material_cost / 1e6,
        
        'shell_installed_M': shell_installed_cost / 1e6,
        'cable_installed_M': cable_installed_cost / 1e6,
        'total_installed_M': total_installed_cost / 1e6,
        
        'cost_per_meter': cost_per_meter,
    }
    
    if verbose:
        print(f"\n{'═'*70}")
        print(f"COST ANALYSIS: {height}m TOWER")
        print(f"{'═'*70}")
        print(f"\n  CONFIGURATION:")
        print(f"    Height: {height}m")
        print(f"    Base diameter: {d_base:.1f}m")
        print(f"    Number of eggs: {len(eggs)}")
        print(f"    Guy levels: {len(guy_indices)}")
        
        print(f"\n  MASS BREAKDOWN:")
        print(f"    Shell (BFRP):     {shell_mass_kg/1000:>10,.0f} tonnes")
        print(f"    Cables (Aramid):  {total_cable_mass_kg/1000:>10,.0f} tonnes")
        print(f"    TOTAL:            {(shell_mass_kg + total_cable_mass_kg)/1000:>10,.0f} tonnes")
        
        print(f"\n  MATERIAL COSTS (USD):")
        print(f"    Shell @ ${COST_BFRP_PER_KG}/kg:    ${shell_material_cost/1e6:>10,.1f} M")
        print(f"    Cables @ ${COST_ARAMID_PER_KG}/kg:  ${cable_material_cost/1e6:>10,.1f} M")
        print(f"    TOTAL MATERIAL:       ${total_material_cost/1e6:>10,.1f} M")
        
        print(f"\n  INSTALLED COSTS (incl. fabrication, transport, installation):")
        print(f"    Shell (×{INSTALL_MULTIPLIER_SHELL}):          ${shell_installed_cost/1e6:>10,.1f} M")
        print(f"    Cables (×{INSTALL_MULTIPLIER_CABLE}):         ${cable_installed_cost/1e6:>10,.1f} M")
        print(f"    TOTAL INSTALLED:      ${total_installed_cost/1e6:>10,.1f} M")
        
        print(f"\n  COST EFFICIENCY:")
        print(f"    Cost per meter:       ${cost_per_meter:>10,.0f} /m")
        print(f"{'═'*70}")
    
    return result


def compare_heights(heights):
    """
    Compare costs across different tower heights.
    """
    print(f"\n{'═'*100}")
    print(f"TOWER COST COMPARISON")
    print(f"{'═'*100}")
    print(f"  Material costs: BFRP @ ${COST_BFRP_PER_KG}/kg, Aramid @ ${COST_ARAMID_PER_KG}/kg")
    print(f"  Install multipliers: Shell ×{INSTALL_MULTIPLIER_SHELL}, Cable ×{INSTALL_MULTIPLIER_CABLE}")
    print(f"{'═'*100}\n")
    
    print(f"{'Height':>8} {'Eggs':>6} {'D_base':>8} {'Mass(t)':>12} {'Material$':>12} "
          f"{'Installed$':>12} {'$/meter':>12}")
    print("-" * 80)
    
    results = []
    for h in heights:
        r = estimate_costs(h, verbose=False)
        results.append(r)
        
        print(f"{h:>8.0f} {r['n_eggs']:>6d} {r['d_base']:>8.1f} "
              f"{r['total_mass_t']:>12,.0f} "
              f"${r['total_material_M']:>10,.1f}M "
              f"${r['total_installed_M']:>10,.1f}M "
              f"${r['cost_per_meter']:>10,.0f}")
    
    print("-" * 80)
    
    # Show scaling factor relative to 1600m
    base = results[0]
    print(f"\n  SCALING RELATIVE TO {base['height']:.0f}m:")
    for r in results:
        height_ratio = r['height'] / base['height']
        cost_ratio = r['total_installed_M'] / base['total_installed_M']
        print(f"    {r['height']:.0f}m: {height_ratio:.2f}× height → "
              f"{cost_ratio:.1f}× cost (${r['total_installed_M']:.0f}M vs ${base['total_installed_M']:.0f}M)")
    
    return results


if __name__ == "__main__":
    # Compare key heights
    heights = [1600, 1800, 2000, 2400, 2800, 3000]
    compare_heights(heights)
    
    # Detailed breakdown for extremes
    print("\n" + "="*80)
    print("DETAILED BREAKDOWNS")
    print("="*80)
    
    estimate_costs(1600)
    estimate_costs(3000)
