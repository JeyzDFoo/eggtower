"""
Output and printing functions for tower analysis.
"""
import numpy as np
from config import TOTAL_HEIGHT, D_BASE, D_TOP, V_REF, N_CABLES_PER_LEVEL, CL_VORTEX, CD, STROUHAL
from materials import CABLE_MATERIALS
from wind_cables import analyze_cable_system, analyze_vortex_shedding


def print_tower_summary(analyzed_eggs):
    """
    Print a summary table of the tower including buckling analysis.
    """
    print("=" * 120)
    print(f"TAPERED EGG TOWER (Height: {TOTAL_HEIGHT}m, Base Diameter: {D_BASE}m, Top Diameter: {D_TOP}m)")
    print("=" * 120)
    
    # Material stress table
    print(f"\n{'─'*60} MATERIAL STRESS {'─'*60}")
    print(f"{'i':<4} {'z_base':<10} {'z_top':<10} {'d(m)':<8} {'h(m)':<8} "
          f"{'Mass(t)':<10} {'Stress(MPa)':<12} {'SF Material':<12}")
    print("-" * 120)
    
    total_mass = 0
    total_volume = 0
    for egg in analyzed_eggs:
        total_mass += egg['mass']
        total_volume += egg['volume']
        
        print(f"{egg['index']:<4} {egg['z_base']:<10.1f} {egg['z_top']:<10.1f} "
              f"{egg['d']:<8.2f} {egg['h']:<8.2f} "
              f"{egg['mass']/1e3:<10.2f} {egg['stress_MPa']:<12.2f} "
              f"{egg['safety_factor']:<12.1f}")
    
    # Buckling analysis table
    print(f"\n{'─'*60} BUCKLING ANALYSIS {'─'*60}")
    print(f"{'i':<4} {'R_crit(m)':<10} {'σ_shell(MPa)':<14} {'σ_euler(MPa)':<14} "
          f"{'σ_actual(MPa)':<14} {'SF Shell':<10} {'SF Euler':<10} {'Governing':<10}")
    print("-" * 120)
    
    min_sf_buckling = np.inf
    critical_egg = None
    
    for egg in analyzed_eggs:
        sf_shell = egg.get('SF_shell_buckling', np.inf)
        sf_euler = egg.get('SF_euler_buckling', np.inf)
        
        if egg['SF_buckling_combined'] < min_sf_buckling:
            min_sf_buckling = egg['SF_buckling_combined']
            critical_egg = egg
        
        print(f"{egg['index']:<4} {egg['R_critical']:<10.2f} "
              f"{egg['sigma_cr_shell_MPa']:<14.1f} {egg['sigma_cr_euler_MPa']:<14.1f} "
              f"{egg['sigma_actual_MPa']:<14.2f} "
              f"{sf_shell:<10.1f} {sf_euler:<10.1f} {egg['governing_mode']:<10}")
    
    print("-" * 120)
    print(f"\nTower Statistics:")
    print(f"  Number of eggs: {len(analyzed_eggs)}")
    print(f"  Total height: {analyzed_eggs[-1]['z_top']:.1f} m")
    print(f"  Total mass: {total_mass/1e3:.1f} tonnes")
    print(f"  Total shell volume: {total_volume:.1f} m³")
    print(f"\n  MATERIAL STRENGTH:")
    print(f"    Bottom egg stress: {analyzed_eggs[0]['stress_MPa']:.2f} MPa")
    print(f"    Safety factor (material): {analyzed_eggs[0]['safety_factor']:.1f}")
    print(f"\n  BUCKLING:")
    print(f"    Critical egg: #{critical_egg['index']} (z = {critical_egg['z_base']:.0f}m)")
    print(f"    Governing mode: {critical_egg['governing_mode']} buckling")
    print(f"    Min safety factor (buckling): {min_sf_buckling:.1f}")
    print(f"\n  Diameter ratio: {analyzed_eggs[1]['d']/analyzed_eggs[0]['d']:.4f}" if len(analyzed_eggs) > 1 else "")


def print_alignment_cable_forces(cable_forces, guy_indices):
    """
    Print alignment cable forces with section-based load distribution.
    
    Alignment cables transfer wind loads from intermediate eggs to guy-anchored
    levels. They carry section-local shear, not full cumulative tower load.
    """
    n_cables_per_side = cable_forces[0].get('n_cables_per_side', 1) if cable_forces else 1
    n_cables_total = cable_forces[0].get('n_cables_total', 2) if cable_forces else 2
    
    print(f"\n{'─'*60} ALIGNMENT CABLE FORCES (Wind: {V_REF} m/s at 10m) {'─'*40}")
    print(f"  Guy levels at eggs: {guy_indices}")
    print(f"  Alignment cables transfer section-local loads to guy-anchored levels")
    print(f"  Distribution: {n_cables_per_side} cables/quarter × 4 = {n_cables_total*2} points")
    print(f"")
    print(f"{'Segment':<15} {'From z(m)':<10} {'To z(m)':<10} {'Length(m)':<10} "
          f"{'F_shear(kN)':<12} {'T_max(MN)':<12} {'Section':<10} {'Guy?':<6}")
    print("-" * 100)
    
    for cable in cable_forces:
        is_guy = "✓" if cable.get('is_guy_level', False) else ""
        section = cable.get('section_idx', 0)
        F_shear = cable.get('F_shear', cable.get('F_horizontal_total', 0))
        print(f"egg {cable['from_egg']} → egg {cable['to_egg']:<4} "
              f"{cable['z1']:<10.1f} {cable['z2']:<10.1f} "
              f"{cable['cable_length']:<10.1f} {F_shear/1000:<12.1f} "
              f"{cable['T_cable_kN']/1000:<12.2f} {section:<10} {is_guy:<6}")
    
    max_cable = max(cable_forces, key=lambda x: x['T_cable'])
    total_shear = sum(c.get('F_shear', 0) for c in cable_forces)
    print("-" * 100)
    print(f"\n  Max alignment cable tension: {max_cable['T_cable_kN']/1000:.2f} MN "
          f"(segment: egg {max_cable['from_egg']} → {max_cable['to_egg']})")
    print(f"  Note: Shear forces are section-local (between guy levels), not cumulative")
    print(f"  Note: Guy cables (not alignment cables) resist full wind load to ground")


def print_cable_forces(cable_forces):
    """
    Print summary of cable forces.
    """
    n_cables_per_side = cable_forces[0].get('n_cables_per_side', 1) if cable_forces else 1
    n_cables_total = cable_forces[0].get('n_cables_total', 2) if cable_forces else 2
    
    # Calculate total cable drag
    total_cable_drag = sum(cf.get('F_cable_drag', 0) for cf in cable_forces)
    total_egg_force = sum(cf.get('F_horizontal_total', 0) - cf.get('F_cable_drag', 0) for cf in cable_forces[:1])  # Ground cable has full egg forces
    
    print(f"\n{'─'*60} CABLE FORCES (Wind: {V_REF} m/s at 10m) {'─'*60}")
    print(f"  Distribution: {n_cables_per_side} cables/quarter × 4 = {n_cables_total*2} attachment points around circumference")
    print(f"  Angular range: 0° (windward) to 90° (perpendicular) with cos(θ) load distribution")
    if total_cable_drag > 0:
        print(f"  Cable wind drag: {total_cable_drag/1e6:.2f} MN total (angular factor 2/π ≈ 0.64 for effective area)")
    print(f"")
    print(f"{'Segment':<15} {'From z(m)':<10} {'To z(m)':<10} {'Length(m)':<10} "
          f"{'F_total(kN)':<12} {'F_cable(kN)':<12} {'T_max(MN)':<12}")
    print("-" * 100)
    
    for cable in cable_forces:
        to_str = 'ground' if cable['to_egg'] == 'ground' else f"egg {cable['to_egg']}"
        # Get cable drag for this segment
        F_cable_drag = cable.get('F_cable_drag', 0)
        print(f"egg {cable['from_egg']} → {to_str:<8} {cable['z1']:<10.1f} {cable['z2']:<10.1f} "
              f"{cable['cable_length']:<10.1f} {cable['F_horizontal_total']/1000:<12.1f} "
              f"{F_cable_drag/1000:<12.1f} {cable['T_cable_kN']/1000:<12.2f}")
    
    max_cable = max(cable_forces, key=lambda x: x['T_cable'])
    print("-" * 100)
    print(f"\n  Maximum cable tension: {max_cable['T_cable_kN']/1000:.2f} MN "
          f"(segment: egg {max_cable['from_egg']} → {max_cable['to_egg']})")
    print(f"  Note: T_max is the tension in the windward cable (at θ=0°)")
    print(f"  Note: F_cable = wind drag on cables at this level (added to F_total)")
    print(f"  Note: Forces include vortex shedding (Cd={CD}, Cl={CL_VORTEX}, combined via SRSS)")


def print_vortex_shedding_analysis(eggs):
    """
    Print vortex shedding analysis for the tower.
    """
    analysis = analyze_vortex_shedding(eggs)
    
    print(f"\n{'═'*100}")
    print(f"  VORTEX SHEDDING ANALYSIS")
    print(f"{'═'*100}")
    print(f"  Strouhal number (St): {STROUHAL}")
    print(f"  Lift coefficient (Cl): {CL_VORTEX}")
    print(f"  Estimated tower natural frequency: {analysis['f_natural']:.3f} Hz")
    print(f"  Period: {1/analysis['f_natural']:.1f} seconds")
    print(f"{'═'*100}")
    
    print(f"\n  {'Egg':<6} {'Height':<10} {'Diameter':<10} {'Wind':<10} {'Shed Freq':<12} {'f_s/f_n':<10} {'Lift Force':<12} {'Risk':<10}")
    print(f"  {'':<6} {'(m)':<10} {'(m)':<10} {'(m/s)':<10} {'(Hz)':<12} {'':<10} {'(kN)':<12} {'':<10}")
    print(f"  {'─'*90}")
    
    for r in analysis['eggs']:
        lock = r['lock_in']
        risk_color = '🔴' if lock['risk_level'] == 'CRITICAL' else ('🟡' if lock['risk_level'] == 'HIGH' else '🟢')
        print(f"  {r['egg_index']:<6} {r['z']:<10.1f} {r['d']:<10.2f} {r['wind_speed']:<10.1f} "
              f"{r['shedding_freq']:<12.3f} {lock['ratio']:<10.2f} {r['F_lift_kN']:<12.1f} {risk_color} {lock['risk_level']:<8}")
    
    print(f"  {'─'*90}")
    print(f"")
    print(f"  SUMMARY:")
    print(f"  {'─'*50}")
    print(f"  Tower natural frequency:    {analysis['f_natural']:.4f} Hz")
    print(f"  Max vortex lift force:      {analysis['max_lift_force_kN']:,.0f} kN")
    
    if analysis['has_lock_in_risk']:
        print(f"")
        print(f"  ⚠️  LOCK-IN RESONANCE DETECTED at eggs: {analysis['critical_eggs']}")
        print(f"      This is CRITICAL - oscillations will be amplified!")
        print(f"      Mitigation options:")
        print(f"        • Add helical strakes to break vortex coherence")
        print(f"        • Install tuned mass dampers")
        print(f"        • Modify geometry to shift shedding frequency")
    else:
        print(f"")
        print(f"  ✓ No lock-in resonance detected")
        print(f"    Shedding frequencies are sufficiently separated from natural frequency")
    
    print(f"")
    print(f"  FORCE COMBINATION:")
    print(f"  {'─'*50}")
    print(f"  Cable forces use SRSS combination: F = √(F_drag² + F_lift²)")
    print(f"  This is conservative - assumes peak lift coincides with full drag")
    print(f"  Lift/Drag ratio = Cl/Cd = {CL_VORTEX}/{CD} = {CL_VORTEX/CD:.2f}")
    print(f"  Force amplification factor: √(1 + (Cl/Cd)²) = {np.sqrt(1 + (CL_VORTEX/CD)**2):.3f}")


def print_cable_material_analysis(cable_forces, material_key=None):
    """
    Print cable sizing for the configured or specified material.
    """
    from config import CABLE_MATERIAL
    
    # Use configured material if not specified
    if material_key is None:
        # Find the key for the configured material
        for key, mat in CABLE_MATERIALS.items():
            if mat.name == CABLE_MATERIAL.name:
                material_key = key
                break
        if material_key is None:
            material_key = 'steel'  # fallback
    
    print(f"\n{'='*120}")
    print("CABLE MATERIAL ANALYSIS")
    print(f"{'='*120}")
    print_laced_cable_analysis(cable_forces, material_key)


def print_laced_cable_analysis(cable_forces, mat_key='steel'):
    """
    Print detailed laced cable analysis with angular distribution.
    Works with any cable material that has fiber_diameter defined.
    """
    mat = CABLE_MATERIALS[mat_key]
    results, total_mass, total_strands = analyze_cable_system(cable_forces, mat_key)
    
    # Check if material has fiber_diameter
    if mat.fiber_diameter is None:
        print(f"  Material {mat.name} does not have fiber_diameter defined.")
        return
    
    d_fiber = mat.fiber_diameter * 1000  # mm
    n_cables_per_side = cable_forces[0].get('n_cables_per_side', 1) if cable_forces else 1
    
    # Material-specific notes
    notes = {
        'steel': '(galvanized for corrosion protection)',
        'dyneema': '(floats on water!)',
        'carbon': '(requires dielectric isolation at anchors)',
        'aramid': '(UV-protected sheath required)',
        'bfrp': '(corrosion-free, non-magnetic)'
    }
    note = notes.get(mat_key, '')
    
    print(f"\n{'═'*100}")
    print(f"  {mat.name.upper()} EGG ALIGNMENT CABLE SYSTEM")
    print(f"{'═'*100}")
    print(f"  Material: {mat.name}")
    print(f"  Tensile strength: {mat.tensile_strength/1e6:.0f} MPa")
    print(f"  Density: {mat.density:.0f} kg/m³ {note}")
    print(f"  Safety factor: {mat.safety_factor:.1f}")
    print(f"")
    print(f"  🔗 LACED CONFIGURATION: {d_fiber:.0f}mm diameter individual ropes")
    print(f"  📐 ANGULAR DISTRIBUTION: {n_cables_per_side} cables per quarter (×4 = {n_cables_per_side*4} around circumference)")
    print(f"{'═'*100}")
    
    print(f"\n  {'Level':<18} {'Diameter':<10} {'Tension':<10} {'Strands':<12} {'Density':<14} {'Length':<10} {'Mass':<10}")
    print(f"  {'':<18} {'(m)':<10} {'(MN)':<10} {'(total)':<12} {'(per m circ)':<14} {'(m)':<10} {'(tonnes)':<10}")
    print(f"  {'─'*90}")
    
    max_strands_at_level = 0
    max_strand_density = 0
    
    for r in results:
        to_str = 'ground' if r['to_egg'] == 'ground' else f"egg {r['to_egg']}"
        total_strands_level = r['total_strands_at_level']
        strand_density = r['strand_density']
        mass_level = r['total_mass_at_level_kg'] / 1000
        
        print(f"  egg {r['from_egg']} → {to_str:<8} {r['diameter']:<10.1f} {r['T_max_kN']/1000:<10.2f} "
              f"{total_strands_level:<12} {strand_density:<14.1f} {r['cable_length']:<10.1f} {mass_level:<10.2f}")
        
        max_strands_at_level = max(max_strands_at_level, total_strands_level)
        max_strand_density = max(max_strand_density, strand_density)
    
    print(f"  {'─'*90}")
    print(f"")
    print(f"  SUMMARY:")
    print(f"  {'─'*50}")
    print(f"  Individual rope diameter:     {d_fiber:.0f} mm")
    print(f"  Cables per quarter circle:    {n_cables_per_side}")
    print(f"  Max strands at one level:     {max_strands_at_level:,}")
    print(f"  Max strand density:           {max_strand_density:.1f} strands/m circumference")
    print(f"  Total strand segments:        {total_strands:,} (all levels)")
    print(f"  Total alignment cable mass:   {total_mass/1000:,.1f} tonnes")
    print(f"")
    
    A_rope = np.pi * (mat.fiber_diameter/2)**2 * mat.fill_factor
    breaking_load = A_rope * mat.tensile_strength / 1000
    working_load = breaking_load / mat.safety_factor
    
    print(f"  ROPE SPECIFICATIONS ({d_fiber:.0f}mm {mat_key.upper()}):")
    print(f"  {'─'*50}")
    print(f"  Breaking load per rope:       {breaking_load:,.0f} kN ({breaking_load/1000:.2f} MN)")
    print(f"  Working load (SF={mat.safety_factor}):        {working_load:,.0f} kN ({working_load/1000:.2f} MN)")
    print(f"  Mass per meter:               {A_rope/mat.fill_factor * mat.density:.1f} kg/m")
    print(f"")
    
    print(f"  STRAND DENSITY INTERPRETATION:")
    print(f"  {'─'*50}")
    print(f"  Strand density = total strands around circumference / πD")
    print(f"  Higher density at lower levels due to cumulative wind load")
    if max_strand_density > 100:
        spacing_mm = 1000 / max_strand_density
        print(f"  At max density: 1 strand every {spacing_mm:.0f}mm around circumference")
    print(f"")
    
    print(f"  PRACTICAL CONSIDERATIONS:")
    print(f"  {'─'*50}")
    if max_strand_density > 50:
        print(f"  ⚠️  High strand density ({max_strand_density:.0f}/m) requires careful attachment design")
        print(f"      Consider larger individual ropes or bundling")
    else:
        print(f"  ✓ Manageable strand density")
    
    print(f"  • Load naturally reduces toward perpendicular (cos distribution)")
    print(f"  • Redundancy: failure of one strand doesn't cause collapse")
    print(f"  • Easier to transport and install than single massive cable")
    print(f"  • Can be tensioned individually for load balancing")
