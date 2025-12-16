#!/usr/bin/env python3
"""
Egg Tower Structural Optimization

A 1600m tall tower made of stacked hollow egg-shaped shells (BFRP).
Each egg is connected to adjacent eggs via guy cables at their maximum diameter points.

Modules:
    config.py       - Configuration constants and material properties
    geometry.py     - Egg geometry and tower building functions
    structural.py   - Stress and buckling analysis
    wind_cables.py  - Wind loading and cable force calculations
    output.py       - Printing and output functions
    visualization.py - Tower diagram plotting
    optimization.py - Diameter optimization algorithms
"""

from config import TOTAL_HEIGHT, D_TOP, N_CABLES_PER_LEVEL, THICKNESS_FACTOR
from geometry import build_tower
from structural import analyze_tower
from wind_cables import calculate_alignment_cable_forces, wind_force_on_egg, analyze_cable_system
from output import print_tower_summary, print_alignment_cable_forces, print_cable_material_analysis, print_vortex_shedding_analysis
from visualization import plot_tapered_tower, plot_wind_analysis
from optimization import binary_search_optimal_diameter
from ring_guys import (calculate_ring_guy_system, print_ring_guy_analysis, 
                       calculate_guy_placement_from_alignment_capacity, print_guy_placement_analysis)
from three_d_visualization import plot_tower_3d, plot_stress_distribution_3d
from active_stabilization import design_active_stabilization
from elevator_system import design_elevator_system
from project_costs import calculate_all_project_costs, print_full_project_costs


def main():
    """Main execution: optimize tower and display results."""
    
    # Run binary search for optimal diameter (minimize total mass)
    optimal = binary_search_optimal_diameter(
        objective='total_mass', 
        d_min=20, 
        d_max=200, 
        tol=1.0
    )
    
    print(f"\n\n{'='*100}")
    print(f"DETAILED ANALYSIS FOR OPTIMAL D_BASE = {optimal['d_base']:.1f}m (thickness factor = {THICKNESS_FACTOR:.2f})")
    print(f"{'='*100}")
    
    # Build the optimal tower (initial analysis without guy forces)
    # Use increased thickness factor to handle guy wire compression
    eggs, ratio = build_tower(
        d_base=optimal['d_base'], 
        d_top=D_TOP, 
        target_height=TOTAL_HEIGHT,
        thickness_factor=THICKNESS_FACTOR
    )
    
    # Initial analysis to get geometry (guy forces added later)
    analyzed_eggs = analyze_tower(eggs)
    
    # Calculate wind forces on each egg
    for egg in analyzed_eggs:
        wind_data = wind_force_on_egg(egg)
        egg['wind_force'] = wind_data['wind_force']
    
    # Calculate guy wire placement based on alignment wire capacity
    # This determines how many eggs the alignment wires can handle
    # before needing a guy wire to transfer the horizontal force
    guy_placement = calculate_guy_placement_from_alignment_capacity(analyzed_eggs)
    print_guy_placement_analysis(guy_placement)
    
    # Use the calculated guy indices for the ring guy system
    guy_indices = guy_placement['guy_indices']
    
    # Ring + ground-anchored guy wire system analysis
    # Pass the calculated guy_indices to use our placement
    guy_system = calculate_ring_guy_system(
        analyzed_eggs, 
        guy_indices=guy_indices  # Use our calculated placement
    )
    
    # RE-ANALYZE with guy wire vertical forces included
    # This gives accurate stress including guy wire compression
    analyzed_eggs = analyze_tower(eggs, guy_system=guy_system)
    
    # Restore wind forces (lost during re-analysis)
    for egg in analyzed_eggs:
        wind_data = wind_force_on_egg(egg)
        egg['wind_force'] = wind_data['wind_force']
    
    # Print tower summary WITH guy forces included
    print_tower_summary(analyzed_eggs)
    
    # Calculate alignment cable forces (section-local loads)
    cable_forces = calculate_alignment_cable_forces(
        analyzed_eggs, 
        guy_indices,
        n_cables_per_side=N_CABLES_PER_LEVEL
    )
    print_alignment_cable_forces(cable_forces, guy_indices)
    
    # Vortex shedding analysis
    print_vortex_shedding_analysis(analyzed_eggs)
    
    # Analyze cable materials
    print_cable_material_analysis(cable_forces)
    
    # Print guy system analysis
    print_ring_guy_analysis(guy_system)
    
    # Plot the tower structure diagram with guy system
    plot_tapered_tower(analyzed_eggs, cable_forces=cable_forces, guy_system=guy_system)
    
    # Plot dedicated wind analysis chart (includes guy cable tensions)
    plot_wind_analysis(analyzed_eggs, cable_forces, guy_system=guy_system)
    
    # Generate 3D visualizations
    guy_system_3d = {
        'anchor_radius': guy_system['anchor_radius'],
        'n_cables': guy_system['n_guys_per_ring'],
        'attachment_heights': [level['z_attach'] for level in guy_system['levels']]
    }
    
    # 3D wireframe tower view (high resolution with more cables)
    plot_tower_3d(analyzed_eggs, guy_system=guy_system_3d, 
                  wireframe_only=False, n_cables_display=40,
                  high_resolution=True, figsize=(20, 28),
                  save_path='egg_tower_3d.png')
    
    # 3D stress distribution visualization
    stresses = [egg['stress_Pa'] for egg in analyzed_eggs]
    plot_stress_distribution_3d(analyzed_eggs, stresses, 
                                save_path='egg_tower_stress_3d.png')
    
    # === COST ANALYSIS ===
    print_cost_analysis(analyzed_eggs, cable_forces, guy_system)
    save_cost_analysis_png(analyzed_eggs, cable_forces, guy_system)


def save_cost_analysis_png(eggs, cable_forces, guy_system, save_path='cost_analysis.png'):
    """
    Generate and save cost analysis as a PNG table image with full project costs.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.table import Table
    from project_costs import calculate_all_project_costs
    
    # === CALCULATE ALL COSTS (same as print function) ===
    COST_BFRP_PER_KG = 10.0
    COST_ARAMID_PER_KG = 30.0
    INSTALL_MULT_SHELL = 3.0
    INSTALL_MULT_CABLE = 2.0
    
    shell_mass_kg = sum(e['mass'] for e in eggs)
    guy_cable_mass_kg = guy_system['total_system_mass_tonnes'] * 1000
    cable_results, alignment_mass_kg, total_strands = analyze_cable_system(cable_forces)
    total_cable_mass_kg = guy_cable_mass_kg + alignment_mass_kg
    total_mass_kg = shell_mass_kg + total_cable_mass_kg
    
    shell_material = shell_mass_kg * COST_BFRP_PER_KG
    cable_material = total_cable_mass_kg * COST_ARAMID_PER_KG
    total_material = shell_material + cable_material
    
    shell_installed = shell_material * INSTALL_MULT_SHELL
    cable_installed = cable_material * INSTALL_MULT_CABLE
    structure_installed = shell_installed + cable_installed
    
    n_eggs = len(eggs)
    n_guy_levels = len(guy_system['levels'])
    tower_height = eggs[-1]['z_base'] + eggs[-1]['h']
    base_diameter = eggs[0]['d']
    
    stab_system = design_active_stabilization(
        n_eggs=n_eggs, n_guy_levels=n_guy_levels,
        n_cables_per_level=9, tower_height=tower_height
    )
    stab_cost = stab_system['total_cost']
    
    elevator_system = design_elevator_system(
        tower_height=tower_height, n_shafts=4,
        passenger_throughput_per_hour=500
    )
    elevator_cost = elevator_system['total_cost']
    elevator_mass_kg = elevator_system['system_mass_tonnes'] * 1000
    total_mass_kg += elevator_mass_kg
    
    # Calculate full project costs
    full_costs = calculate_all_project_costs(
        tower_height=tower_height,
        base_diameter=base_diameter,
        n_eggs=n_eggs,
        n_guy_levels=n_guy_levels,
        shell_mass_tonnes=shell_mass_kg/1000,
        cable_mass_tonnes=total_cable_mass_kg/1000,
        structure_cost=structure_installed,
        stabilization_cost=stab_cost,
        elevator_cost=elevator_cost
    )
    
    materials_systems_total = structure_installed + stab_cost + elevator_cost
    grand_total = full_costs['grand_total_with_contingency']
    
    # === CREATE FIGURE ===
    fig, ax = plt.subplots(figsize=(14, 28))  # Taller for full project costs
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    # Colors
    header_color = '#2C3E50'
    section_color = '#34495E'
    row_color1 = '#ECF0F1'
    row_color2 = '#FFFFFF'
    total_color = '#27AE60'
    text_color = '#2C3E50'
    
    # Title
    ax.text(0.5, 0.97, 'EGG TOWER COST ANALYSIS', fontsize=24, fontweight='bold',
            ha='center', va='top', color=header_color)
    ax.text(0.5, 0.94, f'{tower_height:.0f}m Tower | {n_eggs} Eggs | {n_guy_levels} Guy Levels',
            fontsize=14, ha='center', va='top', color='#7F8C8D')
    
    # Table data
    y_pos = 0.88
    row_height = 0.022
    
    def draw_row(y, label, value, pct=None, color=row_color1, bold=False, indent=0):
        weight = 'bold' if bold else 'normal'
        ax.add_patch(mpatches.Rectangle((0.05, y - row_height/2), 0.9, row_height,
                                         facecolor=color, edgecolor='#BDC3C7', linewidth=0.5))
        ax.text(0.08 + indent*0.02, y, label, fontsize=11, va='center', 
                fontweight=weight, color=text_color)
        ax.text(0.65, y, value, fontsize=11, va='center', ha='right',
                fontweight=weight, color=text_color, family='monospace')
        if pct is not None:
            ax.text(0.92, y, pct, fontsize=11, va='center', ha='right',
                    fontweight=weight, color=text_color, family='monospace')
        return y - row_height
    
    def draw_header(y, title):
        ax.add_patch(mpatches.Rectangle((0.05, y - row_height/2), 0.9, row_height,
                                         facecolor=section_color, edgecolor='#2C3E50', linewidth=1))
        ax.text(0.5, y, title, fontsize=12, fontweight='bold', va='center',
                ha='center', color='white')
        return y - row_height
    
    def draw_total_row(y, label, value, pct=None):
        ax.add_patch(mpatches.Rectangle((0.05, y - row_height/2), 0.9, row_height,
                                         facecolor=total_color, edgecolor='#1E8449', linewidth=1))
        ax.text(0.08, y, label, fontsize=12, fontweight='bold', va='center', color='white')
        ax.text(0.65, y, value, fontsize=12, va='center', ha='right',
                fontweight='bold', color='white', family='monospace')
        if pct is not None:
            ax.text(0.92, y, pct, fontsize=12, va='center', ha='right',
                    fontweight='bold', color='white', family='monospace')
        return y - row_height
    
    # Column headers
    y_pos = draw_header(y_pos, 'COST BREAKDOWN')
    ax.text(0.08, y_pos + row_height, 'Item', fontsize=10, fontweight='bold', va='center')
    ax.text(0.65, y_pos + row_height, 'Cost (USD)', fontsize=10, fontweight='bold', va='center', ha='right')
    ax.text(0.92, y_pos + row_height, '% of Total', fontsize=10, fontweight='bold', va='center', ha='right')
    
    # Structure section
    y_pos -= row_height * 0.5
    y_pos = draw_header(y_pos, 'STRUCTURE')
    y_pos = draw_row(y_pos, 'Shell material (BFRP @ $10/kg)', f'${shell_material/1e6:,.1f}M', 
                     f'{shell_material/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, 'Shell installation (×3.0)', f'${shell_installed/1e6:,.1f}M', 
                     f'{shell_installed/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Cable material (Aramid @ $30/kg)', f'${cable_material/1e6:,.1f}M', 
                     f'{cable_material/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, 'Cable installation (×2.0)', f'${cable_installed/1e6:,.1f}M', 
                     f'{cable_installed/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Subtotal Structure', f'${structure_installed/1e6:,.1f}M', 
                     f'{structure_installed/grand_total*100:.1f}%', '#D5DBDB', bold=True)
    
    # Stabilization section
    y_pos -= row_height * 0.3
    y_pos = draw_header(y_pos, 'ACTIVE STABILIZATION')
    y_pos = draw_row(y_pos, f'Hydraulic actuators ({stab_system["components"]["actuators"]["count"]})', 
                     f'${stab_system["costs"]["actuators"]/1e6:,.2f}M',
                     f'{stab_system["costs"]["actuators"]/grand_total*100:.1f}%', row_color1, indent=1)
    sensors_cost = (stab_system['costs']['accelerometers'] + stab_system['costs']['gps'] + 
                   stab_system['costs']['inclinometers'] + stab_system['costs']['wind_sensors'] + 
                   stab_system['costs']['strain_gauges'])
    y_pos = draw_row(y_pos, 'Sensors & instrumentation', f'${sensors_cost/1e6:,.2f}M',
                     f'{sensors_cost/grand_total*100:.1f}%', row_color2, indent=1)
    controls_cost = stab_system['costs']['controllers'] + stab_system['costs']['software']
    y_pos = draw_row(y_pos, 'Control system & software', f'${controls_cost/1e6:,.2f}M',
                     f'{controls_cost/grand_total*100:.1f}%', row_color1, indent=1)
    infra_cost = (stab_system['costs']['fiber_network'] + stab_system['costs']['power_backup'] + 
                 stab_system['costs']['installation'] + stab_system['costs']['maintenance_year1'])
    y_pos = draw_row(y_pos, 'Power, networking, installation', f'${infra_cost/1e6:,.2f}M',
                     f'{infra_cost/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Subtotal Stabilization', f'${stab_cost/1e6:,.1f}M',
                     f'{stab_cost/grand_total*100:.1f}%', '#D5DBDB', bold=True)
    
    # Elevator section
    y_pos -= row_height * 0.3
    y_pos = draw_header(y_pos, 'SELF-CLIMBING ELEVATOR')
    capsule_cost = (elevator_system['costs']['standard_capsules'] + elevator_system['costs']['express_capsules'] +
                   elevator_system['costs']['service_capsules'] + elevator_system['costs']['emergency_capsules'])
    y_pos = draw_row(y_pos, f'Capsules ({elevator_system["capsules"]["total"]} total, {elevator_system["n_shafts"]} shafts)',
                     f'${capsule_cost/1e6:,.1f}M', f'{capsule_cost/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, f'Guide rails ({elevator_system["infrastructure"]["rail_length_m"]/1000:.0f}km)',
                     f'${elevator_system["costs"]["guide_rails"]/1e6:,.1f}M',
                     f'{elevator_system["costs"]["guide_rails"]/grand_total*100:.1f}%', row_color2, indent=1)
    station_cost = (elevator_system['costs']['stations'] + elevator_system['costs']['power_system'] +
                   elevator_system['costs']['control_system'] + elevator_system['costs']['maintenance_facility'])
    y_pos = draw_row(y_pos, 'Stations, power, controls', f'${station_cost/1e6:,.1f}M',
                     f'{station_cost/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, 'Installation', f'${elevator_system["costs"]["installation"]/1e6:,.1f}M',
                     f'{elevator_system["costs"]["installation"]/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Subtotal Elevator', f'${elevator_cost/1e6:,.1f}M',
                     f'{elevator_cost/grand_total*100:.1f}%', '#D5DBDB', bold=True)
    
    # Materials & Systems subtotal
    y_pos -= row_height * 0.3
    y_pos = draw_row(y_pos, '═══ MATERIALS & SYSTEMS SUBTOTAL', f'${materials_systems_total/1e6:,.1f}M',
                     f'{materials_systems_total/grand_total*100:.1f}%', '#AED6F1', bold=True)
    
    # === FOUNDATION SECTION ===
    fd = full_costs['foundation']
    y_pos -= row_height * 0.5
    y_pos = draw_header(y_pos, 'FOUNDATION')
    cf = fd['central_foundation']
    y_pos = draw_row(y_pos, f'Central foundation ({cf["diameter_m"]:.0f}m × {cf["depth_m"]:.0f}m deep)',
                     f'${cf["cost"]/1e6:,.1f}M', f'{cf["cost"]/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, f'Guy wire anchors ({fd["guy_anchors"]["n_anchors"]} locations)',
                     f'${fd["guy_anchors"]["cost"]/1e6:,.1f}M', 
                     f'{fd["guy_anchors"]["cost"]/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Access roads & site prep',
                     f'${fd["access_roads"]["cost"]/1e6:,.1f}M',
                     f'{fd["access_roads"]["cost"]/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, 'Subtotal Foundation', f'${fd["total_cost"]/1e6:,.1f}M',
                     f'{fd["total_cost"]/grand_total*100:.1f}%', '#D5DBDB', bold=True)
    
    # === CONSTRUCTION SECTION ===
    cn = full_costs['construction']
    y_pos -= row_height * 0.3
    y_pos = draw_header(y_pos, f'CONSTRUCTION ({cn["construction_months"]:.0f} months)')
    y_pos = draw_row(y_pos, 'Cranes & heavy equipment',
                     f'${cn["cranes"]/1e6:,.1f}M', f'{cn["cranes"]/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, 'Scaffolding & temporary works',
                     f'${(cn["scaffolding"]+cn["temporary_works"])/1e6:,.1f}M', 
                     f'{(cn["scaffolding"]+cn["temporary_works"])/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Labor & site facilities',
                     f'${(cn["labor"]+cn["site_facilities"])/1e6:,.1f}M',
                     f'{(cn["labor"]+cn["site_facilities"])/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, 'Transport, QC, weather contingency',
                     f'${(cn["transport"]+cn["quality_control"]+cn["weather_contingency"])/1e6:,.1f}M',
                     f'{(cn["transport"]+cn["quality_control"]+cn["weather_contingency"])/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Subtotal Construction', f'${cn["total_cost"]/1e6:,.1f}M',
                     f'{cn["total_cost"]/grand_total*100:.1f}%', '#D5DBDB', bold=True)
    
    # === ENGINEERING SECTION ===
    eg = full_costs['engineering']
    y_pos -= row_height * 0.3
    y_pos = draw_header(y_pos, 'ENGINEERING & PROJECT MANAGEMENT')
    y_pos = draw_row(y_pos, 'Design & structural engineering',
                     f'${eg["design_engineering"]/1e6:,.1f}M', 
                     f'{eg["design_engineering"]/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, 'Testing (wind tunnel, prototypes)',
                     f'${(eg["wind_tunnel"]+eg["prototype_testing"])/1e6:,.1f}M',
                     f'{(eg["wind_tunnel"]+eg["prototype_testing"])/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Project management & insurance',
                     f'${(eg["project_management"]+eg["insurance"])/1e6:,.1f}M',
                     f'{(eg["project_management"]+eg["insurance"])/grand_total*100:.1f}%', row_color1, indent=1)
    y_pos = draw_row(y_pos, 'Regulatory, legal, contingency',
                     f'${(eg["regulatory"]+eg["legal"]+eg["contingency"])/1e6:,.1f}M',
                     f'{(eg["regulatory"]+eg["legal"]+eg["contingency"])/grand_total*100:.1f}%', row_color2, indent=1)
    y_pos = draw_row(y_pos, 'Subtotal Engineering', f'${eg["total_cost"]/1e6:,.1f}M',
                     f'{eg["total_cost"]/grand_total*100:.1f}%', '#D5DBDB', bold=True)
    
    # === PROJECT CONTINGENCY ===
    y_pos -= row_height * 0.3
    y_pos = draw_row(y_pos, 'PROJECT CONTINGENCY (10%)', f'${full_costs["project_contingency"]/1e6:,.1f}M',
                     f'{full_costs["project_contingency"]/grand_total*100:.1f}%', '#F9E79F', bold=True)
    
    # Grand total
    y_pos -= row_height * 0.5
    y_pos = draw_total_row(y_pos, 'GRAND TOTAL', f'${grand_total/1e6:,.1f}M', '100.0%')
    
    # Key metrics section
    y_pos -= row_height * 1.0
    y_pos = draw_header(y_pos, 'KEY METRICS')
    y_pos = draw_row(y_pos, 'Cost per meter of height', f'${grand_total/tower_height:,.0f}/m', color=row_color1)
    y_pos = draw_row(y_pos, 'Total mass', f'{total_mass_kg/1000:,.0f} tonnes', color=row_color2)
    y_pos = draw_row(y_pos, 'Structure mass', f'{(shell_mass_kg + total_cable_mass_kg)/1000:,.0f} tonnes', color=row_color1)
    y_pos = draw_row(y_pos, 'Elevator mass', f'{elevator_mass_kg/1000:,.0f} tonnes', color=row_color2)
    y_pos = draw_row(y_pos, 'Express travel time', f'{elevator_system["performance"]["travel_time_express_min"]:.1f} min', color=row_color1)
    y_pos = draw_row(y_pos, 'Passenger throughput', f'{elevator_system["performance"]["throughput_per_hour"]:,}/hour', color=row_color2)
    
    # Mass breakdown pie chart
    ax_pie = fig.add_axes([0.55, 0.01, 0.4, 0.10])
    masses = [shell_mass_kg, alignment_mass_kg, guy_cable_mass_kg, elevator_mass_kg]
    labels = ['Shell\n(BFRP)', 'Alignment\nCables', 'Guy\nCables', 'Elevator\nSystem']
    colors = ['#3498DB', '#E74C3C', '#F39C12', '#9B59B6']
    ax_pie.pie(masses, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90,
               textprops={'fontsize': 8})
    ax_pie.set_title('Mass Distribution', fontsize=10, fontweight='bold')
    
    # Cost breakdown pie chart - now shows full project costs
    ax_pie2 = fig.add_axes([0.08, 0.01, 0.4, 0.10])
    costs = [materials_systems_total, fd['total_cost'], cn['total_cost'], 
             eg['total_cost'], full_costs['project_contingency']]
    labels2 = ['Materials &\nSystems', 'Foundation', 'Construction', 'Engineering', 'Contingency']
    colors2 = ['#3498DB', '#E67E22', '#27AE60', '#9B59B6', '#F1C40F']
    ax_pie2.pie(costs, labels=labels2, colors=colors2, autopct='%1.1f%%', startangle=90,
                textprops={'fontsize': 8})
    ax_pie2.set_title('Full Project Cost Distribution', fontsize=10, fontweight='bold')
    
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"✓ Cost analysis saved: '{save_path}'")


def print_cost_analysis(eggs, cable_forces, guy_system):
    """
    Print comprehensive cost analysis table for the tower.
    """
    # Material cost rates (USD per kg)
    COST_BFRP_PER_KG = 10.0       # Mid-range BFRP
    COST_ARAMID_PER_KG = 30.0     # Kevlar 49 fiber
    
    # Installation multipliers
    INSTALL_MULT_SHELL = 3.0
    INSTALL_MULT_CABLE = 2.0
    
    # Calculate masses
    shell_mass_kg = sum(e['mass'] for e in eggs)
    guy_cable_mass_kg = guy_system['total_system_mass_tonnes'] * 1000
    
    # Alignment cable mass from cable forces
    cable_results, alignment_mass_kg, total_strands = analyze_cable_system(cable_forces)
    
    total_cable_mass_kg = guy_cable_mass_kg + alignment_mass_kg
    total_mass_kg = shell_mass_kg + total_cable_mass_kg
    
    # Material costs
    shell_material = shell_mass_kg * COST_BFRP_PER_KG
    cable_material = total_cable_mass_kg * COST_ARAMID_PER_KG
    total_material = shell_material + cable_material
    
    # Installed costs
    shell_installed = shell_material * INSTALL_MULT_SHELL
    cable_installed = cable_material * INSTALL_MULT_CABLE
    structure_installed = shell_installed + cable_installed
    
    # Active stabilization
    n_eggs = len(eggs)
    n_guy_levels = len(guy_system['levels'])
    tower_height = eggs[-1]['z_base'] + eggs[-1]['h']
    
    stab_system = design_active_stabilization(
        n_eggs=n_eggs,
        n_guy_levels=n_guy_levels,
        n_cables_per_level=9,
        tower_height=tower_height
    )
    stab_cost = stab_system['total_cost']
    
    # Self-climbing elevator system
    elevator_system = design_elevator_system(
        tower_height=tower_height,
        n_shafts=4,
        passenger_throughput_per_hour=500
    )
    elevator_cost = elevator_system['total_cost']
    elevator_mass_kg = elevator_system['system_mass_tonnes'] * 1000
    
    # Update total mass to include elevator
    total_mass_kg += elevator_mass_kg
    
    # Grand total
    grand_total = structure_installed + stab_cost + elevator_cost
    
    # Print table
    print(f"\n\n{'═'*90}")
    print(f"COST ANALYSIS SUMMARY")
    print(f"{'═'*90}")
    
    print(f"\n┌{'─'*88}┐")
    print(f"│{'TOWER CONFIGURATION':^88}│")
    print(f"├{'─'*44}┬{'─'*43}┤")
    print(f"│ {'Height':<42} │ {tower_height:>35,.0f} m   │")
    print(f"│ {'Number of eggs':<42} │ {n_eggs:>35,d}     │")
    print(f"│ {'Base diameter':<42} │ {eggs[0]['d']:>35.1f} m   │")
    print(f"│ {'Top diameter':<42} │ {eggs[-1]['d']:>35.1f} m   │")
    print(f"│ {'Guy wire levels':<42} │ {n_guy_levels:>35,d}     │")
    print(f"│ {'Thickness factor':<42} │ {THICKNESS_FACTOR:>35.2f}     │")
    print(f"└{'─'*44}┴{'─'*43}┘")
    
    print(f"\n┌{'─'*88}┐")
    print(f"│{'MASS BREAKDOWN':^88}│")
    print(f"├{'─'*44}┬{'─'*21}┬{'─'*21}┤")
    print(f"│ {'Component':<42} │ {'Mass (tonnes)':>19} │ {'% of Total':>19} │")
    print(f"├{'─'*44}┼{'─'*21}┼{'─'*21}┤")
    print(f"│ {'Shell (BFRP)':<42} │ {shell_mass_kg/1000:>19,.0f} │ {shell_mass_kg/total_mass_kg*100:>18.1f}% │")
    print(f"│ {'Alignment cables (Aramid)':<42} │ {alignment_mass_kg/1000:>19,.0f} │ {alignment_mass_kg/total_mass_kg*100:>18.1f}% │")
    print(f"│ {'Guy cables (Aramid)':<42} │ {guy_cable_mass_kg/1000:>19,.0f} │ {guy_cable_mass_kg/total_mass_kg*100:>18.1f}% │")
    print(f"│ {'Elevator system (rails + capsules)':<42} │ {elevator_mass_kg/1000:>19,.0f} │ {elevator_mass_kg/total_mass_kg*100:>18.1f}% │")
    print(f"├{'─'*44}┼{'─'*21}┼{'─'*21}┤")
    print(f"│ {'TOTAL MASS':<42} │ {total_mass_kg/1000:>19,.0f} │ {100:>18.1f}% │")
    print(f"└{'─'*44}┴{'─'*21}┴{'─'*21}┘")
    
    print(f"\n┌{'─'*88}┐")
    print(f"│{'COST BREAKDOWN':^88}│")
    print(f"├{'─'*44}┬{'─'*21}┬{'─'*21}┤")
    print(f"│ {'Item':<42} │ {'Cost (USD)':>19} │ {'% of Total':>19} │")
    print(f"├{'─'*44}┼{'─'*21}┼{'─'*21}┤")
    print(f"│ {'MATERIALS':^42} │{' ':>21} │{' ':>21} │")
    print(f"│ {'  Shell (BFRP @ $10/kg)':<42} │ ${shell_material/1e6:>17.1f}M │ {shell_material/grand_total*100:>18.1f}% │")
    print(f"│ {'  Cables (Aramid @ $30/kg)':<42} │ ${cable_material/1e6:>17.1f}M │ {cable_material/grand_total*100:>18.1f}% │")
    print(f"│ {'  Subtotal materials':<42} │ ${total_material/1e6:>17.1f}M │ {total_material/grand_total*100:>18.1f}% │")
    print(f"├{'─'*44}┼{'─'*21}┼{'─'*21}┤")
    print(f"│ {'INSTALLATION (fabrication, transport, etc.)':^42} │{' ':>21} │{' ':>21} │")
    print(f"│ {'  Shell installation (×3.0)':<42} │ ${shell_installed/1e6:>17.1f}M │ {shell_installed/grand_total*100:>18.1f}% │")
    print(f"│ {'  Cable installation (×2.0)':<42} │ ${cable_installed/1e6:>17.1f}M │ {cable_installed/grand_total*100:>18.1f}% │")
    print(f"│ {'  Subtotal installed structure':<42} │ ${structure_installed/1e6:>17.1f}M │ {structure_installed/grand_total*100:>18.1f}% │")
    print(f"├{'─'*44}┼{'─'*21}┼{'─'*21}┤")
    print(f"│ {'ACTIVE STABILIZATION SYSTEM':^42} │{' ':>21} │{' ':>21} │")
    print(f"│ {'  Actuators ({:d} hydraulic tensioners)':<42} │ ${stab_system['costs']['actuators']/1e6:>17.2f}M │ {stab_system['costs']['actuators']/grand_total*100:>18.1f}% │".format(stab_system['components']['actuators']['count']))
    print(f"│ {'  Sensors & instrumentation':<42} │ ${(stab_system['costs']['accelerometers']+stab_system['costs']['gps']+stab_system['costs']['inclinometers']+stab_system['costs']['wind_sensors']+stab_system['costs']['strain_gauges'])/1e6:>17.2f}M │ {(stab_system['costs']['accelerometers']+stab_system['costs']['gps']+stab_system['costs']['inclinometers']+stab_system['costs']['wind_sensors']+stab_system['costs']['strain_gauges'])/grand_total*100:>18.1f}% │")
    print(f"│ {'  Control system & software':<42} │ ${(stab_system['costs']['controllers']+stab_system['costs']['software'])/1e6:>17.2f}M │ {(stab_system['costs']['controllers']+stab_system['costs']['software'])/grand_total*100:>18.1f}% │")
    print(f"│ {'  Power, networking, installation':<42} │ ${(stab_system['costs']['fiber_network']+stab_system['costs']['power_backup']+stab_system['costs']['installation']+stab_system['costs']['maintenance_year1'])/1e6:>17.2f}M │ {(stab_system['costs']['fiber_network']+stab_system['costs']['power_backup']+stab_system['costs']['installation']+stab_system['costs']['maintenance_year1'])/grand_total*100:>18.1f}% │")
    print(f"│ {'  Subtotal stabilization':<42} │ ${stab_cost/1e6:>17.1f}M │ {stab_cost/grand_total*100:>18.1f}% │")
    print(f"├{'─'*44}┼{'─'*21}┼{'─'*21}┤")
    print(f"│ {'SELF-CLIMBING ELEVATOR SYSTEM':^42} │{' ':>21} │{' ':>21} │")
    print(f"│ {'  Capsules ({:d} total, {:d} shafts)':<42} │ ${(elevator_system['costs']['standard_capsules']+elevator_system['costs']['express_capsules']+elevator_system['costs']['service_capsules']+elevator_system['costs']['emergency_capsules'])/1e6:>17.1f}M │ {(elevator_system['costs']['standard_capsules']+elevator_system['costs']['express_capsules']+elevator_system['costs']['service_capsules']+elevator_system['costs']['emergency_capsules'])/grand_total*100:>18.1f}% │".format(elevator_system['capsules']['total'], elevator_system['n_shafts']))
    print(f"│ {'  Guide rails ({:.0f}km rack-and-pinion)':<42} │ ${elevator_system['costs']['guide_rails']/1e6:>17.1f}M │ {elevator_system['costs']['guide_rails']/grand_total*100:>18.1f}% │".format(elevator_system['infrastructure']['rail_length_m']/1000))
    print(f"│ {'  Stations, power, controls':<42} │ ${(elevator_system['costs']['stations']+elevator_system['costs']['power_system']+elevator_system['costs']['control_system']+elevator_system['costs']['maintenance_facility'])/1e6:>17.1f}M │ {(elevator_system['costs']['stations']+elevator_system['costs']['power_system']+elevator_system['costs']['control_system']+elevator_system['costs']['maintenance_facility'])/grand_total*100:>18.1f}% │")
    print(f"│ {'  Installation':<42} │ ${elevator_system['costs']['installation']/1e6:>17.1f}M │ {elevator_system['costs']['installation']/grand_total*100:>18.1f}% │")
    print(f"│ {'  Subtotal elevator':<42} │ ${elevator_cost/1e6:>17.1f}M │ {elevator_cost/grand_total*100:>18.1f}% │")
    print(f"├{'─'*44}┼{'─'*21}┼{'─'*21}┤")
    print(f"│ {'GRAND TOTAL':<42} │ ${grand_total/1e6:>17.1f}M │ {100:>18.1f}% │")
    print(f"└{'─'*44}┴{'─'*21}┴{'─'*21}┘")
    
    print(f"\n┌{'─'*88}┐")
    print(f"│{'KEY METRICS':^88}│")
    print(f"├{'─'*44}┬{'─'*43}┤")
    print(f"│ {'Cost per meter of height':<42} │ ${grand_total/tower_height:>32,.0f} /m   │")
    print(f"│ {'Cost per tonne of structure':<42} │ ${structure_installed/(total_mass_kg/1000):>32,.0f} /t   │")
    print(f"│ {'Structure cost':<42} │ ${structure_installed/1e6:>32.1f} M   │")
    print(f"│ {'Stabilization cost':<42} │ ${stab_cost/1e6:>32.1f} M   │")
    print(f"│ {'Elevator cost':<42} │ ${elevator_cost/1e6:>32.1f} M   │")
    print(f"│ {'Systems as % of total':<42} │ {(stab_cost+elevator_cost)/grand_total*100:>32.1f} %   │")
    print(f"│ {'Express travel time (top to bottom)':<42} │ {elevator_system['performance']['travel_time_express_min']:>32.1f} min │")
    print(f"└{'─'*44}┴{'─'*43}┘")
    
    print(f"\n{'═'*90}")
    
    # === FULL PROJECT COSTS (including foundation, construction, engineering) ===
    shell_mass_t = shell_mass_kg / 1000
    cable_mass_t = total_cable_mass_kg / 1000
    
    full_costs = calculate_all_project_costs(
        tower_height=tower_height,
        base_diameter=eggs[0]['d'],
        n_eggs=n_eggs,
        n_guy_levels=n_guy_levels,
        shell_mass_tonnes=shell_mass_t,
        cable_mass_tonnes=cable_mass_t,
        structure_cost=structure_installed,
        stabilization_cost=stab_cost,
        elevator_cost=elevator_cost,
    )
    
    print_full_project_costs(full_costs, tower_height)


if __name__ == '__main__':
    main()
