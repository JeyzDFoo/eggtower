"""
Visualization functions for egg tower diagrams.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from config import TOTAL_HEIGHT, V_REF
from wind_cables import wind_speed_at_height, wind_force_on_egg


def plot_tapered_tower(eggs, cable_forces=None, guy_system=None, save_path='egg_tower_diagram.png',
                       arm_radius_factor=10.0):
    """
    Draw a 2D side-view diagram of the tapered egg tower (to scale).
    
    Parameters:
        eggs: list of egg dictionaries from build_tower()
        cable_forces: list of cable force dicts from calculate_cable_forces()
        guy_system: dict from calculate_ring_guy_system() for ground-anchored guys
        save_path: filename to save the figure
        arm_radius_factor: how far arms extend (as multiple of egg radius)
    """
    if not eggs:
        print("No eggs to plot!")
        return None, None
    
    d_max = eggs[0]['d']
    total_height = eggs[-1]['z_top']
    
    # Determine plot width based on guy system anchor radius
    anchor_radius = guy_system['anchor_radius'] if guy_system else d_max * 1.5
    
    # Calculate figure dimensions
    plot_width = anchor_radius * 2.5
    plot_height = total_height * 1.05
    tower_aspect = plot_height / plot_width
    
    fig_width = 12
    fig_height = fig_width * tower_aspect
    
    max_fig_height = 20
    if fig_height > max_fig_height:
        fig_height = max_fig_height
        fig_width = fig_height / tower_aspect
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Color scheme
    egg_color = '#2E86AB'
    cable_color = '#E94F37'
    arm_color = '#FF8C00'  # Orange for horizontal arms
    guy_color = '#2ECC71'  # Green for ground-anchored guys
    ring_color = '#9B59B6'  # Purple for rings
    annotation_color = '#1B4332'
    ground_color = '#5C4033'
    
    theta = np.linspace(0, 2 * np.pi, 200)
    
    for i, egg in enumerate(eggs):
        d = egg['d']
        h = egg['h']
        z_base = egg['z_base']
        z_center = z_base + h / 2
        a = h / 2
        b = d / 2
        
        # Draw egg outline
        x_outer = b * np.cos(theta)
        z_outer = z_center + a * np.sin(theta)
        ax.plot(x_outer, z_outer, color=egg_color, linewidth=2, zorder=2)
        ax.fill(x_outer, z_outer, color=egg_color, alpha=0.15, zorder=1)
        
        # Label diameter (first, last, and every 5th egg)
        if i == 0 or i == len(eggs) - 1 or (len(eggs) > 10 and i % 5 == 0) or len(eggs) <= 10:
            ax.annotate(f'Ø{d:.1f}m', 
                       xy=(-b, z_center), 
                       xytext=(-b - d_max * 0.12, z_center),
                       fontsize=7, color=annotation_color, va='center', ha='right',
                       arrowprops=dict(arrowstyle='-', color=annotation_color, lw=0.3),
                       zorder=5)
        
        # Alignment cables with horizontal arms to next egg
        if i < len(eggs) - 1:
            next_egg = eggs[i + 1]
            
            # Arms at egg/egg interface (tip of current egg = base of next egg)
            z_interface = z_base + h  # Top of current egg
            
            # Arm extends from tower axis to arm_radius_factor × egg_radius
            r_egg_curr = b  # Current egg radius
            r_egg_next = next_egg['d'] / 2  # Next egg radius
            
            # Arm extends proportionally to egg radius
            r_arm_curr = r_egg_curr * arm_radius_factor
            r_arm_next = r_egg_next * arm_radius_factor
            
            # Draw horizontal arms at this interface (current egg's top)
            ax.plot([r_egg_curr, r_arm_curr], [z_interface, z_interface],
                   color=arm_color, linewidth=2.5, zorder=3, solid_capstyle='round')
            ax.plot([-r_egg_curr, -r_arm_curr], [z_interface, z_interface],
                   color=arm_color, linewidth=2.5, zorder=3, solid_capstyle='round')
            
            # Arm attachment points on arm ends
            ax.plot([r_arm_curr, -r_arm_curr], [z_interface, z_interface],
                   'o', color=arm_color, markersize=3, zorder=4)
            
            # Cable runs from this arm to the arm at the NEXT interface (top of next egg)
            z_next_interface = next_egg['z_base'] + next_egg['h']
            
            # Draw cable from arm end to next arm end
            ax.plot([r_arm_curr, r_arm_next], [z_interface, z_next_interface],
                   color=cable_color, linestyle='-', linewidth=1.2, zorder=3)
            ax.plot([-r_arm_curr, -r_arm_next], [z_interface, z_next_interface],
                   color=cable_color, linestyle='-', linewidth=1.2, zorder=3)
    
    # Top egg arm (at the very top)
    top_egg = eggs[-1]
    z_top_interface = top_egg['z_base'] + top_egg['h']
    r_top = top_egg['d'] / 2
    r_arm_top = r_top * arm_radius_factor
    ax.plot([r_top, r_arm_top], [z_top_interface, z_top_interface],
           color=arm_color, linewidth=2.5, zorder=3, solid_capstyle='round')
    ax.plot([-r_top, -r_arm_top], [z_top_interface, z_top_interface],
           color=arm_color, linewidth=2.5, zorder=3, solid_capstyle='round')
    ax.plot([r_arm_top, -r_arm_top], [z_top_interface, z_top_interface],
           'o', color=arm_color, markersize=3, zorder=4)
    
    # Ground cables from bottom egg arm
    bottom_egg = eggs[0]
    z_bottom_interface = bottom_egg['z_base'] + bottom_egg['h']  # Top of bottom egg
    r_bottom = bottom_egg['d'] / 2
    r_arm_bottom = r_bottom * arm_radius_factor
    pedestal_spread = r_arm_bottom * 1.2  # Ground anchor outside arm reach
    pedestal_depth = -total_height * 0.015
    
    # Draw arm at bottom of first egg (z_base = 0, so arm is at egg base)
    z_base_arm = bottom_egg['z_base']
    ax.plot([r_bottom, r_arm_bottom], [z_base_arm, z_base_arm],
           color=arm_color, linewidth=2.5, zorder=3, solid_capstyle='round')
    ax.plot([-r_bottom, -r_arm_bottom], [z_base_arm, z_base_arm],
           color=arm_color, linewidth=2.5, zorder=3, solid_capstyle='round')
    
    # Ground cables from bottom arm to foundation
    ax.plot([r_arm_bottom, pedestal_spread], [z_base_arm, pedestal_depth],
           color=cable_color, linestyle='-', linewidth=1.5, zorder=3)
    ax.plot([-r_arm_bottom, -pedestal_spread], [z_base_arm, pedestal_depth],
           color=cable_color, linestyle='-', linewidth=1.5, zorder=3)
    
    if cable_forces:
        ground_cable = next((c for c in cable_forces if c['to_egg'] == 'ground'), None)
        if ground_cable:
            T_MN = ground_cable['T_cable_kN'] / 1000
            ax.text(pedestal_spread + d_max * 0.05, pedestal_depth, 
                   f'Ground: {T_MN:.1f} MN',
                   fontsize=7, color=cable_color, va='center', ha='left',
                   fontweight='bold', zorder=5)
    
    ax.plot([r_arm_bottom, -r_arm_bottom], [z_base_arm, z_base_arm],
           'o', color=arm_color, markersize=3, zorder=4)
    
    # Local egg-to-egg anchors
    anchor_color = '#8B4513'
    ax.plot([pedestal_spread, -pedestal_spread], [pedestal_depth, pedestal_depth],
           's', color=anchor_color, markersize=6, zorder=4)
    
    # Draw ring + ground-anchored guy wire system
    if guy_system and 'levels' in guy_system:
        anchor_r = guy_system['anchor_radius']
        
        for level in guy_system['levels']:
            egg_idx = level['egg_index']
            if egg_idx < len(eggs):
                egg = eggs[egg_idx]
                z_ring = egg['z_base'] + egg['h'] / 2  # Ring at egg center
                r_egg = egg['d'] / 2
                
                # Ring extends outward from egg surface
                ring_extension = r_egg * 0.3
                r_ring = r_egg + ring_extension
                
                # Draw ring as horizontal line segments (2D view)
                ax.plot([-r_ring, -r_egg], [z_ring, z_ring], 
                       color=ring_color, linewidth=3, zorder=3)
                ax.plot([r_egg, r_ring], [z_ring, z_ring], 
                       color=ring_color, linewidth=3, zorder=3)
                
                # Ring attachment points
                ax.plot([r_ring, -r_ring], [z_ring, z_ring],
                       'o', color=ring_color, markersize=4, zorder=4)
                
                # Ground-anchored guy wires
                ax.plot([r_ring, anchor_r], [z_ring, pedestal_depth],
                       color=guy_color, linestyle='-', linewidth=1.8, zorder=2)
                ax.plot([-r_ring, -anchor_r], [z_ring, pedestal_depth],
                       color=guy_color, linestyle='-', linewidth=1.8, zorder=2)
                
                # Label guy angle
                angle = level['theta_deg']
                label_x = r_ring + (anchor_r - r_ring) * 0.3
                label_z = z_ring - (z_ring - pedestal_depth) * 0.3
                ax.text(label_x + 5, label_z, f'{angle:.0f}°',
                       fontsize=6, color=guy_color, va='center', ha='left',
                       fontweight='bold', zorder=5)
        
        # Ground anchor points for guy wires
        ax.plot([anchor_r, -anchor_r], [pedestal_depth, pedestal_depth],
               '^', color=guy_color, markersize=8, zorder=4)
        ax.text(anchor_r, pedestal_depth - total_height * 0.02, 
               f'{anchor_r:.0f}m',
               fontsize=7, color=guy_color, va='top', ha='center', zorder=5)
    
    ground_width = max(d_max * 0.8, pedestal_spread * 1.2, anchor_radius * 1.1)
    ax.fill_between([-ground_width, ground_width], [-total_height*0.025, -total_height*0.025], [0, 0],
                    color=ground_color, alpha=0.5, zorder=0)
    ax.plot([-ground_width, ground_width], [0, 0], color=ground_color, linewidth=3, zorder=1)
    
    pedestal_width = bottom_egg['d'] * 0.4
    ax.fill_between([-pedestal_width, pedestal_width], [-total_height*0.02, -total_height*0.02], [0, 0],
                    color='#666666', alpha=0.8, zorder=1)
    
    # Height annotation on far right
    height_x = anchor_radius * 0.15  # Position relative to plot width
    ax.annotate('', xy=(height_x, total_height), xytext=(height_x, 0),
                arrowprops=dict(arrowstyle='<->', color=annotation_color, lw=1.2))
    ax.text(height_x + anchor_radius * 0.02, total_height/2, f'{total_height:.0f}m',
            fontsize=9, color=annotation_color, va='center', ha='left', fontweight='bold')
    
    # Title
    guy_info = f" + {len(guy_system['levels'])} guy levels" if guy_system else ""
    title = f'Egg Tower: {len(eggs)} eggs, {total_height:.0f}m tall{guy_info}'
    ax.set_title(title, fontsize=11, fontweight='bold', pad=10)
    ax.set_xlabel('Width (m)', fontsize=9)
    ax.set_ylabel('Height (m)', fontsize=9)
    
    ax.set_aspect('equal')
    margin = anchor_radius * 0.15
    ax.set_xlim(-anchor_radius - margin, anchor_radius + margin)
    ax.set_ylim(-total_height*0.04, total_height + total_height * 0.04)
    
    ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Legend
    legend_elements = [
        Line2D([0], [0], color=egg_color, linewidth=2, label='Egg shell'),
        Line2D([0], [0], color=arm_color, linewidth=2.5, label=f'Horizontal arms ({arm_radius_factor:.0f}× radius)'),
        Line2D([0], [0], color=cable_color, linestyle='-', linewidth=1.5, label='Alignment cables'),
    ]
    if guy_system:
        legend_elements.extend([
            Line2D([0], [0], color=ring_color, linewidth=3, label='Ring structures'),
            Line2D([0], [0], color=guy_color, linewidth=1.8, label='Ground-anchored guys'),
        ])
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7, 
              framealpha=0.95, edgecolor='gray')
    
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Tower diagram saved: '{save_path}'")
    plt.show()
    
    return fig, ax


def plot_wind_analysis(eggs, cable_forces, guy_system=None, save_path='wind_analysis.png'):
    """
    Create a dedicated wind force analysis chart showing:
    - Wind speed profile vs height
    - Individual wind force on each egg
    - Cumulative horizontal force
    - Alignment cable tensions
    - Guy cable tensions (if provided)
    """
    from wind_cables import wind_force_on_egg, wind_speed_at_height
    
    total_height = eggs[-1]['z_top']
    n_eggs = len(eggs)
    
    # Calculate wind data for all eggs
    wind_data = []
    for egg in eggs:
        wd = wind_force_on_egg(egg)
        wd['z_center'] = egg['z_base'] + egg['h'] / 2
        wd['d'] = egg['d']
        wind_data.append(wd)
    
    # Calculate cumulative force from top down
    cumulative_forces = []
    cum_force = 0
    for i in range(n_eggs - 1, -1, -1):
        cum_force += wind_data[i]['wind_force']
        cumulative_forces.insert(0, cum_force)
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(14, 8), sharey=True)
    
    heights = [wd['z_center'] for wd in wind_data]
    
    # ═══════════════════════════════════════════════════════════════════════════
    # Panel 1: Wind Speed Profile
    # ═══════════════════════════════════════════════════════════════════════════
    ax1 = axes[0]
    z_profile = np.linspace(10, total_height, 100)
    v_profile = [wind_speed_at_height(z) for z in z_profile]
    
    ax1.fill_betweenx(z_profile, 0, v_profile, alpha=0.3, color='#9B59B6')
    ax1.plot(v_profile, z_profile, color='#9B59B6', linewidth=2)
    
    # Mark wind speed at key heights
    for h in [100, 500, 1000, 1600]:
        if h <= total_height:
            v = wind_speed_at_height(h)
            ax1.plot(v, h, 'o', color='#9B59B6', markersize=6)
            ax1.annotate(f'{v:.0f} m/s', xy=(v, h), xytext=(v + 5, h),
                        fontsize=8, va='center')
    
    ax1.set_xlabel('Wind Speed (m/s)', fontsize=10)
    ax1.set_ylabel('Height (m)', fontsize=10)
    ax1.set_title(f'Wind Speed Profile\nv = {V_REF}·(z/10)^0.16', fontsize=11, fontweight='bold')
    ax1.set_xlim(0, max(v_profile) * 1.3)
    ax1.grid(True, alpha=0.3)
    
    # Add reference wind speed annotation
    ax1.axhline(y=10, color='gray', linestyle='--', linewidth=0.5)
    ax1.text(V_REF, 10, f'  v_ref = {V_REF} m/s at 10m', fontsize=8, va='bottom')
    
    # ═══════════════════════════════════════════════════════════════════════════
    # Panel 2: Wind Force per Egg
    # ═══════════════════════════════════════════════════════════════════════════
    ax2 = axes[1]
    forces_kN = [wd['wind_force'] / 1000 for wd in wind_data]
    
    ax2.barh(heights, forces_kN, height=[egg['h'] * 0.8 for egg in eggs], 
             color='#E74C3C', alpha=0.7, edgecolor='#C0392B')
    
    # Label every 5th bar
    for i, (h, f) in enumerate(zip(heights, forces_kN)):
        if i % 5 == 0 or i == n_eggs - 1:
            ax2.text(f + max(forces_kN) * 0.02, h, f'{f:.0f}', 
                    fontsize=7, va='center', color='#C0392B')
    
    ax2.set_xlabel('Wind Force (kN)', fontsize=10)
    ax2.set_title('Force per Egg\nF = ½ρv²CdA', fontsize=11, fontweight='bold')
    ax2.set_xlim(0, max(forces_kN) * 1.2)
    ax2.grid(True, alpha=0.3, axis='x')
    
    # Total force annotation
    total_force = sum(forces_kN)
    ax2.text(0.95, 0.05, f'Total: {total_force/1000:.2f} MN', 
            transform=ax2.transAxes, fontsize=10, fontweight='bold',
            ha='right', va='bottom', color='#C0392B',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    # ═══════════════════════════════════════════════════════════════════════════
    # Panel 3: Wind Forces & Cable Systems
    # ═══════════════════════════════════════════════════════════════════════════
    ax3 = axes[2]
    cum_forces_MN = [f / 1e6 for f in cumulative_forces]
    
    ax3.fill_betweenx(heights, 0, cum_forces_MN, alpha=0.3, color='#27AE60')
    ax3.plot(cum_forces_MN, heights, color='#27AE60', linewidth=2, label='Wind load (cumulative)')
    
    # Add alignment cable tensions if provided
    if cable_forces:
        cable_heights = []
        cable_tensions = []
        for cf in cable_forces:
            if cf['to_egg'] != 'ground':
                cable_heights.append(cf['z1'])
                cable_tensions.append(cf['T_cable_kN'] / 1000)  # MN
        
        ax3.plot(cable_tensions, cable_heights, 'o-', color='#E94F37', 
                linewidth=1.5, markersize=4, label='Alignment cable tension')
        
        # Calculate cable capacity and strand density using analyze_cable_system
        from config import CABLE_MATERIAL
        from wind_cables import analyze_cable_system, size_cable
        
        results, total_mass, total_strands = analyze_cable_system(cable_forces)
        
        cable_capacities = []
        strand_densities = []
        
        for r in results:
            if r['to_egg'] != 'ground':
                # Strand density (strands per meter circumference)
                strand_densities.append(r['strand_density'])
                
                # Calculate capacity for windward cable (for comparison)
                if r['cables_detail']:
                    windward_strands = r['cables_detail'][0]['n_strands']
                else:
                    windward_strands = r['total_strands_at_level'] // 2
                
                A_rope = np.pi * (CABLE_MATERIAL.fiber_diameter / 2)**2 * CABLE_MATERIAL.fill_factor
                working_load_per_rope = A_rope * CABLE_MATERIAL.tensile_strength / CABLE_MATERIAL.safety_factor
                total_capacity = windward_strands * working_load_per_rope / 1e6  # MN
                cable_capacities.append(total_capacity)
        
        if cable_capacities:
            ax3.plot(cable_capacities, cable_heights, 's--', color='#3498DB', 
                    linewidth=1.5, markersize=4, alpha=0.8, label='Alignment cable capacity')
            
            # Add strand density labels
            for i, (cap, h, density) in enumerate(zip(cable_capacities, cable_heights, strand_densities)):
                # Label every few points to avoid clutter
                if i % 3 == 0 or i == len(cable_capacities) - 1:
                    ax3.annotate(f'{density:.0f}/m', (cap, h), xytext=(8, 0), 
                                textcoords='offset points', fontsize=6, color='#3498DB',
                                ha='left', va='center')
    
    # Add guy cable tensions if provided
    if guy_system and 'levels' in guy_system:
        guy_heights = [level['z_attach'] for level in guy_system['levels']]
        guy_tensions = [level['T_per_guy_kN'] / 1000 for level in guy_system['levels']]  # MN
        
        ax3.plot(guy_tensions, guy_heights, 'D-', color='#9B59B6', 
                linewidth=2, markersize=8, label='Guy cable tension (per cable)')
        
        # Annotate guy points
        for h, t in zip(guy_heights, guy_tensions):
            ax3.annotate(f'{t:.1f}', (t, h), xytext=(-10, 5), 
                        textcoords='offset points', fontsize=7, color='#9B59B6',
                        fontweight='bold', ha='right', va='bottom')
    
    ax3.set_xlabel('Force (MN)', fontsize=10)
    ax3.set_title('Wind Forces &\nCable Systems', fontsize=11, fontweight='bold')
    
    # Calculate x-axis limit including guy tensions
    max_val = max(cum_forces_MN)
    if cable_forces and cable_capacities:
        max_val = max(max_val, max(cable_capacities))
    if guy_system and 'levels' in guy_system:
        max_val = max(max_val, max(guy_tensions))
    ax3.set_xlim(0, max_val * 1.15)
    
    ax3.grid(True, alpha=0.3, axis='x')
    ax3.legend(loc='upper right', fontsize=7)
    
    # Max force annotation
    if cable_forces and cable_tensions:
        max_density = max(strand_densities) if strand_densities else 0
        max_guy = max(guy_tensions) if guy_system and 'levels' in guy_system else 0
        annotation_text = f'Alignment: {max(cable_tensions):.1f} MN\n'
        if max_guy > 0:
            annotation_text += f'Guy: {max_guy:.1f} MN\n'
        annotation_text += f'Density: {max_density:.0f}/m'
        ax3.text(0.95, 0.05, annotation_text,
                transform=ax3.transAxes, fontsize=8, fontweight='bold',
                ha='right', va='bottom', color='#333333',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    # Common settings
    for ax in axes:
        ax.set_ylim(0, total_height * 1.02)
    
    from config import CABLE_MATERIAL
    fiber_dia_mm = CABLE_MATERIAL.fiber_diameter * 1000
    cable_name = CABLE_MATERIAL.name.split('(')[0].strip()  # Get short name
    plt.suptitle(f'Wind Load Analysis - {V_REF} m/s ({V_REF * 3.6:.0f} km/h) | {cable_name} {fiber_dia_mm:.0f}mm',
                fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"✓ Wind analysis saved: '{save_path}'")
    plt.show()
    
    return fig, axes
