"""
Visualization diagrams for alternative energy concepts.

Creates physics diagrams for:
1. Vertical Updraft Column (Solar/Wind Chimney)
2. Vortex Shedding (for completeness, though not viable)
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch, Arc, Polygon, FancyBboxPatch
import numpy as np
import os

# Output directory
OUTPUT_DIR = os.path.dirname(__file__)


def draw_updraft_column_physics():
    """
    Draw a physics diagram of the vertical updraft column with valved accelerators.
    
    Shows:
    - Cross-section of stacked eggs forming a chimney
    - Valved ports at each egg waist (leeward suction)
    - Solar heating at base
    - Wind pressure differential (inlet/outlet)
    - Updraft velocity vectors
    - Horizontal turbines
    - Temperature and pressure annotations
    """
    fig, ax = plt.subplots(1, 1, figsize=(16, 20))
    ax.set_xlim(-10, 10)
    ax.set_ylim(-2, 22)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Title
    ax.text(0, 21.5, 'VERTICAL UPDRAFT COLUMN', fontsize=18, fontweight='bold',
            ha='center', va='top')
    ax.text(0, 20.8, 'with Tapered Duct & Energy Conservation', fontsize=12,
            ha='center', va='top', style='italic', color='#006600')
    
    # === DRAW EGG STACK (simplified cross-section) ===
    n_eggs = 6  # Simplified representation
    egg_heights = [2.8, 2.5, 2.2, 1.9, 1.6, 1.3]
    egg_widths = [3.5, 3.0, 2.5, 2.0, 1.6, 1.2]
    shell_color = '#E8D5B7'
    
    y_pos = 0
    egg_centers = []
    waist_positions = []  # For valve annotations
    throat_positions = []  # For injection ports at junctions
    
    for i, (h, w) in enumerate(zip(egg_heights, egg_widths)):
        # Draw egg shell (ellipse outline)
        egg = patches.Ellipse((0, y_pos + h/2), w, h, 
                              facecolor=shell_color, edgecolor='#8B7355',
                              linewidth=2, alpha=0.8)
        ax.add_patch(egg)
        
        # Internal column (hollow center) - 60% of width
        inner_w = w * 0.6
        inner = patches.Ellipse((0, y_pos + h/2), inner_w, h * 0.85,
                                facecolor='#E6F3FF', edgecolor='none', alpha=0.6)
        ax.add_patch(inner)
        
        # Waist position (maximum diameter at ~40% height)
        waist_y = y_pos + h * 0.4
        waist_positions.append((w/2, waist_y, w))
        
        # Throat position (top of egg - junction with next egg)
        throat_y = y_pos + h * 0.92  # Where eggs connect
        if i < len(egg_heights) - 1:  # Not the top egg
            # Width at throat (narrower than waist)
            throat_w = w * 0.35
            throat_positions.append((throat_w, throat_y))
        
        # Draw SUCTION ports at WAIST (leeward side - right side in diagram)
        valve_x = w/2 - 0.1
        valve_rect = patches.Rectangle((valve_x, waist_y - 0.15), 0.25, 0.3,
                                       facecolor='#00AA00', edgecolor='#006600',
                                       linewidth=1.5, zorder=5)
        ax.add_patch(valve_rect)
        
        # Suction arrow into valve (from outside - air drawn IN)
        ax.annotate('', xy=(valve_x + 0.1, waist_y), 
                   xytext=(valve_x + 0.8, waist_y),
                   arrowprops=dict(arrowstyle='->', color='#0066FF',
                                  lw=1.5, mutation_scale=10))
        
        egg_centers.append(y_pos + h/2)
        y_pos += h * 0.92  # Slight overlap
    
    # Draw INJECTION ports at THROAT (windward side - left side, at junctions)
    for throat_w, throat_y in throat_positions:
        # Injection valve at throat (windward = left side)
        inj_x = -throat_w - 0.15
        inj_rect = patches.Rectangle((inj_x, throat_y - 0.15), 0.25, 0.3,
                                     facecolor='#FF6600', edgecolor='#CC4400',
                                     linewidth=1.5, zorder=5)
        ax.add_patch(inj_rect)
        
        # Injection arrow (high pressure wind pushing IN)
        ax.annotate('', xy=(inj_x + 0.15, throat_y), 
                   xytext=(inj_x - 0.6, throat_y),
                   arrowprops=dict(arrowstyle='->', color='#FF6600',
                                  lw=1.5, mutation_scale=10))
    
    # === DRAW HORIZONTAL TURBINES AT THROAT (max velocity) ===
    # Use throat positions (junctions) - velocity is highest where column is narrowest
    turbine_throat_indices = [1, 3]  # Select every other throat for visual clarity
    for i in turbine_throat_indices:
        if i < len(throat_positions):
            throat_w, ty = throat_positions[i]
            
            # Turbine hub
            hub = plt.Circle((0, ty), 0.15, color='#333', zorder=10)
            ax.add_patch(hub)
            
            # Turbine blades (3-blade rotor, viewed from side = line)
            blade_len = 0.6  # Smaller to fit throat
            for angle in [0, 120, 240]:
                rad = np.radians(angle)
                x1 = blade_len * np.cos(rad)
                y1 = ty + blade_len * np.sin(rad) * 0.3  # Foreshortened
                ax.plot([0, x1], [ty, y1], 'k-', linewidth=3, zorder=9)
            
            # Rotation arrow
            arc = Arc((0, ty), 1.4, 0.5, angle=0, theta1=30, theta2=150,
                      color='#0066CC', linewidth=2, linestyle='--')
            ax.add_patch(arc)
    
    # === UPDRAFT ARROWS (inside column) ===
    arrow_x_positions = [-0.4, 0, 0.4]
    for x in arrow_x_positions:
        for y_start in np.linspace(0.5, 14, 9):
            ax.annotate('', xy=(x, y_start + 1.3), xytext=(x, y_start),
                       arrowprops=dict(arrowstyle='->', color='#FF4444',
                                      lw=2.5, mutation_scale=15))
    
    # Central updraft label (uniform velocity - chimney physics)
    ax.text(0, 8, 'UPDRAFT\n33 m/s\n(uniform)', fontsize=11, ha='center', va='center',
            color='#CC0000', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    # === SOLAR HEATING (base) ===
    # Sun rays entering base
    sun_y = -0.5
    for angle in [-30, -15, 0, 15, 30]:
        rad = np.radians(angle)
        x_start = 2.5 * np.tan(rad)
        ax.annotate('', xy=(x_start * 0.3, 0.3), xytext=(x_start, sun_y),
                   arrowprops=dict(arrowstyle='->', color='#FFD700',
                                  lw=2, mutation_scale=12))
    
    # Solar heating zone
    solar_zone = patches.Rectangle((-1.8, 0), 3.6, 1.5, 
                                   facecolor='#FFEB99', alpha=0.4, zorder=1)
    ax.add_patch(solar_zone)
    ax.text(0, 0.3, 'Solar Heating', fontsize=9, ha='center',
            color='#CC8800', fontweight='bold')
    ax.text(0, -0.3, 'ΔP = 5221 Pa', fontsize=10, ha='center', 
            color='#FF6600', fontweight='bold')
    
    # === WIND FLOW (external) - varying with height ===
    # Left side - Windward (high pressure) - show increasing speed with height
    wind_arrows = [
        (-4, 2, 10, 2.0),   # y=2: ~15 m/s (low)
        (-4, 5, 30, 2.5),   # y=5: ~35 m/s (mid)  
        (-4, 8, 50, 3.0),   # y=8: ~50 m/s (upper-mid)
        (-4, 12, 65, 3.0),  # y=12: ~65 m/s (near top)
    ]
    for x_end, y, v, lw in wind_arrows:
        ax.annotate('', xy=(x_end, y), xytext=(-7, y),
                   arrowprops=dict(arrowstyle='->', color='#0066FF',
                                  lw=lw, mutation_scale=15))
        ax.text(-7.3, y, f'{v}', fontsize=7, ha='right', va='center', color='#0066FF')
    
    ax.text(-8.5, 7, 'WIND\n(m/s)', fontsize=10, ha='center', va='center',
            color='#0066FF', fontweight='bold')
    ax.text(-8.5, 5.5, 'v∝z^0.16', fontsize=8, ha='center', color='#0066FF', style='italic')
    
    # === VALVED ACCELERATOR DETAIL BOX ===
    detail_box = FancyBboxPatch((4.5, 9.5), 5.2, 7, boxstyle='round,pad=0.1',
                                facecolor='#E6FFE6', edgecolor='#006600', linewidth=2)
    ax.add_patch(detail_box)
    
    ax.text(7.1, 16, 'VALVED ACCELERATOR', fontsize=10, fontweight='bold',
           ha='center', color='#006600')
    ax.text(7.1, 15.3, '(Multi-Stage Ejector Pump)', fontsize=8,
           ha='center', color='#006600', style='italic')
    
    # Mini diagram showing TWO eggs with junction
    # Lower egg
    mini_egg1 = patches.Ellipse((7.1, 12.5), 2.4, 1.8, facecolor='#F5E6D3',
                               edgecolor='#8B7355', linewidth=1.5)
    ax.add_patch(mini_egg1)
    
    # Upper egg
    mini_egg2 = patches.Ellipse((7.1, 14.2), 2.0, 1.5, facecolor='#F5E6D3',
                               edgecolor='#8B7355', linewidth=1.5)
    ax.add_patch(mini_egg2)
    
    # Column inside (continuous through both eggs)
    mini_col = patches.Rectangle((6.4, 11.8), 1.4, 3.5, facecolor='#CCE5FF',
                                 edgecolor='none', alpha=0.7)
    ax.add_patch(mini_col)
    
    # THROAT (junction) - mark it
    throat_y = 13.35
    ax.plot([6.0, 8.2], [throat_y, throat_y], 'k--', linewidth=1, alpha=0.5)
    ax.text(8.4, throat_y, 'THROAT\n(junction)', fontsize=6, ha='left', va='center', 
            color='#555', style='italic')
    
    # WAIST (max diameter of lower egg)
    waist_y = 12.2
    ax.plot([5.8, 8.4], [waist_y, waist_y], 'b--', linewidth=1, alpha=0.5)
    ax.text(8.6, waist_y, 'WAIST\n(max Ø)', fontsize=6, ha='left', va='center', 
            color='#0066FF', style='italic')
    
    # INJECTION port at THROAT (windward = left side)
    inj_valve = patches.Rectangle((5.5, throat_y - 0.15), 0.35, 0.3, facecolor='#FF6600',
                                  edgecolor='#CC4400', linewidth=2)
    ax.add_patch(inj_valve)
    
    # Injection arrow (wind pushes IN at stagnation point)
    ax.annotate('', xy=(5.9, throat_y), xytext=(4.9, throat_y),
               arrowprops=dict(arrowstyle='->', color='#FF6600', lw=2))
    ax.text(4.7, throat_y + 0.3, 'INJECT', fontsize=7, ha='right', color='#FF6600', fontweight='bold')
    ax.text(4.7, throat_y - 0.3, 'Cp=+0.9', fontsize=7, ha='right', color='#FF6600')
    
    # SUCTION port at WAIST (leeward = right side)
    suct_valve = patches.Rectangle((8.1, waist_y - 0.15), 0.35, 0.3, facecolor='#00CC00',
                              edgecolor='#006600', linewidth=2)
    ax.add_patch(suct_valve)
    
    # Suction arrow (low pressure draws air IN)
    ax.annotate('', xy=(8.15, waist_y), xytext=(9.2, waist_y),
               arrowprops=dict(arrowstyle='<-', color='#0066FF', lw=2))
    ax.text(9.3, waist_y + 0.3, 'SUCTION', fontsize=7, ha='left', color='#0066FF', fontweight='bold')
    ax.text(9.3, waist_y - 0.3, 'Cp=-0.8', fontsize=7, ha='left', color='#0066FF')
    
    # Updraft arrow in mini diagram
    ax.annotate('', xy=(7.1, 14.8), xytext=(7.1, 11.8),
               arrowprops=dict(arrowstyle='->', color='#FF4444', lw=2.5))
    
    # Labels
    valve_labels = [
        '29 suction @ waist (leeward)',
        '28 injection @ throat (windward)',
        'Valve ΔP = 3733 Pa',
        'v_eq: 33.3 m/s (+37%)',
    ]
    y_label = 11.3
    for label in valve_labels:
        ax.text(7.1, y_label, label, fontsize=7, ha='center', color='#004400')
        y_label -= 0.4
    
    # === PHYSICS EQUATIONS BOX ===
    eq_box = FancyBboxPatch((-9, 10), 4.5, 6.5, boxstyle='round,pad=0.1',
                            facecolor='white', edgecolor='#333', linewidth=1)
    ax.add_patch(eq_box)
    
    equations = [
        ('CHIMNEY PHYSICS:', None, True),
        ('', None, False),
        ('Driving ΔP:', '#CC8800', True),
        ('5221 Pa total', '#CC8800', False),
        ('', None, False),
        ('Equilibrium:', '#0066FF', True),
        ('ΔP = K×½ρv²', '#0066FF', False),
        ('v = 33 m/s', '#0066FF', False),
        ('', None, False),
        ('Power budget:', '#006600', True),
        ('2158 kW available', '#006600', False),
        ('1360 kW friction', '#006600', False),
        ('359 kW extracted', '#006600', False),
        ('16.6% efficiency', '#CC0000', False),
    ]
    
    y_eq = 16
    for text, color, is_header in equations:
        if text:
            weight = 'bold' if is_header else 'normal'
            c = '#333' if color is None else color
            ax.text(-6.75, y_eq, text, fontsize=8 if is_header else 7,
                   ha='center', va='center', color=c, fontweight=weight,
                   family='monospace' if not is_header else None)
        y_eq -= 0.42
    
    # === TURBINE SPECS BOX ===
    spec_box = FancyBboxPatch((-9, 3), 4, 5.5, boxstyle='round,pad=0.1',
                              facecolor='#F0FFF0', edgecolor='#228B22', linewidth=1)
    ax.add_patch(spec_box)
    
    specs = [
        'SYSTEM OUTPUT',
        '─────────────',
        '28 x 2m turbines',
        '359 kW total',
        '1102 MWh/year',
        '',
        'CAPEX: $6.4M',
        'Payback: 117 yr',
    ]
    y_spec = 8
    for i, text in enumerate(specs):
        weight = 'bold' if i == 0 else 'normal'
        ax.text(-7, y_spec, text, fontsize=9, ha='center', va='center',
               fontweight=weight, color='#006400')
        y_spec -= 0.6
    
    # === DIMENSION ANNOTATIONS ===
    # Column diameter at base (tapered duct)
    ax.annotate('', xy=(-2.1 * 0.6, 1.4), xytext=(2.1 * 0.6, 1.4),
               arrowprops=dict(arrowstyle='<->', color='#666', lw=1.5))
    ax.text(0, 1.8, 'Tapered: 6m→1.5m', fontsize=8, ha='center', color='#666')
    
    # Height arrow
    ax.annotate('', xy=(4, 0), xytext=(4, 16),
               arrowprops=dict(arrowstyle='<->', color='#666', lw=1.5))
    ax.text(4.3, 8, '1600m', fontsize=9, ha='left', va='center',
           color='#666', rotation=90)
    
    # === LEGEND ===
    legend_y = 0
    ax.add_patch(patches.Rectangle((5, legend_y + 2.4), 0.4, 0.3, 
                                   facecolor='#00AA00', edgecolor='#006600'))
    ax.text(5.6, legend_y + 2.55, 'Suction port (waist, leeward)', fontsize=8, va='center')
    
    ax.add_patch(patches.Rectangle((5, legend_y + 1.8), 0.4, 0.3, 
                                   facecolor='#FF6600', edgecolor='#CC4400'))
    ax.text(5.6, legend_y + 1.95, 'Injection port (throat, windward)', fontsize=8, va='center')
    
    ax.add_patch(patches.Rectangle((5, legend_y + 1.2), 0.4, 0.3,
                                   facecolor='#E6F3FF', edgecolor='#0066CC'))
    ax.text(5.6, legend_y + 1.35, 'Internal updraft column', fontsize=8, va='center')
    
    ax.annotate('', xy=(5.3, legend_y + 0.75), xytext=(5, legend_y + 0.55),
               arrowprops=dict(arrowstyle='->', color='#FF4444', lw=2))
    ax.text(5.6, legend_y + 0.65, 'Updraft flow', fontsize=8, va='center')
    
    ax.annotate('', xy=(5.3, legend_y + 0.15), xytext=(5, legend_y + 0.15),
               arrowprops=dict(arrowstyle='->', color='#0066FF', lw=2))
    ax.text(5.6, legend_y + 0.15, 'External wind / suction', fontsize=8, va='center')
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(OUTPUT_DIR, 'updraft_column_physics.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved: {output_path}")
    return output_path


def draw_combined_energy_diagram():
    """
    Draw a combined diagram comparing all energy generation concepts.
    Shows Option A (Internal Duct) vs Option B (Updraft + Valves) as mutually exclusive.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 10))
    
    # === 1. INTERNAL DUCT TURBINES (Option A) ===
    ax1 = axes[0]
    ax1.set_xlim(-3, 3)
    ax1.set_ylim(-1, 6)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title('OPTION A: Internal Duct\n(Horizontal Flow)', fontsize=12, fontweight='bold',
                  color='#006600')
    
    # Draw egg cross-section
    egg = patches.Ellipse((0, 3), 4, 5, facecolor='#E8D5B7', 
                          edgecolor='#8B7355', linewidth=2)
    ax1.add_patch(egg)
    
    # Throat/duct
    throat = patches.Rectangle((-0.4, 1.5), 0.8, 3, facecolor='#CCE5FF',
                               edgecolor='#0066CC', linewidth=1)
    ax1.add_patch(throat)
    
    # Horizontal wind arrows
    for y in [2.5, 3, 3.5]:
        ax1.annotate('', xy=(0.3, y), xytext=(-2.5, y),
                    arrowprops=dict(arrowstyle='->', color='#0066FF', lw=2))
        ax1.annotate('', xy=(2.5, y), xytext=(0.3, y),
                    arrowprops=dict(arrowstyle='->', color='#0066FF', lw=2))
    
    # Turbine
    hub = plt.Circle((0, 3), 0.15, color='#333', zorder=10)
    ax1.add_patch(hub)
    for angle in [0, 90, 180, 270]:
        rad = np.radians(angle)
        ax1.plot([0, 0.5*np.cos(rad)], [3, 3 + 0.5*np.sin(rad)], 
                'k-', linewidth=3)
    
    ax1.text(0, -0.5, 'Wind flows horizontally\nthrough egg throat', 
             ha='center', fontsize=9, style='italic')
    ax1.text(0, 0.3, '295 kW | 904 MWh/yr', ha='center', fontsize=10,
             fontweight='bold', color='#006600')
    ax1.text(0, -0.1, 'CAPEX $1.34M | 29.6 yr', ha='center', fontsize=8, color='#666')
    
    # === 2. VERTICAL UPDRAFT + VALVED ACCELERATORS (Option B) ===
    ax2 = axes[1]
    ax2.set_xlim(-3, 3)
    ax2.set_ylim(-1, 6)
    ax2.set_aspect('equal')
    ax2.axis('off')
    ax2.set_title('OPTION B: Updraft + Valves\n(Vertical Chimney)', fontsize=12, fontweight='bold',
                  color='#0066CC')
    
    # Draw egg stack (3 simplified)
    y_pos = 0.5
    for h, w in [(2, 3), (1.5, 2.2), (1, 1.5)]:
        egg = patches.Ellipse((0, y_pos + h/2), w, h, facecolor='#E8D5B7',
                              edgecolor='#8B7355', linewidth=2)
        ax2.add_patch(egg)
        inner = patches.Ellipse((0, y_pos + h/2), w*0.6, h*0.8,
                                facecolor='#FFE6E6', edgecolor='none', alpha=0.6)
        ax2.add_patch(inner)
        
        # Valve ports at waist - BOTH sides for dual-port
        waist_y = y_pos + h * 0.4
        # Right (leeward) - suction
        valve_r = patches.Rectangle((w/2 - 0.1, waist_y - 0.08), 0.15, 0.16,
                                    facecolor='#00CC00', edgecolor='#006600', linewidth=1)
        ax2.add_patch(valve_r)
        ax2.annotate('', xy=(w/2, waist_y), xytext=(w/2 + 0.5, waist_y),
                    arrowprops=dict(arrowstyle='<-', color='#0066FF', lw=1))
        # Left (windward) - injection
        valve_l = patches.Rectangle((-w/2 - 0.05, waist_y - 0.08), 0.15, 0.16,
                                    facecolor='#00CC00', edgecolor='#006600', linewidth=1)
        ax2.add_patch(valve_l)
        ax2.annotate('', xy=(-w/2, waist_y), xytext=(-w/2 - 0.5, waist_y),
                    arrowprops=dict(arrowstyle='->', color='#FF6600', lw=1))
        
        y_pos += h * 0.9
    
    # Vertical updraft arrows
    for y in np.linspace(1, 4.5, 5):
        ax2.annotate('', xy=(0, y + 0.6), xytext=(0, y),
                    arrowprops=dict(arrowstyle='->', color='#FF4444', lw=2.5))
    
    # Sun at base
    ax2.text(0, 0.2, '☀️', fontsize=20, ha='center')
    
    # Horizontal turbine
    ax2.plot([-0.6, 0.6], [3, 3], 'k-', linewidth=4)
    hub = plt.Circle((0, 3), 0.1, color='#333', zorder=10)
    ax2.add_patch(hub)
    
    # Valve acceleration note
    ax2.text(2.3, 3, '29 dual\nvalves', fontsize=7, ha='left', color='#006600')
    
    ax2.text(0, -0.5, 'Chimney physics: uniform\n33 m/s throughout',
             ha='center', fontsize=9, style='italic')
    ax2.text(0, 0.3, '359 kW | 1102 MWh/yr', ha='center', fontsize=10,
             fontweight='bold', color='#0066CC')
    ax2.text(0, -0.1, 'CAPEX $6.4M | 117 yr', ha='center', fontsize=8, color='#666')
    
    # === 3. VORTEX SHEDDING ===
    ax3 = axes[2]
    ax3.set_xlim(-3, 3)
    ax3.set_ylim(-1, 6)
    ax3.set_aspect('equal')
    ax3.axis('off')
    ax3.set_title('Vortex Shedding\n(Not Viable)', fontsize=12, fontweight='bold', color='#999')
    
    # Draw oscillating egg
    for offset, alpha in [(0, 0.8), (-0.3, 0.3), (0.3, 0.3)]:
        egg = patches.Ellipse((offset, 3), 3, 4, facecolor='#E8D5B7',
                              edgecolor='#8B7355', linewidth=2, alpha=alpha)
        ax3.add_patch(egg)
    
    # Oscillation arrow
    ax3.annotate('', xy=(-1, 3), xytext=(1, 3),
                arrowprops=dict(arrowstyle='<->', color='#666', lw=2))
    ax3.text(0, 2, '~50mm', ha='center', fontsize=9, color='#666')
    
    # Wind causing vortices
    ax3.annotate('', xy=(-2, 4), xytext=(-3, 4),
                arrowprops=dict(arrowstyle='->', color='#0066FF', lw=2))
    
    # Vortex spirals (simplified)
    for x, y, r in [(1.8, 4.5, 0.3), (2, 3.5, 0.25), (1.9, 2.5, 0.3)]:
        vortex = patches.Arc((x, y), r*2, r*2, angle=0, theta1=0, theta2=270,
                            color='#0066FF', linewidth=1.5, linestyle='--')
        ax3.add_patch(vortex)
    
    ax3.text(0, -0.5, 'Oscillation too small\nfor useful power',
             ha='center', fontsize=9, style='italic', color='#999')
    ax3.text(0, 0.3, '2 kW | 5 MWh/yr', ha='center', fontsize=10,
             fontweight='bold', color='#999')
    ax3.text(0, -0.1, 'NOT VIABLE', ha='center', fontsize=8, color='#CC0000')
    
    # Add "MUTUALLY EXCLUSIVE" annotation between Option A and B
    fig.text(0.37, 0.93, '← MUTUALLY EXCLUSIVE →', ha='center', fontsize=10,
             fontweight='bold', color='#CC6600', style='italic')
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(OUTPUT_DIR, 'energy_concepts_comparison.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved: {output_path}")
    return output_path


def generate_all_diagrams():
    """Generate all alternative energy visualization diagrams."""
    print("\n🎨 Generating Alternative Energy Diagrams...")
    print("=" * 50)
    
    paths = []
    paths.append(draw_updraft_column_physics())
    paths.append(draw_combined_energy_diagram())
    
    print("=" * 50)
    print(f"✅ Generated {len(paths)} diagrams")
    return paths


if __name__ == '__main__':
    generate_all_diagrams()
