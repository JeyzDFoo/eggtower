"""
Visualization for Internal Duct Turbine System.

Generates diagrams showing how the egg shell acts as a wind concentrator
with turbines at the core.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Circle, Wedge, Polygon
from matplotlib.lines import Line2D
import os

# Output path
OUTPUT_DIR = os.path.dirname(__file__)


def draw_internal_duct_diagram(
    egg_diameter: float = 20.0,
    throat_ratio: float = 0.15,
    save_path: str = None,
    figsize: tuple = (14, 10)
):
    """
    Draw a cross-section diagram showing how the internal duct turbine works.
    
    Parameters:
        egg_diameter: Outer diameter of the egg (m)
        throat_ratio: Throat diameter as fraction of internal diameter
        save_path: Path to save the figure (default: internal_duct_diagram.png)
    """
    if save_path is None:
        save_path = os.path.join(OUTPUT_DIR, 'internal_duct_diagram.png')
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # === GEOMETRY ===
    shell_thickness = egg_diameter / 50  # Approximate shell thickness
    d_internal = egg_diameter - 2 * shell_thickness
    d_throat = d_internal * throat_ratio
    
    r_outer = egg_diameter / 2
    r_inner = d_internal / 2
    r_throat = d_throat / 2
    
    # Center of diagram
    cx, cy = 0, 0
    
    # === DRAW EGG SHELL (cross-section as ellipse) ===
    # Outer shell
    theta = np.linspace(0, 2 * np.pi, 100)
    aspect = 1.4  # Egg is taller than wide
    
    x_outer = r_outer * np.cos(theta)
    y_outer = r_outer * aspect * np.sin(theta)
    
    x_inner = r_inner * np.cos(theta)
    y_inner = r_inner * aspect * np.sin(theta)
    
    # Fill shell (between outer and inner)
    ax.fill(x_outer, y_outer, color='#8B7355', alpha=0.6, label='BFRP Shell')
    ax.fill(x_inner, y_inner, color='#E8E8E8', alpha=1.0)  # Interior (light gray)
    
    # Shell outline
    ax.plot(x_outer, y_outer, 'k-', linewidth=2)
    ax.plot(x_inner, y_inner, color='black', linewidth=1, linestyle='--', alpha=0.5)
    
    # === DRAW INLET OPENINGS (windward side - left) ===
    inlet_angle = 60  # degrees, opening arc
    inlet_y_top = r_inner * aspect * 0.5
    inlet_y_bot = -r_inner * aspect * 0.5
    
    # Draw inlet slots on left side of shell
    inlet_width = shell_thickness * 1.5
    for y_pos in [inlet_y_top * 0.6, 0, inlet_y_bot * 0.6]:
        inlet_rect = plt.Rectangle(
            (-r_outer - inlet_width/2, y_pos - r_outer * 0.08),
            inlet_width * 2,
            r_outer * 0.16,
            color='#4A90D9',
            alpha=0.8,
            zorder=5
        )
        ax.add_patch(inlet_rect)
    
    # === DRAW TURBINE AT THROAT ===
    turbine_r = r_throat * 0.8
    
    # Turbine hub
    hub = Circle((cx, cy), turbine_r * 0.25, color='#2C3E50', zorder=10)
    ax.add_patch(hub)
    
    # Turbine blades (3 blades)
    n_blades = 5
    blade_color = '#3498DB'
    for i in range(n_blades):
        angle = i * (360 / n_blades) + 15  # Offset for visual appeal
        angle_rad = np.radians(angle)
        
        # Blade shape (tapered)
        blade_length = turbine_r * 0.9
        blade_width_base = turbine_r * 0.15
        blade_width_tip = turbine_r * 0.05
        
        # Create blade polygon
        cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
        cos_perp, sin_perp = np.cos(angle_rad + np.pi/2), np.sin(angle_rad + np.pi/2)
        
        hub_r = turbine_r * 0.25
        blade_points = [
            (cx + hub_r * cos_a + blade_width_base/2 * cos_perp,
             cy + hub_r * sin_a + blade_width_base/2 * sin_perp),
            (cx + hub_r * cos_a - blade_width_base/2 * cos_perp,
             cy + hub_r * sin_a - blade_width_base/2 * sin_perp),
            (cx + blade_length * cos_a - blade_width_tip/2 * cos_perp,
             cy + blade_length * sin_a - blade_width_tip/2 * sin_perp),
            (cx + blade_length * cos_a + blade_width_tip/2 * cos_perp,
             cy + blade_length * sin_a + blade_width_tip/2 * sin_perp),
        ]
        blade = Polygon(blade_points, color=blade_color, ec='#2980B9', linewidth=1, zorder=9)
        ax.add_patch(blade)
    
    # Rotation arrow around turbine
    rot_arrow = mpatches.FancyArrowPatch(
        (cx + turbine_r * 1.1, cy + turbine_r * 0.3),
        (cx + turbine_r * 0.3, cy + turbine_r * 1.1),
        connectionstyle="arc3,rad=0.4",
        arrowstyle='->,head_length=6,head_width=4',
        color='#E74C3C',
        linewidth=2,
        zorder=15
    )
    ax.add_patch(rot_arrow)
    ax.text(cx + turbine_r * 1.3, cy + turbine_r * 0.8, 'ω', fontsize=14, 
            color='#E74C3C', fontweight='bold')
    
    # === DRAW FLOW ARROWS ===
    arrow_color = '#27AE60'
    
    # Incoming wind (left side)
    for y_offset in [-2, 0, 2]:
        ax.annotate('', 
            xy=(-r_outer * 0.95, y_offset),
            xytext=(-r_outer * 1.8, y_offset),
            arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2.5))
    
    # Wind entering inlets
    for y_pos in [inlet_y_top * 0.6, 0, inlet_y_bot * 0.6]:
        ax.annotate('',
            xy=(-r_inner * 0.7, y_pos),
            xytext=(-r_outer * 0.9, y_pos),
            arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2))
    
    # Convergent flow toward center
    for y_pos in [inlet_y_top * 0.4, 0, inlet_y_bot * 0.4]:
        ax.annotate('',
            xy=(-r_throat * 1.2, y_pos * 0.3),
            xytext=(-r_inner * 0.5, y_pos),
            arrowprops=dict(arrowstyle='->', color='#2ECC71', lw=2.5))
    
    # Flow through turbine (accelerated)
    ax.annotate('',
        xy=(r_throat * 1.5, 0),
        xytext=(-r_throat * 1.5, 0),
        arrowprops=dict(arrowstyle='->', color='#F39C12', lw=4, 
                       connectionstyle='arc3,rad=0'))
    
    # Divergent exhaust flow
    for y_pos in [inlet_y_top * 0.4, 0, inlet_y_bot * 0.4]:
        ax.annotate('',
            xy=(r_inner * 0.6, y_pos),
            xytext=(r_throat * 1.2, y_pos * 0.3),
            arrowprops=dict(arrowstyle='->', color='#E67E22', lw=2))
    
    # Exit flow (right side)
    for y_pos in [inlet_y_top * 0.5, 0, inlet_y_bot * 0.5]:
        ax.annotate('',
            xy=(r_outer * 1.5, y_pos),
            xytext=(r_inner * 0.7, y_pos),
            arrowprops=dict(arrowstyle='->', color='#E67E22', lw=2))
    
    # === DRAW THROAT BOUNDARY (dashed circle) ===
    throat_circle = Circle((cx, cy), r_throat, fill=False, 
                           linestyle=':', color='#9B59B6', linewidth=2)
    ax.add_patch(throat_circle)
    
    # === LABELS ===
    # Title
    ax.text(0, r_outer * aspect + 2, 
            'INTERNAL DUCT TURBINE SYSTEM', 
            ha='center', va='bottom', fontsize=16, fontweight='bold')
    ax.text(0, r_outer * aspect + 0.5, 
            'Egg Shell as Wind Concentrator (Plan View Cross-Section)', 
            ha='center', va='bottom', fontsize=11, style='italic', color='gray')
    
    # Wind direction
    ax.text(-r_outer * 1.8, 3.5, 'WIND', fontsize=12, fontweight='bold', color=arrow_color)
    ax.text(-r_outer * 1.8, 2.5, f'V = 20 m/s', fontsize=10, color=arrow_color)
    
    # Inlet label
    ax.text(-r_outer * 1.1, inlet_y_top * 0.85, 'Inlets\n(25% of\nshell area)', 
            fontsize=9, ha='center', color='#4A90D9', fontweight='bold')
    
    # Convergent section
    ax.text(-r_inner * 0.35, r_inner * aspect * 0.65, 'CONVERGENT\nSECTION', 
            fontsize=9, ha='center', color='#2ECC71', fontweight='bold')
    
    # Throat label
    ax.text(0, -r_throat * 1.8, 'THROAT\n(Turbine Location)', 
            fontsize=10, ha='center', color='#9B59B6', fontweight='bold')
    ax.text(0, -r_throat * 2.8, f'V = {20 * 2.5:.0f} m/s\n(2.5× amplified)', 
            fontsize=9, ha='center', color='#F39C12')
    
    # Divergent section
    ax.text(r_inner * 0.35, r_inner * aspect * 0.65, 'DIVERGENT\nSECTION', 
            fontsize=9, ha='center', color='#E67E22', fontweight='bold')
    
    # Exhaust label
    ax.text(r_outer * 1.4, inlet_y_top * 0.85, 'Exhaust\n(Leeward\nSide)', 
            fontsize=9, ha='center', color='#E67E22', fontweight='bold')
    
    # Shell label
    ax.text(r_outer * 0.7, r_outer * aspect * 0.85, 'BFRP\nShell', 
            fontsize=9, ha='center', color='#5D4E37', fontweight='bold')
    
    # Turbine label
    ax.text(turbine_r * 1.5, -turbine_r * 0.5, 'Turbine\nRotor', 
            fontsize=9, ha='left', color='#2C3E50', fontweight='bold')
    
    # === PHYSICS BOX ===
    physics_text = (
        "PHYSICS:\n"
        "━━━━━━━━━━━━━━━━━━━━━━\n"
        "Continuity: A₁V₁ = A₂V₂\n"
        "V_throat = V_inlet × √(A_inlet/A_throat)\n"
        "\n"
        "Power ∝ V³\n"
        "2.5× velocity → 15.6× power!\n"
        "━━━━━━━━━━━━━━━━━━━━━━\n"
        f"Throat: {throat_ratio*100:.0f}% of diameter\n"
        f"Efficiency: 70% (duct losses)"
    )
    props = dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', 
                 edgecolor='#BDC3C7', alpha=0.95)
    ax.text(r_outer * 1.6, -r_outer * aspect * 0.5, physics_text,
            fontsize=9, fontfamily='monospace', ha='left', va='top',
            bbox=props)
    
    # === LEGEND ===
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#8B7355', 
               markersize=12, label='BFRP Shell'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#4A90D9', 
               markersize=10, label='Inlet Openings'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498DB', 
               markersize=10, label='Turbine Rotor'),
        Line2D([0], [0], color='#27AE60', linewidth=3, label='Incoming Wind'),
        Line2D([0], [0], color='#F39C12', linewidth=3, label='Accelerated Flow'),
        Line2D([0], [0], color='#E67E22', linewidth=3, label='Exhaust Flow'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', 
              framealpha=0.95, fontsize=9)
    
    # === FORMATTING ===
    ax.set_xlim(-r_outer * 2.2, r_outer * 2.8)
    ax.set_ylim(-r_outer * aspect * 1.3, r_outer * aspect * 1.4)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    
    print(f"💨 Saved internal duct diagram to: {save_path}")
    return save_path


def draw_egg_cutaway_3d(save_path: str = None, figsize: tuple = (12, 10)):
    """
    Draw a 3D-style cutaway view of the egg showing the turbine inside.
    """
    if save_path is None:
        save_path = os.path.join(OUTPUT_DIR, 'egg_cutaway_3d.png')
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Egg parameters
    egg_width = 8
    egg_height = 12
    shell_t = 0.3
    
    # Draw egg outline (ellipse)
    theta = np.linspace(0, 2 * np.pi, 100)
    
    # Outer egg shape (slightly egg-shaped, narrower at top)
    def egg_shape(t, w, h, narrow_factor=0.85):
        """Create egg-like shape (narrower at top)"""
        x = w * np.cos(t)
        y = h * np.sin(t)
        # Make top narrower
        narrow = 1 - (1 - narrow_factor) * (np.sin(t) + 1) / 2
        x = x * narrow
        return x, y
    
    x_out, y_out = egg_shape(theta, egg_width/2, egg_height/2, 0.8)
    x_in, y_in = egg_shape(theta, egg_width/2 - shell_t, egg_height/2 - shell_t, 0.8)
    
    # Draw back half of shell (solid)
    ax.fill(x_out, y_out, color='#A0522D', alpha=0.7)
    ax.fill(x_in, y_in, color='#F5DEB3', alpha=0.9)  # Interior visible
    
    # Draw cutaway section (front quarter removed)
    cutaway_mask = (theta > np.pi * 0.1) & (theta < np.pi * 0.9)
    x_cut = x_in.copy()
    y_cut = y_in.copy()
    x_cut[cutaway_mask] = x_out[cutaway_mask]  # Extend to outer edge
    
    # Turbine at center
    turbine_y = 0
    turbine_r = 1.2
    
    # Draw turbine disk (horizontal, seen from angle)
    ellipse_w = turbine_r * 2
    ellipse_h = turbine_r * 0.5  # Foreshortened
    
    turbine_disk = mpatches.Ellipse((0, turbine_y), ellipse_w, ellipse_h,
                                     color='#3498DB', alpha=0.9, zorder=10)
    ax.add_patch(turbine_disk)
    
    # Hub
    hub = mpatches.Ellipse((0, turbine_y), 0.4, 0.15,
                            color='#2C3E50', zorder=11)
    ax.add_patch(hub)
    
    # Flow arrows going through
    for y_off in [-1.5, 0, 1.5]:
        ax.annotate('',
            xy=(3, turbine_y + y_off * 0.3),
            xytext=(-3, turbine_y + y_off * 0.3),
            arrowprops=dict(arrowstyle='->', color='#27AE60', lw=2.5),
            zorder=5)
    
    # Inlet indicators on left
    ax.text(-egg_width/2 - 0.5, 2, '→ INLET', fontsize=10, color='#27AE60', fontweight='bold')
    ax.text(-egg_width/2 - 0.5, -2, '→ INLET', fontsize=10, color='#27AE60', fontweight='bold')
    
    # Exhaust indicator on right  
    ax.text(egg_width/2 + 0.3, 0, 'EXHAUST →', fontsize=10, color='#E67E22', fontweight='bold')
    
    # Labels
    ax.text(0, egg_height/2 + 1, 'EGG SHELL CROSS-SECTION', 
            ha='center', fontsize=14, fontweight='bold')
    ax.text(0, egg_height/2 + 0.3, 'with Internal Duct Turbine', 
            ha='center', fontsize=11, style='italic', color='gray')
    
    ax.text(0, turbine_y - 1.5, 'Turbine Rotor', ha='center', fontsize=10, 
            color='#2C3E50', fontweight='bold')
    
    ax.text(2, 4, 'BFRP\nShell', ha='center', fontsize=10, color='#8B4513')
    
    # Outline
    ax.plot(x_out, y_out, 'k-', linewidth=2)
    
    ax.set_xlim(-7, 7)
    ax.set_ylim(-8, 9)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    
    print(f"🥚 Saved egg cutaway to: {save_path}")
    return save_path


def generate_all_diagrams():
    """Generate all wind energy visualizations."""
    paths = []
    paths.append(draw_internal_duct_diagram())
    paths.append(draw_egg_cutaway_3d())
    return paths


if __name__ == '__main__':
    generate_all_diagrams()
