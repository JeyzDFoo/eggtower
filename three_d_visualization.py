"""
3D Visualization for Egg Tower Structure.

This module provides interactive 3D visualization of the egg tower
using Matplotlib's 3D plotting capabilities.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from config import TOTAL_HEIGHT, N_CABLES_PER_LEVEL


def generate_egg_surface(d, h, z_base, n_theta=50, n_phi=25):
    """
    Generate 3D surface mesh for an egg (ellipsoid).
    
    Parameters:
        d: maximum diameter (m)
        h: height (m)
        z_base: z-coordinate of the egg base
        n_theta: number of points around circumference
        n_phi: number of points from pole to pole
        
    Returns:
        X, Y, Z: mesh coordinates for the egg surface
    """
    a = h / 2  # semi-axis vertical
    b = d / 2  # semi-axis horizontal
    z_center = z_base + h / 2
    
    theta = np.linspace(0, 2 * np.pi, n_theta)
    phi = np.linspace(-np.pi / 2, np.pi / 2, n_phi)
    
    theta, phi = np.meshgrid(theta, phi)
    
    X = b * np.cos(phi) * np.cos(theta)
    Y = b * np.cos(phi) * np.sin(theta)
    Z = z_center + a * np.sin(phi)
    
    return X, Y, Z


def generate_egg_wireframe(d, h, z_base, n_lines=8, n_points=50):
    """
    Generate wireframe lines for an egg (lighter weight than full surface).
    
    Parameters:
        d: maximum diameter (m)
        h: height (m)
        z_base: z-coordinate of the egg base
        n_lines: number of meridian lines
        n_points: points per line
        
    Returns:
        List of (x, y, z) line arrays
    """
    a = h / 2
    b = d / 2
    z_center = z_base + h / 2
    
    lines = []
    
    # Meridian lines (vertical)
    for i in range(n_lines):
        theta = 2 * np.pi * i / n_lines
        phi = np.linspace(-np.pi / 2, np.pi / 2, n_points)
        
        x = b * np.cos(phi) * np.cos(theta)
        y = b * np.cos(phi) * np.sin(theta)
        z = z_center + a * np.sin(phi)
        
        lines.append((x, y, z))
    
    # Parallel lines (horizontal circles)
    for i in range(1, n_lines // 2 + 1):
        phi = np.pi * i / (n_lines + 1) - np.pi / 2
        theta = np.linspace(0, 2 * np.pi, n_points)
        
        x = b * np.cos(phi) * np.cos(theta)
        y = b * np.cos(phi) * np.sin(theta)
        z = np.full_like(theta, z_center + a * np.sin(phi))
        
        lines.append((x, y, z))
    
    return lines


def generate_cable_points(egg1, egg2, n_cables=N_CABLES_PER_LEVEL, arm_radius_factor=10.0):
    """
    Generate cable attachment points between two adjacent eggs with horizontal arms.
    
    Cables attach at arm ends, which extend proportionally to egg radius.
    
    Parameters:
        egg1: lower egg dictionary
        egg2: upper egg dictionary
        n_cables: number of cables around circumference
        arm_radius_factor: arms extend to egg_radius × this factor
        
    Returns:
        cables: List of (start_point, end_point) tuples for each cable
        arms1: List of (egg_surface, arm_end) tuples for lower egg arms
        arms2: List of (egg_surface, arm_end) tuples for upper egg arms
    """
    cables = []
    arms1 = []
    arms2 = []
    
    # Attachment at egg/egg interface (tip of egg1 = base of egg2)
    z1 = egg1['z_base'] + egg1['h']  # Top of lower egg
    z2 = egg2['z_base'] + egg2['h']  # Top of upper egg
    
    r1 = egg1['d'] / 2  # Egg radius at interface
    r2 = egg2['d'] / 2
    
    # Arm extends proportionally to egg radius
    r_arm1 = r1 * arm_radius_factor
    r_arm2 = r2 * arm_radius_factor
    
    for i in range(n_cables):
        theta = 2 * np.pi * i / n_cables
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        
        # Egg surface points
        x1_egg = r1 * cos_t
        y1_egg = r1 * sin_t
        x2_egg = r2 * cos_t
        y2_egg = r2 * sin_t
        
        # Arm end points
        x1_arm = r_arm1 * cos_t
        y1_arm = r_arm1 * sin_t
        x2_arm = r_arm2 * cos_t
        y2_arm = r_arm2 * sin_t
        
        # Cable runs from arm1 end to arm2 end
        cables.append(((x1_arm, y1_arm, z1), (x2_arm, y2_arm, z2)))
        
        # Arms (from egg surface to arm end)
        arms1.append(((x1_egg, y1_egg, z1), (x1_arm, y1_arm, z1)))
        arms2.append(((x2_egg, y2_egg, z2), (x2_arm, y2_arm, z2)))
    
    return cables, arms1, arms2


def plot_tower_3d(eggs, guy_system=None, show_cables=True, show_surface=True, 
                  wireframe_only=False, n_cables_display=8, figsize=(12, 16),
                  elevation=15, azimuth=45, save_path=None, high_resolution=False,
                  arm_radius_factor=10.0):
    """
    Create a 3D visualization of the egg tower with horizontal arms.
    
    Parameters:
        eggs: list of egg dictionaries from build_tower()
        guy_system: dict from calculate_ring_guy_system() for ground-anchored guys
        show_cables: whether to display inter-egg cables and arms
        show_surface: whether to show egg surfaces (if False, uses wireframe)
        wireframe_only: use only wireframe for eggs (faster rendering)
        n_cables_display: number of cables to show per level (for clarity)
        figsize: figure size tuple
        elevation: viewing elevation angle
        azimuth: viewing azimuth angle
        save_path: path to save the figure (optional)
        high_resolution: use higher detail for eggs and cables
        arm_radius_factor: arms extend to egg_radius × this factor
        
    Returns:
        fig, ax: matplotlib figure and 3D axis objects
    """
    if not eggs:
        print("No eggs to plot!")
        return None, None
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    
    # Color scheme
    egg_color = '#2E86AB'
    cable_color = '#E94F37'
    arm_color = '#FF8C00'  # Orange for horizontal arms
    guy_color = '#2ECC71'
    ground_color = '#8B4513'
    
    total_height = eggs[-1]['z_top']
    max_radius = eggs[0]['d'] / 2
    max_arm_radius = max_radius * arm_radius_factor  # For plot sizing
    
    # Determine detail level based on number of eggs and resolution setting
    n_eggs = len(eggs)
    if high_resolution:
        # High resolution settings - very smooth egg surfaces
        n_theta, n_phi = 200, 100
        n_wireframe_lines = 32
    elif n_eggs > 50:
        n_theta, n_phi = 40, 20
        n_wireframe_lines = 8
    elif n_eggs > 20:
        n_theta, n_phi = 60, 30
        n_wireframe_lines = 12
    else:
        n_theta, n_phi = 80, 40
        n_wireframe_lines = 16
    
    # Plot each egg
    for i, egg in enumerate(eggs):
        d = egg['d']
        h = egg['h']
        z_base = egg['z_base']
        
        if wireframe_only or not show_surface:
            # Wireframe representation
            lines = generate_egg_wireframe(d, h, z_base, n_lines=n_wireframe_lines)
            for x, y, z in lines:
                ax.plot(x, y, z, color=egg_color, linewidth=0.5, alpha=0.7)
        else:
            # Surface representation
            X, Y, Z = generate_egg_surface(d, h, z_base, n_theta=n_theta, n_phi=n_phi)
            ax.plot_surface(X, Y, Z, color=egg_color, alpha=0.3, 
                          linewidth=0, antialiased=True)
            # Add wireframe overlay for definition
            ax.plot_wireframe(X, Y, Z, color=egg_color, linewidth=0.2, 
                            alpha=0.5, rstride=5, cstride=5)
        
        # Draw cables and arms to next egg
        if show_cables and i < len(eggs) - 1:
            cables, arms1, arms2 = generate_cable_points(
                egg, eggs[i + 1], 
                n_cables=n_cables_display,
                arm_radius_factor=arm_radius_factor
            )
            
            # Draw cables (from arm end to arm end)
            for (x1, y1, z1), (x2, y2, z2) in cables:
                ax.plot([x1, x2], [y1, y2], [z1, z2], 
                       color=cable_color, linewidth=0.8, alpha=0.7)
            
            # Draw arms at lower level (arms1)
            for (x1, y1, z1), (x2, y2, z2) in arms1:
                ax.plot([x1, x2], [y1, y2], [z1, z2], 
                       color=arm_color, linewidth=1.5, alpha=0.9)
            
            # Draw arms at upper level (arms2) - only for the last egg connection
            if i == len(eggs) - 2:
                for (x1, y1, z1), (x2, y2, z2) in arms2:
                    ax.plot([x1, x2], [y1, y2], [z1, z2], 
                           color=arm_color, linewidth=1.5, alpha=0.9)
    
    # Draw ground-anchored guy cables
    if guy_system is not None:
        anchor_radius = guy_system.get('anchor_radius', max_radius * 3)
        n_guys = guy_system.get('n_cables', 8)
        
        # Guy cable attachment points (typically at several levels)
        guy_levels = guy_system.get('attachment_heights', [total_height * 0.25, 
                                                            total_height * 0.5, 
                                                            total_height * 0.75])
        
        for level_height in guy_levels:
            # Find the egg at this height
            attach_egg = None
            for egg in eggs:
                if egg['z_base'] <= level_height <= egg['z_top']:
                    attach_egg = egg
                    break
            
            if attach_egg is None:
                continue
                
            attach_radius = attach_egg['d'] / 2
            
            for j in range(n_guys):
                theta = 2 * np.pi * j / n_guys
                
                # Attachment point on egg
                x_attach = attach_radius * np.cos(theta)
                y_attach = attach_radius * np.sin(theta)
                z_attach = level_height
                
                # Anchor point on ground
                x_anchor = anchor_radius * np.cos(theta)
                y_anchor = anchor_radius * np.sin(theta)
                z_anchor = 0
                
                ax.plot([x_attach, x_anchor], [y_attach, y_anchor], [z_attach, z_anchor],
                       color=guy_color, linewidth=0.3, alpha=0.8)
                
                # Anchor markers
                ax.scatter([x_anchor], [y_anchor], [z_anchor], 
                          color=guy_color, s=10, marker='^')
    
    # Draw ground plane - size to include arms
    ground_size = max(max_arm_radius * 1.1, 
                     guy_system['anchor_radius'] * 1.2 if guy_system else max_arm_radius * 1.1)
    xx, yy = np.meshgrid(np.linspace(-ground_size, ground_size, 10),
                         np.linspace(-ground_size, ground_size, 10))
    zz = np.zeros_like(xx)
    ax.plot_surface(xx, yy, zz, color=ground_color, alpha=0.2)
    
    # Set labels and title
    ax.set_xlabel('X (m)', fontsize=10)
    ax.set_ylabel('Y (m)', fontsize=10)
    ax.set_zlabel('Height (m)', fontsize=10)
    ax.set_title(f'Egg Tower 3D Visualization\n{n_eggs} eggs, {total_height:.0f}m tall, {arm_radius_factor:.0f}× arm extension', 
                fontsize=12, fontweight='bold')
    
    # Set equal aspect ratio for x and y
    ax.set_box_aspect([1, 1, total_height / (ground_size * 2)])
    
    # Set view angle
    ax.view_init(elev=elevation, azim=azimuth)
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=egg_color, linewidth=2, label='Egg shells'),
    ]
    if show_cables:
        legend_elements.extend([
            Line2D([0], [0], color=arm_color, linewidth=2, label=f'Arms ({arm_radius_factor:.0f}× radius)'),
            Line2D([0], [0], color=cable_color, linewidth=1.5, label='Alignment cables'),
        ])
    if guy_system:
        legend_elements.append(Line2D([0], [0], color=guy_color, linewidth=1.5, 
                                      label='Guy cables'))
    
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"Saved 3D visualization to {save_path}")
    
    return fig, ax


def plot_tower_3d_interactive(eggs, guy_system=None, **kwargs):
    """
    Create an interactive 3D visualization (requires interactive backend).
    
    This is a wrapper around plot_tower_3d that ensures the plot
    can be rotated interactively.
    """
    import matplotlib
    # Switch to interactive backend if available
    try:
        matplotlib.use('TkAgg')
    except:
        pass
    
    fig, ax = plot_tower_3d(eggs, guy_system, **kwargs)
    plt.show()
    return fig, ax


def create_animation_frames(eggs, guy_system=None, n_frames=36, 
                           output_dir='frames', **kwargs):
    """
    Create a series of frames for animation (rotating view).
    
    Parameters:
        eggs: list of egg dictionaries
        guy_system: optional guy cable system dict
        n_frames: number of frames (rotation steps)
        output_dir: directory to save frames
        **kwargs: additional arguments passed to plot_tower_3d
        
    Returns:
        List of saved frame paths
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    frame_paths = []
    
    for i in range(n_frames):
        azimuth = 360 * i / n_frames
        
        fig, ax = plot_tower_3d(eggs, guy_system, azimuth=azimuth, 
                               save_path=None, **kwargs)
        
        frame_path = os.path.join(output_dir, f'frame_{i:03d}.png')
        plt.savefig(frame_path, dpi=100, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close(fig)
        
        frame_paths.append(frame_path)
        print(f"Saved frame {i+1}/{n_frames}")
    
    return frame_paths


def plot_stress_distribution_3d(eggs, stresses, figsize=(14, 10), save_path=None):
    """
    Create a 3D visualization with color-coded stress distribution.
    
    Parameters:
        eggs: list of egg dictionaries
        stresses: list of stress values (Pa) for each egg
        figsize: figure size
        save_path: optional path to save the figure
        
    Returns:
        fig, ax: matplotlib figure and axis
    """
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    
    # Normalize stresses for colormap
    stresses_MPa = np.array(stresses) / 1e6
    norm = Normalize(vmin=0, vmax=max(stresses_MPa) * 1.1)
    cmap = plt.cm.RdYlGn_r  # Red = high stress, Green = low stress
    
    total_height = eggs[-1]['z_top']
    max_radius = eggs[0]['d'] / 2
    
    for i, egg in enumerate(eggs):
        d = egg['d']
        h = egg['h']
        z_base = egg['z_base']
        
        # Get color based on stress
        color = cmap(norm(stresses_MPa[i]))
        
        X, Y, Z = generate_egg_surface(d, h, z_base, n_theta=30, n_phi=15)
        ax.plot_surface(X, Y, Z, color=color, alpha=0.7, linewidth=0)
    
    # Add colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=20, pad=0.1)
    cbar.set_label('Compressive Stress (MPa)', fontsize=10)
    
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Height (m)')
    ax.set_title('Egg Tower Stress Distribution', fontsize=12, fontweight='bold')
    
    ax.set_box_aspect([1, 1, total_height / (max_radius * 4)])
    ax.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved stress visualization to {save_path}")
    
    return fig, ax


def plot_cross_section(eggs, height_fraction=0.5, figsize=(10, 10), save_path=None):
    """
    Plot a horizontal cross-section of the tower at a given height.
    
    Parameters:
        eggs: list of egg dictionaries
        height_fraction: fraction of total height for the cut (0 to 1)
        figsize: figure size
        save_path: optional save path
        
    Returns:
        fig, ax: matplotlib figure and axis
    """
    total_height = eggs[-1]['z_top']
    cut_height = total_height * height_fraction
    
    # Find the egg at this height
    cut_egg = None
    for egg in eggs:
        if egg['z_base'] <= cut_height <= egg['z_top']:
            cut_egg = egg
            break
    
    if cut_egg is None:
        print(f"No egg found at height {cut_height:.1f}m")
        return None, None
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Calculate radius at cut height
    d = cut_egg['d']
    h = cut_egg['h']
    z_base = cut_egg['z_base']
    z_center = z_base + h / 2
    
    a = h / 2  # semi-axis vertical
    b = d / 2  # semi-axis horizontal
    
    # At height z, radius = b * sqrt(1 - ((z - z_center)/a)^2)
    z_rel = (cut_height - z_center) / a
    if abs(z_rel) > 1:
        radius = 0
    else:
        radius = b * np.sqrt(1 - z_rel**2)
    
    t = cut_egg['t_base']  # approximate shell thickness
    
    # Draw outer circle
    theta = np.linspace(0, 2 * np.pi, 100)
    x_outer = radius * np.cos(theta)
    y_outer = radius * np.sin(theta)
    ax.plot(x_outer, y_outer, 'b-', linewidth=2, label='Outer surface')
    
    # Draw inner circle
    x_inner = (radius - t) * np.cos(theta)
    y_inner = (radius - t) * np.sin(theta)
    ax.plot(x_inner, y_inner, 'b--', linewidth=1.5, label='Inner surface')
    
    # Fill the shell
    ax.fill_between(x_outer, y_outer, alpha=0.3, color='blue')
    
    ax.set_xlim(-radius * 1.2, radius * 1.2)
    ax.set_ylim(-radius * 1.2, radius * 1.2)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title(f'Cross-section at {cut_height:.1f}m ({height_fraction*100:.0f}% height)\n'
                f'Egg #{cut_egg["index"]+1}, Radius={radius:.2f}m, Shell={t*1000:.1f}mm')
    ax.legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig, ax


# Example usage and demonstration
if __name__ == "__main__":
    # Import tower building function
    from geometry import build_tower
    
    print("Building egg tower...")
    eggs, ratio = build_tower()
    print(f"Tower built with {len(eggs)} eggs, diameter ratio: {ratio:.4f}")
    
    # Create basic 3D visualization
    print("\nGenerating 3D visualization...")
    
    # Simple guy system for demonstration
    guy_system = {
        'anchor_radius': eggs[0]['d'] * 2.5,
        'n_cables': 8,
        'attachment_heights': [400, 800, 1200]
    }
    
    fig, ax = plot_tower_3d(eggs, guy_system=guy_system, 
                            wireframe_only=True,  # Faster rendering
                            n_cables_display=4,
                            save_path='egg_tower_3d.png')
    
    print("3D visualization complete!")
    plt.show()
