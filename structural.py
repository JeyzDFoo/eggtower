"""
Structural analysis functions for egg tower (stress, buckling).
"""
import numpy as np
from config import ASPECT_RATIO, SHELL_MATERIAL


def shell_buckling_stress(R, t, material=None, knockdown=0.25):
    """
    Calculate critical buckling stress for a curved shell under compression.
    
    Uses ShellMaterial.shell_buckling_stress() method.
    
    Parameters:
        R: radius of curvature (m)
        t: shell thickness (m)
        material: ShellMaterial instance (defaults to SHELL_MATERIAL)
        knockdown: imperfection factor (default 0.25)
    
    Returns:
        Critical buckling stress (Pa)
    """
    if material is None:
        material = SHELL_MATERIAL
    return material.shell_buckling_stress(R, t, knockdown)


def euler_buckling_load(E, I, L_eff):
    """
    Euler critical buckling load for a column.
    
    P_cr = π² × E × I / L_eff²
    
    Parameters:
        E: Young's modulus (Pa)
        I: Second moment of area (m⁴)
        L_eff: Effective length (m)
    
    Returns:
        Critical buckling load (N)
    """
    return (np.pi**2 * E * I) / (L_eff**2)


def egg_buckling_analysis(d, h, t, weight_above, aspect_ratio=ASPECT_RATIO):
    """
    Analyze buckling for a single egg shell.
    
    Returns dict with buckling results.
    """
    r = d / 2
    a = h / 2  # semi-major (vertical)
    b = d / 2  # semi-minor (horizontal)
    
    # Radius of curvature at the pole and equator
    R_pole = b**2 / a
    R_equator = a**2 / b
    R_critical = min(R_pole, R_equator)
    
    # Shell buckling stress
    sigma_cr_shell = shell_buckling_stress(R_critical, t)
    
    # Actual compressive stress from weight above
    base_area = np.pi * (r**2 - (r - t)**2)
    sigma_actual = weight_above / base_area if base_area > 0 else 0
    
    # Second moment of area for thin-walled tube
    I = np.pi * r**3 * t
    
    # Euler buckling (K=1.0 for pinned-pinned with guy cables)
    L_eff = h * 1.0
    P_cr_euler = euler_buckling_load(SHELL_MATERIAL.youngs_modulus, I, L_eff)
    sigma_cr_euler = P_cr_euler / base_area if base_area > 0 else np.inf
    
    # Safety factors
    SF_shell = sigma_cr_shell / sigma_actual if sigma_actual > 0 else np.inf
    SF_euler = sigma_cr_euler / sigma_actual if sigma_actual > 0 else np.inf
    SF_combined = min(SF_shell, SF_euler)
    
    return {
        'R_pole': R_pole,
        'R_equator': R_equator,
        'R_critical': R_critical,
        'sigma_cr_shell_MPa': sigma_cr_shell / 1e6,
        'sigma_cr_euler_MPa': sigma_cr_euler / 1e6,
        'sigma_actual_MPa': sigma_actual / 1e6,
        'SF_shell_buckling': SF_shell,
        'SF_euler_buckling': SF_euler,
        'SF_buckling_combined': SF_combined,
        'governing_mode': 'shell' if SF_shell < SF_euler else 'euler'
    }


def analyze_tower(eggs, guy_system=None):
    """
    Analyze structural properties of the tower.
    Calculate stress at each level considering cumulative weight above
    AND vertical forces from guy wires.
    
    Parameters:
        eggs: list of egg dictionaries from build_tower()
        guy_system: optional dict from calculate_ring_guy_system()
                   If provided, includes guy wire vertical forces in stress
    """
    n = len(eggs)
    g = 9.81
    
    # Pre-calculate guy wire vertical forces at each egg level
    guy_vertical_above = [0.0] * n  # Cumulative guy force from levels above each egg
    
    if guy_system and 'levels' in guy_system:
        for level in guy_system['levels']:
            egg_idx = level['egg_index']
            T_kN = level['T_per_guy_kN']
            T = T_kN * 1000  # Convert to N
            theta = np.radians(level['theta_deg'])
            n_guys = level['n_guys']
            
            # Vertical component = T × cos(θ) × n_guys
            F_vertical = T * np.cos(theta) * n_guys
            
            # This force is added to all eggs BELOW this guy level
            for i in range(egg_idx):
                guy_vertical_above[i] += F_vertical
    
    results = []
    for i, egg in enumerate(eggs):
        # Weight of all eggs above this one
        weight_above = sum(e['mass'] for e in eggs[i+1:]) * g
        
        # Add guy wire vertical forces from levels above
        total_force_above = weight_above + guy_vertical_above[i]
        
        # Stress at base of this egg (including guy forces)
        sigma = total_force_above / egg['base_area'] if egg['base_area'] > 0 else 0
        safety_factor = SHELL_MATERIAL.compressive_strength / sigma if sigma > 0 else np.inf
        
        # Also compute stress from self-weight only (for comparison)
        sigma_self_weight = weight_above / egg['base_area'] if egg['base_area'] > 0 else 0
        
        # Buckling analysis (uses total force)
        buckling = egg_buckling_analysis(
            egg['d'], egg['h'], egg['t_base'], total_force_above, ASPECT_RATIO
        )
        
        results.append({
            **egg,
            'weight_above': weight_above,
            'guy_force_above': guy_vertical_above[i],
            'total_force_above': total_force_above,
            'stress_Pa': sigma,
            'stress_MPa': sigma / 1e6,
            'stress_self_weight_MPa': sigma_self_weight / 1e6,
            'safety_factor': safety_factor,
            **buckling
        })
    
    return results
