"""
Geometry calculations for egg-shaped structural elements.
"""
import numpy as np
from config import ASPECT_RATIO, D_BASE, D_TOP, TOTAL_HEIGHT, SHELL_MATERIAL


def create_egg(d, aspect_ratio=ASPECT_RATIO, thickness_factor=1.0):
    """
    Create a single egg with given diameter.
    
    Parameters:
        d: maximum diameter (m)
        aspect_ratio: height/diameter ratio
        thickness_factor: multiplier for shell thickness (default 1.0)
                         Use >1.0 to increase thickness for higher stress regions
    
    Returns a dictionary with all egg properties.
    """
    h = d * aspect_ratio  # height from diameter
    r = d / 2  # max radius
    
    # Shell thickness (scaled with size, with optional multiplier)
    t_base = (h / 1000) * thickness_factor
    t_tip = (h / 2000) * thickness_factor
    t_avg = (t_base + t_tip) / 2
    
    # Surface area (Knud Thomsen approximation for ellipsoid)
    a = h / 2  # semi-height
    b = d / 2  # semi-width
    p = 1.6075
    surface_area = 4 * np.pi * ((a**p * b**p + a**p * b**p + b**p * b**p) / 3)**(1/p)
    
    # Shell volume and mass
    volume = surface_area * t_avg
    mass = SHELL_MATERIAL.density * volume
    
    # Base cross-sectional area (annular ring)
    r_inner = r - t_base
    base_area = np.pi * (r**2 - r_inner**2)
    
    return {
        'd': d,
        'h': h,
        'r': r,
        't_base': t_base,
        't_avg': t_avg,
        'thickness_factor': thickness_factor,
        'surface_area': surface_area,
        'volume': volume,
        'mass': mass,
        'base_area': base_area
    }


def calculate_diameter_ratio(d_base, d_top, target_height, aspect_ratio=ASPECT_RATIO):
    """
    Find the diameter ratio between consecutive eggs to reach target height.
    Uses binary search to find the ratio that achieves the target height.
    """
    def tower_height(ratio):
        height = 0
        d = d_base
        while d >= d_top:
            h = d * aspect_ratio
            height += h
            d *= ratio
        return height
    
    # Binary search for the right ratio
    low, high = 0.5, 0.99
    while high - low > 1e-6:
        mid = (low + high) / 2
        if tower_height(mid) > target_height:
            high = mid
        else:
            low = mid
    
    return (low + high) / 2


def build_tower(d_base=D_BASE, d_top=D_TOP, target_height=TOTAL_HEIGHT, aspect_ratio=ASPECT_RATIO,
                thickness_factor=1.0):
    """
    Build the tower iteratively from bottom to top.
    
    Parameters:
        d_base: base diameter (m)
        d_top: minimum top diameter (m)
        target_height: target tower height (m)
        aspect_ratio: height/diameter ratio for eggs
        thickness_factor: multiplier for shell thickness (default 1.0)
    
    Returns list of eggs and diameter ratio.
    """
    # Find the diameter ratio that achieves target height
    ratio = calculate_diameter_ratio(d_base, d_top, target_height, aspect_ratio)
    
    eggs = []
    cumulative_height = 0
    d = d_base
    i = 0
    
    while d >= d_top and cumulative_height < target_height:
        egg = create_egg(d, aspect_ratio, thickness_factor=thickness_factor)
        egg['index'] = i
        egg['z_base'] = cumulative_height
        egg['z_top'] = cumulative_height + egg['h']
        
        cumulative_height += egg['h']
        eggs.append(egg)
        
        d *= ratio
        i += 1
    
    return eggs, ratio


def egg_geometry(h, aspect_ratio=5):
    """
    Calculate egg geometry from height (legacy function for compatibility).
    """
    d = h / aspect_ratio
    r = d / 2
    
    # Shell thickness scales with egg size
    t_base = h / 1000
    t_tip = h / 2000
    
    return {
        'h': h,
        'd': d,
        'r': r,
        't_base': t_base,
        't_tip': t_tip,
        't_avg': (t_base + t_tip) / 2
    }


def egg_volume_and_area(h, d, t_avg):
    """
    Calculate egg surface area and shell volume.
    Uses Knud Thomsen approximation for ellipsoid surface area.
    """
    a = h / 2  # semi-axis (height)
    b = d / 2  # semi-axis (width)
    
    # Knud Thomsen approximation
    p = 1.6075
    surface_area = 4 * np.pi * ((a**p * b**p + a**p * b**p + b**p * b**p) / 3)**(1/p)
    
    # Shell volume
    volume = surface_area * t_avg
    
    return {
        'surface_area': surface_area,
        'shell_volume': volume,
        'mass': RHO * volume
    }
