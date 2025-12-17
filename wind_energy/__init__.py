"""
Wind Energy Module for Egg Tower.

This module contains wind energy analysis for the egg tower:

1. Internal Duct Turbines (wind_turbine.py)
   - Uses hollow egg shells as wind concentrators
   - Horizontal turbines in throat of each egg
   
2. Alternative Energy Concepts (alternative_energy.py)
   - Vortex Shedding: Piezoelectric harvesters at junctions
   - Vertical Updraft Column: Solar/wind chimney with horizontal turbines
"""

from .wind_turbine import (
    # Core functions
    load_optimized_tower,
    get_optimized_eggs,
    assess_feasibility,
    print_feasibility_report,
    
    # Internal duct system
    analyze_internal_duct_integration,
    print_internal_duct_analysis,
    calculate_egg_duct_geometry,
    design_internal_duct_system,
    
    # Junction turbines (legacy)
    analyze_tower_turbine_integration,
    print_turbine_analysis,
    
    # Data export/import
    save_wind_energy_analysis,
    load_wind_energy_analysis,
    WIND_ENERGY_DATA_PATH,
    
    # Data classes
    WindTurbine,
    InternalDuctTurbine,
    TURBINE_MODELS,
    INTERNAL_TURBINE_MODELS,
)

__all__ = [
    'load_optimized_tower',
    'get_optimized_eggs',
    'assess_feasibility',
    'print_feasibility_report',
    'analyze_internal_duct_integration',
    'print_internal_duct_analysis',
    'calculate_egg_duct_geometry',
    'design_internal_duct_system',
    'analyze_tower_turbine_integration',
    'print_turbine_analysis',
    'save_wind_energy_analysis',
    'load_wind_energy_analysis',
    'WIND_ENERGY_DATA_PATH',
    'WindTurbine',
    'InternalDuctTurbine',
    'TURBINE_MODELS',
    'INTERNAL_TURBINE_MODELS',
]

# Alternative energy (import separately to avoid circular imports)
from .alternative_energy import (
    analyze_alternative_energy,
    print_alternative_energy_analysis,
    analyze_vortex_energy_harvesting,
    analyze_updraft_column,
    save_alternative_energy_analysis,
)

__all__ += [
    'analyze_alternative_energy',
    'print_alternative_energy_analysis',
    'analyze_vortex_energy_harvesting',
    'analyze_updraft_column',
    'save_alternative_energy_analysis',
]
