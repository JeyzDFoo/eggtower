"""
Wind Energy Module for Egg Tower.

This module contains wind turbine integration analysis for the egg tower,
using the hollow egg shells as wind concentrators (internal duct turbines).
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
