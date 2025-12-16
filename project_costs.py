"""
Foundation, Construction, and Engineering Costs for Egg Tower.

This module estimates the "soft" costs that go beyond materials and systems:
1. Foundation - massive structure to support 12,000+ tonnes
2. Construction - cranes, scaffolding, labor, site work
3. Engineering - design, analysis, testing, certification
4. Project management and contingency
"""
import numpy as np
from dataclasses import dataclass
from typing import Dict


def calculate_foundation_cost(tower_mass_tonnes: float, 
                               tower_height: float,
                               base_diameter: float,
                               guy_anchor_radius: float = 400,
                               n_guy_levels: int = 9,
                               soil_type: str = 'rock') -> Dict:
    """
    Estimate foundation costs for tower and guy wire anchors.
    
    The tower needs:
    1. Central foundation under the base egg
    2. Guy wire anchor foundations (9 anchors at 400m radius)
    
    Parameters:
        tower_mass_tonnes: Total tower mass
        tower_height: Tower height in meters
        base_diameter: Base egg diameter
        guy_anchor_radius: Distance to guy anchors
        n_guy_levels: Number of guy wire levels
        soil_type: 'rock', 'dense_soil', or 'soft_soil'
    """
    
    # Soil bearing capacity and costs vary by type
    soil_params = {
        'rock': {'bearing_MPa': 10, 'excavation_per_m3': 200, 'pile_needed': False},
        'dense_soil': {'bearing_MPa': 0.5, 'excavation_per_m3': 80, 'pile_needed': True},
        'soft_soil': {'bearing_MPa': 0.15, 'excavation_per_m3': 50, 'pile_needed': True},
    }
    soil = soil_params[soil_type]
    
    # === CENTRAL FOUNDATION ===
    # Must support tower weight + guy wire vertical components
    # Guy wires add significant downward force (we calculated ~17x self-weight!)
    total_vertical_load_kN = tower_mass_tonnes * 9.81 * 18  # ~18x for guy forces
    total_vertical_load_MN = total_vertical_load_kN / 1000
    
    # Required foundation area
    required_area_m2 = total_vertical_load_MN / soil['bearing_MPa']
    
    # Foundation is a circular mat with diameter slightly larger than base
    foundation_diameter = max(base_diameter * 1.5, np.sqrt(required_area_m2 * 4 / np.pi))
    foundation_area = np.pi * (foundation_diameter / 2) ** 2
    
    # Foundation depth (rule of thumb: 1/10 to 1/20 of supported height for tall structures)
    foundation_depth = max(10, tower_height / 15)  # At least 10m, up to ~100m for 1600m tower
    
    # Concrete volume
    concrete_volume_m3 = foundation_area * foundation_depth
    
    # Concrete cost ($150-300/m3 for high-strength concrete, placed)
    concrete_cost_per_m3 = 250
    concrete_cost = concrete_volume_m3 * concrete_cost_per_m3
    
    # Reinforcement (typically 100-150 kg steel per m3 concrete)
    rebar_kg_per_m3 = 120
    rebar_mass = concrete_volume_m3 * rebar_kg_per_m3
    rebar_cost_per_kg = 2.0  # Installed
    rebar_cost = rebar_mass * rebar_cost_per_kg
    
    # Excavation
    excavation_volume = concrete_volume_m3 * 1.3  # Over-excavate for working room
    excavation_cost = excavation_volume * soil['excavation_per_m3']
    
    # Piling if needed (for soft soils)
    pile_cost = 0
    if soil['pile_needed']:
        # Deep piles to bedrock or bearing stratum
        pile_depth = 50  # meters
        pile_diameter = 1.5  # meters
        n_piles = int(required_area_m2 / 10)  # ~1 pile per 10 m2
        pile_cost_each = pile_depth * pile_diameter * 3000  # $/m of pile
        pile_cost = n_piles * pile_cost_each
    
    central_foundation_cost = concrete_cost + rebar_cost + excavation_cost + pile_cost
    
    # === GUY WIRE ANCHOR FOUNDATIONS ===
    # 9 anchor points around the tower, each must resist massive tension
    # Peak guy tension is ~60 MN per cable, with 3 cables per anchor point
    max_anchor_load_MN = 60 * 3  # 180 MN per anchor
    
    # Anchors are typically dead-man type (buried concrete block) or rock anchors
    # For rock: drilled and grouted anchors
    # For soil: massive concrete blocks
    
    if soil_type == 'rock':
        # Rock anchors: drilled holes with tensioned cables grouted in
        anchor_depth = 30  # meters into rock
        n_anchor_tendons = 20  # per anchor point
        anchor_cost_each = anchor_depth * n_anchor_tendons * 500  # $500/m per tendon
    else:
        # Dead-man anchors: buried concrete blocks
        # Size based on soil friction and weight
        anchor_volume = max_anchor_load_MN / (soil['bearing_MPa'] * 0.5)  # Conservative
        anchor_cost_each = anchor_volume * concrete_cost_per_m3 * 1.5  # With rebar
    
    n_anchors = 9  # One per guy cable direction
    guy_anchor_cost = n_anchors * anchor_cost_each
    
    # Access roads to anchors (400m radius, 9 locations)
    road_length = guy_anchor_radius * 1.5 * n_anchors  # Approximate
    road_cost_per_m = 500  # Basic construction road
    access_road_cost = road_length * road_cost_per_m
    
    total_foundation_cost = central_foundation_cost + guy_anchor_cost + access_road_cost
    
    return {
        'central_foundation': {
            'diameter_m': foundation_diameter,
            'depth_m': foundation_depth,
            'concrete_m3': concrete_volume_m3,
            'cost': central_foundation_cost,
        },
        'guy_anchors': {
            'n_anchors': n_anchors,
            'cost_each': anchor_cost_each,
            'cost': guy_anchor_cost,
        },
        'access_roads': {
            'length_m': road_length,
            'cost': access_road_cost,
        },
        'total_cost': total_foundation_cost,
    }


def calculate_construction_cost(tower_height: float,
                                  n_eggs: int,
                                  shell_mass_tonnes: float,
                                  cable_mass_tonnes: float) -> Dict:
    """
    Estimate construction costs (labor, equipment, site work).
    
    Construction approach:
    - Eggs fabricated off-site in sections, transported
    - Assembled on-site using climbing cranes
    - Guy cables installed as tower rises
    """
    
    # === CRANE SYSTEM ===
    # Need heavy-lift cranes that can reach 1600m
    # Likely a climbing crane system that rises with the tower
    # Plus ground-based mobile cranes for lower sections
    
    # Climbing crane (custom-built for this project)
    climbing_crane_cost = 25_000_000  # $25M for specialized system
    n_climbing_cranes = 2  # Redundancy
    
    # Mobile cranes for ground work
    mobile_crane_months = 36  # 3 years of crane rental
    mobile_crane_monthly = 150_000
    mobile_crane_cost = mobile_crane_months * mobile_crane_monthly * 4  # 4 cranes
    
    crane_cost = (climbing_crane_cost * n_climbing_cranes) + mobile_crane_cost
    
    # === SCAFFOLDING / ACCESS ===
    # Temporary platforms, safety systems
    # Cost scales with height and number of work levels
    scaffolding_per_level = 200_000
    scaffolding_cost = n_eggs * scaffolding_per_level
    
    # === LABOR ===
    # Estimated construction duration: 4-6 years
    construction_months = 60  # 5 years
    
    # Crew size varies by phase
    avg_workers = 200
    worker_monthly_cost = 8_000  # Fully loaded (salary + benefits + insurance)
    
    # Specialized high-altitude workers (premium)
    specialist_workers = 50
    specialist_monthly = 15_000
    
    labor_cost = (construction_months * 
                  (avg_workers * worker_monthly_cost + 
                   specialist_workers * specialist_monthly))
    
    # === SITE FACILITIES ===
    # Offices, workshops, storage, worker facilities
    site_facilities_cost = 15_000_000
    
    # === TRANSPORTATION ===
    # Moving 12,000+ tonnes of materials to site
    # Assume remote location requiring special logistics
    transport_per_tonne = 500  # $/tonne average
    total_mass = shell_mass_tonnes + cable_mass_tonnes
    transport_cost = total_mass * transport_per_tonne
    
    # === TEMPORARY WORKS ===
    # Formwork, shoring, temporary bracing during construction
    temporary_works = shell_mass_tonnes * 100  # $100 per tonne of shell
    
    # === QUALITY CONTROL ===
    # Testing, inspection, NDT on every joint
    qc_per_egg = 100_000
    qc_cost = n_eggs * qc_per_egg
    
    # === WEATHER DELAYS ===
    # Contingency for weather (high altitude = more wind days)
    weather_contingency = labor_cost * 0.15
    
    total_construction_cost = (crane_cost + scaffolding_cost + labor_cost + 
                               site_facilities_cost + transport_cost + 
                               temporary_works + qc_cost + weather_contingency)
    
    return {
        'cranes': crane_cost,
        'scaffolding': scaffolding_cost,
        'labor': labor_cost,
        'site_facilities': site_facilities_cost,
        'transport': transport_cost,
        'temporary_works': temporary_works,
        'quality_control': qc_cost,
        'weather_contingency': weather_contingency,
        'construction_months': construction_months,
        'total_cost': total_construction_cost,
    }


def calculate_engineering_cost(tower_height: float,
                                project_cost_estimate: float) -> Dict:
    """
    Estimate engineering, design, and project management costs.
    
    For a first-of-its-kind megastructure, engineering is substantial.
    """
    
    # === DESIGN & ENGINEERING ===
    # Structural engineering, wind engineering, geotechnical
    # Rule of thumb: 5-10% of construction cost for complex structures
    design_pct = 0.08
    design_cost = project_cost_estimate * design_pct
    
    # === WIND TUNNEL TESTING ===
    # Physical scale models tested in wind tunnels
    wind_tunnel_cost = 5_000_000
    
    # === PROTOTYPE TESTING ===
    # Full-scale section testing, material qualification
    prototype_cost = 10_000_000
    
    # === CFD / FEA ANALYSIS ===
    # Computational modeling (software, compute time, specialists)
    computational_cost = 3_000_000
    
    # === REGULATORY & CERTIFICATION ===
    # Novel structure = extensive review process
    # Environmental impact, safety certification, permits
    regulatory_cost = 8_000_000
    
    # === PROJECT MANAGEMENT ===
    # 5-year project needs substantial PM
    pm_monthly = 500_000  # PM team fully loaded
    pm_months = 72  # 6 years including pre-construction
    pm_cost = pm_monthly * pm_months
    
    # === INSURANCE ===
    # Construction insurance for megaproject
    # Typically 1-2% of project value
    insurance_cost = project_cost_estimate * 0.015
    
    # === LEGAL & CONTRACTS ===
    legal_cost = 5_000_000
    
    # === CONTINGENCY ===
    # First-of-its-kind = higher risk
    engineering_contingency = (design_cost + wind_tunnel_cost + prototype_cost) * 0.25
    
    total_engineering_cost = (design_cost + wind_tunnel_cost + prototype_cost +
                              computational_cost + regulatory_cost + pm_cost +
                              insurance_cost + legal_cost + engineering_contingency)
    
    return {
        'design_engineering': design_cost,
        'wind_tunnel': wind_tunnel_cost,
        'prototype_testing': prototype_cost,
        'computational': computational_cost,
        'regulatory': regulatory_cost,
        'project_management': pm_cost,
        'insurance': insurance_cost,
        'legal': legal_cost,
        'contingency': engineering_contingency,
        'total_cost': total_engineering_cost,
    }


def calculate_all_project_costs(tower_height: float,
                                  base_diameter: float,
                                  n_eggs: int,
                                  n_guy_levels: int,
                                  shell_mass_tonnes: float,
                                  cable_mass_tonnes: float,
                                  structure_cost: float,
                                  stabilization_cost: float,
                                  elevator_cost: float,
                                  guy_anchor_radius: float = 400) -> Dict:
    """
    Calculate complete project costs including all phases.
    """
    
    tower_mass = shell_mass_tonnes + cable_mass_tonnes
    
    # Materials and systems (already calculated)
    materials_systems = structure_cost + stabilization_cost + elevator_cost
    
    # Foundation
    foundation = calculate_foundation_cost(
        tower_mass_tonnes=tower_mass,
        tower_height=tower_height,
        base_diameter=base_diameter,
        guy_anchor_radius=guy_anchor_radius,
        n_guy_levels=n_guy_levels
    )
    
    # Construction
    construction = calculate_construction_cost(
        tower_height=tower_height,
        n_eggs=n_eggs,
        shell_mass_tonnes=shell_mass_tonnes,
        cable_mass_tonnes=cable_mass_tonnes
    )
    
    # Engineering (based on preliminary estimate)
    preliminary_total = materials_systems + foundation['total_cost'] + construction['total_cost']
    engineering = calculate_engineering_cost(
        tower_height=tower_height,
        project_cost_estimate=preliminary_total
    )
    
    # Grand total
    grand_total = (materials_systems + foundation['total_cost'] + 
                   construction['total_cost'] + engineering['total_cost'])
    
    # Project contingency (additional 10% on top of everything)
    project_contingency = grand_total * 0.10
    grand_total_with_contingency = grand_total + project_contingency
    
    return {
        'materials_systems': {
            'structure': structure_cost,
            'stabilization': stabilization_cost,
            'elevator': elevator_cost,
            'subtotal': materials_systems,
        },
        'foundation': foundation,
        'construction': construction,
        'engineering': engineering,
        'subtotals': {
            'materials_systems': materials_systems,
            'foundation': foundation['total_cost'],
            'construction': construction['total_cost'],
            'engineering': engineering['total_cost'],
        },
        'grand_total': grand_total,
        'project_contingency': project_contingency,
        'grand_total_with_contingency': grand_total_with_contingency,
    }


def print_full_project_costs(costs: Dict, tower_height: float):
    """Print comprehensive project cost breakdown."""
    
    print(f"\n{'═'*90}")
    print(f"COMPLETE PROJECT COST ANALYSIS - {tower_height:.0f}m EGG TOWER")
    print(f"{'═'*90}")
    
    ms = costs['materials_systems']
    print(f"\n┌{'─'*88}┐")
    print(f"│{'MATERIALS & SYSTEMS':^88}│")
    print(f"├{'─'*60}┬{'─'*27}┤")
    print(f"│ {'Structure (shells + cables, installed)':<58} │ ${ms['structure']/1e6:>18.1f}M │")
    print(f"│ {'Active Stabilization System':<58} │ ${ms['stabilization']/1e6:>18.1f}M │")
    print(f"│ {'Elevator System':<58} │ ${ms['elevator']/1e6:>18.1f}M │")
    print(f"├{'─'*60}┼{'─'*27}┤")
    print(f"│ {'Subtotal Materials & Systems':<58} │ ${ms['subtotal']/1e6:>18.1f}M │")
    print(f"└{'─'*60}┴{'─'*27}┘")
    
    fd = costs['foundation']
    print(f"\n┌{'─'*88}┐")
    print(f"│{'FOUNDATION':^88}│")
    print(f"├{'─'*60}┬{'─'*27}┤")
    cf = fd['central_foundation']
    print(f"│ {'Central foundation ({:.0f}m dia × {:.0f}m deep)':<58} │ ${cf['cost']/1e6:>18.1f}M │".format(
        cf['diameter_m'], cf['depth_m']))
    print(f"│ {'Guy wire anchors ({:d} locations)':<58} │ ${fd['guy_anchors']['cost']/1e6:>18.1f}M │".format(
        fd['guy_anchors']['n_anchors']))
    print(f"│ {'Access roads':<58} │ ${fd['access_roads']['cost']/1e6:>18.1f}M │")
    print(f"├{'─'*60}┼{'─'*27}┤")
    print(f"│ {'Subtotal Foundation':<58} │ ${fd['total_cost']/1e6:>18.1f}M │")
    print(f"└{'─'*60}┴{'─'*27}┘")
    
    cn = costs['construction']
    print(f"\n┌{'─'*88}┐")
    print(f"│{'CONSTRUCTION ({:.0f} months = {:.1f} years)':^88}│".format(
        cn['construction_months'], cn['construction_months']/12))
    print(f"├{'─'*60}┬{'─'*27}┤")
    print(f"│ {'Cranes (climbing + mobile)':<58} │ ${cn['cranes']/1e6:>18.1f}M │")
    print(f"│ {'Scaffolding & access':<58} │ ${cn['scaffolding']/1e6:>18.1f}M │")
    print(f"│ {'Labor':<58} │ ${cn['labor']/1e6:>18.1f}M │")
    print(f"│ {'Site facilities':<58} │ ${cn['site_facilities']/1e6:>18.1f}M │")
    print(f"│ {'Transportation':<58} │ ${cn['transport']/1e6:>18.1f}M │")
    print(f"│ {'Temporary works':<58} │ ${cn['temporary_works']/1e6:>18.1f}M │")
    print(f"│ {'Quality control':<58} │ ${cn['quality_control']/1e6:>18.1f}M │")
    print(f"│ {'Weather contingency':<58} │ ${cn['weather_contingency']/1e6:>18.1f}M │")
    print(f"├{'─'*60}┼{'─'*27}┤")
    print(f"│ {'Subtotal Construction':<58} │ ${cn['total_cost']/1e6:>18.1f}M │")
    print(f"└{'─'*60}┴{'─'*27}┘")
    
    en = costs['engineering']
    print(f"\n┌{'─'*88}┐")
    print(f"│{'ENGINEERING & PROJECT MANAGEMENT':^88}│")
    print(f"├{'─'*60}┬{'─'*27}┤")
    print(f"│ {'Design & structural engineering':<58} │ ${en['design_engineering']/1e6:>18.1f}M │")
    print(f"│ {'Wind tunnel testing':<58} │ ${en['wind_tunnel']/1e6:>18.1f}M │")
    print(f"│ {'Prototype & material testing':<58} │ ${en['prototype_testing']/1e6:>18.1f}M │")
    print(f"│ {'Computational analysis (CFD/FEA)':<58} │ ${en['computational']/1e6:>18.1f}M │")
    print(f"│ {'Regulatory & certification':<58} │ ${en['regulatory']/1e6:>18.1f}M │")
    print(f"│ {'Project management':<58} │ ${en['project_management']/1e6:>18.1f}M │")
    print(f"│ {'Insurance':<58} │ ${en['insurance']/1e6:>18.1f}M │")
    print(f"│ {'Legal & contracts':<58} │ ${en['legal']/1e6:>18.1f}M │")
    print(f"│ {'Engineering contingency':<58} │ ${en['contingency']/1e6:>18.1f}M │")
    print(f"├{'─'*60}┼{'─'*27}┤")
    print(f"│ {'Subtotal Engineering':<58} │ ${en['total_cost']/1e6:>18.1f}M │")
    print(f"└{'─'*60}┴{'─'*27}┘")
    
    print(f"\n┌{'─'*88}┐")
    print(f"│{'PROJECT TOTALS':^88}│")
    print(f"├{'─'*60}┬{'─'*27}┤")
    print(f"│ {'Materials & Systems':<58} │ ${costs['subtotals']['materials_systems']/1e6:>18.1f}M │")
    print(f"│ {'Foundation':<58} │ ${costs['subtotals']['foundation']/1e6:>18.1f}M │")
    print(f"│ {'Construction':<58} │ ${costs['subtotals']['construction']/1e6:>18.1f}M │")
    print(f"│ {'Engineering':<58} │ ${costs['subtotals']['engineering']/1e6:>18.1f}M │")
    print(f"├{'─'*60}┼{'─'*27}┤")
    print(f"│ {'SUBTOTAL':<58} │ ${costs['grand_total']/1e6:>18.1f}M │")
    print(f"│ {'Project contingency (10%)':<58} │ ${costs['project_contingency']/1e6:>18.1f}M │")
    print(f"├{'─'*60}┼{'─'*27}┤")
    print(f"│ {'GRAND TOTAL':<58} │ ${costs['grand_total_with_contingency']/1e6:>18.1f}M │")
    print(f"└{'─'*60}┴{'─'*27}┘")
    
    # Key ratios
    gt = costs['grand_total_with_contingency']
    print(f"\n  KEY METRICS:")
    print(f"    Cost per meter of height:  ${gt/tower_height:,.0f}/m")
    print(f"    Materials & Systems:       {costs['subtotals']['materials_systems']/gt*100:.1f}% of total")
    print(f"    Foundation:                {costs['subtotals']['foundation']/gt*100:.1f}% of total")
    print(f"    Construction:              {costs['subtotals']['construction']/gt*100:.1f}% of total")
    print(f"    Engineering:               {costs['subtotals']['engineering']/gt*100:.1f}% of total")
    
    print(f"\n{'═'*90}")


if __name__ == "__main__":
    # Test with 1600m tower parameters
    costs = calculate_all_project_costs(
        tower_height=1600,
        base_diameter=20.4,
        n_eggs=29,
        n_guy_levels=9,
        shell_mass_tonnes=7858,
        cable_mass_tonnes=3417,
        structure_cost=440.8e6,
        stabilization_cost=23.9e6,
        elevator_cost=53.5e6,
    )
    
    print_full_project_costs(costs, tower_height=1600)
