"""
Self-Climbing Robot Capsule Elevator System for Egg Tower.

Instead of traditional cable elevators (problematic at 1600m due to cable mass
and stretch), this system uses self-propelled capsules that climb along 
fixed guide rails mounted to the tower shell.

Advantages:
- No heavy cables that stretch and add load
- Multiple capsules can operate independently on same track
- Can have express and local service zones
- Capsules can be lightweight composite construction
- Regenerative braking recovers energy on descent
- System scales well with height

Technology basis: Similar to ThyssenKrupp MULTI, but with rack-and-pinion
or linear motor propulsion instead of maglev for simplicity and reliability.
"""
import numpy as np
from dataclasses import dataclass
from typing import List, Dict


@dataclass
class ClimbingCapsule:
    """Self-propelled elevator capsule."""
    name: str
    capacity_persons: int         # Passenger capacity
    cabin_mass_kg: float          # Empty cabin mass
    drive_power_kW: float         # Climbing motor power
    max_speed_mps: float          # Maximum speed (m/s)
    cost_usd: float               # Unit cost per capsule


@dataclass
class GuideRail:
    """Vertical guide rail system."""
    name: str
    mass_per_meter_kg: float      # Rail mass per meter
    cost_per_meter_usd: float     # Cost per meter including installation
    max_capsule_mass_kg: float    # Maximum capsule mass supported


# Capsule options
CAPSULES = {
    'standard': ClimbingCapsule(
        name="Tourist Capsule (15 person)",
        capacity_persons=15,
        cabin_mass_kg=1200,       # Lightweight composite cabin with panoramic windows
        drive_power_kW=60,        # ~5 m/s climb rate (leisurely)
        max_speed_mps=5,          # 18 km/h - scenic pace
        cost_usd=400_000          # $400K each - simpler, slower
    ),
    'express': ClimbingCapsule(
        name="Express Capsule (15 person)",
        capacity_persons=15,
        cabin_mass_kg=1400,
        drive_power_kW=100,       # ~8 m/s
        max_speed_mps=8,          # 29 km/h - for staff/return trips
        cost_usd=600_000          # $600K each
    ),
    'service': ClimbingCapsule(
        name="Service/Cargo Capsule",
        capacity_persons=4,
        cabin_mass_kg=1000,
        drive_power_kW=80,        # For cargo
        max_speed_mps=6,
        cost_usd=350_000
    ),
    'emergency': ClimbingCapsule(
        name="Emergency Rescue Capsule",
        capacity_persons=20,
        cabin_mass_kg=1000,
        drive_power_kW=120,       # Fast descent capability
        max_speed_mps=10,
        cost_usd=500_000
    ),
}

# Guide rail options
GUIDE_RAILS = {
    'rack_pinion': GuideRail(
        name="Rack-and-Pinion Rail",
        mass_per_meter_kg=80,     # Steel rack + guide rails
        cost_per_meter_usd=2000,  # Including installation
        max_capsule_mass_kg=5000
    ),
    'linear_motor': GuideRail(
        name="Linear Motor Track",
        mass_per_meter_kg=120,    # Motor stator + rails
        cost_per_meter_usd=5000,  # More expensive but smoother
        max_capsule_mass_kg=8000
    ),
    'friction_drive': GuideRail(
        name="Friction Drive Rail",
        mass_per_meter_kg=60,     # Lighter, simpler
        cost_per_meter_usd=1500,
        max_capsule_mass_kg=3000
    ),
}


def design_elevator_system(tower_height: float, 
                           n_shafts: int = 4,
                           rail_type: str = 'rack_pinion',
                           passenger_throughput_per_hour: int = 2000) -> Dict:
    """
    Design complete self-climbing elevator system.
    
    Parameters:
        tower_height: Total tower height in meters
        n_shafts: Number of independent elevator shafts
        rail_type: Type of guide rail system
        passenger_throughput_per_hour: Target passengers per hour
    
    Returns:
        Dictionary with system design and cost breakdown
    """
    rail = GUIDE_RAILS[rail_type]
    
    # === ZONE CONFIGURATION ===
    # Split tower into zones for efficiency (like supertall buildings)
    # Express zones + local zones
    
    n_zones = max(2, int(tower_height / 400))  # ~400m per zone
    zone_height = tower_height / n_zones
    
    # Sky lobbies at zone transitions
    n_sky_lobbies = n_zones - 1
    
    # === CAPSULE FLEET ===
    # Calculate required capsules based on throughput
    
    # Round trip time (with stops and loading)
    avg_speed = 8  # m/s average including accel/decel
    round_trip_time_s = (2 * tower_height / avg_speed) + (n_zones * 30)  # 30s per stop
    round_trip_time_min = round_trip_time_s / 60
    
    # Capsules needed per shaft
    capsule = CAPSULES['standard']
    trips_per_capsule_per_hour = 60 / round_trip_time_min
    passengers_per_capsule_hour = trips_per_capsule_per_hour * capsule.capacity_persons * 0.7  # 70% fill
    
    capsules_needed = int(np.ceil(passenger_throughput_per_hour / passengers_per_capsule_hour / n_shafts))
    capsules_per_shaft = max(3, capsules_needed)  # Minimum 3 per shaft
    
    total_standard_capsules = capsules_per_shaft * n_shafts
    
    # Add express capsules (1 per 2 shafts)
    n_express = max(2, n_shafts // 2)
    
    # Service capsules (1 per shaft)
    n_service = n_shafts
    
    # Emergency capsules (minimum 4)
    n_emergency = max(4, n_shafts)
    
    # === GUIDE RAILS ===
    # Each shaft needs up + down rails (or bidirectional)
    rail_length_per_shaft = tower_height * 2  # Up and down
    total_rail_length = rail_length_per_shaft * n_shafts
    total_rail_mass = total_rail_length * rail.mass_per_meter_kg
    
    # === STATIONS ===
    # Ground level + sky lobbies + top
    n_stations = 2 + n_sky_lobbies
    station_cost = 500_000  # Per station (platform, controls, safety)
    
    # === POWER SYSTEM ===
    # Peak power when all capsules climbing
    max_climbing_capsules = total_standard_capsules * 0.5  # Half climbing at once
    peak_power_kW = max_climbing_capsules * capsule.drive_power_kW
    
    # Power distribution (cabling, transformers)
    power_system_cost = peak_power_kW * 500  # $500 per kW capacity
    
    # === CONTROL SYSTEM ===
    control_system_cost = 2_000_000  # Dispatch, monitoring, safety
    
    # === MAINTENANCE FACILITIES ===
    # Workshop at base
    maintenance_cost = 1_500_000
    
    # === COST CALCULATION ===
    costs = {
        'standard_capsules': total_standard_capsules * CAPSULES['standard'].cost_usd,
        'express_capsules': n_express * CAPSULES['express'].cost_usd,
        'service_capsules': n_service * CAPSULES['service'].cost_usd,
        'emergency_capsules': n_emergency * CAPSULES['emergency'].cost_usd,
        'guide_rails': total_rail_length * rail.cost_per_meter_usd,
        'stations': n_stations * station_cost,
        'power_system': power_system_cost,
        'control_system': control_system_cost,
        'maintenance_facility': maintenance_cost,
    }
    
    # Installation and commissioning (30% of hardware)
    hardware_total = sum(costs.values())
    costs['installation'] = hardware_total * 0.3
    
    total_cost = sum(costs.values())
    
    # === PERFORMANCE SPECS ===
    travel_time_express = tower_height / CAPSULES['express'].max_speed_mps / 60  # minutes
    travel_time_local = tower_height / capsule.max_speed_mps / 60
    
    # Mass added to tower
    capsule_mass_all = (total_standard_capsules * CAPSULES['standard'].cabin_mass_kg +
                        n_express * CAPSULES['express'].cabin_mass_kg +
                        n_service * CAPSULES['service'].cabin_mass_kg +
                        n_emergency * CAPSULES['emergency'].cabin_mass_kg)
    
    system_mass_tonnes = (total_rail_mass + capsule_mass_all) / 1000
    
    return {
        'tower_height': tower_height,
        'n_shafts': n_shafts,
        'n_zones': n_zones,
        'rail_type': rail.name,
        
        'capsules': {
            'standard': total_standard_capsules,
            'express': n_express,
            'service': n_service,
            'emergency': n_emergency,
            'total': total_standard_capsules + n_express + n_service + n_emergency,
        },
        
        'infrastructure': {
            'rail_length_m': total_rail_length,
            'rail_mass_tonnes': total_rail_mass / 1000,
            'n_stations': n_stations,
            'n_sky_lobbies': n_sky_lobbies,
        },
        
        'performance': {
            'throughput_per_hour': passenger_throughput_per_hour,
            'travel_time_express_min': travel_time_express,
            'travel_time_local_min': travel_time_local,
            'peak_power_kW': peak_power_kW,
        },
        
        'costs': costs,
        'total_cost': total_cost,
        'system_mass_tonnes': system_mass_tonnes,
    }


def print_elevator_report(system: Dict):
    """Print detailed report of elevator system."""
    
    print(f"\n{'═'*80}")
    print(f"SELF-CLIMBING ROBOT CAPSULE ELEVATOR SYSTEM")
    print(f"{'═'*80}")
    
    print(f"\n  SYSTEM OVERVIEW:")
    print(f"    Tower height: {system['tower_height']:.0f}m")
    print(f"    Elevator shafts: {system['n_shafts']}")
    print(f"    Service zones: {system['n_zones']}")
    print(f"    Sky lobbies: {system['infrastructure']['n_sky_lobbies']}")
    print(f"    Guide rail type: {system['rail_type']}")
    
    print(f"\n  CAPSULE FLEET:")
    caps = system['capsules']
    print(f"    Standard (10 person):  {caps['standard']}")
    print(f"    Express (20 person):   {caps['express']}")
    print(f"    Service/Cargo:         {caps['service']}")
    print(f"    Emergency:             {caps['emergency']}")
    print(f"    TOTAL CAPSULES:        {caps['total']}")
    
    print(f"\n  INFRASTRUCTURE:")
    infra = system['infrastructure']
    print(f"    Total rail length: {infra['rail_length_m']:,.0f}m")
    print(f"    Rail mass: {infra['rail_mass_tonnes']:.0f} tonnes")
    print(f"    Stations: {infra['n_stations']}")
    
    print(f"\n  PERFORMANCE:")
    perf = system['performance']
    print(f"    Passenger throughput: {perf['throughput_per_hour']:,}/hour")
    print(f"    Express travel time: {perf['travel_time_express_min']:.1f} min")
    print(f"    Local travel time: {perf['travel_time_local_min']:.1f} min")
    print(f"    Peak power: {perf['peak_power_kW']:.0f} kW")
    
    print(f"\n  COST BREAKDOWN:")
    costs = system['costs']
    for item, cost in costs.items():
        print(f"    {item.replace('_', ' ').title():.<35} ${cost/1e6:>8.2f} M")
    print(f"    {'─'*47}")
    print(f"    {'TOTAL':.<35} ${system['total_cost']/1e6:>8.2f} M")
    
    print(f"\n  SYSTEM MASS: {system['system_mass_tonnes']:.0f} tonnes (added to tower)")
    
    print(f"\n{'═'*80}")


def get_elevator_cost(tower_height: float, n_shafts: int = 4) -> Dict:
    """
    Quick function to get elevator system cost for integration.
    
    Returns dict with total_cost and system details.
    """
    system = design_elevator_system(tower_height, n_shafts)
    return {
        'total_cost': system['total_cost'],
        'system_mass_tonnes': system['system_mass_tonnes'],
        'n_capsules': system['capsules']['total'],
        'n_shafts': n_shafts,
        'system': system,
    }


if __name__ == "__main__":
    # Design for 1600m tower
    print("="*80)
    print("ELEVATOR SYSTEM ANALYSIS")
    print("="*80)
    
    system_1600 = design_elevator_system(
        tower_height=1600,
        n_shafts=4,
        passenger_throughput_per_hour=2000
    )
    print_elevator_report(system_1600)
    
    # Compare with 3000m tower
    system_3000 = design_elevator_system(
        tower_height=3000,
        n_shafts=6,
        passenger_throughput_per_hour=2000
    )
    print_elevator_report(system_3000)
    
    # Summary comparison
    print(f"\n{'═'*80}")
    print("ELEVATOR SYSTEM COMPARISON")
    print(f"{'═'*80}\n")
    
    print(f"{'Height':>8} {'Shafts':>8} {'Capsules':>10} {'Rail(t)':>10} {'Cost':>14}")
    print("-" * 55)
    
    for sys in [system_1600, system_3000]:
        print(f"{sys['tower_height']:>8.0f}m {sys['n_shafts']:>8d} "
              f"{sys['capsules']['total']:>10d} "
              f"{sys['infrastructure']['rail_mass_tonnes']:>10.0f} "
              f"${sys['total_cost']/1e6:>12.1f}M")
