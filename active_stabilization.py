"""
Active Stabilization System for Egg Tower.

Uses guy wires as actuators to counteract wind-induced oscillations.
The system adjusts cable tensions in real-time to dampen sway.

Components:
1. Sensors: Accelerometers, GPS, wind sensors at multiple heights
2. Actuators: Hydraulic or electric tensioners on guy cables
3. Control System: Real-time feedback control to dampen oscillations

Design Philosophy:
- Guy wires already provide static wind resistance
- Adding active tensioning allows dynamic response to gusts
- Opposing cables can push/pull to counteract lateral motion
- System can also dampen vortex-induced vibrations
"""
import numpy as np
from dataclasses import dataclass
from typing import List, Dict


@dataclass
class Actuator:
    """Hydraulic or electric cable tensioner."""
    name: str
    max_force_kN: float          # Maximum tension adjustment
    stroke_m: float              # Maximum cable length change
    response_time_ms: float      # Time to full stroke
    power_kW: float              # Peak power consumption
    cost_usd: float              # Unit cost


@dataclass 
class Sensor:
    """Motion/wind sensor for feedback control."""
    name: str
    measurement: str             # What it measures
    accuracy: str                # Measurement accuracy
    sample_rate_hz: float        # Sampling frequency
    cost_usd: float              # Unit cost


# Component specifications
ACTUATORS = {
    'hydraulic_large': Actuator(
        name="Hydraulic Tensioner (Large)",
        max_force_kN=5000,       # 5 MN capacity
        stroke_m=0.5,            # 50cm stroke
        response_time_ms=200,    # 0.2 second response
        power_kW=150,            # 150 kW peak
        cost_usd=250_000         # $250K each
    ),
    'hydraulic_medium': Actuator(
        name="Hydraulic Tensioner (Medium)", 
        max_force_kN=2000,       # 2 MN capacity
        stroke_m=0.3,            # 30cm stroke
        response_time_ms=150,    # 0.15 second response
        power_kW=75,             # 75 kW peak
        cost_usd=120_000         # $120K each
    ),
    'electric_servo': Actuator(
        name="Electric Servo Tensioner",
        max_force_kN=500,        # 500 kN capacity
        stroke_m=0.2,            # 20cm stroke
        response_time_ms=50,     # 50ms response (faster)
        power_kW=50,             # 50 kW peak
        cost_usd=80_000          # $80K each
    ),
}

SENSORS = {
    'accelerometer': Sensor(
        name="Triaxial Accelerometer",
        measurement="Acceleration (x,y,z)",
        accuracy="±0.001g",
        sample_rate_hz=1000,
        cost_usd=5_000
    ),
    'gps_rtk': Sensor(
        name="RTK GPS Receiver",
        measurement="Position (cm accuracy)",
        accuracy="±2cm horizontal",
        sample_rate_hz=20,
        cost_usd=15_000
    ),
    'inclinometer': Sensor(
        name="MEMS Inclinometer",
        measurement="Tilt angle",
        accuracy="±0.001°",
        sample_rate_hz=100,
        cost_usd=3_000
    ),
    'wind_sensor': Sensor(
        name="Ultrasonic Anemometer",
        measurement="Wind speed/direction",
        accuracy="±0.1 m/s",
        sample_rate_hz=10,
        cost_usd=8_000
    ),
    'strain_gauge': Sensor(
        name="Cable Strain Gauge",
        measurement="Cable tension",
        accuracy="±0.1%",
        sample_rate_hz=500,
        cost_usd=2_000
    ),
}


def design_active_stabilization(n_eggs: int, n_guy_levels: int, 
                                 n_cables_per_level: int = 9,
                                 tower_height: float = 1600) -> Dict:
    """
    Design active stabilization system for the tower.
    
    Parameters:
        n_eggs: Number of eggs in tower
        n_guy_levels: Number of guy wire levels
        n_cables_per_level: Cables per guy level (default 9)
        tower_height: Total tower height in meters
    
    Returns:
        Dictionary with system design and cost breakdown
    """
    
    # === ACTUATORS ===
    # Each guy cable gets a tensioner at the ground anchor
    # Use medium hydraulic for main actuation
    n_actuators = n_guy_levels * n_cables_per_level
    actuator_type = ACTUATORS['hydraulic_medium']
    
    # === SENSORS ===
    # Accelerometers: 2 per egg (orthogonal horizontal axes)
    n_accelerometers = n_eggs * 2
    
    # GPS: 1 per guy level + top
    n_gps = n_guy_levels + 1
    
    # Inclinometers: 1 per egg
    n_inclinometers = n_eggs
    
    # Wind sensors: Every 200m of height
    n_wind_sensors = max(4, int(tower_height / 200))
    
    # Strain gauges: 1 per cable (at actuator)
    n_strain_gauges = n_actuators
    
    # === CONTROL SYSTEM ===
    # Redundant industrial controllers
    n_controllers = 3  # Triple redundancy
    controller_cost = 500_000  # Per controller with software
    
    # Communication network (fiber optic)
    fiber_per_meter = 50  # $/m for armored fiber bundle
    fiber_length = tower_height * 4  # Up/down + redundancy
    fiber_cost = fiber_per_meter * fiber_length
    
    # Power distribution
    power_per_actuator = actuator_type.power_kW
    total_power_kW = n_actuators * power_per_actuator * 0.3  # 30% duty cycle
    
    # UPS/backup power (4 hours at 50% capacity)
    backup_kwh = total_power_kW * 0.5 * 4
    ups_cost_per_kwh = 500  # Battery cost
    ups_cost = backup_kwh * ups_cost_per_kwh
    
    # === COST CALCULATION ===
    costs = {
        'actuators': n_actuators * actuator_type.cost_usd,
        'accelerometers': n_accelerometers * SENSORS['accelerometer'].cost_usd,
        'gps': n_gps * SENSORS['gps_rtk'].cost_usd,
        'inclinometers': n_inclinometers * SENSORS['inclinometer'].cost_usd,
        'wind_sensors': n_wind_sensors * SENSORS['wind_sensor'].cost_usd,
        'strain_gauges': n_strain_gauges * SENSORS['strain_gauge'].cost_usd,
        'controllers': n_controllers * controller_cost,
        'fiber_network': fiber_cost,
        'power_backup': ups_cost,
    }
    
    # Installation and commissioning (50% of hardware)
    hardware_total = sum(costs.values())
    costs['installation'] = hardware_total * 0.5
    
    # Software development and testing
    costs['software'] = 2_000_000  # Custom control algorithms
    
    # Ongoing maintenance (annual, first year included)
    costs['maintenance_year1'] = hardware_total * 0.05
    
    total_cost = sum(costs.values())
    
    # Component counts
    components = {
        'actuators': {
            'count': n_actuators,
            'type': actuator_type.name,
            'unit_cost': actuator_type.cost_usd,
        },
        'accelerometers': {'count': n_accelerometers},
        'gps_receivers': {'count': n_gps},
        'inclinometers': {'count': n_inclinometers},
        'wind_sensors': {'count': n_wind_sensors},
        'strain_gauges': {'count': n_strain_gauges},
        'controllers': {'count': n_controllers},
    }
    
    # Performance specs
    performance = {
        'max_correction_force_MN': n_cables_per_level * actuator_type.max_force_kN / 1000,
        'response_time_ms': actuator_type.response_time_ms,
        'sensor_sample_rate_hz': SENSORS['accelerometer'].sample_rate_hz,
        'total_power_kW': total_power_kW,
        'backup_hours': 4,
    }
    
    return {
        'components': components,
        'costs': costs,
        'total_cost': total_cost,
        'performance': performance,
        'n_guy_levels': n_guy_levels,
        'n_cables_per_level': n_cables_per_level,
    }


def print_stabilization_report(system: Dict, tower_height: float = 1600):
    """Print detailed report of active stabilization system."""
    
    print(f"\n{'═'*80}")
    print(f"ACTIVE STABILIZATION SYSTEM - {tower_height:.0f}m TOWER")
    print(f"{'═'*80}")
    
    print(f"\n  SYSTEM OVERVIEW:")
    print(f"    Guy wire levels: {system['n_guy_levels']}")
    print(f"    Cables per level: {system['n_cables_per_level']}")
    print(f"    Total actuated cables: {system['n_guy_levels'] * system['n_cables_per_level']}")
    
    print(f"\n  ACTUATORS:")
    act = system['components']['actuators']
    print(f"    Type: {act['type']}")
    print(f"    Quantity: {act['count']}")
    print(f"    Max correction force: {system['performance']['max_correction_force_MN']:.1f} MN per level")
    print(f"    Response time: {system['performance']['response_time_ms']:.0f} ms")
    
    print(f"\n  SENSORS:")
    for sensor_type in ['accelerometers', 'gps_receivers', 'inclinometers', 
                         'wind_sensors', 'strain_gauges']:
        if sensor_type in system['components']:
            count = system['components'][sensor_type]['count']
            print(f"    {sensor_type.replace('_', ' ').title()}: {count}")
    
    print(f"\n  CONTROL SYSTEM:")
    print(f"    Redundant controllers: {system['components']['controllers']['count']}")
    print(f"    Sensor sample rate: {system['performance']['sensor_sample_rate_hz']} Hz")
    print(f"    Peak power: {system['performance']['total_power_kW']:.0f} kW")
    print(f"    Backup power: {system['performance']['backup_hours']} hours")
    
    print(f"\n  COST BREAKDOWN:")
    costs = system['costs']
    for item, cost in costs.items():
        print(f"    {item.replace('_', ' ').title():.<30} ${cost/1e6:>8.2f} M")
    print(f"    {'─'*42}")
    print(f"    {'TOTAL':.<30} ${system['total_cost']/1e6:>8.2f} M")
    
    print(f"\n{'═'*80}")


def add_stabilization_to_tower_cost(tower_cost_M: float, n_eggs: int, 
                                     n_guy_levels: int, tower_height: float) -> Dict:
    """
    Add active stabilization cost to tower construction cost.
    
    Returns updated cost breakdown.
    """
    stab = design_active_stabilization(n_eggs, n_guy_levels, 
                                        tower_height=tower_height)
    stab_cost_M = stab['total_cost'] / 1e6
    
    total_with_stab = tower_cost_M + stab_cost_M
    stab_percent = (stab_cost_M / total_with_stab) * 100
    
    return {
        'tower_cost_M': tower_cost_M,
        'stabilization_cost_M': stab_cost_M,
        'total_cost_M': total_with_stab,
        'stabilization_percent': stab_percent,
        'stabilization_system': stab,
    }


if __name__ == "__main__":
    # Example: 1600m tower with 9 guy levels
    print("="*80)
    print("ACTIVE STABILIZATION ANALYSIS")
    print("="*80)
    
    # Design for 1600m tower
    system_1600 = design_active_stabilization(
        n_eggs=28, 
        n_guy_levels=9, 
        n_cables_per_level=9,
        tower_height=1600
    )
    print_stabilization_report(system_1600, tower_height=1600)
    
    # Design for 3000m tower (for comparison)
    system_3000 = design_active_stabilization(
        n_eggs=18,
        n_guy_levels=9,
        n_cables_per_level=9, 
        tower_height=3000
    )
    print_stabilization_report(system_3000, tower_height=3000)
    
    # Add to total tower costs
    print(f"\n{'═'*80}")
    print("TOTAL PROJECT COSTS WITH ACTIVE STABILIZATION")
    print(f"{'═'*80}\n")
    
    # From cost_analysis.py results
    towers = [
        {'height': 1600, 'n_eggs': 28, 'n_guy': 9, 'structure_cost_M': 391.7},
        {'height': 3000, 'n_eggs': 18, 'n_guy': 9, 'structure_cost_M': 13351.6},
    ]
    
    print(f"{'Height':>8} {'Structure':>14} {'Stabilization':>14} {'TOTAL':>14} {'Stab %':>10}")
    print("-" * 65)
    
    for t in towers:
        result = add_stabilization_to_tower_cost(
            t['structure_cost_M'], 
            t['n_eggs'],
            t['n_guy'],
            t['height']
        )
        print(f"{t['height']:>8.0f}m ${result['tower_cost_M']:>12.1f}M "
              f"${result['stabilization_cost_M']:>12.1f}M "
              f"${result['total_cost_M']:>12.1f}M "
              f"{result['stabilization_percent']:>9.1f}%")
    
    print(f"\n{'═'*80}")
