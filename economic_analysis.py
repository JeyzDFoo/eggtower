"""
Economic Analysis for Egg Tower Project.

Revenue projections and financial viability assessment for:
1. Observation/Tourism
2. Telecommunications
3. Hospitality (restaurants, events, hotels)
4. Real estate appreciation effects

Based on comparable structures:
- Burj Khalifa (828m): ~2M visitors/year, $40-150 tickets
- Tokyo Skytree (634m): ~6M visitors/year (first year), $15-30 tickets
- CN Tower (553m): ~1.5M visitors/year, $30-50 tickets
- Shanghai Tower (632m): ~3M visitors/year, $25-45 tickets
"""
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class RevenueScenario:
    """A revenue projection scenario."""
    name: str
    annual_revenue: float  # USD
    description: str


# ============================================================================
# TOURISM & OBSERVATION
# ============================================================================

def calculate_tourism_revenue(
    tower_height: float,
    ticket_price: float = 100,
    visitors_per_year: float = 2_500_000,
    observation_levels: int = 3,
    premium_deck_multiplier: float = 2.0,
    premium_deck_fraction: float = 0.15
) -> Dict:
    """
    Estimate annual tourism/observation revenue.
    
    Parameters:
        tower_height: Tower height in meters
        ticket_price: Base ticket price (USD)
        visitors_per_year: Annual visitor count
        observation_levels: Number of observation decks
        premium_deck_multiplier: Price multiplier for highest deck
        premium_deck_fraction: Fraction of visitors paying premium
    
    Returns:
        Revenue breakdown dictionary
    """
    # Height premium: taller = more valuable experience
    # Burj Khalifa at 828m charges ~$40-150
    # At 1600m (2x height), can command premium pricing
    height_factor = (tower_height / 800) ** 0.5  # Diminishing returns
    
    # Adjusted ticket price based on height
    effective_base_price = ticket_price * height_factor
    premium_price = effective_base_price * premium_deck_multiplier
    
    # Revenue calculation
    standard_visitors = visitors_per_year * (1 - premium_deck_fraction)
    premium_visitors = visitors_per_year * premium_deck_fraction
    
    standard_revenue = standard_visitors * effective_base_price
    premium_revenue = premium_visitors * premium_price
    total_revenue = standard_revenue + premium_revenue
    
    # Gift shop and merchandise (typically 15-25% of ticket revenue)
    merchandise_rate = 0.20
    merchandise_revenue = total_revenue * merchandise_rate
    
    # Photography services (green screen, professional photos)
    photo_revenue_per_visitor = 5  # Average across all visitors
    photo_revenue = visitors_per_year * photo_revenue_per_visitor
    
    return {
        'visitors_per_year': visitors_per_year,
        'effective_base_price': effective_base_price,
        'premium_price': premium_price,
        'height_factor': height_factor,
        'ticket_revenue': total_revenue,
        'merchandise_revenue': merchandise_revenue,
        'photo_revenue': photo_revenue,
        'total_revenue': total_revenue + merchandise_revenue + photo_revenue,
    }


# ============================================================================
# TELECOMMUNICATIONS
# ============================================================================

def calculate_telecom_revenue(
    tower_height: float,
    n_mobile_carriers: int = 5,
    n_broadcasters: int = 8,
    include_emergency: bool = True,
    include_research: bool = True
) -> Dict:
    """
    Estimate annual telecommunications lease revenue.
    
    At 1600m, line-of-sight coverage extends ~140km radius.
    This is extremely valuable real estate for:
    - Mobile carriers (5G/6G base stations)
    - TV/Radio broadcasters
    - Emergency services
    - Weather/atmospheric research
    
    Parameters:
        tower_height: Tower height in meters
        n_mobile_carriers: Number of mobile network tenants
        n_broadcasters: Number of TV/Radio tenants
        include_emergency: Include emergency services lease
        include_research: Include research/weather stations
    """
    # Height premium for coverage area
    # Coverage radius ~ sqrt(2 * R_earth * h) where h is height
    # At 1600m vs 600m, coverage area is ~2.7x larger
    coverage_factor = (tower_height / 600) ** 0.5
    
    # Base lease rates (annual, USD)
    # Based on comparable tower lease rates, scaled for height premium
    mobile_carrier_base = 8_000_000  # $8M/year base for major carrier
    broadcaster_base = 3_000_000     # $3M/year base for TV/radio
    emergency_base = 2_000_000       # $2M/year for emergency services
    research_base = 1_500_000        # $1.5M/year for research stations
    
    # Apply coverage premium
    mobile_lease = mobile_carrier_base * coverage_factor
    broadcaster_lease = broadcaster_base * coverage_factor
    
    # Calculate revenues
    mobile_revenue = n_mobile_carriers * mobile_lease
    broadcast_revenue = n_broadcasters * broadcaster_lease
    emergency_revenue = emergency_base if include_emergency else 0
    research_revenue = research_base * 2 if include_research else 0  # Multiple stations
    
    # Additional: small cell/IoT infrastructure
    iot_revenue = 5_000_000 * coverage_factor
    
    return {
        'coverage_factor': coverage_factor,
        'coverage_radius_km': np.sqrt(2 * 6371 * tower_height / 1000),  # Approximate
        'n_mobile_carriers': n_mobile_carriers,
        'n_broadcasters': n_broadcasters,
        'mobile_lease_each': mobile_lease,
        'broadcaster_lease_each': broadcaster_lease,
        'mobile_revenue': mobile_revenue,
        'broadcast_revenue': broadcast_revenue,
        'emergency_revenue': emergency_revenue,
        'research_revenue': research_revenue,
        'iot_revenue': iot_revenue,
        'total_revenue': mobile_revenue + broadcast_revenue + emergency_revenue + research_revenue + iot_revenue,
    }


# ============================================================================
# HOSPITALITY
# ============================================================================

def calculate_hospitality_revenue(
    tower_height: float,
    n_restaurants: int = 3,
    n_event_spaces: int = 2,
    hotel_rooms: int = 0,  # Optional hotel component
    restaurant_covers_per_day: int = 400,
    avg_spend_per_cover: float = 150,
    events_per_year: int = 200,
    avg_event_revenue: float = 50_000,
    hotel_rate_per_night: float = 500,
    hotel_occupancy: float = 0.75
) -> Dict:
    """
    Estimate annual hospitality revenue.
    
    Includes:
    - Fine dining restaurants at various heights
    - Event/function spaces
    - Optional hotel rooms
    """
    # Height premium for dining experience
    height_premium = (tower_height / 500) ** 0.3
    
    # Restaurant revenue
    operating_days = 365
    restaurant_revenue = (n_restaurants * restaurant_covers_per_day * 
                         operating_days * avg_spend_per_cover * height_premium)
    
    # Event space revenue
    event_revenue = events_per_year * avg_event_revenue * height_premium
    
    # Hotel revenue (if applicable)
    hotel_revenue = 0
    if hotel_rooms > 0:
        room_nights = hotel_rooms * 365 * hotel_occupancy
        hotel_revenue = room_nights * hotel_rate_per_night * height_premium
    
    return {
        'n_restaurants': n_restaurants,
        'n_event_spaces': n_event_spaces,
        'hotel_rooms': hotel_rooms,
        'height_premium': height_premium,
        'restaurant_revenue': restaurant_revenue,
        'event_revenue': event_revenue,
        'hotel_revenue': hotel_revenue,
        'total_revenue': restaurant_revenue + event_revenue + hotel_revenue,
    }


# ============================================================================
# OPERATING COSTS
# ============================================================================

def calculate_operating_costs(
    tower_height: float,
    total_revenue: float,
    n_eggs: int,
    elevator_system_cost: float,
    stabilization_system_cost: float
) -> Dict:
    """
    Estimate annual operating costs.
    
    Includes:
    - Staff and labor
    - Utilities (power, water, HVAC)
    - Maintenance (structure, elevators, systems)
    - Insurance
    - Marketing
    - Property taxes/fees
    """
    # Staffing (scales with height and revenue)
    # Security, operations, hospitality, maintenance, management
    staff_count = int(50 + tower_height / 20)  # ~130 staff for 1600m
    avg_salary = 60_000
    benefits_multiplier = 1.35
    labor_cost = staff_count * avg_salary * benefits_multiplier
    
    # Utilities
    # Power for elevators, lighting, HVAC, systems
    # Estimate based on floor area and height
    power_per_meter = 500  # kWh/year per meter of height
    power_cost_kwh = 0.12
    power_usage = tower_height * power_per_meter * n_eggs / 10
    utility_cost = power_usage * power_cost_kwh
    utility_cost += 2_000_000  # Base utilities (water, waste, etc.)
    
    # Maintenance
    # Structure: ~0.5% of construction cost annually
    # Elevators: ~5% of elevator system cost
    # Active systems: ~3% of stabilization cost
    structure_maintenance = 5_000_000  # Base for unique structure
    elevator_maintenance = elevator_system_cost * 0.05
    systems_maintenance = stabilization_system_cost * 0.03
    maintenance_cost = structure_maintenance + elevator_maintenance + systems_maintenance
    
    # Insurance (property, liability, business interruption)
    # Higher for unique structures
    insurance_cost = total_revenue * 0.03 + 5_000_000  # 3% of revenue + base
    
    # Marketing and sales
    marketing_cost = total_revenue * 0.05  # 5% of revenue
    
    # Property taxes and fees (varies by jurisdiction)
    # Estimate 1% of imputed property value
    property_tax = total_revenue * 2 * 0.01  # ~2x revenue as property value proxy
    
    # Management and overhead
    overhead_cost = total_revenue * 0.03  # 3% of revenue
    
    total_opex = (labor_cost + utility_cost + maintenance_cost + 
                  insurance_cost + marketing_cost + property_tax + overhead_cost)
    
    return {
        'staff_count': staff_count,
        'labor_cost': labor_cost,
        'utility_cost': utility_cost,
        'maintenance_cost': maintenance_cost,
        'insurance_cost': insurance_cost,
        'marketing_cost': marketing_cost,
        'property_tax': property_tax,
        'overhead_cost': overhead_cost,
        'total_opex': total_opex,
        'opex_as_pct_revenue': total_opex / total_revenue * 100 if total_revenue > 0 else 0,
    }


# ============================================================================
# FINANCIAL ANALYSIS
# ============================================================================

def calculate_financial_metrics(
    construction_cost: float,
    annual_revenue: float,
    annual_opex: float,
    discount_rate: float = 0.08,
    analysis_years: int = 30,
    revenue_growth_rate: float = 0.02,
    opex_growth_rate: float = 0.025,
    tax_rate: float = 0.25,
    depreciation_years: int = 40
) -> Dict:
    """
    Calculate key financial metrics for the project.
    
    Parameters:
        construction_cost: Total construction cost (USD)
        annual_revenue: First year annual revenue (USD)
        annual_opex: First year operating costs (USD)
        discount_rate: Discount rate for NPV (default 8%)
        analysis_years: Years of analysis (default 30)
        revenue_growth_rate: Annual revenue growth (default 2%)
        opex_growth_rate: Annual opex growth (default 2.5%)
        tax_rate: Corporate tax rate (default 25%)
        depreciation_years: Straight-line depreciation period
    
    Returns:
        Financial metrics dictionary
    """
    # Annual depreciation
    annual_depreciation = construction_cost / depreciation_years
    
    # Year-by-year cash flows
    cash_flows = []
    cumulative_cash = -construction_cost
    payback_year = None
    
    for year in range(1, analysis_years + 1):
        # Revenue and costs with growth
        year_revenue = annual_revenue * (1 + revenue_growth_rate) ** (year - 1)
        year_opex = annual_opex * (1 + opex_growth_rate) ** (year - 1)
        
        # EBITDA
        ebitda = year_revenue - year_opex
        
        # Taxable income (EBITDA - depreciation)
        taxable_income = ebitda - annual_depreciation
        taxes = max(0, taxable_income * tax_rate)
        
        # Net income
        net_income = ebitda - taxes
        
        # Free cash flow (add back depreciation as non-cash)
        # Simplified: assume minimal CapEx after construction
        fcf = net_income
        
        # Cumulative for payback
        cumulative_cash += fcf
        if payback_year is None and cumulative_cash >= 0:
            # Linear interpolation for fractional year
            prev_cumulative = cumulative_cash - fcf
            fraction = -prev_cumulative / fcf
            payback_year = year - 1 + fraction
        
        cash_flows.append({
            'year': year,
            'revenue': year_revenue,
            'opex': year_opex,
            'ebitda': ebitda,
            'taxes': taxes,
            'net_income': net_income,
            'fcf': fcf,
            'cumulative': cumulative_cash,
        })
    
    # NPV calculation
    npv = -construction_cost
    for i, cf in enumerate(cash_flows):
        npv += cf['fcf'] / (1 + discount_rate) ** (i + 1)
    
    # IRR calculation (iterative)
    def npv_at_rate(rate):
        result = -construction_cost
        for i, cf in enumerate(cash_flows):
            result += cf['fcf'] / (1 + rate) ** (i + 1)
        return result
    
    # Binary search for IRR
    irr_low, irr_high = -0.5, 1.0
    for _ in range(100):
        irr_mid = (irr_low + irr_high) / 2
        if npv_at_rate(irr_mid) > 0:
            irr_low = irr_mid
        else:
            irr_high = irr_mid
    irr = (irr_low + irr_high) / 2
    
    # First year metrics
    first_year = cash_flows[0]
    ebitda_margin = first_year['ebitda'] / first_year['revenue'] * 100
    
    # Year 10 metrics
    year_10 = cash_flows[9] if len(cash_flows) >= 10 else cash_flows[-1]
    
    return {
        'construction_cost': construction_cost,
        'first_year_revenue': annual_revenue,
        'first_year_opex': annual_opex,
        'first_year_ebitda': first_year['ebitda'],
        'first_year_net_income': first_year['net_income'],
        'ebitda_margin_pct': ebitda_margin,
        'payback_years': payback_year,
        'npv': npv,
        'irr': irr,
        'irr_pct': irr * 100,
        'year_10_revenue': year_10['revenue'],
        'year_10_ebitda': year_10['ebitda'],
        'year_10_cumulative': year_10['cumulative'],
        'cash_flows': cash_flows,
    }


# ============================================================================
# SCENARIO ANALYSIS
# ============================================================================

def run_scenario_analysis(
    tower_height: float,
    n_eggs: int,
    construction_cost: float,
    elevator_system_cost: float,
    stabilization_system_cost: float
) -> Dict:
    """
    Run conservative, moderate, and optimistic scenarios.
    """
    scenarios = {}
    
    # Conservative scenario
    tourism_cons = calculate_tourism_revenue(
        tower_height, ticket_price=75, visitors_per_year=1_500_000)
    telecom_cons = calculate_telecom_revenue(
        tower_height, n_mobile_carriers=3, n_broadcasters=5)
    hospitality_cons = calculate_hospitality_revenue(
        tower_height, n_restaurants=2, events_per_year=100, hotel_rooms=0)
    
    revenue_cons = (tourism_cons['total_revenue'] + telecom_cons['total_revenue'] + 
                    hospitality_cons['total_revenue'])
    opex_cons = calculate_operating_costs(
        tower_height, revenue_cons, n_eggs, elevator_system_cost, stabilization_system_cost)
    financial_cons = calculate_financial_metrics(construction_cost, revenue_cons, opex_cons['total_opex'])
    
    scenarios['conservative'] = {
        'tourism': tourism_cons,
        'telecom': telecom_cons,
        'hospitality': hospitality_cons,
        'total_revenue': revenue_cons,
        'opex': opex_cons,
        'financial': financial_cons,
    }
    
    # Moderate scenario
    tourism_mod = calculate_tourism_revenue(
        tower_height, ticket_price=100, visitors_per_year=2_500_000)
    telecom_mod = calculate_telecom_revenue(
        tower_height, n_mobile_carriers=5, n_broadcasters=8)
    hospitality_mod = calculate_hospitality_revenue(
        tower_height, n_restaurants=3, events_per_year=200, hotel_rooms=50)
    
    revenue_mod = (tourism_mod['total_revenue'] + telecom_mod['total_revenue'] + 
                   hospitality_mod['total_revenue'])
    opex_mod = calculate_operating_costs(
        tower_height, revenue_mod, n_eggs, elevator_system_cost, stabilization_system_cost)
    financial_mod = calculate_financial_metrics(construction_cost, revenue_mod, opex_mod['total_opex'])
    
    scenarios['moderate'] = {
        'tourism': tourism_mod,
        'telecom': telecom_mod,
        'hospitality': hospitality_mod,
        'total_revenue': revenue_mod,
        'opex': opex_mod,
        'financial': financial_mod,
    }
    
    # Optimistic scenario
    tourism_opt = calculate_tourism_revenue(
        tower_height, ticket_price=125, visitors_per_year=4_000_000)
    telecom_opt = calculate_telecom_revenue(
        tower_height, n_mobile_carriers=6, n_broadcasters=10)
    hospitality_opt = calculate_hospitality_revenue(
        tower_height, n_restaurants=4, events_per_year=300, hotel_rooms=100)
    
    revenue_opt = (tourism_opt['total_revenue'] + telecom_opt['total_revenue'] + 
                   hospitality_opt['total_revenue'])
    opex_opt = calculate_operating_costs(
        tower_height, revenue_opt, n_eggs, elevator_system_cost, stabilization_system_cost)
    financial_opt = calculate_financial_metrics(construction_cost, revenue_opt, opex_opt['total_opex'])
    
    scenarios['optimistic'] = {
        'tourism': tourism_opt,
        'telecom': telecom_opt,
        'hospitality': hospitality_opt,
        'total_revenue': revenue_opt,
        'opex': opex_opt,
        'financial': financial_opt,
    }
    
    return scenarios


# ============================================================================
# PRINTING FUNCTIONS
# ============================================================================

def print_economic_analysis(scenarios: Dict, construction_cost: float):
    """Print formatted economic analysis results."""
    
    print(f"\n{'═'*100}")
    print(f"ECONOMIC ANALYSIS - EGG TOWER PROJECT")
    print(f"{'═'*100}")
    print(f"Construction Cost: ${construction_cost/1e9:.2f}B")
    print(f"{'═'*100}")
    
    for scenario_name in ['conservative', 'moderate', 'optimistic']:
        s = scenarios[scenario_name]
        fin = s['financial']
        
        print(f"\n{'─'*100}")
        print(f"  {scenario_name.upper()} SCENARIO")
        print(f"{'─'*100}")
        
        # Revenue breakdown
        print(f"\n  ANNUAL REVENUE (Year 1)")
        print(f"    Tourism/Observation:    ${s['tourism']['total_revenue']/1e6:>10,.1f}M  "
              f"({s['tourism']['visitors_per_year']/1e6:.1f}M visitors @ ${s['tourism']['effective_base_price']:.0f})")
        print(f"    Telecommunications:     ${s['telecom']['total_revenue']/1e6:>10,.1f}M  "
              f"({s['telecom']['n_mobile_carriers']} carriers, {s['telecom']['n_broadcasters']} broadcasters)")
        print(f"    Hospitality:            ${s['hospitality']['total_revenue']/1e6:>10,.1f}M  "
              f"({s['hospitality']['n_restaurants']} restaurants, {s['hospitality']['hotel_rooms']} hotel rooms)")
        print(f"    {'─'*50}")
        print(f"    TOTAL REVENUE:          ${s['total_revenue']/1e6:>10,.1f}M")
        
        # Operating costs
        print(f"\n  OPERATING COSTS (Year 1)")
        print(f"    Labor ({s['opex']['staff_count']} staff):       ${s['opex']['labor_cost']/1e6:>10,.1f}M")
        print(f"    Utilities:              ${s['opex']['utility_cost']/1e6:>10,.1f}M")
        print(f"    Maintenance:            ${s['opex']['maintenance_cost']/1e6:>10,.1f}M")
        print(f"    Insurance:              ${s['opex']['insurance_cost']/1e6:>10,.1f}M")
        print(f"    Marketing:              ${s['opex']['marketing_cost']/1e6:>10,.1f}M")
        print(f"    Other:                  ${(s['opex']['property_tax']+s['opex']['overhead_cost'])/1e6:>10,.1f}M")
        print(f"    {'─'*50}")
        print(f"    TOTAL OPEX:             ${s['opex']['total_opex']/1e6:>10,.1f}M  "
              f"({s['opex']['opex_as_pct_revenue']:.1f}% of revenue)")
        
        # Financial metrics
        print(f"\n  FINANCIAL METRICS")
        print(f"    EBITDA (Year 1):        ${fin['first_year_ebitda']/1e6:>10,.1f}M  "
              f"({fin['ebitda_margin_pct']:.1f}% margin)")
        print(f"    Net Income (Year 1):    ${fin['first_year_net_income']/1e6:>10,.1f}M")
        print(f"    Payback Period:         {fin['payback_years']:>10.1f} years")
        print(f"    NPV (8% discount):      ${fin['npv']/1e6:>10,.1f}M")
        print(f"    IRR:                    {fin['irr_pct']:>10.1f}%")
        print(f"    Year 10 Revenue:        ${fin['year_10_revenue']/1e6:>10,.1f}M")
        print(f"    Year 10 Cumulative:     ${fin['year_10_cumulative']/1e6:>10,.1f}M")
    
    # Summary comparison
    print(f"\n{'═'*100}")
    print(f"  SCENARIO COMPARISON SUMMARY")
    print(f"{'═'*100}")
    print(f"\n  {'Metric':<25} {'Conservative':>18} {'Moderate':>18} {'Optimistic':>18}")
    print(f"  {'─'*79}")
    
    cons, mod, opt = scenarios['conservative'], scenarios['moderate'], scenarios['optimistic']
    
    print(f"  {'Annual Revenue':<25} ${cons['total_revenue']/1e6:>14,.0f}M ${mod['total_revenue']/1e6:>14,.0f}M ${opt['total_revenue']/1e6:>14,.0f}M")
    print(f"  {'EBITDA':<25} ${cons['financial']['first_year_ebitda']/1e6:>14,.0f}M ${mod['financial']['first_year_ebitda']/1e6:>14,.0f}M ${opt['financial']['first_year_ebitda']/1e6:>14,.0f}M")
    print(f"  {'EBITDA Margin':<25} {cons['financial']['ebitda_margin_pct']:>14.1f}% {mod['financial']['ebitda_margin_pct']:>14.1f}% {opt['financial']['ebitda_margin_pct']:>14.1f}%")
    print(f"  {'Payback (years)':<25} {cons['financial']['payback_years']:>15.1f} {mod['financial']['payback_years']:>15.1f} {opt['financial']['payback_years']:>15.1f}")
    print(f"  {'NPV':<25} ${cons['financial']['npv']/1e6:>14,.0f}M ${mod['financial']['npv']/1e6:>14,.0f}M ${opt['financial']['npv']/1e6:>14,.0f}M")
    print(f"  {'IRR':<25} {cons['financial']['irr_pct']:>14.1f}% {mod['financial']['irr_pct']:>14.1f}% {opt['financial']['irr_pct']:>14.1f}%")
    
    print(f"\n{'═'*100}")


def save_economic_analysis_png(scenarios: Dict, construction_cost: float, 
                                tower_height: float, save_path: str = 'economic_analysis.png'):
    """Generate and save economic analysis as a PNG chart."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    fig = plt.figure(figsize=(16, 20))
    
    # Title
    fig.suptitle('EGG TOWER ECONOMIC ANALYSIS', fontsize=24, fontweight='bold', y=0.98)
    fig.text(0.5, 0.955, f'{tower_height:.0f}m Tower | Construction Cost: ${construction_cost/1e9:.2f}B',
             ha='center', fontsize=14, color='#7F8C8D')
    
    # Create subplots
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.25, top=0.92, bottom=0.05, left=0.08, right=0.92)
    
    # Colors
    colors = {'conservative': '#3498DB', 'moderate': '#27AE60', 'optimistic': '#E74C3C'}
    
    # 1. Revenue Breakdown (stacked bar)
    ax1 = fig.add_subplot(gs[0, 0])
    scenarios_list = ['conservative', 'moderate', 'optimistic']
    x = np.arange(len(scenarios_list))
    width = 0.6
    
    tourism = [scenarios[s]['tourism']['total_revenue']/1e6 for s in scenarios_list]
    telecom = [scenarios[s]['telecom']['total_revenue']/1e6 for s in scenarios_list]
    hospitality = [scenarios[s]['hospitality']['total_revenue']/1e6 for s in scenarios_list]
    
    ax1.bar(x, tourism, width, label='Tourism', color='#3498DB')
    ax1.bar(x, telecom, width, bottom=tourism, label='Telecom', color='#E74C3C')
    ax1.bar(x, hospitality, width, bottom=np.array(tourism)+np.array(telecom), label='Hospitality', color='#F39C12')
    
    ax1.set_ylabel('Annual Revenue ($M)', fontsize=12)
    ax1.set_title('Revenue by Source', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(['Conservative', 'Moderate', 'Optimistic'])
    ax1.legend(loc='upper left')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add total labels
    totals = [scenarios[s]['total_revenue']/1e6 for s in scenarios_list]
    for i, t in enumerate(totals):
        ax1.text(i, t + 10, f'${t:.0f}M', ha='center', fontsize=11, fontweight='bold')
    
    # 2. Financial Metrics Comparison
    ax2 = fig.add_subplot(gs[0, 1])
    metrics = ['Payback\n(years)', 'IRR\n(%)', 'EBITDA\nMargin (%)']
    x = np.arange(len(metrics))
    width = 0.25
    
    for i, scenario in enumerate(scenarios_list):
        fin = scenarios[scenario]['financial']
        values = [fin['payback_years'], fin['irr_pct'], fin['ebitda_margin_pct']]
        ax2.bar(x + i*width, values, width, label=scenario.capitalize(), color=colors[scenario])
    
    ax2.set_ylabel('Value', fontsize=12)
    ax2.set_title('Key Financial Metrics', fontsize=14, fontweight='bold')
    ax2.set_xticks(x + width)
    ax2.set_xticklabels(metrics)
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. NPV Comparison
    ax3 = fig.add_subplot(gs[1, 0])
    npvs = [scenarios[s]['financial']['npv']/1e9 for s in scenarios_list]
    bars = ax3.bar(scenarios_list, npvs, color=[colors[s] for s in scenarios_list])
    ax3.set_ylabel('NPV ($B)', fontsize=12)
    ax3.set_title('Net Present Value (8% discount, 30 years)', fontsize=14, fontweight='bold')
    ax3.set_xticklabels(['Conservative', 'Moderate', 'Optimistic'])
    ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax3.grid(axis='y', alpha=0.3)
    
    for bar, npv in zip(bars, npvs):
        ax3.text(bar.get_x() + bar.get_width()/2, npv + 0.1, f'${npv:.1f}B', 
                ha='center', fontsize=11, fontweight='bold')
    
    # 4. Cumulative Cash Flow
    ax4 = fig.add_subplot(gs[1, 1])
    for scenario in scenarios_list:
        cf = scenarios[scenario]['financial']['cash_flows']
        years = [c['year'] for c in cf]
        cumulative = [c['cumulative']/1e9 for c in cf]
        ax4.plot(years, cumulative, label=scenario.capitalize(), color=colors[scenario], linewidth=2)
    
    ax4.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax4.set_xlabel('Year', fontsize=12)
    ax4.set_ylabel('Cumulative Cash Flow ($B)', fontsize=12)
    ax4.set_title('Cumulative Cash Flow Over Time', fontsize=14, fontweight='bold')
    ax4.legend()
    ax4.grid(alpha=0.3)
    
    # 5. Operating Cost Breakdown (for moderate scenario)
    ax5 = fig.add_subplot(gs[2, 0])
    mod_opex = scenarios['moderate']['opex']
    opex_labels = ['Labor', 'Utilities', 'Maintenance', 'Insurance', 'Marketing', 'Other']
    opex_values = [
        mod_opex['labor_cost']/1e6,
        mod_opex['utility_cost']/1e6,
        mod_opex['maintenance_cost']/1e6,
        mod_opex['insurance_cost']/1e6,
        mod_opex['marketing_cost']/1e6,
        (mod_opex['property_tax'] + mod_opex['overhead_cost'])/1e6
    ]
    opex_colors = ['#3498DB', '#E74C3C', '#F39C12', '#9B59B6', '#1ABC9C', '#95A5A6']
    
    wedges, texts, autotexts = ax5.pie(opex_values, labels=opex_labels, colors=opex_colors,
                                        autopct='%1.1f%%', startangle=90)
    ax5.set_title(f'Operating Costs Breakdown (Moderate)\nTotal: ${mod_opex["total_opex"]/1e6:.0f}M/year',
                  fontsize=14, fontweight='bold')
    
    # 6. Revenue vs Cost Summary Table
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.axis('off')
    
    # Create summary table
    table_data = [
        ['Metric', 'Conservative', 'Moderate', 'Optimistic'],
        ['Annual Revenue', f'${scenarios["conservative"]["total_revenue"]/1e6:.0f}M',
         f'${scenarios["moderate"]["total_revenue"]/1e6:.0f}M',
         f'${scenarios["optimistic"]["total_revenue"]/1e6:.0f}M'],
        ['Operating Costs', f'${scenarios["conservative"]["opex"]["total_opex"]/1e6:.0f}M',
         f'${scenarios["moderate"]["opex"]["total_opex"]/1e6:.0f}M',
         f'${scenarios["optimistic"]["opex"]["total_opex"]/1e6:.0f}M'],
        ['EBITDA', f'${scenarios["conservative"]["financial"]["first_year_ebitda"]/1e6:.0f}M',
         f'${scenarios["moderate"]["financial"]["first_year_ebitda"]/1e6:.0f}M',
         f'${scenarios["optimistic"]["financial"]["first_year_ebitda"]/1e6:.0f}M'],
        ['Payback Period', f'{scenarios["conservative"]["financial"]["payback_years"]:.1f} years',
         f'{scenarios["moderate"]["financial"]["payback_years"]:.1f} years',
         f'{scenarios["optimistic"]["financial"]["payback_years"]:.1f} years'],
        ['IRR', f'{scenarios["conservative"]["financial"]["irr_pct"]:.1f}%',
         f'{scenarios["moderate"]["financial"]["irr_pct"]:.1f}%',
         f'{scenarios["optimistic"]["financial"]["irr_pct"]:.1f}%'],
        ['30-Year NPV', f'${scenarios["conservative"]["financial"]["npv"]/1e9:.1f}B',
         f'${scenarios["moderate"]["financial"]["npv"]/1e9:.1f}B',
         f'${scenarios["optimistic"]["financial"]["npv"]/1e9:.1f}B'],
    ]
    
    table = ax6.table(cellText=table_data, loc='center', cellLoc='center',
                      colWidths=[0.3, 0.23, 0.23, 0.23])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.8)
    
    # Style header row
    for j in range(4):
        table[(0, j)].set_facecolor('#2C3E50')
        table[(0, j)].set_text_props(color='white', fontweight='bold')
    
    # Style data rows
    for i in range(1, len(table_data)):
        table[(i, 0)].set_facecolor('#ECF0F1')
        table[(i, 0)].set_text_props(fontweight='bold')
        for j in range(1, 4):
            table[(i, j)].set_facecolor('#FFFFFF' if i % 2 == 0 else '#F8F9FA')
    
    ax6.set_title('Financial Summary', fontsize=14, fontweight='bold', pad=20)
    
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nEconomic analysis chart saved to: {save_path}")


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def run_economic_analysis(
    tower_height: float,
    n_eggs: int,
    construction_cost: float,
    elevator_system_cost: float,
    stabilization_system_cost: float,
    save_chart: bool = True
) -> Dict:
    """
    Run complete economic analysis and return results.
    """
    scenarios = run_scenario_analysis(
        tower_height=tower_height,
        n_eggs=n_eggs,
        construction_cost=construction_cost,
        elevator_system_cost=elevator_system_cost,
        stabilization_system_cost=stabilization_system_cost
    )
    
    print_economic_analysis(scenarios, construction_cost)
    
    if save_chart:
        save_economic_analysis_png(scenarios, construction_cost, tower_height)
    
    return scenarios


if __name__ == "__main__":
    # Example usage with typical values from the egg tower analysis
    scenarios = run_economic_analysis(
        tower_height=1600,
        n_eggs=29,
        construction_cost=1_115_000_000,  # $1.115B from cost analysis
        elevator_system_cost=53_500_000,
        stabilization_system_cost=23_900_000
    )
