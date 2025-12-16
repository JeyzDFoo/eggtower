# Egg Tower Structural Optimization - Copilot Instructions

## Project Overview
Optimize a 1600m tall tower made of stacked hollow egg-shaped shells made from Basalt Fiber Reinforced Polymer (BFRP). Each egg is connected to adjacent eggs via guy cables at their maximum diameter points.

## Key Parameters

### Material Properties (BFRP)
- Density: 2000 kg/m³
- Compressive strength: 400 MPa
- Tensile strength: 1500 MPa
- Young's modulus: 40-50 GPa

### Design Constraints
- Total tower height: 1600m (n × h_egg = 1600)
- Egg aspect ratio: h:d between 2:1 and 8:1 (height to max diameter)
- Shell thickness: scales with egg size, typically t = h/800 to h/1200
- Minimum safety factor: 2.0 for compressive stress

## Physics to Model

### Geometry
- Model eggs as ellipsoids with semi-axes: a = h/2 (vertical), b = c = d/2 (horizontal)
- Shell tapers from base to tip (thicker at base)
- Maximum diameter occurs at ~60% of height from base

### Structural Analysis
- **Self-weight stress**: Bottom egg supports (n-1) eggs above it
- **Base stress**: σ = Total_weight / Base_cross_sectional_area
- **Shell volume**: Surface_area × average_thickness
- **Mass per egg**: density × shell_volume

### Key Calculations
```python
# Base cross-sectional area (annular ring)
A_base = π × (r_outer² - r_inner²)

# Compressive stress at base
σ = (n-1) × mass_per_egg × g / A_base

# Safety factor
SF = σ_max / σ
```

## Code Structure Expectations

### Functions should:
- Use SI units (meters, kg, Pascals) consistently
- Return dictionaries with clear keys for results
- Include docstrings explaining geometry assumptions
- Handle edge cases (n=1, very small/large eggs)

### Variable naming:
- `n` = number of eggs
- `h` = individual egg height (m)
- `d` = egg max diameter (m)
- `t_base`, `t_tip` = shell thickness at base/tip (m)
- `sigma` = stress (Pa)
- `rho` = density (kg/m³)
- Use `_MPa` suffix when converting to megapascals

### Optimization goals:
1. **Primary**: Ensure base stress < 400 MPa with SF ≥ 2
2. **Secondary**: Minimize total material volume
3. **Tertiary**: Balance between buildability (fewer eggs) and efficiency

## Wind Loading (Future Enhancement)
- Wind speed increases with height: v(z) = v_ref × (z/10)^0.16
- Dynamic pressure: q = 0.5 × ρ_air × v²
- At 1600m: ~70 m/s wind speed, ~3 kPa pressure
- Guy cables must resist lateral loads

## Expected Outputs
- Parameter sweep tables showing n, h, d, stress, safety factor
- Plots: stress vs n, safety factor vs n, material volume vs n
- Identification of optimal configuration
- Visualization of tower profile (optional)

## Engineering Principles
- Egg shape provides excellent buckling resistance due to double curvature
- Slender structures (high aspect ratio) are more material-efficient but require more eggs
- Trade-off between fewer joints (fewer eggs) vs lower individual weights (more eggs)
- Shell thickness must be sufficient to prevent local buckling

## Testing Approach
- Verify stress calculations against hand calculations for n=2 case
- Check that n=1 gives zero stress on base (no eggs above)
- Ensure all stresses stay below material limits
- Validate that results are physically reasonable (no negative masses, etc.)