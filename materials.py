"""
Material classes for egg tower structural analysis.

Engineering Validation Notes:
============================
Material properties are validated against:
- ISO/ASTM testing standards
- Manufacturer datasheets (Toray, DSM Dyneema, DuPont)
- EN 12385-4, API 9A (steel wire rope)
- ACI 440 guidelines (FRP materials)

Key distinctions:
- Fiber strength vs. Cable/Composite strength (Rule of Mixtures applies)
- Material density vs. Bulk/Effective density (Fill Factor applies)
- Static vs. Dynamic loading considerations for Safety Factors
"""
from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class Material:
    """Base class for structural materials."""
    name: str
    density: float  # kg/m³
    tensile_strength: float  # Pa
    compressive_strength: float  # Pa
    youngs_modulus: float  # Pa
    poisson_ratio: float = 0.3
    
    @property
    def tensile_strength_MPa(self) -> float:
        return self.tensile_strength / 1e6
    
    @property
    def compressive_strength_MPa(self) -> float:
        return self.compressive_strength / 1e6
    
    @property
    def youngs_modulus_GPa(self) -> float:
        return self.youngs_modulus / 1e9


@dataclass
class ShellMaterial(Material):
    """Material for egg shell structures."""
    
    def allowable_compressive_stress(self, safety_factor: float = 2.0) -> float:
        """Return allowable compressive stress with safety factor."""
        return self.compressive_strength / safety_factor
    
    def shell_buckling_stress(self, R: float, t: float, knockdown: float = 0.25) -> float:
        """
        Critical buckling stress for curved shell.
        
        Args:
            R: Radius of curvature (m)
            t: Shell thickness (m)
            knockdown: Knockdown factor for imperfections (default 0.25)
        """
        sigma_cr = (2 * self.youngs_modulus * t) / (R * np.sqrt(3 * (1 - self.poisson_ratio**2)))
        return knockdown * sigma_cr


@dataclass  
class CableMaterial:
    """
    Material for cable/rope elements.
    
    Engineering Notes:
    - tensile_strength: Ultimate tensile strength of the cable material.
      For composites, this is the fiber-dominated longitudinal strength.
    - density: Material density (not bulk density). Effective mass calculations
      must account for fill_factor: ρ_eff = ρ_material × fill_factor
    - safety_factor: Design factor relating working load to breaking load.
      Higher values required for dynamic lifting (5.0) vs static guying (2.5-3.0).
    - fill_factor: Ratio of metallic/fiber area to gross circumscribed area.
      Ranges from 0.46 (fiber core rope) to 0.80 (parallel wire strand).
    """
    name: str
    tensile_strength: float  # Pa (Ultimate Breaking Stress)
    density: float  # kg/m³ (material density, not bulk)
    safety_factor: float = 2.0  # Design factor: MBL / Working Load
    fill_factor: float = 0.8  # A_effective / A_gross (packing efficiency)
    fiber_diameter: Optional[float] = None  # m, for laced configurations
    
    @property
    def tensile_strength_MPa(self) -> float:
        return self.tensile_strength / 1e6
    
    @property
    def allowable_stress(self) -> float:
        """Allowable tensile stress with safety factor."""
        return self.tensile_strength / self.safety_factor
    
    @property
    def allowable_stress_MPa(self) -> float:
        return self.allowable_stress / 1e6
    
    def required_area(self, tension: float) -> float:
        """Required effective cross-sectional area for given tension."""
        return tension / self.allowable_stress
    
    def gross_area(self, tension: float) -> float:
        """Required gross area accounting for fill factor."""
        return self.required_area(tension) / self.fill_factor
    
    def cable_diameter(self, tension: float) -> float:
        """Required cable diameter for given tension."""
        A = self.gross_area(tension)
        return 2 * np.sqrt(A / np.pi)
    
    def cable_mass(self, tension: float, length: float) -> float:
        """Mass of cable for given tension and length."""
        return self.gross_area(tension) * length * self.density
    
    def fiber_count(self, tension: float) -> int:
        """Number of individual fibers needed for laced configuration."""
        if self.fiber_diameter is None:
            raise ValueError("fiber_diameter not set for this material")
        A_fiber = np.pi * (self.fiber_diameter / 2)**2 * self.fill_factor
        return int(np.ceil(self.required_area(tension) / A_fiber))
    
    def fiber_breaking_load(self) -> float:
        """Breaking load per individual fiber (N)."""
        if self.fiber_diameter is None:
            raise ValueError("fiber_diameter not set for this material")
        A_fiber = np.pi * (self.fiber_diameter / 2)**2 * self.fill_factor
        return A_fiber * self.tensile_strength
    
    def fiber_working_load(self) -> float:
        """Working load per fiber with safety factor (N)."""
        return self.fiber_breaking_load() / self.safety_factor
    
    def fiber_mass_per_meter(self) -> float:
        """Mass per meter of individual fiber (kg/m)."""
        if self.fiber_diameter is None:
            raise ValueError("fiber_diameter not set for this material")
        A_fiber = np.pi * (self.fiber_diameter / 2)**2
        return A_fiber * self.density


# =============================================================================
# PRE-DEFINED MATERIALS
# =============================================================================

# -----------------------------------------------------------------------------
# SHELL MATERIAL: BFRP (Basalt Fiber Reinforced Polymer)
# -----------------------------------------------------------------------------
# Basalt fibers extruded from molten volcanic rock at >1400°C.
# Superior alkaline resistance vs E-glass; lower carbon footprint.
# Composite strength governed by Rule of Mixtures: σ_c = σ_f × V_f
# Standard grades: 1000-1200 MPa; Premium grades: up to 1400-1500 MPa
# Fiber volume fraction typically 70-80% for pultruded rods.
# Linearly elastic to failure (brittle) - no yielding behavior.
# -----------------------------------------------------------------------------
BFRP = ShellMaterial(
    name="BFRP (Basalt Fiber Reinforced Polymer)",
    density=2000,  # kg/m³ (validated: 1900-2100 typical)
    tensile_strength=1500e6,  # Pa (premium grade; std is 1100-1300 MPa)
    compressive_strength=400e6,  # Pa
    youngs_modulus=45e9,  # Pa (range: 40-60 GPa, softer than steel/carbon)
    poisson_ratio=0.3
)

# -----------------------------------------------------------------------------
# CABLE MATERIAL: Steel Wire Rope (Grade 1770)
# -----------------------------------------------------------------------------
# Metallurgy: Pearlitic microstructure via patenting and cold-drawing.
# Grade 1770 = "Improved Plow Steel" (IPS) per EN 12385-4 / API 9A.
# Higher grades available (1960 EIPS, 2160 EEIPS) but reduced fatigue life.
# 
# Fill Factor 0.60: Requires compacted/swaged rope construction (e.g., Dyform).
# Standard IWRC ropes achieve only 0.57-0.58; fiber core ropes 0.46-0.50.
# 
# Safety Factor 3.0: Valid for STATIC applications (guy wires, pendants).
# UNSAFE for dynamic lifting - OSHA/ASME mandate SF ≥ 5.0 for hoisting.
# 
# Constructional modulus: ~60 GPa initial, ~110 GPa after pre-stretching
# (differs from solid steel modulus of 200 GPa due to helical wire lay).
# Susceptible to galvanic corrosion - requires zinc galvanization outdoors.
# -----------------------------------------------------------------------------
STEEL_CABLE = CableMaterial(
    name="Steel Wire Rope (Grade 1770)",
    tensile_strength=1770e6,  # Pa (Grade 1770 per EN 12385-4)
    density=7850,  # kg/m³ (carbon steel; bulk ρ = 7850 × 0.6 = 4710)
    safety_factor=3.0,  # Static loading only; use 5.0+ for lifting
    fill_factor=0.60,  # Compacted strand construction required
    fiber_diameter=0.100  # 100mm heavy structural wire rope
)

# -----------------------------------------------------------------------------
# CABLE MATERIAL: BFRP Cable
# -----------------------------------------------------------------------------
# Pultruded basalt fiber rods or bundled strands.
# 1500 MPa: Premium/aerospace grade (std commercial is 1100-1300 MPa).
# Requires fiber volume fraction ~75-80% with near-perfect alignment.
# 
# Fill Factor 0.70: Represents bundled rods with interstitial potting,
# or solid rod where 70% is the internal fiber content by weight.
# 
# Safety Factor 2.5: Working load = 40% UTS (1/2.5 = 0.4).
# Below creep rupture threshold (50-60% UTS for BFRP).
# Provides buffer for environmental degradation (ACI 440 C_E factors).
# 
# Advantages: Corrosion-free, non-magnetic, dielectric.
# Disadvantages: Lower stiffness than steel (50-60 GPa vs 200 GPa);
# in hybrid systems, steel carries load until yield.
# Best used in pre-stressed applications where strain is locked in.
# -----------------------------------------------------------------------------
BFRP_CABLE = CableMaterial(
    name="BFRP Cable",
    tensile_strength=1500e6,  # Pa (premium pultruded grade)
    density=2000,  # kg/m³ (≈25% of steel weight)
    safety_factor=2.5,  # Appropriate for permanent anchors/tendons
    fill_factor=0.70  # Bundled rods or high-fiber-content solid rod
)

# -----------------------------------------------------------------------------
# CABLE MATERIAL: Carbon Fiber (Toray T700S)
# -----------------------------------------------------------------------------
# T700S fiber properties: 4900 MPa tensile, 230 GPa modulus.
# 
# 3500 MPa specification: High-efficiency Parallel Wire Strand (PWS)
# or dry fiber bundle. Standard pultruded rods achieve ~2500-2900 MPa
# with 60% fiber volume fraction.
# 
# Best balance of high stiffness (230 GPa, comparable to steel) and
# high specific strength - preferred for long-span bridge stay cables.
# 
# Safety Factor 2.5: Per EMPA/Japanese CFCC standards (2.2-2.5 typical).
# Exceptional fatigue resistance; main failure mode is anchorage stress.
# Fibers are weak in shear - complex wedge/potting anchorages required.
# 
# CRITICAL: Carbon fiber is ELECTRICALLY CONDUCTIVE.
# Galvanic corrosion risk with aluminum/steel fittings in marine env.
# Requires dielectric isolation (glass fiber sleeves, ceramic coatings).
# 
# Specific Strength: 223 km breaking length (9.7× steel efficiency).
# -----------------------------------------------------------------------------
CARBON_FIBER_CABLE = CableMaterial(
    name="Carbon Fiber Cable (T700S)",
    tensile_strength=3500e6,  # Pa (high-efficiency PWS; std rod ~2800 MPa)
    density=1600,  # kg/m³ (60% fiber/40% resin composite)
    safety_factor=2.5,  # Per international bridge cable standards
    fill_factor=0.65  # Parallel wire strand or pultruded rod
)

# -----------------------------------------------------------------------------
# CABLE MATERIAL: Dyneema SK78 (UHMWPE)
# -----------------------------------------------------------------------------
# Ultra-High Molecular Weight Polyethylene - highest specific strength
# of any industrial fiber. SK78 grade: improved creep vs SK75.
# 
# 3600 MPa: Validated fiber/material strength (38-40 cN/dtex).
# Rope efficiency: 12-strand braid realizes 80-90% of fiber strength.
# 50mm SK78 rope typically breaks at 160-190 tonnes (1600-1900 kN).
# 
# Density 975 kg/m³: THE CABLE FLOATS IN SEAWATER (SG < 1.0).
# Eliminates propeller fouling risk; self-weight penalty negligible.
# 
# ⚠️  CRITICAL WARNING: CREEP RUPTURE RISK ⚠️
# Safety Factor 2.0 = Working load at 50% MBL.
# At this load level, measurable creep occurs - UNSAFE for permanent
# static loads. Rope may fail within weeks/months, especially if T > 25°C.
# 
# SF 2.0 VALID ONLY FOR: Short-duration engineered lifts where weight
# is critical and load duration is minimal.
# For permanent mooring/guying: Use SF 5.0-8.0 minimum.
# For lifting slings: OSHA mandates SF 5.0.
# 
# Fill Factor 0.80: Requires Parallel Core construction with tight jacket
# or heavily pre-stretched/compacted braid. Stiff, difficult to splice.
# 
# Specific Strength: 376 km breaking length (16.4× steel efficiency).
# -----------------------------------------------------------------------------
DYNEEMA = CableMaterial(
    name="Dyneema SK78 (UHMWPE)",
    tensile_strength=3600e6,  # Pa (fiber strength; rope lower due to braid)
    density=975,  # kg/m³ (floats in seawater - SG 0.97)
    safety_factor=2.0,  # ⚠️ SHORT-DURATION USE ONLY - creep risk!
    fill_factor=0.80,  # Parallel core or compacted braid
    fiber_diameter=0.050  # 50mm nominal rope diameter
)

# -----------------------------------------------------------------------------
# CABLE MATERIAL: Aramid (Kevlar 49)
# -----------------------------------------------------------------------------
# Kevlar 49: Structural/aerospace grade (higher modulus than K29 ballistic).
# DuPont datasheet: 2760-2800 MPa typical; 3000 MPa in optimized strands.
# 
# Fill Factor 0.75: Indicates Parafil (Parallel Filament) construction.
# Parallel dry fibers in polymeric sheath; 25-30% air void content.
# Parallel construction eliminates inter-fiber abrasion and maximizes stiffness.
# 
# Safety Factor 2.5: Standard for dielectric antenna guy wires (ADSS cables).
# Sufficient tension to prevent slack while within static fatigue limits.
# Unlike Dyneema, aramid does NOT suffer significant creep at room temperature.
# 
# Vulnerabilities:
# - UV Degradation: Rapid strength loss in sunlight; sheathing required.
# - Axial Compression Fatigue: Internal kinking if rope goes slack and
#   snap-loaded. Must maintain positive tension at all times.
# 
# Advantages: Dielectric (non-conductive), high modulus, no creep.
# Specific Strength: 212 km breaking length (9.2× steel efficiency).
# -----------------------------------------------------------------------------
ARAMID_CABLE = CableMaterial(
    name="Aramid (Kevlar 49)",
    tensile_strength=3000e6,  # Pa (optimized Parafil strand)
    density=1440,  # kg/m³
    safety_factor=2.5,  # Standard for static dielectric applications
    fill_factor=0.75,  # Parafil parallel filament construction
    fiber_diameter=0.075  # 75mm Parafil rope (standard heavy-duty size)
)


# =============================================================================
# MATERIAL REGISTRY AND DEFAULTS
# =============================================================================
# Comparative Performance (Specific Strength / Breaking Length):
# ┌────────────────┬──────────┬──────────┬───────────────┬──────────────┐
# │ Material       │ σ (MPa)  │ ρ (kg/m³)│ Break Len (km)│ vs Steel     │
# ├────────────────┼──────────┼──────────┼───────────────┼──────────────┤
# │ Steel Wire     │ 1770     │ 7850     │ 23.0          │ 1.0×         │
# │ BFRP           │ 1500     │ 2000     │ 76.5          │ 3.3×         │
# │ Aramid K49     │ 3000     │ 1440     │ 212.4         │ 9.2×         │
# │ Carbon T700    │ 3500     │ 1600     │ 223.0         │ 9.7×         │
# │ Dyneema SK78   │ 3600     │ 975      │ 376.5         │ 16.4×        │
# └────────────────┴──────────┴──────────┴───────────────┴──────────────┘
# 
# Breaking Length = σ / (ρ × g) [km] - theoretical max cable length
# supporting its own weight before failure.
# =============================================================================

CABLE_MATERIALS = {
    'steel': STEEL_CABLE,
    'bfrp': BFRP_CABLE,
    'carbon': CARBON_FIBER_CABLE,
    'dyneema': DYNEEMA,
    'aramid': ARAMID_CABLE
}

# Default materials
DEFAULT_SHELL_MATERIAL = BFRP
DEFAULT_CABLE_MATERIAL = DYNEEMA
