"""
====================================================================================
  BIOPCL BIOPRINTED MANDIBLE -- MATHEMATICAL MODEL & BIOMECHANICAL SIMULATION
  Anatomical Construct Analysis | Gibson-Ashby Porosity | FEA Nodal Mapping
====================================================================================

MATERIAL SPECIFICATION (User Required Final BioPCL):
  - Density                       : 1200 kg/m^3 (1.2 g/cm^3)
  - Young's Modulus               : 350 MPa (0.35 GPa)
  - Poisson's Ratio               : 0.30
  - Tensile Yield Strength        : 15 MPa
  - Tensile Ultimate Strength     : 23 MPa
  - Compressive Yield Strength    : 15 MPa
  - Compressive Ultimate Strength : 23 MPa

COMPUTATIONAL FRAMEWORK:
  1. FEA Data Transformation: Map 55,226 ANSYS simulation nodes (ABS baseline)
     to BioPCL stress/strain/safety factor space via linear elastic scaling.
  2. Gibson-Ashby Bioprinted Scaffold Mechanics: Porosity-dependent modulus
     E(P) and yield strength sigma_y(P).
  3. Frost's Mechanostat Evaluation: Tissue-level microstrain classification
     (disuse, homeostasis, osteoinduction, overload).
  4. 7-Zone Anatomical Mandible Construct Analysis under Physiological
     Occlusal Loads (bilateral molar bite, 600 N total force).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# SECTION 1 -- MATERIAL DATABASE & GIBSON-ASHBY POROSITY MODEL
# ==============================================================================

class BioPCLMaterial:
    """
    BioPCL (Polycaprolactone) Bioprinted Biomaterial formulation.
    Includes solid properties and Gibson-Ashby cellular scaffold models.

    References:
      - Zein et al. (2002), Biomaterials 23(4), 1169-1185
      - Hutmacher (2000), Biomaterials 21(24), 2529-2543
      - Gibson & Ashby (1997), Cellular Solids, Cambridge Univ. Press
    """
    def __init__(self,
                 density_kg_m3=1200.0,
                 E_solid_mpa=350.0,
                 nu=0.30,
                 sigma_yield_t=15.0,
                 sigma_ult_t=23.0,
                 sigma_yield_c=15.0,
                 sigma_ult_c=23.0):
        self.rho_solid = density_kg_m3          # kg/m^3
        self.E_solid = E_solid_mpa              # MPa
        self.nu = nu
        self.sigma_yt = sigma_yield_t           # MPa (tensile yield)
        self.sigma_ut = sigma_ult_t             # MPa (tensile ultimate)
        self.sigma_yc = sigma_yield_c           # MPa (compressive yield)
        self.sigma_uc = sigma_ult_c             # MPa (compressive ultimate)
        self.G_solid = self.E_solid / (2.0 * (1.0 + self.nu))  # Shear modulus

    def effective_properties(self, porosity, n_E=2.0, n_sigma=1.5):
        """
        Gibson-Ashby model for bioprinted porous scaffolds:
          E(P) / E_0     = (1 - P)^n_E       (n_E ~ 2.0 for open-cell PCL)
          sigma_y(P)/sy0 = (1 - P)^n_sigma    (n_sigma ~ 1.5)
          rho(P)         = rho_0 * (1 - P)
        """
        P = np.clip(porosity, 0.0, 0.85)
        rel_density = 1.0 - P
        E_eff = self.E_solid * (rel_density ** n_E)
        sigma_yt_eff = self.sigma_yt * (rel_density ** n_sigma)
        sigma_ut_eff = self.sigma_ut * (rel_density ** n_sigma)
        rho_eff = self.rho_solid * rel_density
        return {
            "porosity": P,
            "rel_density": rel_density,
            "E_eff_MPa": E_eff,
            "sigma_yield_MPa": sigma_yt_eff,
            "sigma_ult_MPa": sigma_ut_eff,
            "density_kg_m3": rho_eff
        }

    def tissue_microstrain(self, scaffold_stress_mpa, porosity=0.0,
                            bone_volume_fraction=0.30):
        """
        Frost's Mechanostat Theory — Tissue-Level Microstrain Estimation.

        The Frost mechanostat applies to BONE TISSUE growing into the
        scaffold, not to the polymer itself. We use a Voigt (upper-bound)
        composite modulus to estimate the strain experienced by
        regenerating bone at the scaffold-tissue interface:

          E_composite = BVF * E_bone + (1-BVF) * E_scaffold(P)

        where BVF = bone volume fraction (typically 0.20-0.40 at 12 weeks
        post-implantation for PCL scaffolds, per Cipitria et al. 2012).

        Microstrain = (sigma / E_composite) * 1e6

        Frost classification thresholds (for bone tissue):
          < 50 ue       : Disuse / bone resorption
          50 - 1500 ue  : Physiological homeostasis
          1500 - 3000 ue: Mild overload / osteoinductive (bone formation)
          > 3000 ue     : Pathologic overload / microdamage

        Reference: Frost (2003), J Musculoskelet Neuronal Interact 3(2)
        """
        props = self.effective_properties(porosity)
        E_scaffold = props["E_eff_MPa"]
        E_bone = 15000.0  # MPa — representative cortical bone modulus
        E_composite = bone_volume_fraction * E_bone + (1.0 - bone_volume_fraction) * E_scaffold
        microstrain = (scaffold_stress_mpa / E_composite) * 1e6
        return microstrain

    def classify_mechanostat(self, microstrain_val):
        """Return Frost mechanostat classification string."""
        if microstrain_val < 50:
            return "Disuse/Resorption"
        elif microstrain_val < 1500:
            return "Homeostasis"
        elif microstrain_val < 3000:
            return "Osteoinductive"
        else:
            return "Pathologic Overload"


BIOPCL = BioPCLMaterial()

# ==============================================================================
# SECTION 2 -- FEA DATA PARSER & BIOPCL TRANSFORMATION
# ==============================================================================

def load_and_transform_fea_data(data_dir="data", E_abs_mpa=2000.0):
    """
    Reads ANSYS export files and transforms from ABS baseline to BioPCL.

    The ANSYS tabular_stress_data.csv contains tab-separated:
      Column 1: Node Number (integer)
      Column 2: Equivalent (von-Mises) Stress in Pascals

    For linear static structural analysis with prescribed force BCs:
      - Stress field sigma(x) is invariant under material change
        (sigma depends on geometry + loads, not on E when Poisson's ratio
        change is small: nu_ABS=0.35 vs nu_BioPCL=0.30 -> <1.5% deviation)
      - Displacement field scales as: u_BioPCL = u_ABS * (E_ABS / E_BioPCL)
    """
    stress_file = os.path.join(data_dir, "tabular_stress_data.csv")
    vol_file = os.path.join(data_dir, "volume.txt")

    # 1. Read nodal von-Mises stresses (Pa -> MPa)
    node_ids = []
    stresses_mpa = []
    with open(stress_file, "r") as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue  # skip header
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    node_ids.append(int(parts[0]))
                    stresses_mpa.append(float(parts[1]) / 1e6)
                except ValueError:
                    continue

    stresses_mpa = np.array(stresses_mpa)
    node_ids = np.array(node_ids)

    # 2. Read volume from ANSYS export
    volume_m3 = 2.427e-5  # default fallback from volume.txt
    if os.path.exists(vol_file):
        with open(vol_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                for p in parts:
                    try:
                        val = float(p)
                        if 1e-6 < val < 1e-3:
                            volume_m3 = val
                    except ValueError:
                        pass

    # 3. Deformation values from Total_deformation.xls (ANSYS ABS baseline)
    max_def_abs_m = 2.4829e-5   # Maximum total deformation (m)
    avg_def_abs_m = 5.9981e-6   # Average total deformation (m)

    # 4. Scale deformations to BioPCL (linear elasticity: u ~ 1/E)
    modulus_ratio = E_abs_mpa / BIOPCL.E_solid   # = 2000/350 = 5.714
    max_def_biopcl_mm = (max_def_abs_m * 1000.0) * modulus_ratio
    avg_def_biopcl_mm = (avg_def_abs_m * 1000.0) * modulus_ratio

    # 5. Nodal Safety Factors for Solid BioPCL (no porosity)
    sf_yield = BIOPCL.sigma_yt / np.maximum(stresses_mpa, 1e-6)
    sf_ult = BIOPCL.sigma_ut / np.maximum(stresses_mpa, 1e-6)

    safe_mask = stresses_mpa <= BIOPCL.sigma_yt
    yield_mask = (stresses_mpa > BIOPCL.sigma_yt) & (stresses_mpa <= BIOPCL.sigma_ut)
    fail_mask = stresses_mpa > BIOPCL.sigma_ut

    # 6. Construct mass
    mass_solid_g = (volume_m3 * BIOPCL.rho_solid) * 1000.0  # grams

    return {
        "node_ids": node_ids,
        "stresses_mpa": stresses_mpa,
        "volume_m3": volume_m3,
        "volume_cm3": volume_m3 * 1e6,
        "mass_solid_g": mass_solid_g,
        "max_def_abs_mm": max_def_abs_m * 1000.0,
        "avg_def_abs_mm": avg_def_abs_m * 1000.0,
        "max_def_biopcl_mm": max_def_biopcl_mm,
        "avg_def_biopcl_mm": avg_def_biopcl_mm,
        "modulus_ratio": modulus_ratio,
        "sf_yield": sf_yield,
        "sf_ult": sf_ult,
        "counts": {
            "total": len(stresses_mpa),
            "safe": int(np.sum(safe_mask)),
            "yield": int(np.sum(yield_mask)),
            "exceed": int(np.sum(fail_mask))
        },
        "stats": {
            "max": float(np.max(stresses_mpa)),
            "min": float(np.min(stresses_mpa)),
            "mean": float(np.mean(stresses_mpa)),
            "median": float(np.median(stresses_mpa)),
            "std": float(np.std(stresses_mpa)),
            "p95": float(np.percentile(stresses_mpa, 95.0)),
            "p99": float(np.percentile(stresses_mpa, 99.0))
        }
    }


# ==============================================================================
# SECTION 3 -- 7-ZONE ANATOMICAL BIOPCL MANDIBLE CONSTRUCT MODEL
# ==============================================================================

class BioPCLZoneConstruct:
    """
    Represents one anatomical zone of the bioprinted mandible construct.

    The cross-section properties (A, Ixx) are derived from the original
    mandible geometry anatomical measurements (Korioth & Hannam 1994).
    The design porosity varies per zone for optimized mechanical + biological
    performance.
    """
    def __init__(self, name, arc_pos, cross_area_mm2, Ixx_mm4,
                 design_porosity=0.20, color='#00d4ff'):
        self.name = name
        self.arc_pos = arc_pos     # normalized [0=condyle, 1=symphysis]
        self.A = cross_area_mm2    # cross-section area (mm^2)
        self.Ixx = Ixx_mm4         # second moment of area (mm^4)
        self.porosity = design_porosity
        self.color = color
        self.props = BIOPCL.effective_properties(design_porosity)
        # c = distance from neutral axis to extreme fiber
        # For mandibular cross-sections (roughly elliptical), approximate
        # height h from I_xx ~ A*h^2/12 -> h = sqrt(12*I/A), c = h/2
        self.c_neutral = np.sqrt(3.0 * Ixx_mm4 / cross_area_mm2)  # mm

    def compute_response(self, M_nmm, V_n, F_axial_n):
        """
        Compute stress state at this zone under given internal forces.

        Parameters:
          M_nmm    : Bending moment at this zone (N.mm)
          V_n      : Transverse shear force (N)
          F_axial_n: Axial force component (N)

        Returns dict with von Mises stress, safety factors, microstrain.
        """
        E_eff = self.props["E_eff_MPa"]
        sigma_y = self.props["sigma_yield_MPa"]
        sigma_u = self.props["sigma_ult_MPa"]

        # Bending stress: sigma_b = M*c / I
        sigma_bend = (M_nmm * self.c_neutral) / self.Ixx

        # Axial stress: sigma_a = F / A
        sigma_axial = abs(F_axial_n) / self.A

        # Total normal stress (bending + axial, same direction)
        sigma_total = sigma_bend + sigma_axial

        # Shear stress (Jourawski/Timoshenko approximation)
        # tau = V*Q / (I*b), where Q ~ 0.4*I/c and b ~ A/(2c)
        # This simplifies to tau = 0.8 * V / A (rectangular section approx.)
        Q_approx = self.Ixx * 0.40 / self.c_neutral
        b_width = self.A / (2.0 * self.c_neutral)
        tau_shear = (V_n * Q_approx) / (self.Ixx * b_width)

        # Von Mises equivalent stress
        sigma_vm = np.sqrt(sigma_total**2 + 3.0 * tau_shear**2)

        # Safety factors
        sf_yield = sigma_y / max(sigma_vm, 1e-6)
        sf_ult = sigma_u / max(sigma_vm, 1e-6)

        # Microstrain at scaffold-tissue interface (using composite modulus)
        microstrain = BIOPCL.tissue_microstrain(sigma_vm, self.porosity)

        return {
            "zone": self.name,
            "porosity": self.porosity,
            "E_eff_MPa": E_eff,
            "sigma_vm_MPa": sigma_vm,
            "sigma_bend_MPa": sigma_bend,
            "sigma_axial_MPa": sigma_axial,
            "tau_shear_MPa": tau_shear,
            "SF_yield": sf_yield,
            "SF_ult": sf_ult,
            "microstrain": microstrain,
            "mechanostat": BIOPCL.classify_mechanostat(microstrain)
        }


ZONES_BIOPCL = [
    BioPCLZoneConstruct("Condylar Head",    0.00,  95.0,  280.0,
                         design_porosity=0.15, color='#00d4ff'),
    BioPCLZoneConstruct("Condylar Neck",    0.12,  80.0,  210.0,
                         design_porosity=0.10, color='#0099cc'),
    BioPCLZoneConstruct("Ramus",            0.30, 320.0, 1850.0,
                         design_porosity=0.30, color='#ff6b35'),
    BioPCLZoneConstruct("Mandibular Angle", 0.45, 380.0, 2200.0,
                         design_porosity=0.25, color='#ffcc00'),
    BioPCLZoneConstruct("Molar Region",     0.62, 290.0, 1600.0,
                         design_porosity=0.20, color='#cc44ff'),
    BioPCLZoneConstruct("Premolar Region",  0.80, 260.0, 1400.0,
                         design_porosity=0.30, color='#00ffaa'),
    BioPCLZoneConstruct("Symphysis",        1.00, 220.0, 1100.0,
                         design_porosity=0.20, color='#39ff14'),
]

# Physiological load case: Bilateral Molar Bite (600 N total)
# Bending moments computed from curved beam model (Korioth & Hannam 1994)
# M(s) at each zone, V(s) shear, and axial force from muscle resultants
# These are derived from the force/moment equilibrium along the mandible arc
# under 600 N bilateral molar occlusion with masseter + temporalis activation.
ZONE_LOADS = {
    "Condylar Head":    {"M": 4800.0, "V": 180.0, "F_ax": 120.0},
    "Condylar Neck":    {"M": 6200.0, "V": 220.0, "F_ax": 150.0},
    "Ramus":            {"M": 9800.0, "V": 310.0, "F_ax": 200.0},
    "Mandibular Angle": {"M":12500.0, "V": 380.0, "F_ax": 250.0},
    "Molar Region":     {"M":15200.0, "V": 420.0, "F_ax": 300.0},
    "Premolar Region":  {"M":11800.0, "V": 350.0, "F_ax": 220.0},
    "Symphysis":        {"M": 8500.0, "V": 280.0, "F_ax": 180.0},
}


# ==============================================================================
# SECTION 4 -- PLOTTING & PUBLICATION SUITE
# ==============================================================================

plt.style.use('dark_background')
DARK_BG  = '#050a0f'
PANEL_BG = '#0b1420'
GRID_C   = '#0f2035'
TEXT_C   = '#c8ddf0'
DIM_C    = '#5a7fa0'
ACCENT   = '#00d4ff'

def apply_dark_style(ax, xlabel='', ylabel='', title=''):
    ax.set_facecolor(PANEL_BG)
    ax.grid(True, color=GRID_C, linewidth=0.6, alpha=0.8)
    ax.tick_params(colors=DIM_C, labelsize=9)
    for spine in ax.spines.values():
        spine.set_edgecolor('#1a3a5c')
    if xlabel: ax.set_xlabel(xlabel, color=DIM_C, fontsize=10)
    if ylabel: ax.set_ylabel(ylabel, color=DIM_C, fontsize=10)
    if title:  ax.set_title(title, color=ACCENT, fontsize=11, fontweight='bold', pad=10)


def generate_all_figures(fea_data, zone_results):
    """Generate all publication-quality figures."""

    # ==================================================================
    # FIGURE 1: Nodal Stress Distribution & Safety Classification
    # ==================================================================
    print("  [1/4] Nodal Stress Distribution & Safety Classification...")
    fig1 = plt.figure(figsize=(18, 7), facecolor=DARK_BG)
    gs1 = gridspec.GridSpec(1, 2, figure=fig1, wspace=0.30,
                            left=0.06, right=0.96, top=0.88, bottom=0.12)

    # Left: Histogram with BioPCL yield/ultimate thresholds
    ax1 = fig1.add_subplot(gs1[0, 0])
    stresses = fea_data["stresses_mpa"]
    ax1.hist(stresses, bins=80, range=(0, 45), color='#0099cc',
             edgecolor='white', linewidth=0.3, alpha=0.85)
    ax1.axvline(x=15.0, color='#39ff14', linestyle='--', linewidth=2.0,
                label=f'Yield Strength ({BIOPCL.sigma_yt} MPa)')
    ax1.axvline(x=23.0, color='#ff3366', linestyle='--', linewidth=2.0,
                label=f'Ultimate Strength ({BIOPCL.sigma_ut} MPa)')
    ax1.axvline(x=fea_data["stats"]["mean"], color='#ffcc00', linestyle='-',
                linewidth=1.5, label=f'Mean ({fea_data["stats"]["mean"]:.2f} MPa)')
    ax1.axvline(x=fea_data["stats"]["p95"], color='#cc44ff', linestyle=':',
                linewidth=1.5, label=f'95th pct ({fea_data["stats"]["p95"]:.2f} MPa)')
    ax1.legend(fontsize=8, labelcolor=TEXT_C, loc='upper right')
    apply_dark_style(ax1, xlabel='Von Mises Stress (MPa)', ylabel='Node Count',
                     title='ANSYS Nodal Stress Distribution (BioPCL Thresholds)')

    # Right: Safety classification bar chart
    ax2 = fig1.add_subplot(gs1[0, 1])
    categories = ['Safe Zone\n(<= 15 MPa)', 'Yield Zone\n(15-23 MPa)',
                  'Stress Conc.\n(> 23 MPa)']
    counts = [fea_data["counts"]["safe"], fea_data["counts"]["yield"],
              fea_data["counts"]["exceed"]]
    tot = fea_data["counts"]["total"]
    pcts = [c / tot * 100.0 for c in counts]
    colors = ['#39ff14', '#ffcc00', '#ff3366']

    bars = ax2.bar(categories, pcts, color=colors, edgecolor='white',
                   linewidth=0.8, width=0.50)
    for bar, pct, count in zip(bars, pcts, counts):
        ax2.text(bar.get_x() + bar.get_width()/2, pct + 2.0,
                 f'{pct:.2f}%\n({count:,} nodes)',
                 ha='center', va='bottom', color=TEXT_C, fontsize=9,
                 fontweight='bold')
    ax2.set_ylim(0, 115)
    apply_dark_style(ax2, ylabel='Percentage of Nodes (%)',
                     title='Structural Integrity Classification (BioPCL)')

    fig1.suptitle('BIOPCL BIOPRINTED MANDIBLE -- FEA NODAL STRESS ANALYSIS',
                  color=ACCENT, fontsize=13, fontweight='bold', y=0.97)
    plt.savefig("fig1_biopcl_nodal_stress.png", dpi=180,
                bbox_inches='tight', facecolor=DARK_BG)
    plt.close(fig1)

    # ==================================================================
    # FIGURE 2: Gibson-Ashby Porosity Sweep
    # ==================================================================
    print("  [2/4] Gibson-Ashby Bioprinted Porosity Analysis...")
    fig2 = plt.figure(figsize=(18, 7), facecolor=DARK_BG)
    gs2 = gridspec.GridSpec(1, 2, figure=fig2, wspace=0.35,
                            left=0.07, right=0.95, top=0.88, bottom=0.12)

    P_sweep = np.linspace(0.0, 0.60, 100)
    E_sweep = np.array([BIOPCL.effective_properties(p)["E_eff_MPa"]
                        for p in P_sweep])
    sy_sweep = np.array([BIOPCL.effective_properties(p)["sigma_yield_MPa"]
                         for p in P_sweep])
    su_sweep = np.array([BIOPCL.effective_properties(p)["sigma_ult_MPa"]
                         for p in P_sweep])
    mass_sweep = np.array([
        (fea_data["volume_cm3"] * BIOPCL.effective_properties(p)["density_kg_m3"]) / 1000.0
        for p in P_sweep
    ])

    # Left: Modulus & Yield Strength vs Porosity
    ax3 = fig2.add_subplot(gs2[0, 0])
    ax3t = ax3.twinx()
    ax3.plot(P_sweep * 100, E_sweep, color='#00d4ff', linewidth=2.5,
             label='Effective Modulus E(P)')
    ax3t.plot(P_sweep * 100, sy_sweep, color='#ff6b35', linewidth=2.5,
              linestyle='--', label='Yield Strength sy(P)')
    ax3t.plot(P_sweep * 100, su_sweep, color='#ff3366', linewidth=2.0,
              linestyle=':', label='Ultimate Strength su(P)')
    ax3.set_ylabel("Young's Modulus E (MPa)", color='#00d4ff', fontsize=10)
    ax3t.set_ylabel('Strength (MPa)', color='#ff6b35', fontsize=10)
    ax3.tick_params(axis='y', colors='#00d4ff')
    ax3t.tick_params(axis='y', colors='#ff6b35')
    # Mark zone design porosities
    for z in ZONES_BIOPCL:
        zp = z.porosity * 100
        ax3.axvline(x=zp, color=z.color, linewidth=0.8, alpha=0.5, linestyle=':')
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3t.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, fontsize=8,
               labelcolor=TEXT_C, loc='upper right')
    apply_dark_style(ax3, xlabel='Bioprinted Scaffold Porosity P (%)',
                     title='Gibson-Ashby Mechanical Properties vs. Porosity')

    # Right: Construct mass vs porosity
    ax4 = fig2.add_subplot(gs2[0, 1])
    ax4.plot(P_sweep * 100, mass_sweep, color='#39ff14', linewidth=2.5)
    ax4.fill_between(P_sweep * 100, mass_sweep, alpha=0.15, color='#39ff14')
    for z in ZONES_BIOPCL:
        zp = z.porosity
        zm = (fea_data["volume_cm3"] * BIOPCL.effective_properties(zp)["density_kg_m3"]) / 1000.0
        ax4.plot(zp * 100, zm, 'o', color=z.color, markersize=8, zorder=5)
        ax4.annotate(z.name.split()[0], xy=(zp * 100, zm),
                     xytext=(zp * 100 + 2, zm + 0.3),
                     color=z.color, fontsize=7,
                     arrowprops=dict(arrowstyle='->', color=z.color, lw=0.6))
    apply_dark_style(ax4, xlabel='Bioprinted Scaffold Porosity P (%)',
                     ylabel='Construct Mass (grams)',
                     title='BioPCL Construct Mass vs. Scaffold Porosity')

    fig2.suptitle('GIBSON-ASHBY BIOPRINTED SCAFFOLD OPTIMIZATION',
                  color=ACCENT, fontsize=13, fontweight='bold', y=0.97)
    plt.savefig("fig2_biopcl_gibson_ashby.png", dpi=180,
                bbox_inches='tight', facecolor=DARK_BG)
    plt.close(fig2)

    # ==================================================================
    # FIGURE 3: 7-Zone Anatomical Construct Response
    # ==================================================================
    print("  [3/4] 7-Zone Anatomical Construct Response...")
    fig3 = plt.figure(figsize=(20, 14), facecolor=DARK_BG)
    gs3 = gridspec.GridSpec(2, 2, figure=fig3, hspace=0.40, wspace=0.30,
                            left=0.07, right=0.96, top=0.92, bottom=0.06)

    zone_names = [r["zone"] for r in zone_results]
    zone_short = [n.split()[0] for n in zone_names]
    zone_colors = [z.color for z in ZONES_BIOPCL]
    x_idx = np.arange(len(zone_names))

    # (0,0) Von Mises Stress per zone with yield line
    ax5 = fig3.add_subplot(gs3[0, 0])
    vm_vals = [r["sigma_vm_MPa"] for r in zone_results]
    bars5 = ax5.bar(x_idx, vm_vals, color=zone_colors, edgecolor='white',
                    linewidth=0.5, width=0.55)
    ax5.axhline(y=BIOPCL.sigma_yt, color='#39ff14', linestyle='--',
                linewidth=2.0, label=f'BioPCL Yield ({BIOPCL.sigma_yt} MPa)')
    ax5.axhline(y=BIOPCL.sigma_ut, color='#ff3366', linestyle='--',
                linewidth=1.5, alpha=0.7, label=f'BioPCL Ultimate ({BIOPCL.sigma_ut} MPa)')
    for bar, val in zip(bars5, vm_vals):
        ax5.text(bar.get_x() + bar.get_width()/2, val + 0.3,
                 f'{val:.2f}', ha='center', va='bottom',
                 color=TEXT_C, fontsize=8, fontweight='bold')
    ax5.set_xticks(x_idx)
    ax5.set_xticklabels(zone_short, fontsize=8, color=TEXT_C, rotation=30)
    ax5.legend(fontsize=8, labelcolor=TEXT_C)
    apply_dark_style(ax5, ylabel='Von Mises Stress (MPa)',
                     title='Von Mises Stress per Anatomical Zone')

    # (0,1) Safety Factor per zone
    ax6 = fig3.add_subplot(gs3[0, 1])
    sf_vals = [r["SF_yield"] for r in zone_results]
    sf_colors = ['#39ff14' if s > 2.0 else ('#ffcc00' if s > 1.0 else '#ff3366')
                 for s in sf_vals]
    bars6 = ax6.bar(x_idx, sf_vals, color=sf_colors, edgecolor='white',
                    linewidth=0.5, width=0.55)
    ax6.axhline(y=1.0, color='#ff3366', linestyle='--', linewidth=1.5,
                label='SF = 1.0 (Failure)')
    ax6.axhline(y=2.0, color='#39ff14', linestyle='--', linewidth=1.0,
                alpha=0.6, label='SF = 2.0 (Design Target)')
    for bar, val in zip(bars6, sf_vals):
        ax6.text(bar.get_x() + bar.get_width()/2, val + 0.05,
                 f'{val:.2f}', ha='center', va='bottom',
                 color=TEXT_C, fontsize=8, fontweight='bold')
    ax6.set_xticks(x_idx)
    ax6.set_xticklabels(zone_short, fontsize=8, color=TEXT_C, rotation=30)
    ax6.legend(fontsize=8, labelcolor=TEXT_C)
    apply_dark_style(ax6, ylabel='Safety Factor (Yield)',
                     title='Yield Safety Factor per Zone')

    # (1,0) Stress component breakdown
    ax7 = fig3.add_subplot(gs3[1, 0])
    sb = [r["sigma_bend_MPa"] for r in zone_results]
    sa = [r["sigma_axial_MPa"] for r in zone_results]
    ts = [r["tau_shear_MPa"] for r in zone_results]
    w = 0.25
    ax7.bar(x_idx - w, sb, w, label='Bending', color='#00d4ff', alpha=0.85)
    ax7.bar(x_idx,     sa, w, label='Axial',   color='#39ff14', alpha=0.85)
    ax7.bar(x_idx + w, ts, w, label='Shear',   color='#ff6b35', alpha=0.85)
    ax7.set_xticks(x_idx)
    ax7.set_xticklabels(zone_short, fontsize=8, color=TEXT_C, rotation=30)
    ax7.legend(fontsize=8, labelcolor=TEXT_C)
    apply_dark_style(ax7, ylabel='Stress (MPa)',
                     title='Stress Component Breakdown per Zone')

    # (1,1) Microstrain & Mechanostat classification
    ax8 = fig3.add_subplot(gs3[1, 1])
    ue_vals = [r["microstrain"] for r in zone_results]
    ue_colors = []
    for ue in ue_vals:
        if ue < 50:
            ue_colors.append('#5a7fa0')
        elif ue < 1500:
            ue_colors.append('#39ff14')
        elif ue < 3000:
            ue_colors.append('#ffcc00')
        else:
            ue_colors.append('#ff3366')
    bars8 = ax8.barh(zone_short, ue_vals, color=ue_colors, edgecolor='white',
                     linewidth=0.5, height=0.55)
    ax8.axvline(x=50, color='#5a7fa0', linestyle=':', linewidth=1.0,
                label='Disuse (<50 ue)')
    ax8.axvline(x=1500, color='#ffcc00', linestyle='--', linewidth=1.5,
                label='Osteoinductive (1500 ue)')
    ax8.axvline(x=3000, color='#ff3366', linestyle='--', linewidth=1.5,
                label='Overload (3000 ue)')
    for bar, val, mech in zip(bars8, ue_vals, [r["mechanostat"] for r in zone_results]):
        ax8.text(val + 30, bar.get_y() + bar.get_height()/2,
                 f'{val:.0f} ue ({mech})', va='center',
                 color=TEXT_C, fontsize=8)
    ax8.legend(fontsize=7, labelcolor=TEXT_C, loc='lower right')
    apply_dark_style(ax8, xlabel='Tissue Microstrain (ue)',
                     title='Frost Mechanostat Classification per Zone')

    fig3.suptitle('BIOPCL MANDIBLE -- 7-ZONE ANATOMICAL BIOMECHANICAL RESPONSE\n'
                  '(Bilateral Molar Bite, 600 N Total)',
                  color=ACCENT, fontsize=13, fontweight='bold', y=0.98)
    plt.savefig("fig3_biopcl_zone_analysis.png", dpi=180,
                bbox_inches='tight', facecolor=DARK_BG)
    plt.close(fig3)

    # ==================================================================
    # FIGURE 4: Deformation Comparison & Material Summary
    # ==================================================================
    print("  [4/4] Deformation Comparison & Material Properties...")
    fig4 = plt.figure(figsize=(18, 7), facecolor=DARK_BG)
    gs4 = gridspec.GridSpec(1, 2, figure=fig4, wspace=0.35,
                            left=0.07, right=0.95, top=0.88, bottom=0.15)

    # Left: Deformation comparison ABS vs BioPCL
    ax9 = fig4.add_subplot(gs4[0, 0])
    labels_def = ['Maximum\nDeformation', 'Average\nDeformation']
    abs_def = [fea_data["max_def_abs_mm"] * 1000,
               fea_data["avg_def_abs_mm"] * 1000]   # um
    pcl_def = [fea_data["max_def_biopcl_mm"] * 1000,
               fea_data["avg_def_biopcl_mm"] * 1000] # um
    x_def = np.arange(len(labels_def))
    w_def = 0.30
    b1 = ax9.bar(x_def - w_def/2, abs_def, w_def, label='ABS (E=2000 MPa)',
                 color='#5a7fa0', edgecolor='white')
    b2 = ax9.bar(x_def + w_def/2, pcl_def, w_def, label='BioPCL (E=350 MPa)',
                 color='#00d4ff', edgecolor='white')
    ax9.axhline(y=150.0, color='#ff3366', linestyle='--', linewidth=1.5,
                label='Osseointegration Limit (150 um)')
    for bar, val in zip(b1, abs_def):
        ax9.text(bar.get_x() + bar.get_width()/2, val + 2,
                 f'{val:.1f} um', ha='center', va='bottom',
                 color=DIM_C, fontsize=8)
    for bar, val in zip(b2, pcl_def):
        ax9.text(bar.get_x() + bar.get_width()/2, val + 2,
                 f'{val:.1f} um', ha='center', va='bottom',
                 color=TEXT_C, fontsize=9, fontweight='bold')
    ax9.set_xticks(x_def)
    ax9.set_xticklabels(labels_def, fontsize=10, color=TEXT_C)
    ax9.legend(fontsize=8, labelcolor=TEXT_C)
    apply_dark_style(ax9, ylabel='Deformation (um)',
                     title='Deformation Comparison: ABS vs BioPCL')

    # Right: Material properties summary table
    ax10 = fig4.add_subplot(gs4[0, 1])
    ax10.axis('off')
    table_data = [
        ['Property', 'BioPCL Value'],
        ['Density', f'{BIOPCL.rho_solid:.0f} kg/m^3'],
        ["Young's Modulus (E)", f'{BIOPCL.E_solid:.0f} MPa'],
        ["Poisson's Ratio (nu)", f'{BIOPCL.nu:.2f}'],
        ['Tensile Yield Strength', f'{BIOPCL.sigma_yt:.0f} MPa'],
        ['Tensile Ultimate Strength', f'{BIOPCL.sigma_ut:.0f} MPa'],
        ['Compressive Yield Strength', f'{BIOPCL.sigma_yc:.0f} MPa'],
        ['Compressive Ultimate Strength', f'{BIOPCL.sigma_uc:.0f} MPa'],
        ['Shear Modulus (G)', f'{BIOPCL.G_solid:.1f} MPa'],
        ['Construct Volume', f'{fea_data["volume_cm3"]:.2f} cm^3'],
        ['Solid Construct Mass', f'{fea_data["mass_solid_g"]:.2f} g'],
    ]
    table = ax10.table(cellText=table_data, loc='center',
                       cellLoc='left', edges='horizontal')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.8)
    for key, cell in table.get_celld().items():
        cell.set_edgecolor('#1a3a5c')
        cell.set_facecolor(PANEL_BG)
        cell.set_text_props(color=TEXT_C)
        if key[0] == 0:  # header row
            cell.set_text_props(color=ACCENT, fontweight='bold')
            cell.set_facecolor('#0d1a2a')
    ax10.set_title('BioPCL Material Property Card', color=ACCENT,
                   fontsize=12, fontweight='bold', pad=15)

    fig4.suptitle('DEFORMATION ANALYSIS & MATERIAL SPECIFICATION',
                  color=ACCENT, fontsize=13, fontweight='bold', y=0.97)
    plt.savefig("fig4_biopcl_deformation_material.png", dpi=180,
                bbox_inches='tight', facecolor=DARK_BG)
    plt.close(fig4)


# ==============================================================================
# SECTION 5 -- MAIN EXECUTION & CONSOLE REPORT
# ==============================================================================

def main():
    # Handle Windows console encoding
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        except Exception:
            pass

    print("\n" + "=" * 85)
    print("    BIOPCL BIOPRINTED MANDIBLE -- COMPUTATIONAL BIOMECHANICAL ANALYSIS")
    print("=" * 85 + "\n")

    # ---- Load & transform FEA data ----
    fea_data = load_and_transform_fea_data()

    print("--- 1. BIOPCL MATERIAL & CONSTRUCT SPECIFICATIONS ---------------------")
    print(f"  Material        : Bioprinted BioPCL (Polycaprolactone)")
    print(f"  Density         : {BIOPCL.rho_solid:.0f} kg/m^3 ({BIOPCL.rho_solid/1000:.2f} g/cm^3)")
    print(f"  Young's Modulus : {BIOPCL.E_solid:.0f} MPa (0.35 GPa)")
    print(f"  Poisson's Ratio : {BIOPCL.nu:.2f}")
    print(f"  Yield Strength  : {BIOPCL.sigma_yt:.0f} MPa (tensile & compressive)")
    print(f"  Ultimate Str.   : {BIOPCL.sigma_ut:.0f} MPa (tensile & compressive)")
    print(f"  Shear Modulus   : {BIOPCL.G_solid:.1f} MPa")
    print(f"  Construct Vol.  : {fea_data['volume_cm3']:.2f} cm^3")
    print(f"  Solid Mass      : {fea_data['mass_solid_g']:.2f} g\n")

    print(f"--- 2. FEA NODAL STRESS ANALYSIS ({fea_data['counts']['total']:,} Nodes) --------")
    s = fea_data['stats']
    print(f"  Max Von Mises   : {s['max']:.2f} MPa")
    print(f"  99th Percentile : {s['p99']:.2f} MPa")
    print(f"  95th Percentile : {s['p95']:.2f} MPa")
    print(f"  Mean Stress     : {s['mean']:.2f} MPa")
    print(f"  Median Stress   : {s['median']:.2f} MPa")
    print(f"  Std. Deviation  : {s['std']:.2f} MPa")
    print(f"  Max Deform. ABS : {fea_data['max_def_abs_mm']:.4f} mm ({fea_data['max_def_abs_mm']*1000:.1f} um)")
    print(f"  Max Deform. PCL : {fea_data['max_def_biopcl_mm']:.4f} mm ({fea_data['max_def_biopcl_mm']*1000:.1f} um)")
    print(f"  Avg Deform. PCL : {fea_data['avg_def_biopcl_mm']:.4f} mm ({fea_data['avg_def_biopcl_mm']*1000:.1f} um)")
    print(f"  Scaling Factor  : E_ABS/E_PCL = {fea_data['modulus_ratio']:.4f}\n")

    tot = fea_data['counts']['total']
    c = fea_data['counts']
    print("--- 3. CLINICAL SAFETY CLASSIFICATION -----------------------------------")
    print(f"  [SAFE]   <= {BIOPCL.sigma_yt} MPa : {c['safe']:>6,} nodes ({c['safe']/tot*100:.2f}%)")
    print(f"  [YIELD]  {BIOPCL.sigma_yt}-{BIOPCL.sigma_ut} MPa : {c['yield']:>6,} nodes ({c['yield']/tot*100:.2f}%)")
    print(f"  [EXCEED] > {BIOPCL.sigma_ut} MPa  : {c['exceed']:>6,} nodes ({c['exceed']/tot*100:.2f}%)\n")

    # ---- Gibson-Ashby scaffold optimization table ----
    print("--- 4. GIBSON-ASHBY SCAFFOLD OPTIMIZATION TABLE -------------------------")
    print(f"  {'P(%)':>5}  {'RelRho':>7}  {'E_eff(MPa)':>10}  {'sy(MPa)':>8}  {'Mass(g)':>8}  {'Median ue':>10}  {'Mechanostat':>20}")
    print("  " + "-" * 82)
    for p_pct in [0, 10, 15, 20, 25, 30, 40, 50]:
        p = p_pct / 100.0
        props = BIOPCL.effective_properties(p)
        mass_g = (fea_data['volume_cm3'] * props['density_kg_m3']) / 1000.0
        # Use MEDIAN stress & composite modulus for mechanostat
        ue = BIOPCL.tissue_microstrain(fea_data['stats']['median'], p)
        status = BIOPCL.classify_mechanostat(ue)
        print(f"  {p_pct:>5.0f}  {props['rel_density']:>7.2f}  {props['E_eff_MPa']:>10.1f}  "
              f"{props['sigma_yield_MPa']:>8.2f}  {mass_g:>8.2f}  {ue:>10.0f}  {status:>20}")
    print()

    # ---- 7-Zone anatomical analysis ----
    print("--- 5. 7-ZONE ANATOMICAL RESPONSE (Bilateral Molar Bite, 600 N) ---------")
    print(f"  {'Zone':<20} {'sigma_VM':>8} {'sigma_b':>8} {'tau':>8} {'SF_y':>6} {'SF_u':>6} {'ue':>8} {'Mechanostat':>18}")
    print("  " + "-" * 98)

    zone_results = []
    for z in ZONES_BIOPCL:
        loads = ZONE_LOADS[z.name]
        res = z.compute_response(loads["M"], loads["V"], loads["F_ax"])
        zone_results.append(res)
        print(f"  {res['zone']:<20} {res['sigma_vm_MPa']:>8.2f} {res['sigma_bend_MPa']:>8.2f} "
              f"{res['tau_shear_MPa']:>8.2f} {res['SF_yield']:>6.2f} {res['SF_ult']:>6.2f} "
              f"{res['microstrain']:>8.0f} {res['mechanostat']:>18}")

    max_vm_zone = max(zone_results, key=lambda r: r['sigma_vm_MPa'])
    print(f"\n  Critical zone: {max_vm_zone['zone']} "
          f"(sigma_VM = {max_vm_zone['sigma_vm_MPa']:.2f} MPa, "
          f"SF_yield = {max_vm_zone['SF_yield']:.2f})\n")

    # ---- Generate figures ----
    print("--- 6. GENERATING PUBLICATION FIGURES ------------------------------------")
    generate_all_figures(fea_data, zone_results)

    print("\n" + "=" * 85)
    print("  [DONE] All analyses complete. 4 figures saved to working directory.")
    print("=" * 85 + "\n")


if __name__ == "__main__":
    main()
