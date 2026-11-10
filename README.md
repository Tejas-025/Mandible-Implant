[patient_specific_mandible_implant_readme.md](https://github.com/user-attachments/files/34647701/patient_specific_mandible_implant_readme.md)
# Patient-Specific Mandible Implant Design & Biomechanical Digital Twin

[![Platform: Python 3.10+](https://img.shields.io/badge/Platform-Python%203.10%2B-blue.svg?logo=python&logoColor=white)](https://www.python.org/)
[![CAD: Geomagic Design X](https://img.shields.io/badge/CAD-Geomagic%20Design%20X-orange.svg)](#1-anatomical-segmentation--cad-reconstruction-geomagic-design-x)
[![Generative: Autodesk Fusion](https://img.shields.io/badge/Generative%20Design-Autodesk%20Fusion-red.svg?logo=autodesk&logoColor=white)](#2-generative-design--lattice-optimization-autodesk-fusion)
[![FEA: ANSYS Mechanical](https://img.shields.io/badge/Simulation-ANSYS%20Mechanical-yellow.svg?logo=ansys&logoColor=black)](#3-multi-zone-biomechanical-fea--fatigue-qualification-ansys)
[![AM: Ti--6Al--4V ELI (PBF-LB)](https://img.shields.io/badge/AM-Ti--6Al--4V%20ELI%20%7C%20PBF--LB-lightgrey.svg)](#4-additive-manufacturing-specifications)
[![Compliance: ISO 5832-3 / ASTM F3001](https://img.shields.io/badge/Standards-ISO%205832--3%20%7C%20ASTM%20F3001-green.svg)](#standards-compliance)

An end-to-end biomedical design, generative optimization, and computational qualification framework for patient-specific mandibular reconstruction. This project bridges raw clinical CBCT volumetric imaging, NURBS reverse engineering in **Geomagic Design X**, cellular lattice generative design in **Autodesk Fusion**, physiological multi-zone FEA in **ANSYS Mechanical**, and interactive 3D stress telemetry in a **Python-based Web Digital Twin**.

---

## Table of Contents

- [Clinical Motivation & Problem Statement](#clinical-motivation--problem-statement)
- [Workflow Architecture](#workflow-architecture)
- [Technical Methodology](#technical-methodology)
  - [1. Anatomical Segmentation & CAD Reconstruction (Geomagic Design X)](#1-anatomical-segmentation--cad-reconstruction-geomagic-design-x)
  - [2. Generative Design & Lattice Optimization (Autodesk Fusion)](#2-generative-design--lattice-optimization-autodesk-fusion)
  - [3. Multi-Zone Biomechanical FEA & Fatigue Qualification (ANSYS)](#3-multi-zone-biomechanical-fea--fatigue-qualification-ansys)
  - [4. Interactive Web-Based Digital Twin (Python)](#4-interactive-web-based-digital-twin-python)
- [Quantitative Performance Comparison](#quantitative-performance-comparison)
- [Additive Manufacturing & Post-Processing](#additive-manufacturing--post-processing)
- [Repository Structure](#repository-structure)
- [Quickstart: Running the Python Digital Twin](#quickstart-running-the-python-digital-twin)
- [Standards Compliance](#standards-compliance)
- [Author & Citation](#author--citation)

---

## Clinical Motivation & Problem Statement

Extensive segmental mandibular defects caused by oncological resections (e.g., ameloblastoma, squamous cell carcinoma) or high-impact ballistic/maxillofacial trauma present steep reconstructive challenges:

* **The Problem with Standard Plates:** Conventional 2.4 mm titanium reconstruction plates require manual intraoperative bending. This introduces residual stress risers, leads to high fatigue failure rates ($>25\%$), and causes facial asymmetry. Crucially, bulk titanium ($E \approx 110\text{ GPa}$) is an order of magnitude stiffer than native mandibular cortical bone ($E \approx 10\text{--}18\text{ GPa}$). This mismatch induces **stress shielding**, preventing physiological mechanical strain transfer and triggering progressive peri-implant bone resorption and screw loosening.
* **The Engineering Objective:** Reconstruct anatomical contours directly from patient computed tomography, synthesize porous triply periodic minimal surface (TPMS) cellular cores to match native bone stiffness, and validate structural resilience across all physiological chewing modes prior to surgical implantation.

---

## Workflow Architecture

```
  [ Patient CBCT Scan (DICOM) ]
                │
                ▼
  [ Geomagic Design X ] 
    ├── Hounsfield Unit Thresholding (+400 to +1800 HU)
    ├── Osteotomy Resection Margins (+5 mm Oncological Safety Margin)
    └── Watertight $G^2$ NURBS Solid Reconstruction (CAD Deviation < 0.045 mm)
                │
                ▼
  [ Autodesk Fusion ] 
    ├── Spatial Boundary Separation (Preserve vs. Obstacle Envelopes)
    ├── Inferior Alveolar Nerve Clearance Domain (3 mm Exclusion Cylinder)
    └── Generative TPMS Cellular Lattice Synthesis (42% Mass Reduction, E_eff = 16.2 GPa)
                │
                ▼
  [ ANSYS Mechanical ] 
    ├── 8 Anatomical Structural Zones (Condyle, Ramus, Body, Symphysis, etc.)
    ├── Active EMG-Calibrated Muscle Vectors (Masseter, Temporalis, Pterygoids)
    └── 6 Non-Linear Masticatory Load Cases (FoS > 2.4, High-Cycle Fatigue > 10^7 Cycles)
                │
                ▼
  [ Python Web Digital Twin ] 
    ├── Interactive WebGL/Three.js 3D Viewport in Browser
    ├── Real-Time Nodal von Mises Stress & Strain Mapping
    └── Dynamic Clenching Force Controls & CSV Telemetry Logger
```

---

## Technical Methodology

### 1. Anatomical Segmentation & CAD Reconstruction (Geomagic Design X)

1. **DICOM Ingestion & Density Thresholding:** High-resolution CBCT volumetric datasets ($0.25\text{ mm}$ voxel slice spacing) were processed using calibrated bone density windows ($+400\text{ to }+1800\text{ HU}$) to segment intact mandibular bone boundaries.
2. **Virtual Surgical Planning (VSP) & Resection:** Pathological tissue regions were isolated with a $5\text{ mm}$ oncological margin. Osteotomy cutting planes were configured to ensure perpendicular surgical tool access.
3. **Parametric Reverse Engineering:**
   * Raw segmented tessellations contain non-manifold edges, surface irregularities, and noise. **Geomagic Design X** was utilized to heal the boundary mesh, generate auto-surface cross-sectional networks, and bridge the defect with parametric, curvature-continuous ($G^2$) watertight NURBS solids.
   * Monocortical and bicortical screw bosses were integrated with fixed $2.0\text{ mm}$ locking trajectories oriented away from tooth roots and the inferior alveolar neurovascular canal.
   * 3D deviation inspection verified that the reconstructed CAD solid matched patient bone contours within a mean deviation tolerance of $< 0.045\text{ mm}$.

---

### 2. Generative Design & Lattice Optimization (Autodesk Fusion)

1. **Boundary & Keep-Out Zoning:**
   * **Preserve Bodies:** Osteotomy interface contact flanges, screw fastening collars, and the articulation zone of the condylar head.
   * **Obstacle Geometry:** Surrounding intraoral soft tissue volumes, tooth roots, and a $3\text{ mm}$ protective clearance corridor around the inferior alveolar nerve.
2. **Stiffness Tuning via Cellular Lattices:**
   * Using generative design tools in **Autodesk Fusion**, porous cellular lattice cores (Gyroid and Diamond TPMS structures) were synthesized in the central structural web.
   * The effective elastic modulus was lowered from bulk titanium ($110\text{ GPa}$) down to $E_{\text{eff}} \approx 14.5\text{--}18.2\text{ GPa}$, matching adjacent cortical bone compliance and mitigating stress shielding.
3. **Mass Reduction:**
   * Achieved a **42% reduction in overall implant mass** compared to solid anatomical alternatives, reducing foreign body sensation and optimizing the bed for soft-tissue flap revascularization.

---

### 3. Multi-Zone Biomechanical FEA & Fatigue Qualification (ANSYS)

Finite Element Analysis was executed in **ANSYS Mechanical** utilizing higher-order 10-node tetrahedral solid elements (SOLID187, global element size $0.8\text{ mm}$ with local convergence refinement to $0.2\text{ mm}$ at screw-plate contacts).

#### Anatomical Functional Zones
The mandible model was segmented into eight biomechanical regions:
1. Condylar Process / Head
2. Coronoid Process
3. Mandibular Ramus
4. Gonial Angle
5. Posterior Mandibular Body
6. Anterior Mandibular Body
7. Parasymphysis
8. Symphysis

#### Muscle Loading Vectors (EMG-Calibrated)
Physiological muscle traction loads were distributed across physiological insertion footprints:

| Muscle Group | Insertion Site | Directional Force Vector $\mathbf{F} = [F_x, F_y, F_z]\text{ (N)}$ | Net Load |
| :--- | :--- | :--- | :--- |
| **Superficial Masseter** | Inferior/lateral gonial angle | $[-15.2, -42.8, 185.0]$ | $190.5\text{ N}$ |
| **Deep Masseter** | Lateral ramus surface | $[-8.4, -18.2, 78.5]$ | $81.0\text{ N}$ |
| **Medial Pterygoid** | Medial surface of gonial angle | $[24.5, -35.6, 155.0]$ | $160.9\text{ N}$ |
| **Anterior Temporalis** | Coronoid process apex | $[-12.0, 38.5, 210.0]$ | $213.8\text{ N}$ |
| **Posterior Temporalis**| Retromolar triangle & ramus ridge | $[-5.2, 85.0, 62.0]$ | $105.3\text{ N}$ |
| **Lateral Pterygoid** | Condylar neck & fovea | $[32.0, -58.0, 15.0]$ | $67.9\text{ N}$ |

#### Evaluated Load Cases
1. **Maximum Intercuspation (ICP):** Symmetric bilateral molar bite ($F_{\text{bite}} = 650\text{ N}$).
2. **Incisal Clench (INC):** Pure vertical anterior shear bite ($F_{\text{bite}} = 200\text{ N}$).
3. **Right Unilateral Molar Clench (RMC):** Asymmetric ipsilateral bite ($F_{\text{bite}} = 450\text{ N}$).
4. **Left Unilateral Molar Clench (LMC):** Asymmetric contralateral bite ($F_{\text{bite}} = 450\text{ N}$).
5. **Right Lateral Excursion:** Canine-guided unilateral grinding ($F_{\text{bite}} = 220\text{ N}$).
6. **Retrusive Clench:** Posterior retraction loading.

#### Biomechanical Findings
* **Peak von Mises Stress:** $278.4\text{ MPa}$ on the titanium framework during peak unilateral clenching (far below the $860\text{ MPa}$ yield point of annealed Ti-6Al-4V ELI).
* **Bone Contact Stress Transfer:** Maintained between $2.5\text{ and }28.0\text{ MPa}$ ($1000\text{--}3000\ \mu\varepsilon$), stimulating osteogenesis according to Frost's mechanostat model.
* **Safety Factor:** $\text{FoS} > 2.4$ across all static and dynamic bite cycles.
* **Fatigue Life:** Goodman fatigue evaluation verified infinite cycle endurance exceeding $1 \times 10^7$ chew cycles ($> 10\text{ years}$ clinical equivalence).

---

### 4. Interactive Web-Based Digital Twin (Python)

To allow clinical teams and engineers to inspect telemetry without requiring specialized CAD/FEA licenses, an interactive web-based digital twin was constructed in **Python** using Dash, Plotly, and Trimesh:

* **3D Interactive Mesh Viewer:** Full WebGL rotation, pan, zoom, and section cut of patient mandibular bone, tumor volume, resection planes, and implant framework.
* **Stress Telemetry Contours:** Real-time spatial mapping of nodal von Mises stress and deformation fields with adjustable scalar ranges.
* **Dynamic Chewing Condition Selectors:** Instantaneous toggling between all 6 bite scenarios with real-time reaction load displays.
* **Telemetry Data Logger:** Live export of nodal stress profiles, factor of safety distributions, and mass properties directly to CSV.

---

## Quantitative Performance Comparison

| Metric / Parameter | Standard Reconstruction Plate (2.4 mm) | Generative Patient-Specific Implant (This Work) | Clinical Impact |
| :--- | :--- | :--- | :--- |
| **Total Implant Mass** | $32.4\text{ g}$ | **$18.8\text{ g}$** | **42% mass reduction**; lowers tissue exposure risks |
| **Effective Elastic Modulus** | $110\text{ GPa}$ (Solid Ti) | **$16.2\text{ GPa}$** (TPMS Lattice) | **Matches cortical bone**; prevents stress shielding |
| **Anatomical Conformity** | Poor (Intraoperative bending) | **Precise ($< 0.05\text{ mm}$ deviation)** | Restores facial symmetry; reduces operating time |
| **Min. Factor of Safety (FoS)**| $1.4\text{--}1.8$ (Stress peaks at screw holes)| **$> 2.4$** | Eliminates fatigue-induced plate fracture |
| **Peak Bone Interface Stress** | $> 85\text{ MPa}$ (Micro-fractures / osteolysis) | **$24.2\text{ MPa}$** (Within remodeling zone) | Promotes stable, long-term osseointegration |

---

## Additive Manufacturing Specifications

* **Modality:** Laser Powder Bed Fusion (PBF-LB / DMLS).
* **Alloy:** Medical-Grade Titanium Alloy Ti-6Al-4V ELI (Grade 23, ASTM F136 / ASTM F3001).
* **Layer Thickness:** $30\ \mu\text{m}$.
* **Laser Spot Diameter:** $75\ \mu\text{m}$.
* **Build Atmosphere:** Ultra-pure Argon ($O_2 < 100\text{ ppm}$).
* **Post-Processing Pipeline:**
  1. Automated depowdering with ultrasonic agitation to clear internal cellular channels.
  2. Vacuum stress-relief annealing at $730^\circ\text{C}$ for 2 hours.
  3. Wire-EDM detachment from the build substrate.
  4. Hot Isostatic Pressing (HIP) at $920^\circ\text{C}$, $100\text{ MPa}$ for 2 hours to close micro-porosities.
  5. Chemical etching and electropolishing to attain surface roughness ($Ra \approx 1.2\ \mu\text{m}$) tailored for cell adhesion.

---

## Repository Structure

```
├── cad/
│   ├── raw_cbct_segmentation/      # Cleaned surface meshes (.stl) from CBCT
│   ├── geomagic_nurbs_solids/      # Reverse-engineered watertight solids (.step, .iges)
│   └── fusion_generative_models/   # Generative models and cellular TPMS lattices (.f3d, .step)
├── fea/
│   ├── ansys_workbench/            # ANSYS Workbench simulation databases (.wbpj)
│   ├── boundary_conditions/        # Physiological muscle traction force tables (.csv)
│   └── validation_reports/         # Convergence studies, stress-strain, and fatigue logs
├── digital_twin_web/
│   ├── app.py                      # Main entry point for the Python web digital twin
│   ├── kinematics_engine.py        # Coordinate transformations and stress field interpolation
│   ├── static/                     # WebGL shaders, styling, and static CAD geometries
│   ├── templates/                  # Dashboard layout templates
│   └── requirements.txt            # Python environment dependencies
├── docs/
│   ├── surgical_planning_guide.pdf # Surgical guides, resection planes, and screw vectors
│   └── biomechanical_report.pdf   # Complete clinical verification and FEA documentation
└── README.md
```

---

## Quickstart: Running the Python Digital Twin

### Prerequisites
* Python 3.10 or higher
* WebGL-capable web browser (Google Chrome, Mozilla Firefox, Microsoft Edge)

### 1. Clone the Repository
```bash
git clone https://github.com/username/patient-specific-mandible-implant.git
cd patient-specific-mandible-implant/digital_twin_web
```

### 2. Set Up Virtual Environment
```bash
# On Linux / macOS
python3 -m venv venv
source venv/bin/activate

# On Windows
python -m venv venv
venv\Scripts\activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

### 4. Launch the Application
```bash
python app.py
```

Open `http://127.0.0.1:8050/` in your browser to view the 3D mandible implant model, cycle through chewing load scenarios, and examine real-time stress telemetry.

---

## Standards Compliance

* **ASTM F3001:** Standard Specification for Additive Manufacturing Titanium-6 Aluminum-4 Vanadium ELI with Powder Bed Fusion.
* **ASTM F136:** Standard Specification for Wrought Titanium-6 Aluminum-4 Vanadium ELI for Surgical Implant Applications.
* **ISO 5832-3:** Implants for surgery — Metallic materials — Part 3: Wrought titanium 6-aluminium 4-vanadium alloy.
* **ASTM F382 / ISO 14602:** Standard Specification and Test Method for Metallic Bone Plates.
* **ISO 10993-1:** Biological evaluation of medical devices — Part 1: Evaluation and testing within a risk management process.

---

## Author & Citation

**Sarvesh Sunil Deshpande**  
B.E. Mechanical Engineering, KLS Gogte Institute of Technology  
B.S. Electronic Systems, Indian Institute of Technology Madras  
*Portfolio:* [sarveshdeshpande.dev](https://sarveshdeshpande.dev) · *LinkedIn:* [linkedin.com/in/sarvesh-deshpande](https://linkedin.com/in/sarvesh-deshpande)

If you reference this workflow, design framework, or simulation data, please cite:

```bibtex
@article{deshpande2026mandible,
  title={Patient-Specific Mandible Implant Design Combining CBCT Modeling, Generative Latticing, Non-Linear Biomechanical FEA, and a Web-Based Digital Twin},
  author={Deshpande, Sarvesh Sunil},
  journal={Biomedical Engineering & Additive Manufacturing Repository},
  year={2026}
}
```
