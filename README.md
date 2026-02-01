V.A.L.L.Y. Project: Universal Computational Framework (v1.7)
Redefining Agility in Drug Discovery and Materials Science

🔬 Overview
V.A.L.L.Y. (Vibrational Analysis for Ligand Likelihood Yielding) is an ultra-fast computational engine designed to identify critical dynamic sites in complex structures. 
By integrating Anisotropic Network Models (ANM) with a proprietary Triple-Factor Validation logic, the framework bridges the gap between theoretical physics and experimental crystallography.

🚀 Key Features (v1.7 Universal Edition)Triple-Factor Validation System (R2):Physics Engine: 
Calculates intrinsic dynamics using the Hessian Matrix.Informed Heuristics: Automatically discriminates between catalytic and allosteric sites using geometric exclusion.Experimental Ground Truth: Real-time correlation ($r$) with crystallographic B-factors.Universal Core: Support for both Viral Proteases (SARS-CoV-2, Dengue) and Mineral/Crystalline Structures.Extreme Performance: Analysis completed in seconds on standard consumer hardware (Core i3, 8GB RAM).

# Clone the repository
git clone https://github.com/Lionnava/vally-scan.git

# Install dependencies
pip install prody numpy pandas matplotlib scipy

Running the Universal Engine:

from vally_scan_v1_7_universal import vally_universal_engine

# For Viral Research (e.g., SARS-CoV-2)
vally_universal_engine('6LU7.pdb', active_site_residues=[41, 145], mode='viral')

# For Mineral Research
vally_universal_engine('quartz.pdb', mode='mineral')

📊 Validated ResultsThe framework has been cross-validated against high-priority viral targets:| Target | PDB ID | Pearson Correlation ($r$) | Key Discovery || :--- | :--- | :--- | :--- || SARS-CoV-2 Mpro | 6LU7 | 0.67 | Identified Allosteric Cluster (277-279) || Dengue Protease | 2FOM | 0.43 | High-flexibility region (43-45) |

🎓 About the Author
Ing. Lionell E. Nava Ramos Researcher at UPTMA - Centro de Investigación en Informática (CII). Accepted Presenter at the APS Global Physics Summit 2026, Denver, CO.
