# Enzyme Clusters and Reaction Simulation in CoQ Biosynthesis

This repository contains the computational framework and analysis tools for the research described in the paper "Complete Enzyme Clustering Enhances Substrate Channeling in Coenzyme Q Biosynthesis" (https://www.biorxiv.org/content/10.1101/2025.05.24.655883v1) 

## 📖 Overview

This project uses coarse-grained molecular dynamics simulations to investigate how transient assemblies of sequential metabolic enzymes enhances reaction flux in the coenzyme Q (CoQ) biosynthetic pathway. The simulations model the COQ metabolon using experimentally measured protein-protein interaction strengths and demonstrate that complete enzyme clustering enables substrate channeling, dramatically enhancing CoQ production efficiency.

## 📁 Repository Structure
```
CoQ_metabolon/
├── example_simulation/                     # Example LAMMPS simulation setup
│   ├── metabolon_final.in                  # Main LAMMPS input file
│   ├── system.data                         # LAMMPS data file with initial particle coordinates
│   └── generate_system_data.py             # Script to generate system.data
│
├── analysis/                               # Analysis tools and figure generation
│   ├── utils.py                            # Core analysis functions
│   ├── fix_indexing.ipynb                  # Converts VMD 1-based indexing to 0-based indexing
│   ├── paper_figures.ipynb                 # Jupyter notebook for paper figure generation
│   ├── cluster_analysis_code/              # VMD/TCL scripts for enzyme cluster detection
│   └── rdf_code/                           # VMD/TCL scripts for Radial distribution function analysis
│
├── environments.yml                        # Conda environment specification file
└── README.md                               # Project documentation (this file)

```

## ⚙️ Simulation Model

### Coarse-Grained Representation
- **Enzymes (COQ3-7, COQ9)**: Hard spheres (radius 10 Å) with adhesive interaction sites (radius 3.5 Å) representing active sites
- **Substrates/Products**: Hard spheres (radius 1 Å) representing metabolic intermediates
- **Crowders**: Hard spheres (radius 10 Å) without active sites

### Interaction Parameters
- **Protein-protein interactions (εenz-enz)**: Based on experimentally measured Kd constants
- **Active site-substrate interactions (εpatch-lig)**: 4kT for cognate pairs, 2kT for non-cognate
- **Enzyme-substrate interactions (εenz-lig)**: 2kT (non-specific binding)

### Key Simulation Parameters
- **Temperature**: 310 K (physiological)
- **Integration timestep**: 25 fs
- **Equilibration**: 625 ns
- **Production run**: 7.5-15 μs depending on analysis
- **System size**: 800 Å cubic box

## 📊 Simulation and Analysis Pipeline

- Create the Conda environment using the provided environment.yml
- Execute LAMMPS simulation with metabolon.in and system.data input files. 
- After the simulation finishes, open VMD and run "source run_all.tcl" for cluster analysis and rdf analysis
- VMD uses 1-based indexing. To convert to 0-based indexing (for Python analysis), run fix_indexing.ipynb
- Use appropriate functions in analysis/utils.py for analysis

## Software versions
- LAMMPS version 12Dec2018 (stable) was used for all simulations, except for substrate replenishing simulations (LAMMPS version 3Mar2020 linked to python 3.12)
- VMD version 1.9.4a57

## Citation

If you use this code in your research, please cite:
```bibtex
@article {Wang2025.05.24.655883,
	author = {Wang, Dianzhuo and Gottinger, Andrea and Jeong, Jio and Nicoll, Callum R. and Liu, Junlang and Kadav, Tereza and Cecchini, Domiziana and Malatesta, Marco and Heck, Albert J.R. and Mattevi, Andrea and Shakhnovich, Eugene},
	title = {Complete Enzyme Clustering Enhances Coenzyme Q Biosynthesis via Substrate Channeling},
	elocation-id = {2025.05.24.655883},
	year = {2025},
	doi = {10.1101/2025.05.24.655883},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/05/28/2025.05.24.655883},
	eprint = {https://www.biorxiv.org/content/early/2025/05/28/2025.05.24.655883.full.pdf},
	journal = {bioRxiv}
    }
```
## Contact

For questions about the code, please contact:
- **Dianzhuo(John) Wang**: johnwang@g.harvard.edu
- **Jio Jeong**: jiojeong@g.harvard.edu 

## License

This project is licensed under the MIT License 