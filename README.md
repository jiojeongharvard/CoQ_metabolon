# Enzyme Clusters and Reaction Simulation in CoQ Biosynthesis

This repository contains the computational framework and analysis tools for the research described in the paper "Complete Enzyme Clustering Enhances Substrate Channeling in Coenzyme Q Biosynthesis" (https://www.biorxiv.org/content/10.1101/2025.05.24.655883v1) 

## 📖 Overview

This project uses coarse-grained molecular dynamics simulations to investigate how transient assemblies of sequential metabolic enzymes enhances reaction flux in the coenzyme Q (CoQ) biosynthetic pathway. The simulations model the COQ metabolon using experimentally measured protein-protein interaction strengths and demonstrate that complete enzyme clustering enables substrate channeling, dramatically enhancing CoQ production efficiency.

## 📁 Repository Structure
```
CoQ_metabolon/
├── example_simulation/          # Example LAMMPS simulation setup
│   ├── metabolon_final.in       # Main LAMMPS input file
│   ├── system.data              # LAMMPS data file with initial particle coordinates
│   └── generate_system_data/   # Script to generate system.data
├── analysis/                    # Analysis tools and figure generation
│   ├── utils.py                 # Core analysis functions 
│   ├── paper_figures.ipynb     # Jupyter notebook used to produce paper figures
│   ├── cluster_analysis_code/  # VMD/TCL scripts for enzyme cluster detection  using distance-based criteria
│   └── rdf_code/                # Radial distribution function analysis tools
├── environments.yml			# Needed for conda environment generation
└── README.md                    # This file
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

- Create conda environment using environment.yml
- Run LAMMPS simulation with metabolon.in and system.data. 
- Once simulation finishes, run "source run_all.tcl" in VMD for cluster analysis and rdf analysis
- Run fix.indexing.ipynb on result to update VMD's 1-indexing to 0-indexing before analysis
- Use appropriate functions in utils.py for analysis

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