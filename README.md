# COQ Metabolon: Complete Enzyme Clustering Enhances Substrate Channeling in Coenzyme Q Biosynthesis

This repository contains the computational framework and analysis tools for the research described in the paper "Complete Enzyme Clustering Enhances Substrate Channeling in Coenzyme Q Biosynthesis" by Wang et al.

## Overview

This project uses coarse-grained molecular dynamics simulations to investigate how metabolon formation (transient assemblies of sequential metabolic enzymes) enhances metabolic flux in the coenzyme Q (CoQ) biosynthetic pathway. The simulations model the COQ metabolon using experimentally measured protein-protein interaction strengths and demonstrate that complete enzyme clustering enables substrate channeling, dramatically enhancing CoQ production efficiency.

## Key Findings

- **Enzyme Clustering Enhances Flux**: Complete metabolon clustering dramatically improves metabolic flux through substrate channeling
- **Phase Transition Behavior**: The COQ metabolon operates at the critical region of a phase-like transition where small changes in protein-protein interaction strength can affect both enzyme clustering and metabolic output
- **Completeness Over Fine Structure**: Metabolon completeness (containing all required enzymes) rather than precise spatial arrangement is essential for efficient substrate channeling
- **Network Optimization**: The experimentally observed COQ interaction network is evolutionarily selected to promote formation of complete enzyme clusters

## 📁 Repository Structure
CoQ_metabolon/
├── example_simulation/          # Example LAMMPS simulation setup
│   ├── metabolon_final.in       # Main LAMMPS input file
│   ├── system.data              # System data file with particle coordinates
│   └── generate_system_data/   # Scripts to generate initial configurations
├── analysis/                    # Analysis tools and figure generation
│   ├── utils.py                 # Core analysis functions (Python)
│   ├── paper_figures.ipynb     # Jupyter notebook to reproduce paper figures
│   ├── cluster_analysis_code/  # VMD/TCL scripts for cluster detection
│   └── rdf_code/                # Radial distribution function analysis tools
└── README.md                    # This file


## Simulation Model

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

## Analysis Tools

### Core Functions (`analysis/utils.py`)
- `readlogfile()`: Parse LAMMPS output files
- `average_runs()`: Average results across multiple simulation runs
- `readcluster()`: Analyze cluster formation and composition
- `read_complete_clusters()`: Track complete metabolon formation
- `average_counts_atom_type()`: Monitor substrate/product concentrations

### Cluster Analysis (`cluster_analysis_code/`)
- Identifies enzyme clusters using distance-based criteria
- Tracks cluster size, composition, and completeness over time
- Analyzes substrate channeling efficiency

## Data Format

### Input Files
- `system.data`: LAMMPS data file with initial particle positions and topology
- `metabolon_final.in`: LAMMPS input script with simulation parameters

### Output Files
- `output.log`: LAMMPS log file with thermodynamic data and reaction counts
- `cluster_results_skip_first_frame/`: Cluster analysis results
- Trajectory files for visualization and further analysis

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

For questions about the code or methodology, please contact:
- **Dianzhuo(John) Wang: johnwang@g.harvard.edu
- **Jio Jeong**: jiojeong@g.harvard.edu 

## License

This project is licensed under the MIT License - see the LICENSE file for details.