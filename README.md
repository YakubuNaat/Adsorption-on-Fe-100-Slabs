# CO₂ Adsorption on Transition-Metal-Doped Fe(100) Slabs

## Overview
This repository contains Python scripts for constructing transition-metal-doped iron slabs,
adsorbing CO₂ molecules, and performing DFT calculations to extract features, energies and generate
CSV datasets for machine learning applications.

The workflow is designed for computational catalysis and surface science studies.

## Workflow
1. Generate Fe(100) slabs doped with single transition metal atoms
2. Create isolated single atoms for reference calculations
3. Adsorb CO₂ on doped slabs
4. Perform DFT calculations using GPAW
5. Extract energies and export results to CSV

## Files Description
- `singleatom.py` – Generates isolated transition metal atoms
- `cifs.py` – Builds Fe(100) slabs doped with transition metals
- `adsorption.py` – Adsorbs CO₂ on the doped slabs
- `dftcomplete.py` – Runs DFT calculations and outputs CSV files

## Tools & Libraries
- ASE
- GPAW
- NumPy
- Python

## Applications
- CO₂ adsorption energy prediction
- Catalysis screening
- Machine learning dataset generation

## Author
Abdulai Yakubu  
Undergraduate Chemist | Computational Chemistry | Quantum Computing
