import os
import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from gpaw import GPAW, PW, FermiDirac

# Define directories
cif_directory = '/home/yabdulai/lustre/naat/structures'
output_directory = '/home/yabdulai/lustre/naat/dftnaat/direct/new'
optimized_directory = '/home/yabdulai/lustre/naat/dftnaat/direct/optimized_structures'
gpw_directory = '/home/yabdulai/lustre/naat/dftnaat/direct/gpw_files'

# Ensure directories exist
for dir_path in [output_directory, optimized_directory, gpw_directory]:
    os.makedirs(dir_path, exist_ok=True)

# Initialize result storage
results = []
missing_pristine_rows = []
missing_fe_atom_rows = []
pristine_energy = None
fe_atom_energy = None

# Get CIF files
atom_cifs = {file.replace('_atom.cif', ''): file for file in os.listdir(cif_directory) if file.endswith('_atom.cif')}
doped_cifs = {file.split('Fe_dopedwith_')[1].replace('.cif', ''): file for file in os.listdir(cif_directory) if file.startswith('Fe_dopedwith_')}

# Function to modify the slab by adsorbing CO2
def modify_slab_with_CO2(slab):
    CO2 = Atoms('CO2', positions=[
        [0.0, 0.0, 0.0],      # Carbon
        [0.0, 0.76, 0.58],    # Oxygen
        [0.0, -0.76, 0.58]    # Oxygen
    ])

    atom_number = 26
    adsorption_site = slab[atom_number].position
    adsorption_height = 1.8  

    CO2.translate(adsorption_site + np.array([0, 0, adsorption_height]))

    if CO2.positions[:, 2].min() < slab.positions[:, 2].max():
        z_shift = slab.positions[:, 2].max() - CO2.positions[:, 2].min() + 1.0
        CO2.translate([0, 0, z_shift])

    slab += CO2
    return slab

# Process each doped slab
for dopant, slab_file in doped_cifs.items():
    try:
        slab_atoms = read(os.path.join(cif_directory, slab_file))
        slab_atoms.pbc = [True, True, False]

        # Define surface atoms and apply constraints
        surface_threshold = 1.5
        min_z = min(slab_atoms.positions[:, 2])
        surface_atoms = [atom.index for atom in slab_atoms if atom.position[2] > min_z + surface_threshold]
        mask = [atom.index not in surface_atoms for atom in slab_atoms]
        slab_atoms.set_constraint(FixAtoms(mask=mask))

        # Set up and run slab optimization
        slab_calc = GPAW(mode=PW(50),
                         xc='PBE',
                         kpts=(2, 2, 1),
                         occupations=FermiDirac(0.03),
                         txt=os.path.join(output_directory, f'gpaw_output_{slab_file}.txt'),
                         spinpol=True)

        slab_atoms.calc = slab_calc
        relaxer = BFGS(slab_atoms)
        relaxer.run(fmax=5)

        optimized_cif_path = os.path.join(optimized_directory, slab_file)
        write(optimized_cif_path, slab_atoms)

        slab_gpw_path = os.path.join(gpw_directory, f'Fe_dopedwith_{dopant}.gpw')
        slab_calc.write(slab_gpw_path)

        # Extract properties
        total_energy = slab_atoms.get_potential_energy()
        homo, lumo = slab_calc.get_homo_lumo()
        efermi = slab_calc.get_fermi_level()
        bandgap = lumo - homo
        chemical_potential = (lumo + homo) / 2
        entropy_st, free_energy = np.nan, np.nan

        with open(os.path.join(output_directory, f'gpaw_output_{slab_file}.txt'), 'r') as file:
            for line in file:
                if "Entropy (-ST):" in line:
                    entropy_st = float(line.split()[-1])
                elif "Free energy:" in line:
                    free_energy = float(line.split()[-1])

        entropy = entropy_st / 0.03

        # Check if this is pristine Fe slab
        if slab_file.endswith('_Fe.cif'):
            pristine_energy = total_energy
            for row in missing_pristine_rows:
                row['Pristine Energy (eV)'] = pristine_energy

        # Load optimized slab for CO2 adsorption
        slab_atoms = read(optimized_cif_path)
        mask = [atom.index for atom in slab_atoms]
        constraint = FixAtoms(mask=mask)
        slab_atoms.set_constraint(constraint)
        co2_slab_atoms = modify_slab_with_CO2(slab_atoms)
        co2_slab_atoms.pbc = [True, True, False]

        # Set up and run CO2 slab optimization
        co2_slab_calc = GPAW(mode=PW(50),
                             xc='PBE',
                             kpts=(2, 2, 1),
                             occupations=FermiDirac(0.03),
                             txt=os.path.join(output_directory, f'gpaw_output_CO2_{slab_file}.txt'),
                             spinpol=True)

        co2_slab_atoms.calc = co2_slab_calc
        relaxer = BFGS(co2_slab_atoms)
        relaxer.run(fmax=5)
        optimized_cif_path = os.path.join(optimized_directory, f'{co2_slab_file}')
        write(optimized_cif_path, co2_slab_atoms)

        co2_slab_energy = co2_slab_atoms.get_potential_energy()
        co2_slab_gpw_path = os.path.join(gpw_directory, f'co2_Fe_dopedwith_{dopant}.gpw')
        co2_slab_calc.write(co2_slab_gpw_path)

        # Process individual dopant atom energy
        atom_file = atom_cifs.get(dopant, None)
        atom_energy = None

        if atom_file:
            atom_atoms = read(os.path.join(cif_directory, atom_file))
            atom_atoms.pbc = [True, True, True]
            atom_calc = GPAW(mode=PW(50),
                             xc='PBE',
                             kpts=(1, 1, 1),
                             occupations=FermiDirac(0.03),
                             txt=os.path.join(output_directory, f'gpaw_output_{atom_file}.txt'))

            atom_atoms.calc = atom_calc
            atom_energy = atom_atoms.get_potential_energy()
            atom_gpw_path = os.path.join(gpw_directory, f'atom_{dopant}.gpw')
            atom_calc.write(atom_gpw_path)

        # Store Fe atom energy if applicable
        if atom_file and atom_file.startswith('Fe_atom'):
            fe_atom_energy = atom_energy
            for row in missing_fe_atom_rows:
                row['Fe Atom Energy (eV)'] = fe_atom_energy

        # Store results
        results.append({
            'CIF file': slab_file,
            'Atom CIF file': atom_file if atom_file else np.nan,
            'Total Energy (eV)': total_energy,
            'Dopant Energy (eV)': atom_energy if atom_energy else np.nan,
            'CO2 Slab Energy (eV)': co2_slab_energy,
            'HOMO (eV)': homo,
            'LUMO (eV)': lumo,
            'Fermi Level (eV)': efermi,
            'Band Gap (eV)': bandgap,
            'Chemical Potential (eV)': chemical_potential,
            'Entropy (eV/K)': entropy,
            'Free Energy (eV)': free_energy,
            'Pristine Energy (eV)': pristine_energy if pristine_energy is not None else np.nan
        })

    except Exception as e:
        print(f"Error processing {slab_file}: {e}")

# Save results to CSV
df_results = pd.DataFrame(results)
df_results.to_csv(os.path.join(output_directory, 'DFT_results.csv'), index=False)

print("Processing complete. Results saved.")

