from ase import Atoms
from ase.io import write
import os

# List of elements
elements = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
            'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Pb', 'Bi','Sn', 'Sb']


vacuum = 10.0


output_dir = './structs'
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist


for element in elements:

    atom = Atoms(element)
    
    
    atom.set_cell([vacuum, vacuum, vacuum], scale_atoms=True)
    atom.center()


    filename = os.path.join(output_dir, f'{element}_atom.cif')
    
    # Save the structure to a CIF file
    write(filename, atom, format='cif')

print("done")
