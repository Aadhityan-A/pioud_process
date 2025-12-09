import numpy as np
import os
import re
from ase.units import Bohr, Ha

# To extract the atomic numbers from pw.x input file.
def extract_atomic_numbers(input_file="pw_1.in"):
    """
    Extracts the atomic number for each element defined in the ATOMIC_SPECIES
    section of a Quantum Espresso pw.in input file.

    Args:
        input_file (str): The path to the Quantum Espresso pw.in input file.

    Returns:
        dict: A dictionary where keys are the atomic symbols and values are
              their corresponding atomic numbers. Returns an empty dictionary
              if the ATOMIC_SPECIES section is not found or if there's an error.
    """
    atomic_numbers = {}
    array_atomic_number= np.array([])
    try:
        with open(input_file, 'r') as f:
            in_atomic_species = False
            for line in f:
                line = line.strip()
                if line.startswith("ATOMIC_POSITIONS"):
                    in_atomic_species = True
                    continue
                elif line.startswith("/"):
                    in_atomic_species = False
                    continue

                if in_atomic_species and line and not line.startswith("!"):
                    parts = line.split()
                    if len(parts) >= 3:
                        symbol = parts[0]
                        # You'll need a way to map the symbol to the atomic number.
                        # This can be done using a dictionary or a library.
                        atomic_number = get_atomic_number(symbol)
                        if atomic_number is not None:
                            atomic_numbers[symbol] = atomic_number
                            array_atomic_number = np.append(array_atomic_number,atomic_number)
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    return atomic_numbers,array_atomic_number

def get_atomic_number(symbol):
    """
    Returns the atomic number for a given atomic symbol.

    Args:
        symbol (str): The atomic symbol (e.g., 'Si', 'O', 'H').

    Returns:
        int or None: The atomic number if the symbol is recognized, None otherwise.
    """
    atomic_data = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
        'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
        'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
        'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
        'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
        'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
        'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
        'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
        'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
        'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
        'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
        'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105,
        'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
        'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
    }
    return atomic_data.get(symbol)


def read_position_dat(filename='positions_cen.dat', output_in_angs=True):
    positions = np.loadtxt(filename)    # Bohr
    if output_in_angs:
        positions = positions * Bohr

    return positions

def read_velocities_dat(filename='velocities_cen.dat'):
    velocities = np.loadtxt(filename)
    return velocities

def read_stress_dat(filename='stress_cen.dat', output_in_ev_per_angs=True):
    stress = np.loadtxt(filename)    # Ha/Bohr
    if output_in_ev_per_angs:
        stress = stress * Ha / Bohr

    return stress
    
def read_forces_dat(filename='forces_cen.dat', output_in_ev_per_angs=True):
    forces = np.loadtxt(filename)    # Ha/Bohr
    if output_in_ev_per_angs:
        forces = forces * Ha / Bohr

    return forces

def read_energies(filename='pimd.out', output_in_ev=True):
    energies = np.loadtxt(filename, usecols=4, skiprows=20)    # Ha 
    #Current version skiprows=20 check.
    if output_in_ev:
        energies = energies * Ha
    # with open(filename, 'r') as f:
    #     lines = f.readlines()

    # # Filter only the lines that start with spaces and contain actual data
    # data_lines = []
    # for line in lines:
    #     # Check if line starts with spaces and contains numbers
    #     if line.startswith('1)block'):  # Lines with actual data start with multiple spaces
    #         data_lines.append(line.strip())

    # # Process each line manually to create a list of lists
    # data = []
    # for line in data_lines:
    #     # Split the line and convert each element to float
    #     values = [float(x) for x in line.split()]
    #     data.append(values)

    # # Convert to numpy array
    # data_array = np.array(data)

    # # Extract the 4th column (index 3)
    # energies = data_array[:, 4] * Ha

    return energies

# To extract the cell paramerters from pw.x input file.
"""This will extarct cell parameters directly. If it i swritten in the pw_1.in file. 
Else it will try to extarct from the ibrav, alat like values. 
Warning: It seems, currenetly this is not fully functional. 
But the direct extraction seems to be work fine."""

def extract_cell_parameters(qe_input_file):
    """
    Extract cell parameters from a Quantum ESPRESSO input file.

    Parameters:
    -----------
    qe_input_file : str
        Path to the Quantum ESPRESSO input file (pw.in) or input file content as string

    Returns:
    --------
    numpy.ndarray
        3x3 array containing the cell parameters (lattice vectors) in Angstroms
    """
    try:
        if os.path.isfile(qe_input_file):
            with open(qe_input_file, 'r') as f:
                content = f.read()
        else:
            # Assume input is the file content string
            content = qe_input_file
    except Exception as e:
        raise FileNotFoundError(f"Error reading input file: {e}")

    # Extract ibrav parameter
    ibrav_match = re.search(r'ibrav\s*=\s*(-?\d+)', content, re.IGNORECASE)
    if not ibrav_match:
        # Try CELL_PARAMETERS block if ibrav is not specified
        return extract_from_cell_parameters_block(content)

    ibrav = int(ibrav_match.group(1))

    if ibrav == 0:
        # For ibrav=0, cell parameters must be specified explicitly
        return extract_from_cell_parameters_block(content)

    # Extract lattice parameters
    params = extract_lattice_parameters(content)
    a = params['a']
    b = params['b']
    c = params['c']
    cosab = params['cosab']
    cosac = params['cosac']
    cosbc = params['cosbc']

    return generate_cell_from_ibrav(ibrav, a, b, c, cosab, cosac, cosbc)

def extract_lattice_parameters(content):
    """Extract all lattice parameters from QE input"""
    params = {
        'a': 0.0, 'b': 0.0, 'c': 0.0,
        'cosab': 0.0, 'cosac': 0.0, 'cosbc': 0.0
    }

    # Extract a or celldm(1)
    a_match = re.search(r'a\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    celldm1_match = re.search(r'celldm\s*\(\s*1\s*\)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)

    if a_match:
        params['a'] = float(a_match.group(1))  # Already in Angstrom
    elif celldm1_match:
        params['a'] = float(celldm1_match.group(1)) * 0.529177  # Convert from Bohr to Angstrom
    else:
        raise ValueError("Could not find lattice parameter 'a' or 'celldm(1)'")

    # Extract b/a ratio or directly b
    b_match = re.search(r'b\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    celldm2_match = re.search(r'celldm\s*\(\s*2\s*\)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)

    if b_match:
        params['b'] = float(b_match.group(1))
    elif celldm2_match:
        params['b'] = float(celldm2_match.group(1)) * params['a']
    else:
        params['b'] = params['a']  # Default: b = a

    # Extract c/a ratio or directly c
    c_match = re.search(r'c\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    celldm3_match = re.search(r'celldm\s*\(\s*3\s*\)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)

    if c_match:
        params['c'] = float(c_match.group(1))
    elif celldm3_match:
        params['c'] = float(celldm3_match.group(1)) * params['a']
    else:
        params['c'] = params['a']  # Default: c = a

    # Extract cosines of angles
    cosab_match = re.search(r'cosab\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    celldm4_match = re.search(r'celldm\s*\(\s*4\s*\)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    if cosab_match:
        params['cosab'] = float(cosab_match.group(1))
    elif celldm4_match:
        params['cosab'] = float(celldm4_match.group(1))

    cosac_match = re.search(r'cosac\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    celldm5_match = re.search(r'celldm\s*\(\s*5\s*\)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    if cosac_match:
        params['cosac'] = float(cosac_match.group(1))
    elif celldm5_match:
        params['cosac'] = float(celldm5_match.group(1))

    cosbc_match = re.search(r'cosbc\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    celldm6_match = re.search(r'celldm\s*\(\s*6\s*\)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    if cosbc_match:
        params['cosbc'] = float(cosbc_match.group(1))
    elif celldm6_match:
        params['cosbc'] = float(celldm6_match.group(1))

    return params

def extract_from_cell_parameters_block(content):
    """Extract cell parameters from the CELL_PARAMETERS block"""
    cell_block_match = re.search(r'CELL_PARAMETERS\s*(?:\{([^}]*)\})?(\s\S]*?)(?=\n\s*\n|\Z)', content, re.IGNORECASE)
    if not cell_block_match:
        return None

    # Extract unit (default is alat)
    unit = cell_block_match.group(1).strip() if cell_block_match.group(1) else 'alat'

    # Extract the cell parameters
    cell_block = cell_block_match.group(2).strip().split('\n')
    cell_lines = []

    for line in cell_block:
        line = line.strip()
        if not line or line.startswith('!') or line.startswith('#'):
            continue

        try:
            values = [float(val) for val in line.split()[:3]]
            if len(values) == 3:
                cell_lines.append(values)
        except ValueError:
            continue

    if len(cell_lines) != 3:
        return None

    cell_params = np.array(cell_lines)

    # Adjust for different units
    if unit.lower() == 'alat':
        alat = extract_alat(content)
        cell_params *= alat
    elif unit.lower() == 'bohr':
        cell_params *= 0.529177  # Convert to Angstrom

    return cell_params

def generate_cell_from_ibrav(ibrav, a, b, c, cosab, cosac, cosbc):
    """
    Generate the cell matrix based on ibrav and lattice parameters.

    All ibrav types from Quantum ESPRESSO are implemented.
    Returns cell parameters in Angstroms.
    """
    # Helper calculations for trigonometric values
    if abs(cosab) > 1.0 or abs(cosac) > 1.0 or abs(cosbc) > 1.0:
        raise ValueError("Cosine values must be between -1 and 1")

    sinab = np.sqrt(1 - cosab**2) if abs(cosab) < 1.0 else 0
    sinac = np.sqrt(1 - cosac**2) if abs(cosac) < 1.0 else 0
    sinbc = np.sqrt(1 - cosbc**2) if abs(cosbc) < 1.0 else 0

    # Handle all ibrav types
    if ibrav == 1:  # Simple cubic
        return np.array([
            [a, 0.0, 0.0],
            [0.0, a, 0.0],
            [0.0, 0.0, a]
        ])

    elif ibrav == 2:  # Face-centered cubic
        return np.array([
            [0.0, a/2, a/2],
            [a/2, 0.0, a/2],
            [a/2, a/2, 0.0]
        ])

    elif ibrav == 3:  # Body-centered cubic
        return np.array([
            [-a/2, a/2, a/2],
            [a/2, -a/2, a/2],
            [a/2, a/2, -a/2]
        ])

    elif ibrav == 4:  # Hexagonal
        return np.array([
            [a, 0.0, 0.0],
            [-a/2, a*np.sqrt(3)/2, 0.0],
            [0.0, 0.0, c]
        ])

    elif ibrav == 5:  # Trigonal or Rhombohedral
        tx = np.sqrt((1-cosab)/2)
        ty = np.sqrt((1-cosab)/6)
        tz = np.sqrt((1+2*cosab)/3)

        return np.array([
            [a, 0.0, 0.0],
            [-a/2, a*np.sqrt(3)/2, 0.0],
            [a*tx, a*ty, a*tz]
        ])

    elif ibrav == -5:  # Trigonal or Rhombohedral (alternative)
        tx = np.sqrt((1-cosab)/2)
        ty = np.sqrt((1-cosab)/6)
        tz = np.sqrt((1+2*cosab)/3)
        u = tz - 2*np.sqrt(2)*ty
        v = tz + np.sqrt(2)*ty

        return np.array([
            [a*tx, a*ty, a*tz],
            [-a*tx, a*ty, a*tz],
            [0.0, -a*v, a*u]
        ])

    elif ibrav == 6:  # Simple tetragonal
        return np.array([
            [a, 0.0, 0.0],
            [0.0, a, 0.0],
            [0.0, 0.0, c]
        ])

    elif ibrav == 7:  # Body-centered tetragonal
        return np.array([
            [a/2, -a/2, c/2],
            [a/2, a/2, c/2],
            [-a/2, -a/2, c/2]
        ])

    elif ibrav == 8:  # Simple orthorhombic
        return np.array([
            [a, 0.0, 0.0],
            [0.0, b, 0.0],
            [0.0, 0.0, c]
        ])

    elif ibrav == 9:  # Base-centered orthorhombic (C-type)
        return np.array([
            [a/2, b/2, 0.0],
            [-a/2, b/2, 0.0],
            [0.0, 0.0, c]
        ])

    elif ibrav == -9:  # Alternative base-centered orthorhombic (A-type)
        return np.array([
            [a, 0.0, 0.0],
            [0.0, b/2, c/2],
            [0.0, -b/2, c/2]
        ])

    elif ibrav == 10:  # Face-centered orthorhombic
        return np.array([
            [a/2, 0.0, c/2],
            [a/2, b/2, 0.0],
            [0.0, b/2, c/2]
        ])

    elif ibrav == 11:  # Body-centered orthorhombic
        return np.array([
            [a/2, b/2, c/2],
            [-a/2, b/2, c/2],
            [-a/2, -b/2, c/2]
        ])

    elif ibrav == 12:  # Simple monoclinic (unique c-axis)
        # gamma is the angle between a and b vectors
        gamma = np.arccos(cosab)
        return np.array([
            [a, 0.0, 0.0],
            [b*cosab, b*sinab, 0.0],
            [0.0, 0.0, c]
        ])

    elif ibrav == -12:  # Simple monoclinic (unique b-axis)
        # alpha is the angle between b and c vectors
        alpha = np.arccos(cosbc)
        return np.array([
            [a, 0.0, 0.0],
            [0.0, b, 0.0],
            [0.0, c*cosbc, c*sinbc]
        ])

    elif ibrav == 13:  # Monoclinic (unique c-axis, base centered)
        gamma = np.arccos(cosab)
        return np.array([
            [a/2, -b*sinab/2, 0.0],
            [a/2, b*sinab/2, 0.0],
            [0.0, 0.0, c]
        ])

    elif ibrav == -13:  # Monoclinic (unique b-axis, base centered)
        alpha = np.arccos(cosbc)
        return np.array([
            [a/2, 0.0, -c*sinbc/2],
            [a/2, 0.0, c*sinbc/2],
            [0.0, b, 0.0]
        ])

    elif ibrav == 14:  # Triclinic
        # We need all three cosines to calculate the cell
        # Calculate all angles in radians
        gamma = np.arccos(cosab)
        beta = np.arccos(cosac)
        alpha = np.arccos(cosbc)

        # Helper calculations (volume term)
        omega = np.sqrt(1.0 + 2.0*cosbc*cosac*cosab - cosbc**2 - cosac**2 - cosab**2)

        # First vector is along x-axis
        v1 = np.array([a, 0.0, 0.0])

        # Second vector is in xy-plane
        v2 = np.array([b*cosab, b*sinab, 0.0])

        # Third vector has all components
        v3x = c*cosac
        v3y = c*(cosbc - cosac*cosab)/sinab
        v3z = c*omega/sinab
        v3 = np.array([v3x, v3y, v3z])

        return np.vstack((v1, v2, v3))

    else:
        raise ValueError(f"ibrav = {ibrav} is not a valid Quantum ESPRESSO lattice type")

def extract_alat(content):
    """Extract the lattice parameter 'alat'"""
    a_match = re.search(r'a\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    if a_match:
        return float(a_match.group(1))

    celldm1_match = re.search(r'celldm\s*\(\s*1\s*\)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content, re.IGNORECASE)
    if celldm1_match:
        return float(celldm1_match.group(1)) * 0.529177  # Bohr to Angstrom

    return 1.0


def extract_cell_from_qe(filename='pw_1.in'):
    """This will extarct cell parameters directly. If it i swritten in the pw_1.in file. Else it woill try to extarct from the 
    ibrav, alat like values. Warning: It seems, currenetly this is not fully functional. But the direct extraction  seems to be work fine."""
    cell_vectors = []
    reading_cell = False
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        for line in lines:
            # Look for CELL_PARAMETERS block
            if 'CELL_PARAMETERS' in line.upper():
                reading_cell = True
                # Check if units are specified
                if 'alat' in line.lower():
                    print("Note: Cell parameters are in units of alat")
                elif 'angstrom' in line.lower():
                    print("Note: Cell parameters are in Angstroms")
                elif 'bohr' in line.lower():
                    print("Note: Cell parameters are in Bohr")
                continue
                
            # If we're in the cell parameters block, read the vectors
            if reading_cell and len(cell_vectors) < 3:
                # Try to convert line to vector
                try:
                    vector = [float(x) for x in line.split()]
                    if len(vector) == 3:
                        cell_vectors.append(vector)
                except ValueError:
                    continue
                    
        if cell_vectors:
            cell_array = np.array(cell_vectors)
            print("\nCell vectors as NumPy array:")
            print(cell_array)
            return cell_array
        else:
            cell_array = extract_cell_parameters(filename)
            print("No cell parameters found in the file!")
            return cell_array
            
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found!")
        return None
    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return None