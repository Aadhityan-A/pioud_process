# pioud_process

A Python package for generating input and processing output of PIOUD (A Path Integral Molecular Dynamicscode integrated with Quantum ESPRESSO).

## Features

- **GUI Interface**: A user-friendly graphical interface for:
  - Generating input files for PIOUD from and .cif crystal structure files.
  - Creating submission scripts for HPC clusters (PBS/SLURM).
  - Visualizing and setting up simulation parameters.
- **Analysis Tools**:
  - Time averaging of properties.
  - Radial Distribution Function (RDF) analysis.
  - SOAP (Smooth Overlap of Atomic Positions) analysis.
  - MBTR (Many-Body Tensor Representation) analysis.
  - Autocorrelation function calculation.

## Installation

To install the package, clone the repository and install it using pip:

```bash
git clone https://github.com/Aadhityan-A/pioud_process.git
cd pioud_process
pip install .
```

For development mode, use:
    
```bash
pip install -e .
```

## Usage

### Graphical User Interface (GUI)

After installation, you can launch the GUI using the command line:

```bash
pioud-process
```

The GUI allows you to:
1.  **File Input**: Upload CIF files and set calculation prefixes.
2.  **PIMD Parameters**: Configure Path Integral Molecular Dynamics settings (beads, blocks, steps, temperature, etc.).
3.  **QE Parameters**: Set up Quantum ESPRESSO SCF parameters (cutoffs, k-points, smearing, etc.).
4.  **HPCC Submit**: Generate submission scripts for PBS or SLURM schedulers.

### Library Usage

You can also use the analysis modules directly in your Python scripts:

```python
from pioud_process.time_average import calculate_qe_time_averaged_properties

# Example usage
# results = calculate_qe_time_averaged_properties(path_to_data, start_step, end_step)
```

## Dependencies

- numpy
- ase
- scipy
- dscribe
- skmatter
- tkinter (usually included with Python)

## License

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info".

See the [LICENSE](LICENSE) file for details.
