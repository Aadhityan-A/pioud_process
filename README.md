# PIOUD SOAP

A Python package for PIOUD (Path Integral Quantum Dynamics) input generator and output analyzer.

## Installation

You can install the package directly from PyPI:

```bash
pip install pioud_process
```

Or you can install from source:

```bash
git clone https://github.com/yourusername/pioud_soap.git
cd pioud_process
pip install -e .
```

## Usage

See the jupyter notebooks inside the examples folder.

### GUI
- To generate Input
- - open gui_pioud_input.ipynb under the examples folder.
- To process output
```
python gui/gui_pioud_process.py
```

## Required Input Files

The package expects the following files in the specified path to process output:

- `pw_1.in`: Quantum ESPRESSO input file
- `positions.dat`: Atomic positions
- `forces.dat`: Forces on atoms
- `pimd.out`: Energy output file

## Dependencies

- numpy
- ase (Atomic Simulation Environment)
- scipy
- dscribe
- scikit-matter

## License

MIT