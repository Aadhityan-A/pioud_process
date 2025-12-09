import os
import numpy as np
from ase.io import write
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from . import mbtr_analysis
from . import pioud_extract
from . import soap_analysis


def pioud_process(path, centroid=False, descriptor=None, des_parameters=None, list_of_folder=False, start_idx=0, end_idx=None):
    """
    Process PIOUD output files and optionally perform descriptor analysis.

    Args:
        path: Path to the directory containing PIOUD output files, or list of paths
        centroid: If True, use centroid positions and forces
        descriptor: Type of descriptor analysis ('soap', 'mbtr', or None)
        des_parameters: Parameters for descriptor analysis
        list_of_folder: If True, path is a list of folders to process
        start_idx: Starting index for slicing data arrays (default: 0)
        end_idx: Ending index for slicing data arrays (default: None, meaning end of array)
    """
    ### Extraction of values from pioud.x output ###
    if list_of_folder == False:
        # Provide the file name inside the function to change default file name.
        # The default file names: positions-positions_cen.dat; forces-forces_cen.dat
        # energies-pimd.out; array_of_atomic_number-pw_1.in
        # path='./'
        cell_array = pioud_extract.extract_cell_from_qe(filename=os.path.join(path, "pw_1.in"))

        if centroid:
            positions = pioud_extract.read_position_dat(filename=os.path.join(path, 'positions_cen.dat'))[start_idx:end_idx]
            forces = pioud_extract.read_forces_dat(filename=os.path.join(path, 'forces_cen.dat'))[start_idx:end_idx]
            energies = pioud_extract.read_energies(filename=os.path.join(path, 'pimd.out'))[start_idx:end_idx]
        else:
            positions = pioud_extract.read_position_dat(filename=os.path.join(path, 'positions.dat'))[start_idx:end_idx]
            forces = pioud_extract.read_forces_dat(filename=os.path.join(path, 'forces.dat'))[start_idx:end_idx]
            energies = pioud_extract.read_energies(filename=os.path.join(path, 'pimd.out'))[start_idx:end_idx]

        atomic_numbers_dict, array_of_atomic_number = pioud_extract.extract_atomic_numbers(input_file=os.path.join(path, "pw_1.in"))

        atomic_numbers = array_of_atomic_number
        total_atoms = len(array_of_atomic_number)

        if not (len(positions) == len(forces) == len(energies)):
            print("positions length", len(positions))
            print("forces length", len(forces))
            print("energies length", len(energies))
            raise Exception("Number of rows in positions, forces and energies array are not equal")

        no_of_configs = len(energies)
        positions = positions.reshape(no_of_configs, total_atoms, 3)
        forces = forces.reshape(no_of_configs, total_atoms, 3)

        ase_configs_list = []
        for config in range(no_of_configs):
            atoms = Atoms(positions=positions[config], numbers=atomic_numbers, cell=cell_array, pbc=True)
            sp_calc = SinglePointCalculator(energy=energies[config], forces=forces[config], atoms=atoms)
            atoms.set_calculator(sp_calc)
            ase_configs_list.append(atoms)

        write(os.path.join(path, 'ase_configs.xyz'), ase_configs_list)
        # write first config as xsf and cif for sanity check of the structure
        write(os.path.join(path, 'ase_config_1.xsf'), ase_configs_list[0])
        write(os.path.join(path, 'ase_config_1.cif'), ase_configs_list[0])

    else:
        ase_configs_list = []
        for each_path in path:
            cell_array = pioud_extract.extract_cell_from_qe(filename=os.path.join(each_path, "pw_1.in"))
            if centroid:
                positions = pioud_extract.read_position_dat(filename=os.path.join(each_path, 'positions_cen.dat'))[start_idx:end_idx]
                forces = pioud_extract.read_forces_dat(filename=os.path.join(each_path, 'forces_cen.dat'))[start_idx:end_idx]
                energies = pioud_extract.read_energies(filename=os.path.join(each_path, 'pimd.out'))[start_idx:end_idx]
            else:
                positions = pioud_extract.read_position_dat(filename=os.path.join(each_path, 'positions.dat'))[start_idx:end_idx]
                forces = pioud_extract.read_forces_dat(filename=os.path.join(each_path, 'forces.dat'))[start_idx:end_idx]
                energies = pioud_extract.read_energies(filename=os.path.join(each_path, 'pimd.out'))[start_idx:end_idx]

            atomic_numbers_dict, array_of_atomic_number = pioud_extract.extract_atomic_numbers(input_file=os.path.join(each_path, "pw_1.in"))

            atomic_numbers = array_of_atomic_number
            total_atoms = len(array_of_atomic_number)

            if not (len(positions) == len(forces) == len(energies)):
                print("positions length", len(positions))
                print("forces length", len(forces))
                print("energies length", len(energies))
                raise Exception("Number of rows in positions, forces and energies array are not equal")

            no_of_configs = len(energies)
            positions = positions.reshape(no_of_configs, total_atoms, 3)
            forces = forces.reshape(no_of_configs, total_atoms, 3)

            for config in range(no_of_configs):
                atoms = Atoms(positions=positions[config], numbers=atomic_numbers, cell=cell_array, pbc=True)
                sp_calc = SinglePointCalculator(energy=energies[config], forces=forces[config], atoms=atoms)
                atoms.set_calculator(sp_calc)
                ase_configs_list.append(atoms)

        write(os.path.join(path[0], 'ase_configs.xyz'), ase_configs_list)
        # write first config as xsf and cif for sanity check of the structure
        write(os.path.join(path[0], 'ase_config_1.xsf'), ase_configs_list[0])
        write(os.path.join(path[0], 'ase_config_1.cif'), ase_configs_list[0])

    ### Perform SOAP MBTR analysis ###
    if descriptor == None:
        output_file_name_no_analysis = 'ase_configs.xyz'
        return output_file_name_no_analysis
    elif descriptor == 'soap':
        if list_of_folder:
            output_file_name_after_soap_analysis = soap_analysis.soap_analysis_from_pioudop(path[0], des_parameters)
            print("Output filename:", output_file_name_after_soap_analysis)
            return output_file_name_after_soap_analysis
        else:
            output_file_name_after_soap_analysis = soap_analysis.soap_analysis_from_pioudop(path, des_parameters)
            print("Output filename:", output_file_name_after_soap_analysis)
            return output_file_name_after_soap_analysis
    elif descriptor == 'mbtr':
        if list_of_folder:
            output_file_name_after_mbtr_analysis = mbtr_analysis.mbtr_analysis_from_pioudop(path[0], des_parameters)
            print("Output filename:", output_file_name_after_mbtr_analysis)
            return output_file_name_after_mbtr_analysis
        else:
            output_file_name_after_mbtr_analysis = mbtr_analysis.mbtr_analysis_from_pioudop(path, des_parameters)
            print("Output filename:", output_file_name_after_mbtr_analysis)
            return output_file_name_after_mbtr_analysis