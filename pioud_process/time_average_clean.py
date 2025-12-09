import numpy as np
import os
from . import pioud_extract


def calculate_qe_time_averaged_properties(output_folder_path, start_row, end_row):
    """
    Calculate time-averaged properties (positions, velocities, forces, and stress) 
    from PIOUD molecular dynamics simulation output files.
    
    This function reads data from Quantum Espresso format output files where each row 
    contains x,y,z coordinates/forces/velocities for all atoms at a particular timestep,
    and computes time averages over the specified range.
    
    Parameters:
    -----------
    output_folder_path : str
        Path to the folder containing the PIOUD output files
        (positions_cen.dat, velocities_cen.dat, forces_cen.dat, stress_cen.dat)
    start_row : int
        Starting row number (0-indexed) for time averaging
    end_row : int
        Ending row number (0-indexed, exclusive) for time averaging
        
    Returns:
    --------
    dict
        Dictionary containing time-averaged properties with keys:
        'positions', 'velocities', 'forces', 'stress'
        
    Notes:
    ------
    - Input files are expected to be in Quantum Espresso format
    - Each row contains data for all atoms at one timestep
    - Output files are saved with '_timeaveraged.dat' suffix
    - Units are preserved from input files (typically Angstrom for positions, 
      eV/Angstrom for forces)
    """
    
    # Validate input parameters
    if start_row < 0:
        raise ValueError("start_row must be non-negative")
    if end_row <= start_row:
        raise ValueError("end_row must be greater than start_row")
    
    # Dictionary to store results
    results = {}
    
    # File mappings for input and output
    file_mappings = {
        'positions': {
            'input': 'positions_cen.dat',
            'output': 'positions_cen_timeaveraged.dat',
            'reader': pioud_extract.read_position_dat
        },
        'velocities': {
            'input': 'velocities_cen.dat', 
            'output': 'velocities_cen_timeaveraged.dat',
            'reader': pioud_extract.read_velocities_dat
        },
        'forces': {
            'input': 'forces_cen.dat',
            'output': 'forces_cen_timeaveraged.dat', 
            'reader': pioud_extract.read_forces_dat
        },
        'stress': {
            'input': 'stress_cen.dat',
            'output': 'stress_cen_timeaveraged.dat',
            'reader': pioud_extract.read_stress_dat
        }
    }
    
    print(f"Calculating time averages from row {start_row} to {end_row-1}")
    print(f"Processing data from: {output_folder_path}")
    
    # Process each property
    for property_name, file_info in file_mappings.items():
        input_file = os.path.join(output_folder_path, file_info['input'])
        output_file = os.path.join(output_folder_path, file_info['output'])
        
        # Check if input file exists
        if not os.path.exists(input_file):
            print(f"Warning: {input_file} not found. Skipping {property_name}.")
            continue
            
        try:
            # Read data using the appropriate reader function
            print(f"Reading {property_name} from {file_info['input']}...")
            
            if property_name == 'positions':
                data = file_info['reader'](filename=input_file, output_in_angs=True)
            elif property_name == 'forces':
                data = file_info['reader'](filename=input_file, output_in_ev_per_angs=True)
            elif property_name == 'stress':
                data = file_info['reader'](filename=input_file, output_in_ev_per_angs=True)
            else:  # velocities
                data = file_info['reader'](filename=input_file)
            
            # Validate row indices
            if end_row > len(data):
                print(f"Warning: end_row ({end_row}) exceeds data length ({len(data)}) for {property_name}")
                print(f"Using available data up to row {len(data)-1}")
                actual_end_row = len(data)
            else:
                actual_end_row = end_row
                
            if start_row >= len(data):
                raise ValueError(f"start_row ({start_row}) exceeds data length ({len(data)}) for {property_name}")
            
            # Extract the specified time range
            data_slice = data[start_row:actual_end_row]
            
            # Calculate time average
            time_averaged_data = np.mean(data_slice, axis=0)
            
            # Store result
            results[property_name] = time_averaged_data
            
            # Save time-averaged data to file
            save_time_averaged_data(time_averaged_data, output_file, property_name)
            
            print(f"Time-averaged {property_name} saved to {file_info['output']}")
            print(f"  Shape: {time_averaged_data.shape}")
            print(f"  Averaged over {actual_end_row - start_row} timesteps")
            
        except Exception as e:
            print(f"Error processing {property_name}: {str(e)}")
            continue
    
    # Print summary
    print("\nTime averaging completed successfully!")
    print(f"Processed properties: {list(results.keys())}")
    print(f"Results saved in: {output_folder_path}")
    
    return results


def save_time_averaged_data(data, output_file, property_name):
    """
    Save time-averaged data to file in Quantum Espresso format.
    
    Parameters:
    -----------
    data : numpy.ndarray
        Time-averaged data array
    output_file : str
        Path to output file
    property_name : str
        Name of the property being saved (for header comment)
    """
    
    # Create header comment
    header = f"# Time-averaged {property_name} data\n"
    header += f"# Generated by PIOUD time averaging module\n"
    header += f"# Data shape: {data.shape}\n"
    
    # Handle different data shapes
    if len(data.shape) == 1:
        # 1D data (e.g., stress tensor components)
        formatted_data = data.reshape(1, -1)
    else:
        # 2D data (positions, forces, velocities)
        formatted_data = data
    
    # Save data with appropriate formatting
    try:
        # Use scientific notation with sufficient precision
        np.savetxt(output_file, formatted_data, 
                  fmt='%15.8e', 
                  header=header.strip(),
                  comments='')
        
    except Exception as e:
        print(f"Error saving {output_file}: {str(e)}")
        raise


def validate_time_averaging_inputs(output_folder_path, start_row, end_row):
    """
    Validate inputs for time averaging calculation.
    
    Parameters:
    -----------
    output_folder_path : str
        Path to output folder
    start_row : int
        Starting row index
    end_row : int
        Ending row index
        
    Returns:
    --------
    bool
        True if inputs are valid, raises exception otherwise
    """
    
    # Check if folder exists
    if not os.path.exists(output_folder_path):
        raise FileNotFoundError(f"Output folder not found: {output_folder_path}")
    
    if not os.path.isdir(output_folder_path):
        raise NotADirectoryError(f"Path is not a directory: {output_folder_path}")
    
    # Check row indices
    if not isinstance(start_row, int) or not isinstance(end_row, int):
        raise TypeError("start_row and end_row must be integers")
    
    if start_row < 0:
        raise ValueError("start_row must be non-negative")
    
    if end_row <= start_row:
        raise ValueError("end_row must be greater than start_row")
    
    # Check if at least one required file exists
    required_files = ['positions_cen.dat', 'velocities_cen.dat', 'forces_cen.dat', 'stress_cen.dat']
    found_files = []
    
    for filename in required_files:
        filepath = os.path.join(output_folder_path, filename)
        if os.path.exists(filepath):
            found_files.append(filename)
    
    if not found_files:
        raise FileNotFoundError(f"No required data files found in {output_folder_path}. "
                              f"Expected files: {required_files}")
    
    print(f"Found data files: {found_files}")
    return True


def get_data_info(output_folder_path):
    """
    Get information about available data files and their dimensions.
    
    Parameters:
    -----------
    output_folder_path : str
        Path to output folder
        
    Returns:
    --------
    dict
        Dictionary with file information
    """
    
    file_info = {}
    file_mappings = {
        'positions': {'file': 'positions_cen.dat', 'reader': pioud_extract.read_position_dat},
        'velocities': {'file': 'velocities_cen.dat', 'reader': pioud_extract.read_velocities_dat},
        'forces': {'file': 'forces_cen.dat', 'reader': pioud_extract.read_forces_dat},
        'stress': {'file': 'stress_cen.dat', 'reader': pioud_extract.read_stress_dat}
    }
    
    for property_name, mapping in file_mappings.items():
        filepath = os.path.join(output_folder_path, mapping['file'])
        
        if os.path.exists(filepath):
            try:
                # Read data to get shape information
                if property_name == 'positions':
                    data = mapping['reader'](filename=filepath, output_in_angs=True)
                elif property_name in ['forces', 'stress']:
                    data = mapping['reader'](filename=filepath, output_in_ev_per_angs=True)
                else:
                    data = mapping['reader'](filename=filepath)
                
                file_info[property_name] = {
                    'file': mapping['file'],
                    'exists': True,
                    'shape': data.shape,
                    'timesteps': data.shape[0],
                    'size_mb': os.path.getsize(filepath) / (1024*1024)
                }
                
            except Exception as e:
                file_info[property_name] = {
                    'file': mapping['file'],
                    'exists': True,
                    'error': str(e)
                }
        else:
            file_info[property_name] = {
                'file': mapping['file'],
                'exists': False
            }
    
    return file_info


def quick_time_average(output_folder_path, fraction_start=0.5, fraction_end=1.0):
    """
    Perform time averaging over a fraction of the total simulation time.
    
    Parameters:
    -----------
    output_folder_path : str
        Path to output folder
    fraction_start : float
        Starting fraction of simulation (0.0 to 1.0)
    fraction_end : float  
        Ending fraction of simulation (0.0 to 1.0)
        
    Returns:
    --------
    dict
        Time-averaged properties
    """
    
    # Get data info to determine total timesteps
    info = get_data_info(output_folder_path)
    
    # Find a file that exists to get timestep count
    timesteps = None
    for prop_info in info.values():
        if prop_info.get('exists', False) and 'timesteps' in prop_info:
            timesteps = prop_info['timesteps']
            break
    
    if timesteps is None:
        raise ValueError("Could not determine number of timesteps from data files")
    
    # Calculate row indices
    start_row = int(timesteps * fraction_start)
    end_row = int(timesteps * fraction_end)
    
    print(f"Quick time averaging over last {fraction_end-fraction_start:.1%} of simulation")
    print(f"Total timesteps: {timesteps}")
    print(f"Using rows {start_row} to {end_row-1}")
    
    return calculate_qe_time_averaged_properties(output_folder_path, start_row, end_row)