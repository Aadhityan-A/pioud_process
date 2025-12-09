import numpy as np
from ase.units import Bohr, Ha

"""
Based on
https://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation
"""

def autocorrelation(x, max_lag=None):
    """
    Calculate the autocorrelation of the input array.
    
    Parameters:
    -----------
    x : numpy.ndarray
        Input time series data
    max_lag : int, optional
        Maximum lag to calculate the autocorrelation for.
        Default is None, which calculates for all possible lags.
        
    Returns:
    --------
    acf : numpy.ndarray
        Autocorrelation function from lag 0 to max_lag
    """
    # Ensure the input is a numpy array and subtract the mean
    x = np.array(x)
    x_mean = np.mean(x)
    x = x - x_mean
    
    # Calculate the variance
    var = np.var(x)
    
    n = len(x)
    if max_lag is None:
        max_lag = n - 1
    else:
        max_lag = min(max_lag, n - 1)
    
    # Calculate the autocorrelation using numpy's correlate function
    acf = np.correlate(x, x, mode='full')[n-1:n+max_lag+1] / (n * var)
    
    return acf

if __name__ == "__main__":

    # Read the energy values
    # extracted_energy = read_energies(filename="/home/aadhityan/Documents/GitHub/pioud_soap/test_files_1/pimd.out")
    extracted_energy = read_energies(filename="/home/aadhityan/Documents/GitHub/pioud_process/test_files/pimd.out")
    print("extracted_energy",extracted_energy)

    # Calculate the autocorrelation of the energy values
    max_lag = 1000  # Calculate autocorrelation up to lag 100 for better visualization
    energy_acf = autocorrelation(extracted_energy, max_lag=max_lag)
    print("\nAutocorrelation of energy values:")
    print(energy_acf[:10])  # Print first 10 values for clarity

    #plot the autocorrelation
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 6))
        plt.plot(energy_acf)
        plt.title('Autocorrelation of Energy Values')
        plt.xlabel('Time/Lag')
        plt.ylabel('Autocorrelation')
        plt.grid(True)
        plt.savefig('energy_autocorrelation.png')
        plt.close()
        print("\nAutocorrelation plot saved as 'energy_autocorrelation.png'")
    except ImportError:
        print("\nMatplotlib not available. Skipping plotting.")

