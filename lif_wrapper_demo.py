import ctypes
import numpy as np
import matplotlib.pyplot as plt

def get_spike_trace_array_lif(resistance_lif, capacitance_lif, voltage_threshold_lif, current_array_lif, dt_lif):
    if not current_array_lif.flags['C_CONTIGUOUS']:
        current_array_lif = current_array_lif.copy(order='C')
    _get_spike_trace_array_lif = ctypes.cdll.LoadLibrary('./lif_dll.dll').get_spike_trace_array_lif
    _get_spike_trace_array_lif.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, np.ctypeslib.ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_double, np.ctypeslib.ndpointer(ctypes.c_int)]
    _get_spike_trace_array_lif.restype = None
    current_length_lif = current_array_lif.shape[0]
    spike_trace_array_lif = np.zeros(current_length_lif, dtype=np.int)
    _get_spike_trace_array_lif(resistance_lif, capacitance_lif, voltage_threshold_lif, current_array_lif, current_length_lif, dt_lif, spike_trace_array_lif)
    return spike_trace_array_lif


if __name__ == '__main__':
    # Read the transduction current data
    time_array_lif, current_array_lif = np.genfromtxt('../test/csvs/test_current.csv', delimiter=',').T
    time_array_lif *= 1e-3 # Convert ms to sec
    current_array_lif *= 1e-3 # Convert mA to A
    # Input LIF parameters
    resistance_lif = 5e8 # ohm
    capacitance_lif = 1e-11 # F
    voltage_threshold_lif = 29e-3 # V
    dt_lif = 1e-4 # sec
    # Calculate spike trace
    spike_trace_array_lif = get_spike_trace_array_lif(resistance_lif, capacitance_lif, voltage_threshold_lif, current_array_lif, dt_lif)
    # Plot the result
    fig, axs = plt.subplots()    
    axs.plot(spike_trace_array_lif, '-')
    axs.set_ylim(-2, 2)