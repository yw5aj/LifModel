import numpy as np
import matplotlib.pyplot as plt
import ctypes

# Load dll for LIF model
_get_spike_trace_array_lesniak = ctypes.cdll.LoadLibrary('./lesniak_dll.dll'
    ).get_spike_trace_array_lesniak
_get_spike_trace_array_lesniak.argtypes = [ctypes.c_double, ctypes.c_double,
    ctypes.c_double, np.ctypeslib.ndpointer(np.uintp, ndim=1), 
    ctypes.c_int, ctypes.c_double, np.ctypeslib.ndpointer(ctypes.c_int),
    ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int)]
_get_spike_trace_array_lesniak.restype = None

# SI units
DT = 1E-5 # sec
RESISTANCE_LIF = 10 * 5e8 # ohm
CAPACITANCE_LIF = 1e-11 # F
VOLTAGE_THRESHOLD_LIF = 30e-3 # V

def current_array_to_spike_array(current_array, mcnc_grouping):
    # Initialize output array
    spike_array = np.zeros(current_array.shape[0], dtype=np.int)
    # Call C function
    current_array_pp = (current_array.__array_interface__['data'][0]\
        + np.arange(current_array.shape[0])*current_array.strides[0]
        ).astype(np.uintp)
    _get_spike_trace_array_lesniak(RESISTANCE_LIF, CAPACITANCE_LIF,
        VOLTAGE_THRESHOLD_LIF, current_array_pp, current_array.shape[0], DT, 
        spike_array, mcnc_grouping.size, mcnc_grouping)
    return spike_array

if __name__ == '__main__':
    mcnc_grouping = np.array([8, 5, 3, 1])
    current_array = np.loadtxt('../test/csvs/trans_current.csv', delimiter=','
        )[:, 1:].copy()
    a = current_array_to_spike_array(current_array, mcnc_grouping)