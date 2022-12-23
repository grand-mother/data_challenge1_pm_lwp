import numpy as np
from XDU_electronic_chain.efield2voltage import efield2voltage
from XDU_electronic_chain.antenna_effective_length import CEL
from XDU_electronic_chain.galactic_noise import gala
from XDU_electronic_chain.LNA import LNA_get
from XDU_electronic_chain.filters import filter_get
from XDU_electronic_chain.voltage2adc import voltage2adc
import electronic_chain.trace_functions as ectf
from misc_functions import time_passed

# The main function defining the Efield -> ADC conversion pipeline
def manual_pipeline(traces_t, phi, theta, du_count, original_trace_length, sampling_time=0.5, return_all_traces=True):
    # !!!!!!!!!!!!!!!!!!!!!!! Make antenna model pluggable!!!!!!!!!

    all_traces_t = []
    all_traces_f = []

    time_passed(True)

    # Compute the antenna effective length
    Lce_complex, antennas11_complex_short = CEL(theta, phi, traces_t.shape[-1], sampling_time, unit=1)

    # Convert Efield to Voltage
    traces_t, traces_f = efield2voltage(traces_t, sampling_time, Lce_complex)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))
        all_traces_f.append(np.copy(traces_f))

    # Generate galactic noise
    galactic_v_time, v_complex_double = gala(du_count, sampling_rate=sampling_time, lst=18)

    # Apply galactic noise
    traces_t, traces_f = ectf.add_traces_randomized(traces_t, [galactic_v_time, v_complex_double], traces_f, sampling_time)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))
        all_traces_f.append(np.copy(traces_f))

    # Read LNA coefficient
    # ToDo: Not optimal for trace by trace, as it reads the LNA for every trace
    LNA_coefficient = LNA_get(antennas11_complex_short, sampling_time)
    # Apply LNA coefficient
    traces_t, traces_f = ectf.multiply_traces(traces_f, LNA_coefficient, sampling_time)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))
        all_traces_f.append(np.copy(traces_f))

    # Read cable and filter coefficient
    cable_coefficient, filter_coefficient = filter_get(sampling_time)

    # Apply cable coefficient
    traces_t, traces_f = ectf.multiply_traces(traces_f, cable_coefficient, sampling_time)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))
        all_traces_f.append(np.copy(traces_f))

    # Apply filter coefficient
    traces_t, traces_f = ectf.multiply_traces(traces_f, filter_coefficient, sampling_time)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))
        all_traces_f.append(np.copy(traces_f))

    # Restore length of traces to the original one (XDU chain expands them to 2000)
    traces_t = ectf.restore_traces_length(traces_t, traces_t.shape[-1], original_trace_length)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))

    # Convert voltage to ADC counts
    traces_t = voltage2adc(traces_t)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))

    if return_all_traces:
        return all_traces_t, all_traces_f
    else:
        return traces_t
