import numpy as np
from PM_functions.efield2voltage import efield2voltage_pm
from XDU_electronic_chain.voltage2adc import voltage2adc
from electronic_chain.trace_functions import adjust_traces
import electronic_chain.trace_functions as ectf
import math
from misc_functions import time_passed

# The main function defining the Efield -> ADC conversion pipeline
def manual_pipeline(traces_t, phi, theta, du_count, original_trace_length, sampling_time=0.5, galactic_noise_coefficient=None, LNA_coefficient=None, cable_coefficient=None, filter_coefficient=None, return_all_traces=True):
    all_traces_t = []
    all_traces_f = []

    time_passed(True)

    # Adjust the traces lenghts - some pipelines need it
    traces_t = adjust_traces(traces_t, sampling_time)


    # Stuff needed for converting efield to voltage
    fs = 1 / sampling_time * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    freqs_tot = np.fft.rfftfreq(N, d=sampling_time * 1e-9) / 1e6  # MHz
    # Convert Efield to Voltage
    traces_t, traces_f = efield2voltage_pm(traces_t, theta, phi, freqs_tot, sampling_time)
    # If storing the traces from all stages is requested, append the traces to list of traces
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))
        all_traces_f.append(np.copy(traces_f))

    # Apply galactic noise
    if galactic_noise_coefficient is not None:
        traces_t, traces_f = ectf.add_traces_randomized(traces_t, galactic_noise_coefficient, traces_f, sampling_time)
        if return_all_traces:
            all_traces_t.append(np.copy(traces_t))
            all_traces_f.append(np.copy(traces_f))

    # Apply LNA coefficient
    if LNA_coefficient is not None:
        traces_t, traces_f = ectf.multiply_traces(traces_f, LNA_coefficient, sampling_time, irfft=True)
        if return_all_traces:
            all_traces_t.append(np.copy(traces_t))
            all_traces_f.append(np.copy(traces_f))

    # Apply cable coefficient
    if cable_coefficient is not None:
        traces_t, traces_f = ectf.multiply_traces(traces_f, cable_coefficient, sampling_time, irfft=True)
        if return_all_traces:
            all_traces_t.append(np.copy(traces_t))
            all_traces_f.append(np.copy(traces_f))

    # Apply filter coefficient
    if filter_coefficient is not None:
        traces_t, traces_f = ectf.multiply_traces(traces_f, filter_coefficient, sampling_time, irfft=True)
        if return_all_traces:
            all_traces_t.append(np.copy(traces_t))
            all_traces_f.append(np.copy(traces_f))

    # Restore length of traces to the original one
    traces_t = ectf.restore_traces_length(traces_t, traces_t.shape[-1], original_trace_length)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))

    # Convert voltage to ADC counts
    traces_t = voltage2adc(traces_t)
    if return_all_traces:
        all_traces_t.append(np.copy(traces_t))

    # If requested, return the traces from all stages, both in time and frequency
    if return_all_traces:
        return all_traces_t, all_traces_f
    # By default, return only the final traces in time
    else:
        return traces_t
