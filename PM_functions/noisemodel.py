"""
Handels noise models: galactic and a gaussian noise which can be used to mimic electronic noise
"""

import numpy as np
import math
import h5py
from PM_functions.functions import *
import PM_functions.config as PM_config
from functools import lru_cache
import electronic_chain.ec_config as ec_config


@lru_cache(maxsize=1)
def read_galactic_noise(GALAshowFile):
    GALAshow = h5py.File(GALAshowFile, 'r')
    GALAvoltage = np.transpose(GALAshow['v_amplitude'])
    GALAfreq = GALAshow['freq_all']

    return GALAvoltage, GALAfreq

def generate_galacticnoise(du_count, sampling_time=0.5, lst=18, **kwargs):
    trace_len = kwargs["traces_t"].shape[-1]
    return real_generate_galacticnoise(du_count, trace_len, sampling_time, lst)

@lru_cache(maxsize=1)
def real_generate_galacticnoise(du_count, trace_len, sampling_time=0.5, lst=18):
    """
    Calculates galactic noise in Volatge for 3 antenna arms / antenna
    
    :param no_an: int, number of antennas
    :param hr_lst: int, lst hour e.g lst =18
    :param sampling_time: float, time binning of the data
    :param N: int, length of time trace array
    :return: galactic noise in frequency and time domain for each antenna
    """

    GALAshowFile = PM_config.PM_files_path + "/30_250galactic.mat"

    GALAvoltage, GALAfreq = read_galactic_noise(GALAshowFile)
    f_start = GALAfreq[0]
    f_num = len(GALAfreq)
    f_end = GALAfreq[f_num-1]
    V_amplitude = GALAvoltage[:, :, lst-1]
    GALAfreq = np.array(GALAfreq).flatten()

    fs = 1 / sampling_time * 1000  # sampling frequency, MHZ
    # N = math.ceil(fs)
    N = trace_len
    freqs = np.fft.rfftfreq(N, d=sampling_time * 1e-9) / 1e6  # MHz

    amp_interp = get_interpolation(freqs, GALAfreq, V_amplitude.T, f_start, f_end) #(3,n)

    ramp = np.random.normal(loc=0, scale=amp_interp[np.newaxis,...], size=(du_count,*amp_interp.shape))
    rphi = 2 * np.pi * np.random.random(size=(du_count,*amp_interp.shape))

    V_interp_fr = 2 * len(freqs) * ramp * np.exp(-1j * rphi)
    V_interp_t = np.fft.irfft(V_interp_fr, n=N)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"galactic_noise": [V_interp_t, V_interp_fr]}
    # Outside pipeline return - raw values
    else:
        return V_interp_t, V_interp_fr
