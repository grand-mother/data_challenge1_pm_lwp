import numpy as np
import math
from scipy import interpolate
import h5py
from PM_functions.functions import *
import PM_functions.config as PM_config
from functools import lru_cache
import electronic_chain.ec_config as ec_config

"""
 handels noise models: galactic and a gaussian noise which can be used to mimic elecronic noise
 """

@lru_cache(maxsize=1)
def read_galactic_noise(GALAshowFile):
    GALAshow = h5py.File(GALAshowFile, 'r')
    # GALApsd_dbm = np.transpose(GALAshow['psd_narrow_huatu'])
    # GALApower_dbm = np.transpose(GALAshow['p_narrow_huatu'])
    GALAvoltage = np.transpose(GALAshow['v_amplitude'])
    # GALAvoltage = np.moveaxis(np.transpose(GALAshow['v_amplitude']), 1, 2)
    # GALApower_mag = np.transpose(GALAshow['p_narrow'])
    GALAfreq = GALAshow['freq_all']

    return GALAvoltage, GALAfreq

def generate_galacticnoise(du_count, sampling_time=0.5, lst=18, **kwargs):
    return real_generate_galacticnoise(du_count, sampling_time, lst)

@lru_cache(maxsize=1)
def real_generate_galacticnoise(du_count, sampling_time=0.5, lst=18):
    """
    Calculates galactic noise in Volatge for 3 antenna arms / antenna
    
    :param no_an: int, number of antennas
    :param hr_lst: int, lst hour e.g lst =18
    :param sampling_time: float, time binning of the data
    :param N: int, length of time trace array
    :return: galactic noise in frequency and time domain for each antenna
    """


    GALAshowFile = PM_config.PM_files_path + "/30_250galactic.mat"
    # GALAshow = h5py.File(GALAshowFile, 'r')
    # GALApsd_dbm = np.transpose(GALAshow['psd_narrow_huatu'])
    # GALApower_dbm = np.transpose(GALAshow['p_narrow_huatu'])
    # GALAvoltage = np.transpose(GALAshow['v_amplitude'])
    # GALApower_mag = np.transpose(GALAshow['p_narrow'])
    # GALAfreq = GALAshow['freq_all']

    GALAvoltage, GALAfreq = read_galactic_noise(GALAshowFile)

    f_start = GALAfreq[0]
    f_num = len(GALAfreq)
    #f_end = GALAfreq[f_num]
    f_end = GALAfreq[f_num-1]
    V_amplitude = GALAvoltage[:, :, lst-1]
    GALAfreq = np.array(GALAfreq).flatten()
    R = 50
    unit_uv = 1e6
    #V_amplitude = 2 * np.sqrt(GALApower_mag * R) * unit_uv

    fs = 1 / sampling_time * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    freqs = np.fft.rfftfreq(N, d=sampling_time * 1e-9) / 1e6  # MHz

    V_interp_fr = np.zeros([du_count,3,len(freqs)],dtype ='complex')
    V_interp_t = np.zeros([du_count,3,N])

    amp_interp = get_interpolation(freqs, GALAfreq, V_amplitude.T, f_start, f_end) #(3,n)


    for i in np.arange(du_count):
        for j in np.arange(3):
            for k, val in enumerate(freqs):
            
                randomamp = np.random.normal(loc=0, scale=amp_interp[j,k])
                phi = 2 * np.pi * np.random.random()
                V_interp_fr[i,j,k] = 2 * len(freqs) * randomamp * np.exp(-1j * phi) #2freq_len normalization assuming np.fft,double check!

                
        
        V_interp_t[i] = np.fft.irfft(V_interp_fr[i],n=N)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"galactic_noise": [V_interp_t, V_interp_fr]}
    # Outside pipeline return - raw values
    else:
        return V_interp_t, V_interp_fr


def add_gaussiannoise(no_an,N, vrms):
    """ generate normal random noise
    Parameters:
    -----------
        no_an: int
                number of antennas
        N: int
            length of the time trace
        vrms: float
            noise rms, e.g vrmsnoise= 15 meu-V, get from hardware measurement
    Returns:
    ----------
        numpy array
        noise  trace in time domain

    """
    noisetrace_t = np.zeros([no_an,3,N])
    for i in np.arange(no_an):
        for j in np.arange(3):
            noisetrace_t[i,j] = np.random.normal(0, vrms, size=N)


    return noisetrace_t
