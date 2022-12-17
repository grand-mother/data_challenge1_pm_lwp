import h5py
from electronic_chain.XDU_electronic_chain.functions import *
import electronic_chain.XDU_electronic_chain.config as config
import random

# =========================================galacticnoise get=============================================
def gala(lst, N, f0, f1, M):
    # This Python file uses the following encoding: utf-8


    # = == == == == This program is used as a subroutine to complete the calculation and expansion of galactic noise == == == == =
    #  ----------------------input - ---------------------------------- %
    # lst：Select the galactic noise LST at the LST moment
    # N is the extended length
    # f0 is the frequency resolution, f1 is the frequency point of the unilateral spectrum
    # % ----------------------output - ---------------------------------- %
    # v_complex_double, galactic_v_time

    GALAshowFile = config.XDU_files_path+"/30_250galactic.mat"
    GALAshow = h5py.File(GALAshowFile, 'r')
    GALApsd_dbm = np.transpose(GALAshow['psd_narrow_huatu'])
    GALApower_dbm = np.transpose(GALAshow['p_narrow_huatu'])
    GALAvoltage = np.transpose(GALAshow['v_amplitude'])
    GALApower_mag = np.transpose(GALAshow['p_narrow'])
    GALAfreq = GALAshow['freq_all']

    f_start = GALAfreq[0]
    f_num = len(GALAfreq)
    #f_end = GALAfreq[f_num]
    f_end = GALAfreq[f_num-1]


    v_complex_double = np.zeros((M, N, 3), dtype=complex)
    galactic_v_time = np.zeros((M, N, 3), dtype=float)
    galactic_v_m_single = np.zeros((M, int(N / 2) + 1, 3), dtype=float)
    galactic_v_p_single = np.zeros((M, int(N / 2) + 1, 3), dtype=float)
    unit_uv = 1e6
    V_amplitude = GALAvoltage[:, :, lst-1]

    aa = np.zeros((M, f_num, 3), dtype=float)
    phase = np.zeros((M, f_num, 3), dtype=float)
    v_complex = np.zeros((M, f_num, 3), dtype=complex)
    for mm in range(M): # antennas
        for ff in range(f_num): # freq
            for pp in range(3): # ports
                aa[mm, ff, pp] = np.random.normal(loc=0, scale=V_amplitude[ff, pp]) # Generates a normal distribution with 0 as the mean and V_amplitude[ff, pp] as the standard deviation
                phase[mm, ff, pp] = 2 * np.pi * random.random()  # phase of random Gauss noise
                v_complex[mm, ff, pp] = abs(aa[mm, ff, pp] * N / 2) * np.exp(1j * phase[mm, ff, pp])

    for kk in range(M):
        for port in range(3):
            [f, v_complex_double[kk, :, port]] = expan(N, f0, f_start, f_end, v_complex[kk, :, port])
            # print(v_complex_double[k, :, port])
        [galactic_v_time[kk], galactic_v_m_single[kk], galactic_v_p_single[kk]] = ifftget(v_complex_double[kk], N, f1, 2)

    return v_complex_double, galactic_v_time
