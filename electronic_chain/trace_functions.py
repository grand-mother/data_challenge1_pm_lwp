import numpy as np
import math
import random
from electronic_chain.XDU_electronic_chain.functions import ifftget
from electronic_chain.XDU_electronic_chain.efield2voltage import efield2voltage

def adjust_traces(ex, ey, ez, Ts):
    """Adjust the traces length and positioning"""
    # This Python file uses the following encoding: utf-8

    # This program, as a subroutine, completes the reading of the shower time domain signal, ===========================================================
    # intercepts part of the signal or extends the length of the signal, and generates parameters according to the time domain requirements == == == == =
    #  ----------------------input - ---------------------------------- %
    # filename:path
    # Ts:time interval
    # % show_flag:flag of showing picture

    # % ----------------------output - ---------------------------------- %
    # % t:time sequence, unit:ns
    # % E_shower_cut:Corresponding to the triple polarization component of the time series，unit:uv
    # % fs % sampling frequency, unit:MHZ
    # % f0 % base frequency, Frequency resolution，unit:MHZ
    # % f  % frequency sequence，unit:MHz
    # % f1 % Unilateral spectrum frequency sequence，unit:MHz

    # = == == == == == == == == == =Time-frequency parameter generation == == == == == == == == == == == == == == == == == == == == ==
    # = == == =In order to make the frequency resolution 1MHz, the number of sampling points = sampling frequency == == == =
    fs = 1 / Ts * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    f0 = fs / N  # base frequency, Frequency resolution
    f = np.arange(0, N) * f0  # frequency sequence
    # Take only half, pay attention to odd and even numbers, odd numbers: the floor(N / 2 + 1) is conjugated to floor(N / 2 + 1) + 1;
    # Even number: floor(N / 2 + 1)-1 and floor(N / 2 + 1) + 1 are conjugated;

    # = == == == Change the original signal length to be the same as N======================
    ex_cut = np.zeros((N))
    ey_cut = np.zeros((N))
    ez_cut = np.zeros((N))

    lx = len(ex)
    if N <= lx:
        # ============================In order to avoid not getting the peak value, judge whether the peak value is within N == == == == == == == =
        posx = np.argmax(ex)
        posy = np.argmax(ey)
        posz = np.argmax(ez)
        hang = max(posx, posy, posz)
        if hang >= N:
            ex_cut[0: N - 500] = ex[hang - (N - 500): hang]
            ex_cut[N - 500: N] = ex[hang: hang + 500]

            ey_cut[0: N - 500] = ey[hang - (N - 500): hang]
            ey_cut[N - 500: N] = ey[hang: hang + 500]

            ez_cut[0: N - 500] = ez[hang - (N - 500): hang]
            ez_cut[N - 500: N] = ez[hang: hang + 500]
        else:
            ex_cut[0: N] = ex[0: N]
            ey_cut[0: N] = ey[0: N]
            ez_cut[0: N] = ez[0: N]
    else:
        ex_cut[0: lx] = ex[0:]
        ey_cut[0: lx] = ey[0:]
        ez_cut[0: lx] = ez[0:]

    return np.array(ex_cut), np.array(ey_cut), np.array(ez_cut)

def apply_noise_to_trace(V_t, V_f_complex, noise_model):
    noise_model_t, noise_model_f = noise_model
    # ======Galaxy noise power spectrum density, power, etc.=====================
    # ===========Voltage with added noise=======================================
    Voc_noise_t = np.zeros_like(V_t)
    Voc_noise_complex = np.zeros_like(V_t, dtype=complex)
    for p in range(Voc_noise_complex.shape[1]):
        Voc_noise_t[:, p] = V_t[:, p] + noise_model_t[random.randint(a=0, b=noise_model_t.shape[0]-1), :, p]
        Voc_noise_complex[:, p] = V_f_complex[:, p] + noise_model_f[random.randint(a=0, b=noise_model_f.shape[0]-1), :, p]

    # [Voc_noise_t_ifft, Voc_noise_m_single, Voc_noise_p_single] = ifftget(Voc_noise_complex, N, f1, 2)
    return Voc_noise_t, Voc_noise_complex

def apply_electronic_chain_to_trace(V_f_complex, electronic_chain, return_all=False, convert_to_adc=True):
    v_f_list = [V_f_complex]
    v_t_list = [None]
    adc_list = []
    # Loop through all the parts of the electronic chain dictionary
    for (key,part) in electronic_chain.items():
        v_f_list.append(np.zeros_like(V_f_complex, dtype=complex))
        # v_f_list.append(np.zeros_like(V_f_complex))
        v_prev = v_f_list[-2]
        v_new = v_f_list[-1]
        for p in range(V_f_complex.shape[1]):
            v_new[:, p] = part[:,p] * v_prev[:, p] + 0
        if return_all or key == list(electronic_chain.keys())[-1]:

            # ToDo: This should be gotten from the trace parameters
            Ts = 0.5
            fs = 1 / Ts * 1000  # sampling frequency, MHZ
            N = math.ceil(fs)
            f0 = fs / N  # base frequency, Frequency resolution
            f = np.arange(0, N) * f0  # frequency sequence
            f1 = f[0:int(N / 2) + 1]

            [V_t, _, _] = ifftget(v_new, N, f1, 2)
            v_t_list.append(V_t)

            # ====================Voltage after ADC======================================
            if convert_to_adc:
                Length_AD = 14  # Significant bit of the value, plus a sign bit in addition to this
                Range_AD = 1.8 * 1e6  # Vpp,unit:uv
                delta = Range_AD / 2 / (2 ** (Length_AD - 1))
                adc_list.append(np.sign(v_t_list[-1]) * np.floor(abs(v_t_list[-1]) / delta) * delta)

    if return_all:
        return v_t_list, v_f_list, adc_list
    else:
        return v_t_list[-1], v_f_list[-1], adc_list[-1]


# The main function defining the Efield -> ADC conversion for a single trace
def efield_to_adc(ex, ey, ez, phi, theta, antenna_model, noise_model, electronic_chain, efield2voltage_func=efield2voltage, return_voltages=False):
    # Convert Efield to Voltage
    Voc_shower_t, Voc_shower_complex = efield2voltage_func(ex, ey, ez, phi, theta, 0.5, antenna_model)
    # Apply noise to Voltage
    Voc_noise_t, Voc_noise_complex = apply_noise_to_trace(Voc_shower_t, Voc_shower_complex, noise_model)
    # Apply electronic chain to Voltage, resulting in ADC counts
    v_t, v_f, adc = apply_electronic_chain_to_trace(Voc_shower_complex, electronic_chain, return_all=True)
    adc = adc[-1]

    if return_voltages:
        return adc[:,0], adc[:,1], adc[:,2], [Voc_shower_t, Voc_noise_t, *v_t], [Voc_shower_complex, Voc_noise_complex, *v_f]
    else:
        return adc[:,0], adc[:,1], adc[:,2]
