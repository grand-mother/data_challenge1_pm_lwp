import numpy as np
import math
import random
import electronic_chain.ec_config as ec_config
from XDU_electronic_chain import ifftgetn, fftgetn

def read_traces(tree, **kwargs):
    """Reads and combines traces from the tree
    Assumes that get_entry/get_event etc. was already executed on the tree"""

    # Name of tree instead of instance given
    if type(tree) is str:
        tree = kwargs[tree]

    # Read the traces
    ex, ey, ez = np.array(tree.trace_x), np.array(tree.trace_y), np.array(tree.trace_z)

    # Stack the traces
    traces_t = np.stack([ex, ey, ez], axis=-2)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"traces_t": traces_t, "trace_length": traces_t.shape[-1], "du_count": tree.du_count}
    # Outside pipeline return - raw values
    else:
        return traces_t

def read_sampling_time(tree, **kwargs):
    """Reads the sampling time from the tree
    Assumes that get_entry/get_event etc. was already executed on the tree"""

    # Name of tree instead of instance given
    if type(tree) is str:
        tree = kwargs[tree]

    # Read the traces
    sampling_time = tree.t_bin_size

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"sampling_time": sampling_time}
    # Outside pipeline return - raw values
    else:
        return sampling_time

def read_angles(tree, **kwargs):
    """Reads zenith and azimuth angles from the tree and converts them to the ZHAireS way
    Assumes that get_entry/get_event etc. was already executed on the tree"""

    # Name of tree instead of instance given
    if type(tree) is str:
        tree = kwargs[tree]

    # Read the traces
    e_theta = 180 - tree.shower_zenith
    e_phi = tree.shower_azimuth - 180
    if e_theta < 0: e_theta = 360 - e_theta
    if e_phi < 0: e_phi = 360 + e_phi

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"e_theta": e_theta, "e_phi": e_phi}
    # Outside pipeline return - raw values
    else:
        return e_theta, e_phi

# def adjust_traces(ex, ey, ez, Ts):
def adjust_traces(traces_t, sampling_time, **kwargs):
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

    original_trace_length = traces_t.shape[-1]

    # = == == == == == == == == == =Time-frequency parameter generation == == == == == == == == == == == == == == == == == == == == ==
    # = == == =In order to make the frequency resolution 1MHz, the number of sampling points = sampling frequency == == == =
    Ts = sampling_time
    fs = 1 / Ts * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    f0 = fs / N  # base frequency, Frequency resolution
    f = np.arange(0, N) * f0  # frequency sequence
    # Take only half, pay attention to odd and even numbers, odd numbers: the floor(N / 2 + 1) is conjugated to floor(N / 2 + 1) + 1;
    # Even number: floor(N / 2 + 1)-1 and floor(N / 2 + 1) + 1 are conjugated;

    ex, ey, ez = traces_t[...,0,:], traces_t[...,1,:], traces_t[...,2,:]

    # = == == == Change the original signal length to be the same as N======================
    ex_cut = np.zeros((*ex.shape[:-1],N))
    ey_cut = np.zeros((*ex.shape[:-1],N))
    ez_cut = np.zeros((*ex.shape[:-1],N))

    lx = ex.shape[-1]

    if N <= lx:
        # ============================In order to avoid not getting the peak value, judge whether the peak value is within N == == == == == == == =
        posx = np.argmax(ex)
        posy = np.argmax(ey)
        posz = np.argmax(ez)
        hang = max(posx, posy, posz)
        if hang >= N:
            ex_cut[..., 0: N - 500] = ex[..., hang - (N - 500): hang]
            ex_cut[..., N - 500: N] = ex[..., hang: hang + 500]

            ey_cut[..., 0: N - 500] = ey[..., hang - (N - 500): hang]
            ey_cut[..., N - 500: N] = ey[..., hang: hang + 500]

            ez_cut[..., 0: N - 500] = ez[..., hang - (N - 500): hang]
            ez_cut[..., N - 500: N] = ez[..., hang: hang + 500]
        else:
            ex_cut[..., 0: N] = ex[..., 0: N]
            ey_cut[..., 0: N] = ey[..., 0: N]
            ez_cut[..., 0: N] = ez[..., 0: N]
    else:
        ex_cut[..., 0: lx] = ex[..., 0:]
        ey_cut[..., 0: lx] = ey[..., 0:]
        ez_cut[..., 0: lx] = ez[..., 0:]

    traces_t = np.stack([ex_cut, ey_cut, ez_cut], axis=-2)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"traces_t": traces_t, "trace_length": traces_t.shape[-1], "original_trace_length": original_trace_length}
    # Outside pipeline return - raw values
    else:
        # return np.array(ex_cut), np.array(ey_cut), np.array(ez_cut)
        return traces_t

def multiply_traces(traces_f, multiplier, sampling_time = 0.5, **kwargs):
    """Multiply traces in frequency domain, return in time and frequency domain"""
    # Frequency part
    traces_f *= multiplier

    # Time part
    Ts = sampling_time
    fs = 1 / Ts * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    f0 = fs / N  # base frequency, Frequency resolution
    f = np.arange(0, N) * f0  # frequency sequence
    f1 = f[0:int(N / 2) + 1]

    # [V_t, _, _] = ifftget(v_new, N, f1, 2)
    [traces_t, _, _] = ifftgetn(traces_f, N, f1, 2)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"traces_t": traces_t, "traces_f": traces_f}
    # Outside pipeline return - raw values
    else:
        return traces_t, traces_f

def add_traces(traces_t, addend, traces_f = None, sampling_time = 0.5, **kwargs):
    """Add addend to traces"""
    """addend is either a time array, or a list [addend_time, addend_frequency]"""

    if type(addend) == list:
        addend_t, addend_f = addend
    else:
        addend_t = addend
        addend_f = None

    # Add the noise to data
    traces_t += addend_t

    # If addend in frequency domain was supplied, add it
    if addend_f is not None:
        traces_f += addend_f
    # Only time addend was supplied
    else:
        # Calculate the frequency trace with fft
        fs = 1 / sampling_time * 1000  # sampling frequency, MHZ
        N = math.ceil(fs)
        f0 = fs / N  # base frequency, Frequency resolution
        f = np.arange(0, N) * f0  # frequency sequence
        f1 = f[0:int(N / 2) + 1]

        [traces_f, _, _] = fftgetn(traces_t, N, f1)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"traces_t": traces_t, "traces_f": traces_f}
    # Outside pipeline return - raw values
    else:
        return traces_t, traces_f


def add_traces_randomized(traces_t, addend, traces_f = None, sampling_time = 0.5, **kwargs):
    """Add addend to traces, but randomly select/shuffle the addend"""
    """addend is either a time array, or a list [addend_time, addend_frequency]"""

    if type(addend) == list:
        addend_t, addend_f = addend
    else:
        addend_t = addend
        addend_f = None

    # For bulk of traces, permutation of the addend array
    if traces_t.shape == addend_t.shape:
        # Generate the randomised indices
        ind = np.random.permutation(np.arange(addend_t.shape[0]))
    # For single trace, select random value from the addend array
    else:
        ind = random.randint(a=0, b=addend_t.shape[0] - 1)

    # Add the noise to data
    traces_t += addend_t[ind]

    # If addend in frequency domain was supplied, add it
    if addend_f is not None:
        traces_f += addend_f[ind]
    # Only time addend was supplied
    else:
        # Calculate the frequency trace with fft
        fs = 1 / sampling_time * 1000  # sampling frequency, MHZ
        N = math.ceil(fs)
        f0 = fs / N  # base frequency, Frequency resolution
        f = np.arange(0, N) * f0  # frequency sequence
        f1 = f[0:int(N / 2) + 1]

        [traces_f, _, _] = fftgetn(traces_t, N, f1)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"traces_t": traces_t, "traces_f": traces_f}
    # Outside pipeline return - raw values
    else:
        return traces_t, traces_f

def restore_traces_length(traces_t, trace_length, original_trace_length, **kwargs):
    """Interpolate traces to the original length (which is changed in the XDU chain)"""

    # X coordinates of original_trace_length number of points through the (current) trace_length
    # x = np.tile(np.arange(original_trace_length)*(trace_length-1)/(original_trace_length-1), (traces_t.shape[0], 1))
    x = np.arange(original_trace_length)*(trace_length-1)/(original_trace_length-1)
    xp = np.arange(trace_length)

    # print(multiInterp2(x, xp, traces_t[:,0]).shape)

    # Slow, but works ;)
    if traces_t.ndim>2:
        tx = np.array([np.interp(x, xp, traces_t[i,0]) for i in range(traces_t.shape[0])])
        ty = np.array([np.interp(x, xp, traces_t[i,1]) for i in range(traces_t.shape[0])])
        tz = np.array([np.interp(x, xp, traces_t[i,2]) for i in range(traces_t.shape[0])])
    else:
        tx = np.array(np.interp(x, xp, traces_t[0]))
        ty = np.array(np.interp(x, xp, traces_t[0]))
        tz = np.array(np.interp(x, xp, traces_t[0]))

    traces_t = np.stack([tx, ty, tz], axis=-2)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"traces_t": traces_t}
    # Outside pipeline return - raw values
    else:
        return traces_t


def apply_noise_to_trace(V_t, V_f_complex, noise_model):
    noise_model_t, noise_model_f = noise_model

    # For bulk of traces, best to shuffle the noise array
    if V_t.shape == noise_model_t.shape:
        # Shuffle the noise arrays
        np.random.shuffle(noise_model_t)
        np.random.shuffle(noise_model_f)

        # Add the noise to data
        Voc_noise_t = V_t + noise_model_t
        Voc_noise_complex = V_f_complex + noise_model_f
    # For single trace, select random value from the noise array
    else:
        Voc_noise_t = V_t + noise_model_t[random.randint(a=0, b=noise_model_t.shape[0]-1)]
        Voc_noise_complex = V_t + noise_model_f[random.randint(a=0, b=noise_model_f.shape[0]-1)]

    return Voc_noise_t, Voc_noise_complex

def apply_noise_to_trace_old(V_t, V_f_complex, noise_model):
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
        # v_f_list.append(np.zeros_like(V_f_complex, dtype=complex))
        # v_f_list.append(np.zeros_like(V_f_complex))
        # v_prev = v_f_list[-2]
        # v_new = v_f_list[-1]
        # print(v_new.shape)
        # exit()
        # for p in range(V_f_complex.shape[1]):
        #     v_new[:, p] = part[:,p] * v_prev[:, p] + 0

        # v_new[:] = part * v_prev

        # Append the last voltage convolved with the new electronic chain element to the list of results
        v_f_list.append(v_f_list[-1] * part)

        if return_all or key == list(electronic_chain.keys())[-1]:

            # ToDo: This should be gotten from the trace parameters
            Ts = 0.5
            fs = 1 / Ts * 1000  # sampling frequency, MHZ
            N = math.ceil(fs)
            f0 = fs / N  # base frequency, Frequency resolution
            f = np.arange(0, N) * f0  # frequency sequence
            f1 = f[0:int(N / 2) + 1]

            # [V_t, _, _] = ifftget(v_new, N, f1, 2)
            [V_t, _, _] = ifftgetn(v_f_list[-1], N, f1, 2)
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
