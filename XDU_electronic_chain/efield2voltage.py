import electronic_chain.ec_config as ec_config
from XDU_electronic_chain.antenna_effective_length import CEL
from XDU_electronic_chain.functions import *

# def efield2voltage(ex, ey, ez, phi, theta, dt, antenna_model, antenna_model_calculation_function=CEL):
def efield2voltage(traces_t, sampling_time, Lce_complex, **kwargs):
    # ToDo: This should be gotten from the trace parameters
    Ts = sampling_time
    fs = 1 / Ts * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    f0 = fs / N  # base frequency, Frequency resolution
    f = np.arange(0, N) * f0  # frequency sequence
    f1 = f[0:int(N / 2) + 1]


    # [Lce_complex, antennas11_complex_short] = antenna_model_calculation_function(theta, phi, len(ex), f0, 1.0)
    # [Lce_complex, antennas11_complex_short] = antenna_model_calculation_function(e_theta, e_phi, traces_t.shape[1], f0, 1.0)
    # np.save("Lce_complex", Lce_complex)
    # np.save("antennas11_complex_short", antennas11_complex_short)
    # exit()
    # Lce_complex = np.load("Lce_complex.npy")

    # ======Open circuit voltage of air shower=================
    # Lcehang = Lce_complex.shape[0]
    # Lcelie = Lce_complex.shape[2]


    # Edata = np.stack([ex, ey, ez], axis=1)
    #
    # [E_shower_fft, E_shower_fft_m_single, E_shower_fft_p_single] = fftget(Edata, N, f1)  # Frequency domain signal

    # Edata1 = np.stack([ex, ey, ez], axis=-2)
    Edata1 = traces_t

    [E_shower_fft, E_shower_fft_m_single, E_shower_fft_p_single] = fftgetn(Edata1, N, f1)  # Frequency domain signal

    # E_shower_fft = np.moveaxis(E_shower_fft1, 0, 1)

    # Voc_shower_complex = np.zeros((Lcehang, Lcelie), dtype=complex)
    # # Frequency domain signal
    # for p in range(Lcelie):
    #     Voc_shower_complex[:, p] = Lce_complex[:, 0, p] * E_shower_fft[:, 0] + Lce_complex[:, 1, p] * E_shower_fft[:, 1] + Lce_complex[:, 2, p] * E_shower_fft[:, 2] + 0

    Lce_complex = np.moveaxis(Lce_complex, 0, 2)


    # Voc_shower_complex = Lce_complex[:, 0] * E_shower_fft[:, 0, np.newaxis] + Lce_complex[:, 1] * E_shower_fft[:, 1, np.newaxis] + Lce_complex[:, 2] * E_shower_fft[:, 2, np.newaxis]
    # print(E_shower_fft.shape, Lce_complex.shape)
    # print(Lce_complex[0].shape, E_shower_fft[..., 0, :].shape)
    # print(Lce_complex[0].shape, E_shower_fft[0, 0, :].shape, (Lce_complex[0]*E_shower_fft[0, 0, :]).shape)
    # a = Lce_complex[0] * E_shower_fft[..., 0, np.newaxis, :]
    # b = Lce_complex[0] * E_shower_fft[0, 0, :]
    # print(a[0][np.where(a[0]!=0)], b[np.where(b!=0)], np.all(a[0]==b))
    # Voc_shower_complex = Lce_complex[0] * E_shower_fft[..., 0, :] + Lce_complex[1] * E_shower_fft[..., 1, :] + Lce_complex[2] * E_shower_fft[..., 2, :]
    # Complicated broadcasting:
    # Lce_complex is [3,3,trace_length], E_shower_fft is [3,trace_length] in case of a single trace, or [number_of_traces, 3, trace_length]
    # Lce_complex[0] * E_shower_fft[..., 0, :] is (3, trace_length)*(trace_length)=(3, trace_length), can be done for 1 trace
    # for multiple traces it is (3, trace_length)*(number_of_traces,trace_length) which can't be mutiplied
    # However, Lce_complex[0] * E_shower_fft[..., 0, np.newaxis, :] is (3, trace_length)*(number_of_traces, 1, trace_length)
    # = (number_of_traces, 3, trace_length) - works both for single trace and for multiple traces
    Voc_shower_complex = Lce_complex[0] * E_shower_fft[..., 0, np.newaxis, :] + Lce_complex[1] * E_shower_fft[..., 1, np.newaxis, :] + Lce_complex[2] * E_shower_fft[..., 2, np.newaxis, :]

    # print(Voc_shower_complex[np.where(Voc_shower_complex!=0)])
    # print(Voc_shower_complex1[np.where(Voc_shower_complex1!=0)])
    # print(np.all(Voc_shower_complex==Voc_shower_complex1))
    # exit()

    [Voc_shower_t, Voc_shower_m_single, Voc_shower_p_single] = ifftgetn(Voc_shower_complex, N, f1, 2)

    # Voc_shower_m_single1 = np.moveaxis(Voc_shower_m_single1, 0, 1)

    # print("sss", Voc_shower_t1.shape)

    # Voc_shower_complex = np.moveaxis(Voc_shower_complex, 0, 1)

    # time domain signal
    # [Voc_shower_t, Voc_shower_m_single, Voc_shower_p_single] = ifftget(Voc_shower_complex, N, f1, 2)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"traces_t": Voc_shower_t, "traces_f": Voc_shower_complex}
    # Outside pipeline return - raw values
    else:
        return Voc_shower_t, Voc_shower_complex

def efield2voltage_old(ex, ey, ez, phi, theta, dt, antenna_model, antenna_model_calculation_function=CEL):
    # ToDo: This should be gotten from the trace parameters
    Ts = dt
    fs = 1 / Ts * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    f0 = fs / N  # base frequency, Frequency resolution
    f = np.arange(0, N) * f0  # frequency sequence
    f1 = f[0:int(N / 2) + 1]


    [Lce_complex, antennas11_complex_short] = antenna_model_calculation_function(theta, phi, len(ex), f0, 1.0)
    # np.save("Lce_complex", Lce_complex)
    # np.save("antennas11_complex_short", antennas11_complex_short)
    # exit()
    # Lce_complex1 = np.load("Lce_complex.npy")
    # print(Lce_complex[np.where(Lce_complex!=0)],Lce_complex1[np.where(Lce_complex!=0)])
    # print(np.all(Lce_complex==Lce_complex1))
    # exit()

    # ======Open circuit voltage of air shower=================
    Lcehang = Lce_complex.shape[0]
    Lcelie = Lce_complex.shape[2]


    Edata = np.stack([ex, ey, ez], axis=1)

    [E_shower_fft, E_shower_fft_m_single, E_shower_fft_p_single] = fftget(Edata, N, f1)  # Frequency domain signal

    Voc_shower_complex = np.zeros((Lcehang, Lcelie), dtype=complex)
    # # Frequency domain signal
    for p in range(Lcelie):
        Voc_shower_complex[:, p] = Lce_complex[:, 0, p] * E_shower_fft[:, 0] + Lce_complex[:, 1, p] * E_shower_fft[:, 1] + Lce_complex[:, 2, p] * E_shower_fft[:, 2] + 0


    # time domain signal
    [Voc_shower_t, Voc_shower_m_single, Voc_shower_p_single] = ifftget(Voc_shower_complex, N, f1, 2)

    return Voc_shower_t, Voc_shower_complex
