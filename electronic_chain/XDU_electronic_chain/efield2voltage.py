from electronic_chain.XDU_electronic_chain.antenna_effective_length import CEL
from electronic_chain.XDU_electronic_chain.functions import *

def efield2voltage(ex, ey, ez, phi, theta, dt, antenna_model, antenna_model_calculation_function=CEL):
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
    # Lce_complex = np.load("Lce_complex.npy")

    # ======Open circuit voltage of air shower=================
    Lcehang = Lce_complex.shape[0]
    Lcelie = Lce_complex.shape[2]


    Edata = np.stack([ex, ey, ez], axis=1)
    # print(N, f1, Edata.shape)
    # exit()

    [E_shower_fft, E_shower_fft_m_single, E_shower_fft_p_single] = fftget(Edata, N, f1)  # Frequency domain signal

    Voc_shower_complex = np.zeros((Lcehang, Lcelie), dtype=complex)
    # Frequency domain signal
    for p in range(Lcelie):
        Voc_shower_complex[:, p] = Lce_complex[:, 0, p] * E_shower_fft[:, 0] + Lce_complex[:, 1, p] * E_shower_fft[:, 1] + Lce_complex[:, 2, p] * E_shower_fft[:, 2] + 0
    # time domain signal
    [Voc_shower_t, Voc_shower_m_single, Voc_shower_p_single] = ifftget(Voc_shower_complex, N, f1, 2)

    return Voc_shower_t, Voc_shower_complex
