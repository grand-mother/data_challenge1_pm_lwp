import PM_functions.readantennamodel as an

# def efield2voltage_pm(ex, exy, ez, Zenith, Azimuth, N, dt, freqs):
def efield2voltage_pm(ex, ey, ez, Zenith, Azimuth, dt, antenna_model):
    """""Voltage calculation from E field traces - by Pragati Mitra
    Parameters:
        Etrace_cut : time trace of E field : expected shape (3,N)
        Zenith: shower zenith in degree
        Azimuth: shower azimuth in degree
        N: chopped length of the trace with frequency resolution 1 MHz (defin code- fs = 1 / dt * 1000  # sampling frequency, MHZ
            N = math.ceil(fs))
        freqs: frequencies (real part)
    Returns:
        Voltage trace (time domain)
    """

    fs = 1 / dt * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    freqs = np.fft.rfftfreq(N, d=dt * 1e-9) / 1e6

    ## read antenna response function Leff from files in theta, phi direction

    print("0.1")
    table_ewarm_new = an.get_tabulated('./antennamodel/GP300Antenna_EWarm_leff.npy')
    print("0.2")
    table_snarm_new = an.get_tabulated('./antennamodel/GP300Antenna_SNarm_leff.npy')
    print("0.3")
    table_zarm_new = an.get_tabulated('./antennamodel/GP300Antenna_Zarm_leff.npy')
    print("1")
    # interpolated L_eff in 3 arms for given zenith,azimuth

    lt1, lp1 = an.get_interp(table_ewarm_new, Zenith, Azimuth, N, dt * 1e-9)
    lt2, lp2 = an.get_interp(table_snarm_new, Zenith, Azimuth, N, dt * 1e-9)
    lt3, lp3 = an.get_interp(table_zarm_new, Zenith, Azimuth, N, dt * 1e-9)
    lt = np.array([lt1, lt2, lt3]).T
    lp = np.array([lp1, lp2, lp3]).T
    print("2")
    def rotationmatrix(theta, phi):
        """
        Rotation Matrix
        """

        t, p = np.deg2rad(theta), np.deg2rad(phi)
        ct, st = np.cos(t), np.sin(t)
        cp, sp = np.cos(p), np.sin(p)
        rot = np.zeros([3, 3])
        rot[0, 0], rot[0, 1], rot[0, 2] = ct * cp, -sp, st * cp
        rot[1, 0], rot[1, 1], rot[1, 2] = ct * sp, cp, st * sp
        rot[2, 0], rot[2, 1], rot[2, 2] = -st, 0, ct
        return rot

    def xyz_thetaphi(vec_in, theta, phi):
        """ cartesian x,yz to theta, phi coordinate
        """

        rot = rotationmatrix(theta, phi)
        rot_in = rot.T  # transpose is inverse
        vec_out = np.matmul(rot_in, vec_in.T)
        return vec_out[:2]  # not interested in E_r, neglibible

    Etrace_cut = np.array([ex, ey, ez])

    # Fold antenna response in and get Voc

    #################   elctric field trace theta phi and FFT ###########################################
    E_tp_t = np.zeros([N, 2])  # xyz to theta-phi

    print("3")
    for k in range(N):
        E_tp_t[k] = xyz_thetaphi(Etrace_cut[:, k], Zenith, Azimuth)

    E_tp_fft = np.array(np.fft.rfft(E_tp_t.T), dtype='complex')  # fft , shape(2,N)
    # voltage calculation
    # ======Open circuit voltage of air shower=================

    Voc_shower_complex = np.zeros([len(freqs), 3], dtype=complex)

    # Frequency domain signal after folding antenna response
    for p in range(3):
        Voc_shower_complex[:, p] = lt[:, p] * E_tp_fft[0] + lp[:, p] * E_tp_fft[1]
        Voc_shower_t = np.fft.irfft(Voc_shower_complex.T, n=N)  # (3,N)

    return np.moveaxis(Voc_shower_t, 0, 1), Voc_shower_complex
