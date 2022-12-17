from functions import *
from scipy import interpolate
from numpy.ma import log10, abs

# ===========================================LNAparameter get===========================================
def LNA_get(antennas11_complex_short, N, f0, unit):
    # This Python file uses the following encoding: utf-8

    # from complex_expansion import expan

    # = == == == == This program is used as a subroutine to complete the calculation and expansion of the LNA partial pressure coefficient == == == == =
    #  ----------------------input - ---------------------------------- %
    # antennas11_complex_short is the program after the interpolation of the antenna standing wave test result
    # LNA address, s2p file (delete the previous string of s2p file in advance, put the test results of the three ports in the LNASparameter folder, and name them 1 2 3 in turn, corresponding to xyz)
    # The unit of test frequency is hz
    # N is the extended length
    # If unit is 0, the test data is in the form of real and imaginary parts, and 1 is in the form of db and phase.
    # f0 is the frequency resolution,
    # % ----------------------output - ---------------------------------- %
    # rou1 rou2 rou3

    z0 = 50 # characteristic impedance

    antenna_Gama_complex = np.zeros((N, 3), dtype=complex)
    for p in range(3):
        # Antenna related parameter calculation
        antennas11_short = antennas11_complex_short
        f0 = 1
        f_start = 30
        f_end = 250
        [f, antenna_Gama_complex[:, p]] = expan(N, f0, f_start, f_end, antennas11_short[:, p])

    Zin_antenna = z0 * (1 + antenna_Gama_complex) / (1 - antenna_Gama_complex)

    LNA_Gama_complex = np.zeros((N, 3), dtype=complex)  # 3 ports
    LNA_s21_complex = np.zeros((N, 3), dtype=complex)
    for p in range(3):
        #  LNA parameter
        str_p = str(p + 1)
        LNA_Address = ".//LNASparameter//" + str_p + ".s2p"
        freq = np.loadtxt(LNA_Address, usecols=0) / 1e6  # Hz to MHz
        if unit == 0:
            res11 = np.loadtxt(LNA_Address, usecols=1)
            ims11 = np.loadtxt(LNA_Address, usecols=2)
            res21 = np.loadtxt(LNA_Address, usecols=3)
            ims21 = np.loadtxt(LNA_Address, usecols=4)
            dbs21 = 20 * log10(abs(res21 + 1j * ims21))

        elif unit == 1:
            dbs11 = np.loadtxt(LNA_Address, usecols=1)
            degs11 = np.loadtxt(LNA_Address, usecols=2)
            mags11 = 10 ** (dbs11 / 20)
            res11 = mags11 * np.cos(degs11 / 180 * math.pi)
            ims11 = mags11 * np.sin(degs11 / 180 * math.pi)

            dbs21 = np.loadtxt(LNA_Address, usecols=3)
            degs21 = np.loadtxt(LNA_Address, usecols=4)
            mags21 = 10 ** (dbs21 / 20)
            res21 = mags21 * np.cos(degs21 / 180 * math.pi)
            ims21 = mags21 * np.sin(degs21 / 180 * math.pi)

        if p == 0:
            dBs21 = np.zeros((3, len(freq)))
        dBs21[p] = dbs21

        # Interpolate to (30:1:250)MHz
        freqnew = np.arange(30, 251, 1)
        f_res11 = interpolate.interp1d(freq, res11, kind="cubic")
        res11new = f_res11(freqnew)
        f_ims11 = interpolate.interp1d(freq, ims11, kind="cubic")
        ims11new = f_ims11(freqnew)
        s11_complex = res11new + 1j * ims11new
        [f, LNA_Gama_complex[:, p]] = expan(N, f0, f_start, f_end, s11_complex)

        f_res21 = interpolate.interp1d(freq, res21, kind="cubic")
        res21new = f_res21(freqnew)
        f_ims21 = interpolate.interp1d(freq, ims21, kind="cubic")
        ims21new = f_ims21(freqnew)
        s21_complex = res21new + 1j * ims21new
        [f, LNA_s21_complex[:, p]] = expan(N, f0, f_start, f_end, s21_complex)

    Zin_LNA = z0 * (1 + LNA_Gama_complex) / (1 - LNA_Gama_complex)

    # Partial pressure coefficient
    rou1 = Zin_LNA / (Zin_antenna + Zin_LNA)
    rou2 = (1 + LNA_Gama_complex) / (1 - antenna_Gama_complex * LNA_Gama_complex)
    rou3 = LNA_s21_complex / (1 + LNA_Gama_complex)

    return rou1, rou2, rou3
