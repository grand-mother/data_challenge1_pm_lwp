from functions import *
from scipy import interpolate

# ===============================================Filterparameter get====================================
def filter_get(N, f0, unit):
    # This Python file uses the following encoding: utf-8


    # = == == == == This program is used as a subroutine to complete the calculation and expansion of the S parameters of the cable and filter == == == == =
    #  ----------------------input - ---------------------------------- %
    # filter :path of s2p file
    # The filter test data is in the Filterparameter folder (delete the previous string of the s2p file in advance and name it 1, because the three port filters are the same)
    # The cable test data is in the cableparameter folder (delete the previous string of the s2p file in advance and name it cable because the three ports are the same)
    # The unit of test frequency is hz
    # N is the extended length
    # If unit is 0, the test data is in the form of real and imaginary parts, and 1 is in the form of db and phase.
    # f0 is the frequency resolution,
    # % ----------------------output - ---------------------------------- %
    # cable_coefficient
    # filter_coefficient

    Gain_VGA = -1.5  # dB
    r_balun = 630 * 2 / 650
    # test filter without VGA
    # Gain_VGA = 0  # dB
    # r_balun = 1
    cable_Gama_complex = np.zeros((N, 3), dtype=complex)  # 3 ports
    cable_s21_complex = np.zeros((N, 3), dtype=complex)
    for p in range(3):
        #  cable参数
        # str_p = str(p + 1)
        # cable_Address = ".//data//cableparameter//" + str_p + ".s2p"
        cable_Address = ".//cableparameter//cable.s2p"
        freq = np.loadtxt(cable_Address, usecols=0) / 1e6  # HZ to MHz
        if unit == 0:
            res11 = np.loadtxt(cable_Address, usecols=1)
            ims11 = np.loadtxt(cable_Address, usecols=2)
            res21 = np.loadtxt(cable_Address, usecols=3)
            ims21 = np.loadtxt(cable_Address, usecols=4)
            dbs21 = 20 * log10(abs(res21 + 1j * ims21))
            dbs11 = 20 * log10(abs(res11 + 1j * ims11))

        elif unit == 1:
            dbs11 = np.loadtxt(cable_Address, usecols=1)
            degs11 = np.loadtxt(cable_Address, usecols=2)
            mags11 = 10 ** (dbs11 / 20)
            res11 = mags11 * np.cos(degs11 / 180 * math.pi)
            ims11 = mags11 * np.sin(degs11 / 180 * math.pi)

            dbs21 = np.loadtxt(cable_Address, usecols=3)
            degs21 = np.loadtxt(cable_Address, usecols=4)
            mags21 = 10 ** (dbs21 / 20)
            res21 = mags21 * np.cos(degs21 / 180 * math.pi)
            ims21 = mags21 * np.sin(degs21 / 180 * math.pi)

        if p == 0:
            dBs21_cable = np.zeros((3, len(freq)))
            dBs11_cable = np.zeros((3, len(freq)))
        dBs21_cable[p] = dbs21
        dBs11_cable[p] = dbs11

        f_start = 30
        f_end = 250
        #  Interpolate to (30:1:250)MHz
        freqnew = np.arange(30, 251, 1)
        f_res11 = interpolate.interp1d(freq, res11, kind="cubic")
        res11new = f_res11(freqnew)
        f_ims11 = interpolate.interp1d(freq, ims11, kind="cubic")
        ims11new = f_ims11(freqnew)
        s11_complex = res11new + 1j * ims11new
        [f, cable_Gama_complex[:, p]] = expan(N, f0, f_start, f_end, s11_complex)

        f_res21 = interpolate.interp1d(freq, res21, kind="cubic")
        res21new = f_res21(freqnew)
        f_ims21 = interpolate.interp1d(freq, ims21, kind="cubic")
        ims21new = f_ims21(freqnew)
        s21_complex = res21new + 1j * ims21new
        [f, cable_s21_complex[:, p]] = expan(N, f0, f_start, f_end, s21_complex)

    cable_coefficient = (1 + cable_Gama_complex) * cable_s21_complex

    #  =================filter=========================================
    filter_Gama_complex = np.zeros((N, 3), dtype=complex)  # 3 ports
    filter_s21_complex = np.zeros((N, 3), dtype=complex)
    for p in range(3):
        #  filter parameter
        # str_p = str(p + 1)
        # filter_Address = ".//data//filterparameter//" + str_p + ".s2p"
        filter_Address = ".//filterparameter//1.s2p"
        freq = np.loadtxt(filter_Address, usecols=0) / 1e6  # HZ to MHz
        if unit == 0:
            res11 = np.loadtxt(filter_Address, usecols=1)
            ims11 = np.loadtxt(filter_Address, usecols=2)
            res21_filter = np.loadtxt(filter_Address, usecols=3)
            ims21_filter = np.loadtxt(filter_Address, usecols=4)
            dbs21 = 20 * log10(abs(res21 + 1j * ims21))
            dbs11 = 20 * log10(abs(res11 + 1j * ims11))
            dbs21_add_VGA = dbs21 + Gain_VGA + 20 * log10(r_balun)
            mags21 = 10 ** (dbs21_add_VGA / 20)
            degs21 = np.angle(res21_filter + ims21_filter)  # Phase radians
            res21 = mags21 * np.cos(degs21)
            ims21 = mags21 * np.sin(degs21)
        elif unit == 1:
            dbs11 = np.loadtxt(filter_Address, usecols=1)
            degs11 = np.loadtxt(filter_Address, usecols=2)
            mags11 = 10 ** (dbs11 / 20)
            res11 = mags11 * np.cos(degs11 / 180 * math.pi)
            ims11 = mags11 * np.sin(degs11 / 180 * math.pi)

            dbs21 = np.loadtxt(filter_Address, usecols=3)
            dbs21_add_VGA = dbs21 + Gain_VGA + 20 * log10(r_balun)
            degs21 = np.loadtxt(filter_Address, usecols=4)
            mags21 = 10 ** (dbs21_add_VGA / 20)
            res21 = mags21 * np.cos(degs21 / 180 * math.pi)
            ims21 = mags21 * np.sin(degs21 / 180 * math.pi)

        # Filter S parameter display
        if p == 0:
            dBs21 = np.zeros((3, len(freq)))
            dBs11 = np.zeros((3, len(freq)))
            dBs21_add_VGA = np.zeros((3, len(freq)))
        dBs21[p] = dbs21
        dBs11[p] = dbs11
        dBs21_add_VGA[p] = dbs21_add_VGA
        #  Interpolate to (30:1:250)MHz
        freqnew = np.arange(30, 251, 1)
        f_res11 = interpolate.interp1d(freq, res11, kind="cubic")
        res11new = f_res11(freqnew)
        f_ims11 = interpolate.interp1d(freq, ims11, kind="cubic")
        ims11new = f_ims11(freqnew)
        s11_complex = res11new + 1j * ims11new
        [f, filter_Gama_complex[:, p]] = expan(N, f0, f_start, f_end, s11_complex)

        f_res21 = interpolate.interp1d(freq, res21, kind="cubic")
        res21new = f_res21(freqnew)
        f_ims21 = interpolate.interp1d(freq, ims21, kind="cubic")
        ims21new = f_ims21(freqnew)
        s21_complex = res21new + 1j * ims21new
        [f, filter_s21_complex[:, p]] = expan(N, f0, f_start, f_end, s21_complex)

    filter_coefficient = (1 + filter_Gama_complex) * filter_s21_complex

    return cable_coefficient, filter_coefficient
