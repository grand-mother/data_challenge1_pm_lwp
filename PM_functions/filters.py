import electronic_chain.ec_config as ec_config
from PM_functions.functions import *
import PM_functions.config as PM_config
from scipy import interpolate
from numpy.ma import log10, abs


# ===============================================Filterparameter get====================================
def filter_get(freqs_tot, unit):
    # This Python file uses the following encoding: utf-8


    # = == == == == This program is used as a subroutine to complete the calculation and expansion of the S parameters of the cable and filter == == == == =
    #  ----------------------input - ---------------------------------- %
    # filter :path of s2p file
    # The filter test data is in the Filterparameter folder (delete the previous string of the s2p file in advance and name it 1, because the three port filters are the same)
    # The cable test data is in the cableparameter folder (delete the previous string of the s2p file in advance and name it cable because the three ports are the same)
    # The unit of test frequency is hz
    # freqs_tot : array of new frequency to be interpolated
    # If unit is 0, the test data is in the form of real and imaginary parts, and 1 is in the form of db and phase.

    # % ----------------------output - ---------------------------------- %
    # cable_coefficient
    # filter_coefficient

    Gain_VGA = -1.5  # dB
    r_balun = 630 * 2 / 650
    # test filter without VGA
    # Gain_VGA = 0  # dB
    # r_balun = 1
    cable_Gama_complex = np.zeros([3,len(freqs_tot)], dtype=complex)  # 3 ports
    cable_s21_complex = np.zeros([3,len(freqs_tot)], dtype=complex)
    res11_all = []
    ims11_all = []
    res21_all = []
    ims21_all = []
    for p in range(3):
        #  cable参数
        # str_p = str(p + 1)
        cable_Address = PM_config.PM_files_path+"/cableparameter/cable.s2p"
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

        res11_all.append(res11)
        ims11_all.append(ims11)
        res21_all.append(res21)
        ims21_all.append(ims21)
    res11_3 = np.array([res11_all[0], res11_all[1], res11_all[2]])  # 3 channel (3,n)
    ims11_3 = np.array([ims11_all[0], ims11_all[1], ims11_all[2]])  
    res21_3 = np.array([res21_all[0], res21_all[1], res21_all[2]])  
    ims21_3 = np.array([ims21_all[0], ims21_all[1], ims21_all[2]])  

    f_start = 30
    f_end = 250

    # PM: New Interpolation: gets rid of the expan term should return len(1001) (previous 2000)

    res11_new = get_interpolation(freqs_tot, freq, res11_3, f_start, f_end)

    ims11_new = get_interpolation(freqs_tot, freq, ims11_3, f_start, f_end)
    for p in np.arange(3):
        cable_Gama_complex[p] = res11_new[p] + 1j * ims11_new[p]

   # [f, cable_Gama_complex[:, p]] = expan(N, f0, f_start, f_end, s11_complex)

    res21_new = get_interpolation(freqs_tot, freq, res21_3, f_start, f_end)

    ims21_new = get_interpolation(freqs_tot, freq, ims21_3, f_start, f_end)
    for p in np.arange(3):
        cable_s21_complex[p] = res21_new[p] + 1j * ims21_new[p]

    #[f, cable_s21_complex[:, p]] = expan(N, f0, f_start, f_end, s21_complex)

    cable_coefficient = (1 + cable_Gama_complex) * cable_s21_complex

    #  =================filter=========================================

    filter_Gama_complex = np.zeros([3,len(freqs_tot)], dtype=complex)  # 3 ports
    filter_s21_complex = np.zeros([3,len(freqs_tot)], dtype=complex)
    res11_all = []
    ims11_all = []
    res21_all = []
    ims21_all = []
    for p in range(3):
        #  filter parameter
        # str_p = str(p + 1)
        filter_Address = PM_config.PM_files_path+"/filterparameter/1.s2p"
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

        res11_all.append(res11)
        ims11_all.append(ims11)
        res21_all.append(res21)
        ims21_all.append(ims21)
    res11_3 = np.array([res11_all[0], res11_all[1], res11_all[2]])  # 3 channel (3,n)
    ims11_3 = np.array([ims11_all[0], ims11_all[1], ims11_all[2]])  
    res21_3 = np.array([res21_all[0], res21_all[1], res21_all[2]])  
    ims21_3 = np.array([ims21_all[0], ims21_all[1], ims21_all[2]])  

    # PM: New Interpolation: gets rid of the expan term should return len(1001) (previous 2000)

    res11_new = get_interpolation(freqs_tot, freq, res11_3, f_start, f_end)

    ims11_new = get_interpolation(freqs_tot, freq, ims11_3, f_start, f_end)
    for p in np.arange(3):
        filter_Gama_complex[p] = res11_new[p] + 1j * ims11_new[p]

    res21_new = get_interpolation(freqs_tot, freq, res21_3, f_start, f_end)

    ims21_new = get_interpolation(freqs_tot, freq, ims21_3, f_start, f_end)
    for p in np.arange(3):
        filter_s21_complex[p] = res21_new[p] + 1j * ims21_new[p]



    filter_coefficient = (1 + filter_Gama_complex) * filter_s21_complex

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"cable_coefficient": cable_coefficient, "filter_coefficient": filter_coefficient}
    # Outside pipeline return - raw values
    else:
        return cable_coefficient, filter_coefficient


    # return cable_coefficient, filter_coefficient
