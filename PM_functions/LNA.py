import electronic_chain.ec_config as ec_config
from PM_functions.functions import *
import PM_functions.config as PM_config
from scipy import interpolate
from numpy.ma import log10, abs

# ===========================================LNAparameter get===========================================
def LNA_get(freqs_tot, unit):
    # This Python file uses the following encoding: utf-8

    # = == == == == This program is used as a subroutine to complete the calculation and expansion of the LNA partial pressure coefficient,it also reads antenna VSWR == == == == =
    #  ----------------------input - ---------------------------------- %
    # freqs_tot : array of new frequencies to be interpolated (MHz)
    # antennas11_complex_short is the program after the interpolation of the antenna standing wave test result
    # LNA address, s2p file (delete the previous string of s2p file in advance, put the test results of the three ports in the LNASparameter folder, and name them 1 2 3 in turn, corresponding to xyz)
    # The unit of test frequency is hz

    # If unit is 0, the test data is in the form of real and imaginary parts, and 1 is in the form of db and phase.
    # % ----------------------output - ---------------------------------- %
    # rou1 rou2 rou3

    #====================== Antenna VSWR part

    s11_complex = np.zeros([3,len(freqs_tot)], dtype=complex)  # 3 ports
    re_all = []
    im_all = []

    for p in range(3):
        str_p = str(p + 1)
        filename = PM_config.PM_files_path + "/antennaVSWR//" + str_p + ".s1p"
        freq = np.loadtxt(filename, usecols=0) / 1e6  # HZ to MHz
        f_start = freq[0]
        f_end = freq[-1] # can later be input 30-250 MHz
        
        if unit == 0:
            re = np.loadtxt(filename, usecols=1)
            im = np.loadtxt(filename, usecols=2)
            db = 20 * log10(abs(re + 1j * im))
        elif unit == 1:
            db = np.loadtxt(filename, usecols=1)
            deg = np.loadtxt(filename, usecols=2)
            mag = 10 ** (db / 20)
            re = mag * np.cos(deg / 180 * math.pi)
            im = mag * np.sin(deg / 180 * math.pi)
        if p == 0:
            dB = np.zeros((3, len(freq)))
        dB[p] = db
        re_all.append(re)
        im_all.append(im)
    re_3 = np.array([re_all[0],re_all[1],re_all[2]])  #3 channel (3,n)
    im_3 = np.array([im_all[0],im_all[1],im_all[2]])  #3channel (3,n)
    #PM: New Interpolation: gets rid of the expan term should return len(1001) (previous 2000)

    re_new = get_interpolation(freqs_tot,freq,re_3,f_start,f_end)

    im_new = get_interpolation(freqs_tot,freq,im_3,f_start,f_end)
    for p in np.arange(3):
        s11_complex[p] = re_new[p] + 1j * im_new[p]

    z0 = 50 # characteristic impedance

    antenna_Gama_complex = s11_complex 
    '''
    for p in range(3):
        # Antenna related parameter calculation
        antennas11_short = s11_complex #antennas11_complex_short
        #f0 = 1
        f_start = 30
        f_end = 250
        [f, antenna_Gama_complex[:, p]] = expan(N, f0, f_start, f_end, antennas11_short[:, p])
    '''
    Zin_antenna = z0 * (1 + antenna_Gama_complex) / (1 - antenna_Gama_complex)

    # =========================== LNA part ========================
    LNA_Gama_complex = np.zeros([3,len(freqs_tot)], dtype=complex)  # 3 ports
    LNA_s21_complex = np.zeros([3,len(freqs_tot)], dtype=complex)
    res11_all =[]
    ims11_all =[]
    res21_all= []
    ims21_all = []
    for p in range(3):
        #  LNA parameter
        str_p = str(p + 1)
        LNA_Address = PM_config.PM_files_path+"/LNASparameter/" + str_p + ".s2p"
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
        res11_all.append(res11)
        ims11_all.append(ims11)
        res21_all.append(res21)
        ims21_all.append(ims21)
    res11_3 = np.array([res11_all[0], res11_all[1], res11_all[2]])  # 3 channel (3,n)
    ims11_3 = np.array([ims11_all[0], ims11_all[1], ims11_all[2]])  # 3channel (3,n)
    res21_3 = np.array([res21_all[0], res21_all[1], res21_all[2]])  # 3 channel (3,n)
    ims21_3 = np.array([ims21_all[0], ims21_all[1], ims21_all[2]])  # 3channel (3,n)


    # PM: New Interpolation: gets rid of the expan term should return len(1001) (previous 2000)

    res11_new = get_interpolation(freqs_tot, freq, res11_3, f_start, f_end)

    ims11_new = get_interpolation(freqs_tot, freq, ims11_3, f_start, f_end)
    for p in np.arange(3):
        LNA_Gama_complex[p] = res11_new[p] + 1j * ims11_new[p]
    #LNA_Gama_complex = s11_complex

    res21_new = get_interpolation(freqs_tot, freq, res21_3, f_start, f_end)

    ims21_new = get_interpolation(freqs_tot, freq, ims21_3, f_start, f_end)
    for p in np.arange(3):
        LNA_s21_complex[p] = res21_new[p] + 1j * ims21_new[p]


    #LNA_s21_complex = s21_complex

    Zin_LNA = z0 * (1 + LNA_Gama_complex) / (1 - LNA_Gama_complex)

    # Partial pressure coefficient
    rou1 = Zin_LNA / (Zin_antenna + Zin_LNA)
    rou2 = (1 + LNA_Gama_complex) / (1 - antenna_Gama_complex * LNA_Gama_complex)
    rou3 = LNA_s21_complex / (1 + LNA_Gama_complex)

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"LNA_coefficient": rou1*rou2*rou3}
    # Outside pipeline return - raw values
    else:
        return rou1*rou2*rou3

    # return rou1, rou2, rou3
