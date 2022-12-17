import h5py
from scipy import interpolate
from electronic_chain.XDU_electronic_chain.functions import *
import electronic_chain.XDU_electronic_chain.config as config
from functools import lru_cache

# This decorator stores the function result in memory, so the function is called only once, and next calls just read the result from memory
@lru_cache(maxsize=1)
def CEL(e_theta, e_phi, N, f0, unit):
    # This Python file uses the following encoding: utf-8

    # from complex_expansion import expan

    # = == == == == This program is used as a subroutine to complete the calculation and expansion of the 30-250MHz complex equivalent length == == == == =
    #  ----------------------input- ---------------------------------- %
    # filename address, S1P file (delete the previous string of s1p file in advance, put the test results of the three ports in the antennaVSWR folder, and name them 1 2 3 in turn)
    # N is the extended length
    # e_theta, e_phi is the direction of incidence
    # If unit is 0, the test data is in the form of real and imaginary parts, and 1 is in the form of db and phase.
    # f0 is the frequency resolution,
    # % ----------------------output - ---------------------------------- %
    # f frequency sequence, the default unit is MHz
    # Lce_complex_expansion is the equivalent length of a specific incident direction
    # s11_complex is the antenna test data

    # Complex electric field 30-250MHz
    REfile = config.XDU_files_path+"/Complex_RE.mat"
    RE = h5py.File(REfile, 'r')
    RE_zb = np.transpose(RE['data_rE_ALL'])
    re_complex = RE_zb.view('complex')
    f_radiation = np.transpose(RE['f_radiation'])  # mhz
    effective = max(f_radiation.shape[0], f_radiation.shape[1])

    e_radiation = inter(re_complex, e_theta, e_phi)

    # 测试s1p
    s11_complex = np.zeros((effective, 3), dtype=complex)  # 3 ports
    for p in range(3):
        str_p = str(p + 1)
        filename = config.XDU_files_path+"/antennaVSWR/" + str_p + ".s1p"
        freq = np.loadtxt(filename, usecols=0) / 1e6  # HZ to MHz
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

        # Interpolation is a data of 30-250mhz interval 1mhz
        freqnew = np.arange(30, 251, 1)
        f_re = interpolate.interp1d(freq, re, kind="cubic")
        renew = f_re(freqnew)
        f_im = interpolate.interp1d(freq, im, kind="cubic")
        imnew = f_im(freqnew)
        s11_complex[:, p] = renew + 1j * imnew

    print("3_3")

    # %Reduced current
    z0 = 50
    a1 = 1
    I_complex = 1 / math.sqrt(z0) * (1 - s11_complex) * a1

    # %Denominator
    eta = 120 * math.pi
    c = 3 * 1e8
    f_unit = 1e6
    lamda = c / (f_radiation * f_unit)  # m
    lamda = np.transpose(lamda)
    k = 2 * math.pi / lamda
    fenmu = 1j * (I_complex / (2 * lamda) * eta)

    # Equivalent length
    # Extend the frequency band
    f1 = f_radiation[0][0]
    f2 = f_radiation[0][-1]

    print("3_4")

    Lce_complex_short = np.zeros((effective, 3, 3), dtype=complex)
    Lce_complex_expansion = np.zeros((N, 3, 3), dtype=complex)
    for i in range(3):  # Polarization i = 1, 2, 3 respectively represent xyz polarization
        for p in range(3):
            # Xyz polarization of a single port
            Lce_complex_short[:, i, p] = e_radiation[:, p, i] / fenmu[:, p]
            [f, Lce_complex_expansion[:, i, p]] = expan(N, f0, f1, f2, Lce_complex_short[:, i, p])

    print("3_5")

    return Lce_complex_expansion, s11_complex


