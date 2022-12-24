import numpy as np
from scipy.fftpack import fft
from scipy.fftpack import ifft
import math
from scipy import interpolate
from numpy.ma import log10, abs


# ================================================FFT get=============================================
def fftget(data_ori, N, f1):
    # This Python file uses the following encoding: utf-8

    # = == == == == This program is used as a subroutine to complete the FFT of data and generate parameters according to requirements == == == == =
    #  ----------------------input- ---------------------------------- %
    # % data_ori:time domain data, matrix form
    # % show_flag:flag of showing picture
    # % N:number of FFT points
    # % f1:Unilateral frequency
    # % ----------------------output - ---------------------------------- %
    # % data_fft:Frequency domain complex data
    # % data_fft_m_single:Frequency domain amplitude unilateral spectrum
    # % data_fft:Frequency domain phase

    lienum = data_ori.shape[1]
    data_fft = np.zeros((N, lienum), dtype=complex)
    data_fft_m = np.zeros((int(N), lienum))
    # data_fft_m_single = np.zeros((int(N/2), lienum))
    # data_fft_p = np.zeros((int(N), lienum))
    # data_fft_p_single = np.zeros((int(N/2), lienum))

    for i in range(lienum):
        data_fft[:, i] = fft(data_ori[:, i])

        data_fft_m[:, i] = abs(data_fft[:, i]) * 2 / N  # Amplitude
        data_fft_m[0] = data_fft_m[0] / 2

        data_fft_m_single = data_fft_m[0: len(f1)]  # unilateral

        data_fft_p = np.angle(data_fft, deg=True)  # phase
        data_fft_p = np.mod(data_fft_p, 2 * 180)
        # data_fft_p_deg = np.rad2deg(data_fft_p)
        data_fft_p_single = data_fft_p[0: len(f1)]

    return np.array(data_fft), np.array(data_fft_m_single), np.array(data_fft_p_single)


# =====================================IFFT get=================================================
def ifftget(data_ori, N, f1, true):
    # This Python file uses the following encoding: utf-8


    # %= == == == == This program is used as a subroutine to complete the Fourier change of data and generate parameters according to requirements == == == == =
    # % ----------------------input - ---------------------------------- %
    # % data_ori:Frequency domain data, complex numbers
    # % true  1 indicates that the complex number is synthesized, that is, the amplitude is the real amplitude. 2 indicates that the complex number is obtained after Fourier transform;
    # % N:number of FFT points
    # % t:time sequence
    # ns
    # % ----------------------output - ---------------------------------- %
    # % data_ifft :time domain data

    lienum = data_ori.shape[1]

    # %= == == == == == == == == == == == == == == First draw the spectrum phase == == == == == ==
    data_ori_m = np.zeros((int(N), lienum))
    data_ori_p = np.zeros((int(N), lienum))
    if true == 1:
        for i in range(lienum):
            data_ori_m[:, i] = abs(data_ori[:, i])  # Amplitude
            data_ori_m_single = data_ori_m[0: len(f1)]  # unilateral

            data_ori_p[:, i] = np.angle(data_ori[:, i], deg=True)  # phase
            data_ori_p[:, i] = np.mod(data_ori_p[:, i], 2 * 180)
            data_ori_p_single = data_ori_p[0: len(f1)]

    elif true == 2:
        for i in range(lienum):
            data_ori_m[:, i] = abs(data_ori[:, i]) * 2 / N
            data_ori_m[0] = data_ori_m[0] / 2

            data_ori_m_single = data_ori_m[0: len(f1)]  # double to single

            data_ori_p = np.angle(data_ori, deg=True)  # phase
            data_ori_p = np.mod(data_ori_p, 2 * 180)  # (-pi,pi) to (0,2pi)
            data_ori_p_single = data_ori_p[0: len(f1)]

    # % % Time domain
    data_ifft = np.zeros((N, lienum))
    for i in range(lienum):
        data_ifft[:, i] = ifft(data_ori[:, i]).real

    return np.array(data_ifft), np.array(data_ori_m_single), np.array(data_ori_p_single)
#==================================interpolation=======================================
def inter(data_complex_five,e_theta,e_phi):
    # This Python file uses the following encoding: utf-8

    # =================This subroutine is an interpolation procedure for a five-dimensional function in a specific theta phi direction======================
    # data_complex_five is the original data, 5 dimensions
    # e_theta,e_phi is the incident direction

    #    Four adjacent points
    down_theta = math.floor(e_theta)
    up_theta = math.ceil(e_theta)
    down_phi = math.floor(e_phi)
    up_phi = math.ceil(e_phi)

    a = abs(round(e_theta) - e_theta)
    b = abs(round(e_phi) - e_phi)

    numf=data_complex_five.shape[0]
    #    interpolation
    data_complex = np.zeros((181, 361), dtype=complex)
    data_new = np.zeros((numf, 3, 3), dtype=complex)
    for i in range(numf):
        for j in range(3):
            for k in range(3):
                data_complex[:, :] = data_complex_five[i, j, k, :, :]
                L1 = data_complex[down_theta, down_phi]
                L2 = data_complex[up_theta, down_phi]
                L3 = data_complex[down_theta, up_phi]
                L4 = data_complex[up_theta, up_phi]
                rt1 = (e_theta - down_theta) / 1.
                rt0 = 1.0 - rt1
                rp1 = (e_phi - down_phi) / 1.
                rp0 = 1.0 - rp1

                data_new[i, j, k] = rt0 * rp0 * L1 + rt1 * rp0 * L2 + rt0 * rp1 * L3 + rt1 * rp1 * L4
    return np.array(data_new)
# ==========================================complex expansion========================
def expan(N, f0, f1, f2, data):
    # This Python file uses the following encoding: utf-8

    # = == == == == This procedure is used as a subroutine to complete the expansion of the spectrum == == == == =
    # % N is the number of frequency points, that is, the spectrum that needs to be expanded
    # % f0 is the frequency step, MHz
    # % f1 is the starting frequency of the spectrum to be expanded, f2 is the cutoff frequency of the spectrum to be expanded
    # % The program only considers that the length of the expanded data is less than floor(N / 2), such as N = 10, the length of the expanded data <= 5; N = 9, the length of the expanded data <= 4
    # data 1 dimension

    f = np.arange(0, N) * f0  # Frequency sequence
    effective = len(data)
    delta_start = abs(f - f1)  # Difference from f1
    delta_end = abs(f - f2)  # Difference from f2
    f_hang_start = np.where(delta_start == min(delta_start))  # The row with the smallest difference
    f_hang_start = f_hang_start[0][0]
    f_hang_end = np.where(delta_end == min(delta_end))
    f_hang_end = f_hang_end[0][0]
    data_expansion = np.zeros((N), dtype=complex)
    if f_hang_start == 0:
        data_expansion[0] = data[0]
        add = np.arange(f_hang_end + 1, N - effective + 1, 1)
        duichen = np.arange(N - 1, N - effective + 1 - 1, -1)
        data_expansion[add] = 0
        data_expansion[f_hang_start: f_hang_end + 1] = data
        data_expansion[duichen] = data[1:].conjugate()
    else:
        a1 = np.arange(0, f_hang_start - 1 + 1, 1).tolist()
        a2 = np.arange(f_hang_end + 1, N - f_hang_start - effective + 1, 1).tolist()
        a3 = np.arange(N - f_hang_start + 1, N, 1).tolist()  # Need to make up 0;
        add = a1 + a2 + a3
        add = np.array(add)
        duichen = np.arange(N - f_hang_start, N - f_hang_start - effective, -1)
        data_expansion[add] = 0
        data_expansion[f_hang_start: f_hang_end + 1] = data[:]
        data_expansion[duichen] = data.conjugate()

    return f, data_expansion


#======================= PM: New Interpolation
def get_interpolation(freq_interp,freq_data,data,f_start,f_stop):
    """
    interpolates between the given frequency range, the frequency content below and above the measured
    range (normally 30-250 MHz) remain 0.

    :param freq_interp: array
                         required frequency to be interpolated (MHz)
    :param freq_data: array
                    frequency from data (MHz)

    :param data: array 3xn; n = len(freq_data)
                various measured hardware parameters (like LNA,filter,) or noise for 3 channels
    :param f_start: start frequency (MHZ) from data, normally 30 MHZ
    :param f_stop: highest frequency (MHz) from dat, normall 250 mHz
    :return:
       interpolated array 3xn, n= len(freq_interp)
    """
    amp_interp = np.zeros([3,len(freq_interp)]) #(3,n)
    interpolated1 = interpolate.interp1d(freq_data, data[0], kind='cubic')
    interpolated2 = interpolate.interp1d(freq_data, data[1], kind='cubic')
    interpolated3 = interpolate.interp1d(freq_data, data[2], kind='cubic')

    for idx, val in enumerate(freq_interp):
        if (val >= f_start) and (val <= f_stop):  # interpolating between min max frequency
            amp_interp[0,idx] = interpolated1(val)
            amp_interp[1,idx] = interpolated2(val)
            amp_interp[2,idx] = interpolated3(val)
    return amp_interp