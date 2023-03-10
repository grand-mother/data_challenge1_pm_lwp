import electronic_chain.ec_config as ec_config
import PM_functions.readantennamodel as an
import PM_functions.config as PM_config
import numpy as np
from scipy.spatial.transform import Rotation as R
from functools import lru_cache
from misc_functions import time_passed

@lru_cache(maxsize=1)
def read_antenna_files(old=False):
    time_passed(True)
    if old:
        suffix = ".npy"
        table_ewarm_new = an.get_tabulated(PM_config.PM_files_path + f"/GP300Antenna_EWarm_leff{suffix}")
        table_snarm_new = an.get_tabulated(PM_config.PM_files_path + f"/GP300Antenna_SNarm_leff{suffix}")
        table_zarm_new = an.get_tabulated(PM_config.PM_files_path + f"/GP300Antenna_Zarm_leff{suffix}")
    else:
        suffix = "_reshaped_float32.npy"
        # suffix = "_reshaped_compressed_float32.npz"
        table_ewarm_new = an.get_tabulated_mod(PM_config.PM_files_path + f"/GP300Antenna_EWarm_leff{suffix}")
        table_snarm_new = an.get_tabulated_mod(PM_config.PM_files_path + f"/GP300Antenna_SNarm_leff{suffix}")
        table_zarm_new = an.get_tabulated_mod(PM_config.PM_files_path + f"/GP300Antenna_Zarm_leff{suffix}")
    print(time_passed())

    return table_ewarm_new, table_snarm_new, table_zarm_new

def efield2voltage_pm(traces_t, e_theta, e_phi, freqs, sampling_time=0.5, old=False, **kwargs):
    """""Voltage calculation from E field traces - by Pragati Mitra
    Parameters:
        traces_t : time traces of E field for x,y,z, (:, 3, :)
        Zenith: shower zenith in degree
        Azimuth: shower azimuth in degree
        freqs: frequencies (real part)
    Returns:
        Voltage trace (time domain)
    """

    Zenith, Azimuth = e_theta, e_phi

    ## read antenna response function Leff from files in theta, phi direction
    table_ewarm_new, table_snarm_new, table_zarm_new = read_antenna_files(old=old)

    time_passed(True)
    # interpolated L_eff in 3 arms for given zenith,azimuth

    N = traces_t.shape[-1]
    dt = sampling_time

    # The interpolations below take almost all the time in this function
    if old:
        lt1, lp1 = an.get_interp(table_snarm_new, Zenith, Azimuth, N, dt * 1e-9)
        lt2, lp2 = an.get_interp(table_ewarm_new, Zenith, Azimuth, N, dt * 1e-9)
        lt3, lp3 = an.get_interp(table_zarm_new, Zenith, Azimuth, N, dt * 1e-9)
    else:
        lt1, lp1 = an.get_interp_mod(table_snarm_new, Zenith, Azimuth, N, dt * 1e-9)
        lt2, lp2 = an.get_interp_mod(table_ewarm_new, Zenith, Azimuth, N, dt * 1e-9)
        lt3, lp3 = an.get_interp_mod(table_zarm_new, Zenith, Azimuth, N, dt * 1e-9)
    lt = np.array([lt1, lt2, lt3]).T
    lp = np.array([lp1, lp2, lp3]).T

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

    # E_traces_fft = np.array(np.fft.rfft(np.moveaxis(traces_t, -2, -1)), dtype='complex')  # fft , shape(2,N)
    E_traces_fft = np.array(np.fft.rfft(traces_t), dtype='complex')  # fft , shape(2,N)

    # Fold antenna response in and get Voc
    # ToDo: Do the same stuff for 2D traces as for 3D traces array
    if traces_t.ndim<3:
        E_tp_fft = np.zeros([len(freqs), 2], dtype='complex') # xyz to theta-phi
        for k in range(len(freqs)):
            E_tp_fft[k] = xyz_thetaphi(E_traces_fft[:, k], Zenith, Azimuth)
    else:
        rm = R.from_matrix(rotationmatrix(Zenith, Azimuth).T)
        # rotation can be applied only on (x,3) arrays, so we have to do a lot of reshaping back and forth
        ee = np.moveaxis(E_traces_fft, 1, 2)
        E_tp_fft = rm.apply(ee.reshape(-1, 3)).reshape((traces_t.shape[0], -1, 3))


    Voc_shower_complex = np.moveaxis(lt * E_tp_fft[..., 0, np.newaxis] + lp * E_tp_fft[..., 1, np.newaxis], -2, -1)
    Voc_shower_t = np.fft.irfft(Voc_shower_complex, n=N)  # (3,N)

    print(time_passed())

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        # return {"traces_t": Voc_shower_t, "traces_f": Voc_shower_complex}
        return {"traces_t": Voc_shower_t, "traces_f": Voc_shower_complex}
    # Outside pipeline return - raw values
    else:
        # return Voc_shower_t, Voc_shower_complex
        return Voc_shower_t, Voc_shower_complex

