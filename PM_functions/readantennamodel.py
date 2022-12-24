import numpy as np

def get_tabulated(filename):
    # mmap_mode reads only the necessary parts of the file from HD when needed
    # Speed reduction (cached by the OS) nommap->mmap 0.81 -> 0.27
    # Using arrays stored as float32 and avoiding dtype conversion ->0.02 (using float64 and avoiding dtype conversion does not change much? Not sure if I tested properly)
    f, R, X, theta, phi, lefft, leffp, phaset, phasep = np.load(filename, mmap_mode="r")
    # f, R, X, theta, phi, lefft, leffp, phaset, phasep = np.load(filename)
    n_f = f.shape[0]
    n_theta = len(np.unique(theta[0, :]))
    n_phi = R.shape[1] // n_theta
    shape = (n_f, n_phi, n_theta)

    dtype = "f4"
    f = f[:, 0].astype(dtype) * 1.0e6  # MHz --> Hz

    theta = theta[0, :n_theta].astype(dtype)  # deg
    phi = phi[0, ::n_theta].astype(dtype)  # deg

    # Those are not needed, so don't read them from the mmaped file
    # R = R.reshape(shape).astype(dtype)  # Ohm
    # X = X.reshape(shape).astype(dtype)  # Ohm
    lefft = lefft.reshape(shape).astype(dtype)  # m
    leffp = leffp.reshape(shape).astype(dtype)  # m
    # RK TODO: Make sure going from rad to deg does not affect calculations somewhere else.
    phaset = phaset.reshape(shape).astype(dtype)  # deg
    phasep = phasep.reshape(shape).astype(dtype)  # deg
    mydict=  {'frequency':f,
    'theta': theta,
    'phi': phi,
    # 'resistance':R,
    # 'reactance':X,
    'leff_theta':lefft,
    'phase_theta':phaset,
    'leff_phi':leffp,
    'phase_phi': phasep}
        
    return mydict


def get_tabulated_mod(filename):
    """As get_tabulated, but for files with reshaped arrays"""
    # mmap_mode reads only the necessary parts of the file from HD when needed
    # Readout uncached, cached by OS/HDD, rest of the voltage calc:
    # Original (float64): ~48 s, 0.7, ~2.1
    # Original, mmap: ~26.5 s, 0.25, ~2
    # Below, no dtype conversion in the script:
    # Reshaped, float64: ~48 s, ~0.5, ~2.4
    # Reshaped, float64, mmap: ~5.7 s, 0.025, ~2
    # Reshaped, float32: ~22.9 s, ~0.33, ~2.1
    # Reshaped, float32, mmap: ~2.9 s, ~0.017, ~2.0
    # Compressed, reshaped, float32: ~6.5 s
    # Speed reduction (cached by the OS) nommap->mmap 0.81 -> 0.27
    # Using arrays stored as float32 and avoiding dtype conversion ->0.02 (using float64 and avoiding dtype conversion does not change much? Not sure if I tested properly)
    # compressed is unaddected by mmap, always 3.87 s
    f, R, X, theta, phi, lefft, leffp, phaset, phasep = np.load(filename, mmap_mode="r")
    # f, R, X, theta, phi, lefft, leffp, phaset, phasep = np.load(filename)
    n_f = f.shape[1]
    # n_theta = len(np.unique(theta[0, :]))
    n_theta = len(np.unique(theta[:, 0]))
    # print(n_theta)
    # exit()
    # n_theta = 181
    n_phi = R.shape[0] // n_theta
    shape = (n_phi, n_theta, n_f)

    # dtype = "f4"
    dtype = np.float32
    # f = f[0, :].astype(dtype) * 1.0e6  # MHz --> Hz
    f = f[0, :] * 1.0e6  # MHz --> Hz
    theta = theta[:n_theta, 0]#.astype(dtype)  # deg
    phi = phi[::n_theta, 0]#.astype(dtype)  # deg
    # theta = np.arange(181).astype(dtype)
    # phi = np.arange(361).astype(dtype)
    # print(theta, phi)
    # exit()

    # Those are not needed, so don't read them from the mmaped file
    # R = R.reshape(shape).astype(dtype)  # Ohm
    # X = X.reshape(shape).astype(dtype)  # Ohm
    lefft = lefft.reshape(shape)  # m
    leffp = leffp.reshape(shape)  # m
    # RK TODO: Make sure going from rad to deg does not affect calculations somewhere else.
    phaset = phaset.reshape(shape)  # deg
    phasep = phasep.reshape(shape)  # deg
    mydict = {'frequency': f,
              'theta': theta,
              'phi': phi,
              # 'resistance':R,
              # 'reactance':X,
              'leff_theta': lefft,
              'phase_theta': phaset,
              'leff_phi': leffp,
              'phase_phi': phasep}

    return mydict


def get_interp(table,zenith,azimuth,dlength,dt):
    dtheta = table['theta'][1] -  table['theta'][0] # deg
    rt1 = (zenith - table['theta'][0]) / dtheta
    it0 = int(np.floor(rt1) % table['theta'].size)
    it1 = it0 + 1
    if it1 == table['theta'].size:  # Prevent overflow
        it1, rt1 = it0, 0
    else:
        rt1 -= np.floor(rt1)
    rt0 = 1 - rt1

    dphi = table['phi'][1] - table['phi'][0]  # deg
    rp1 = (azimuth - table['phi'][0]) / dphi
    ip0 = int(np.floor(rp1) % table['phi'].size)
    ip1 = ip0 + 1
    if ip1 == table['phi'].size:  # Results are periodic along phi
        ip1 = 0
    rp1 -= np.floor(rp1)
    rp0 = 1 - rp1

    def fftfreq(n, dt):
            
            return np.fft.rfftfreq(n, dt)




    def interp(v):
        fp = (
                rp0 * rt0 * v[:, ip0, it0]
                + rp1 * rt0 * v[:, ip1, it0]
                + rp0 * rt1 * v[:, ip0, it1]
                + rp1 * rt1 * v[:, ip1, it1]
            )
        return np.interp(x, xp, fp, left=0, right=0)

    x = fftfreq(dlength,dt)
    xp = table['frequency']  # frequency [Hz]

    # The interpolations below took almost all the time in this function
    ltr = interp(table['leff_theta'])  # LWP. m
    lta = interp(np.deg2rad(table['phase_theta']))  # LWP. rad
    lpr = interp(table['leff_phi'])  # LWP. m
    lpa = interp(np.deg2rad(table['phase_phi']))  # LWP. rad

    lt = ltr * np.exp(1j * lta)
    lp = lpr * np.exp(1j * lpa)

    return lt,lp


def get_interp_mod(table, zenith, azimuth, dlength, dt):
    dtheta = table['theta'][1] - table['theta'][0]  # deg
    rt1 = (zenith - table['theta'][0]) / dtheta
    it0 = int(np.floor(rt1) % table['theta'].size)
    it1 = it0 + 1
    if it1 == table['theta'].size:  # Prevent overflow
        it1, rt1 = it0, 0
    else:
        rt1 -= np.floor(rt1)
    rt0 = 1 - rt1

    dphi = table['phi'][1] - table['phi'][0]  # deg
    rp1 = (azimuth - table['phi'][0]) / dphi
    ip0 = int(np.floor(rp1) % table['phi'].size)
    ip1 = ip0 + 1
    if ip1 == table['phi'].size:  # Results are periodic along phi
        ip1 = 0
    rp1 -= np.floor(rp1)
    rp0 = 1 - rp1

    def fftfreq(n, dt):

        return np.fft.rfftfreq(n, dt)

    def interp(v):
        fp = (
                rp0 * rt0 * v[ip0, it0, :]
                + rp1 * rt0 * v[ip1, it0, :]
                + rp0 * rt1 * v[ip0, it1, :]
                + rp1 * rt1 * v[ip1, it1, :]
        )
        return np.interp(x, xp, fp, left=0, right=0)

    x = fftfreq(dlength, dt)
    xp = table['frequency']  # frequency [Hz]

    # The interpolations below took almost all the time in this function
    ltr = interp(table['leff_theta'])  # LWP. m
    lta = interp(np.deg2rad(table['phase_theta']))  # LWP. rad
    lpr = interp(table['leff_phi'])  # LWP. m
    lpa = interp(np.deg2rad(table['phase_phi']))  # LWP. rad

    lt = ltr * np.exp(1j * lta)
    lp = lpr * np.exp(1j * lpa)

    return lt, lp





    

