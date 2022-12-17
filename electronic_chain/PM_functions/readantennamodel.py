import numpy as np
import matplotlib.pyplot as plt




def get_tabulated(filename):
    print("loading")
    f, R, X, theta, phi, lefft, leffp, phaset, phasep = np.load(filename)
    print("loaded")
    n_f = f.shape[0]
    n_theta = len(np.unique(theta[0, :]))
    n_phi = int(R.shape[1] / n_theta)
    shape = (n_f, n_phi, n_theta)

    dtype = "f4"
    f = f[:, 0].astype(dtype) * 1.0e6  # MHz --> Hz
    theta = theta[0, :n_theta].astype(dtype)  # deg
    phi = phi[0, ::n_theta].astype(dtype)  # deg
    R = R.reshape(shape).astype(dtype)  # Ohm
    X = X.reshape(shape).astype(dtype)  # Ohm
    lefft = lefft.reshape(shape).astype(dtype)  # m
    leffp = leffp.reshape(shape).astype(dtype)  # m
    # RK TODO: Make sure going from rad to deg does not affect calculations somewhere else.
    phaset = phaset.reshape(shape).astype(dtype)  # deg
    phasep = phasep.reshape(shape).astype(dtype)  # deg
    mydict=  {'frequency':f,
    'theta': theta,
    'phi': phi,
    'resistance':R,
    'reactance':X,
    'leff_theta':lefft,
    'phase_theta':phaset,
    'leff_phi':leffp,
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
    
    ltr = interp(table['leff_theta'])  # LWP. m
    lta = interp(np.deg2rad(table['phase_theta']))  # LWP. rad
    lpr = interp(table['leff_phi'])  # LWP. m
    lpa = interp(np.deg2rad(table['phase_phi']))  # LWP. rad
    lt = ltr * np.exp(1j * lta)
    lp = lpr * np.exp(1j * lpa)
    return lt,lp
    
    
    

