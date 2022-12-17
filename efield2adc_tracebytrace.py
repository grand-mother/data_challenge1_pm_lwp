#!/usr/bin/python

from grand.io.root_trees import *
from electronic_chain.trace_functions import *

from electronic_chain.XDU_electronic_chain.LNA import LNA_get
from electronic_chain.XDU_electronic_chain.antenna_effective_length import CEL
from electronic_chain.XDU_electronic_chain.galactic_noise import gala
from electronic_chain.XDU_electronic_chain.filters import filter_get


def main():
    # Open the ROOT tree with simulation shower data
    tshower = ShowerEventSimdataTree("../data/Stshp_LWP_S23d_Proton_5.6_87.9_290.0_1/Coarse2.root")
    # Open the ROOT tree with traces
    tefield = EfieldEventTree("../data/Stshp_LWP_S23d_Proton_5.6_87.9_290.0_1/Coarse2.root")
    # Open the ROOT tree with simulation run info
    trunefieldsimdata = EfieldRunSimdataTree("../data/Stshp_LWP_S23d_Proton_5.6_87.9_290.0_1/Coarse2.root")
    trunefieldsimdata.get_entry(0)
    # In this file antennas are misordered. This is the dictionary that gives ROOT file trace number for given trace file number
    # trace_file_to_root = {95: 62, 94: 95, 93: 42, 92: 19, 91: 21, 90: 63, 89: 14, 88: 41, 87: 0, 86: 45, 85: 54, 84: 22, 83: 20, 82: 7, 81: 69, 80: 77, 79: 68, 78: 81, 77: 30, 76: 76, 75: 17, 74: 3, 73: 18, 72: 32, 71: 53, 70: 55, 69: 4, 68: 90, 67: 82, 66: 40, 65: 12, 64: 93, 63: 13, 62: 38, 61: 44, 60: 84, 59: 15, 58: 70, 57: 58, 56: 75, 55: 78, 54: 23, 53: 8, 52: 34, 51: 29, 50: 9, 49: 37, 48: 59, 47: 79, 46: 16, 45: 83, 44: 86, 43: 47, 42: 36, 41: 92, 40: 51, 39: 66, 38: 49, 37: 31, 36: 57, 35: 6, 34: 50, 33: 61, 32: 56, 31: 89, 30: 48, 29: 65, 28: 67, 27: 64, 26: 71, 25: 28, 24: 52, 23: 27, 22: 25, 21: 91, 20: 35, 19: 85, 18: 33, 17: 87, 16: 24, 15: 88, 14: 46, 13: 39, 12: 73, 11: 72, 10: 5, 9: 94, 8: 26, 7: 43, 6: 11, 5: 60, 4: 80, 3: 74, 2: 10, 1: 1, 0: 2}

    tefield1 = EfieldEventTree("../result/myefield.root")

    for i in range(tshower.get_entries()):
        tefield.get_entry(i)
        # A bug in root_trees? Get entry should not be necessary
        tshower.get_entry(i)
        tefield1.copy_contents(tefield)

        # e_theta = float(list1[5])
        # e_phi = float(list1[6])
        e_theta = 180-tshower.shower_azimuth
        e_phi = 180+tshower.shower_zenith

        print("theta is:", e_theta, "degree", tshower.shower_azimuth)
        print("phi is:", e_phi, "degree", tshower.shower_zenith)

        # savetxt = 1  # savetxt=0 Don't save output file ,savetxt=1 save output file
        show_flag = 0

        # ===================================Arbitrary input file first generates input parameters needed by subroutine================================================
        Ts = trunefieldsimdata.t_bin_size

        fs = 1 / Ts * 1000  # sampling frequency, MHZ
        N = math.ceil(fs)
        f0 = fs / N  # base frequency, Frequency resolution
        f = np.arange(0, N) * f0 # frequency sequence
        f1 = f[0:int(N / 2) + 1]

        # =======Equivalent length================
        [Lce_complex, antennas11_complex_short] = CEL(e_theta, e_phi, N, f0, 1)
        # np.save("Lce_complex", Lce_complex)
        # np.save("antennas11_complex_short", antennas11_complex_short)
        # exit()
        # Lce_complex = np.load("Lce_complex.npy")
        # antennas11_complex_short = np.load("antennas11_complex_short.npy")

        Lcehang = Lce_complex.shape[0]
        Lcelie = Lce_complex.shape[2]
        print("4")
        # ======Galaxy noise power spectrum density, power, etc.=====================
        lst = 18
        [galactic_v_complex_double, galactic_v_time] = gala(lst, N, f0, f1, tefield.du_count)
        # np.save("galactic_v_complex_double", galactic_v_complex_double)
        # np.save("galactic_v_time", galactic_v_time)
        # galactic_v_complex_double = np.load("galactic_v_complex_double.npy")
        # galactic_v_time = np.load("galactic_v_time.npy")
        print("5")
        # =================== LNA=====================================================
        [rou1_complex, rou2_complex, rou3_complex] = LNA_get(antennas11_complex_short, N, f0, 1)
        print("6")
        # =======================  cable  filter VGA balun=============================================
        [cable_coefficient, filter_coefficient] = filter_get(N, f0, 1)

        print("starting main loop")
        tefield1.trace_x.clear()
        tefield1.trace_y.clear()
        tefield1.trace_z.clear()
        for trace_num in range(len(tefield.trace_x)):
            # This is for Coarse2.root only, with manually restructured order of antennas
            # Works with Stshp_LWP_S23d_Proton_5.6_87.9_290.0_1 and Coarse2.root that I've put there
            # Reorderring was needed for sane comparison with the original electronic chain results
            # ex = tefield.trace_x[trace_file_to_root[trace_num]]
            # ey = tefield.trace_y[trace_file_to_root[trace_num]]
            # ez = tefield.trace_z[trace_file_to_root[trace_num]]
            # This is assuming the ROOT source file has a correct order of antennas
            ex = tefield.trace_x[trace_num]
            ey = tefield.trace_y[trace_num]
            ez = tefield.trace_z[trace_num]

            antenna_model = None
            noise_model = [galactic_v_time, galactic_v_complex_double]
            electronic_chain = {"LNA":rou1_complex * rou2_complex * rou3_complex, "cable":cable_coefficient, "filter":filter_coefficient}

            ex, ey, ez = adjust_traces(ex, ey, ez, Ts)

            adc0, adc1, adc2, v_t, v_f = efield_to_adc(ex, ey, ez, e_phi, e_theta, antenna_model, noise_model, electronic_chain, return_voltages=True)

            tefield1.trace_x.append(v_t[-1][:,0].astype(np.float32).tolist())
            tefield1.trace_y.append(v_t[-1][:,1].astype(np.float32).tolist())
            tefield1.trace_z.append(v_t[-1][:,2].astype(np.float32).tolist())

        tefield1.fill()
        tefield1.write()


if __name__ == '__main__':
	main()
