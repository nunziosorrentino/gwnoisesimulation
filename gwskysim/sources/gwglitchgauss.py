import numpy as np

from gwskysim.sources.gwnoisesource import GWNoiseSource


class GWGlitchGauss(GWNoiseSource): #gaussian glitch

    def __init__(self, name, t0,tau,SNR=None):
        GWNoiseSource.__init__(self, name, t0,tau,SNR)


    def __call__(self, times_array):

        times_array = times_array[times_array <= self.t0 + 5 * self.tau]
        times_array = times_array[times_array >= self.t0 - 5 * self.tau]
        if np.array(times_array).size==0:
            return None          

        if self.SNR is None:

            #print("no noise--->random amplitude")

            h0 = np.random.uniform(10 ** (-22), 5 * 10 ** (-21))

            return h0 * np.exp(-((times_array - self.t0) ** 2) / (2 * self.tau ** 2)),times_array[0]



        else:

            h0=10**-20
            aux = h0 * np.exp(-((times_array - self.t0) ** 2) / (2 * (self.tau) ** 2))
            snr_test = GWNoiseSource.SNR(aux, times_array, self.f_asd_tab, self.asd_tab)

            AUX= (self.SNR/ snr_test) * aux

            #AUX2 = GWNoiseSource.gluer(AUX, times_array_reduced, self.t0, times_array)

            return AUX,times_array[0]



