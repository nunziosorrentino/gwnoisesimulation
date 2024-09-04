import numpy as np

from gwskysim.sources.gwnoisesource import GWNoiseSource


class GWGlitchWhistleLike(GWNoiseSource): #whistlelike glitch

    def __init__(self,name,t0,tau, f0,SNR=None):
        GWNoiseSource.__init__(self, name, t0 ,tau,SNR)
        self.f0 = f0

    def __call__(self, times_array):
        times_array = times_array[times_array <= self.t0 + 8 * self.tau]
        times_array = times_array[times_array >= self.t0 - 8 * self.tau]
        if np.array(times_array).size==0:
            return None         

        if self.SNR is None:
        # if self.SNR is None and self.detector is None :
            #print("no asd curve provided--> random amplitude")
            h0 = np.random.uniform(10 ** (-22), 5 * 10 ** (-21))
            phi =2.* np.pi * self.f0*(times_array-self.t0)*(1.-3*self.tau*(times_array-self.t0)**2.)
            return  h0*np.sin(phi)*np.exp(-(times_array-self.t0)**2/(2*self.tau)),times_array[0]

        else:
            #print("computing with the given snr")
            h0 = 10 ** (-20)

            phi = 2. * np.pi * self.f0 * (times_array - self.t0) * (1. - 3 * self.tau * (times_array - self.t0) ** 2.)
            aux = h0 * np.sin(phi) * np.exp(-(times_array - self.t0) ** 2 / (2 * self.tau))

            snr_test = GWNoiseSource.SNR(aux, times_array, self.f_asd_tab, self.asd_tab)

            AUX= (self.SNR / snr_test) * aux

            return AUX,times_array[0]

