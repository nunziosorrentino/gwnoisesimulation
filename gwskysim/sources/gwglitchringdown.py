import numpy as np

from gwskysim.sources.gwnoisesource import GWNoiseSource


class GWGlitchRingdown(GWNoiseSource): #ringdown

    def __init__(self,name,t0,tau, f0,SNR=None):
        GWNoiseSource.__init__(self, name, t0 ,tau,SNR)
        self.f0 = f0

    def __call__(self, times_array):
        times_array = times_array[times_array <= self.t0 + 5 * self.tau]
        times_array = times_array[times_array >= self.t0 - 5 * self.tau]
        if np.array(times_array).size==0:
            return None         

        if self.SNR is None:
        #if self.SNR is None and self.detector is None :
            h0 = np.random.uniform(10 ** (-21), 9.9 * 10 ** (-21))

            Q = self.tau * ((2 ** 0.5) * np.pi * self.f0)
            #aux = h0 * np.exp(-((times_array - self.t0) ** 2) / (2 * self.tau ** 2))
            l=len(times_array)
            #print(l)

            AuX2 = np.zeros(l)

            for i in range(l):
                #print(times_array[i])

                if times_array[i]> self.t0:

                    aux = h0 * np.exp(-((times_array[i] - self.t0) ** 2) / (2 * self.tau ** 2))
                    AuX2[i]=np.sin(2 * np.pi * self.f0 * (times_array[i] - self.t0)) * aux


            return AuX2,times_array[0]

        else:
            h0=10**-21
            Q = self.tau * ((2 ** 0.5) * np.pi * self.f0)
            # aux = h0 * np.exp(-((times_array - self.t0) ** 2) / (2 * self.tau ** 2))
            l = len(times_array)
            #print(l)

            AuX2 = np.zeros(l)

            for i in range(l):
                # print(times_array[i])

                if times_array[i] > self.t0:

                    aux = h0 * np.exp(-((times_array[i] - self.t0) ** 2) / (2 * self.tau ** 2))
                    AuX2[i] = np.sin(2 * np.pi * self.f0 * (times_array[i] - self.t0)) * aux
            snr_test = GWNoiseSource.SNR(AuX2, times_array, self.f_asd_tab, self.asd_tab)
            return (self.SNR / snr_test) * AuX2,times_array[0]
